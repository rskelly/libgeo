/*
 * polygonize_ff.cpp
 *
 *  Created on: Nov 16, 2018
 *      Author: rob
 */

#include <thread>
#include <condition_variable>
#include <tuple>

#include <gdal.h>
#include <gdal_alg.h>
#include <ogr_feature.h>
#include <ogrsf_frmts.h>

#include <geos/geom/PrecisionModel.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/CoordinateSequenceFactory.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/MultiPolygon.h>
#include <geos/geom/Point.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/operation/union/CascadedPolygonUnion.h>

#include "raster.hpp"

using namespace geo::raster;
using namespace geos::geom;
using namespace geos::operation::geounion;

using namespace geo::raster;

int __fid;
std::mutex __mtx;
std::mutex __ogrmtx;

class Poly {
public:
	std::vector<Polygon*> geoms;
	int id;
	int minx;
	int miny;
	int maxx;
	int maxy;

	Poly() :
		id(0),
		minx(99999999), miny(99999999),
		maxx(-99999999), maxy(-99999999) {}

	void update(int id, int minx, int miny, int maxx, int maxy, Polygon* geom) {
		this->id = id;
		if(minx < this->minx) this->minx = minx;
		if(miny < this->miny) this->miny = miny;
		if(maxx > this->maxx) this->maxx = maxx;
		if(maxy > this->maxy) this->maxy = maxy;
		geoms.push_back(geom);
	}

	void generate(OGRLayer* layer, GeometryFactory::unique_ptr& fact, GEOSContextHandle_t& gctx) {
		int fid;
		{
			std::lock_guard<std::mutex> lk(__mtx);
			fid = ++__fid;

		}

		Geometry* geom = geos::operation::geounion::CascadedPolygonUnion::CascadedPolygonUnion::Union(&geoms);
		if(geom->getGeometryTypeId() != GEOS_MULTIPOLYGON) {
			std::vector<Geometry*>* gs = new std::vector<Geometry*>();
			gs->push_back(geom->clone());
			geom = fact->createMultiPolygon(gs);
		}

		OGRGeometry* ogeom = OGRGeometryFactory::createFromGEOS(gctx, (GEOSGeom) geom);
		if(!geom)
			g_runerr("Null geometry.");
		OGRFeature feat(layer->GetLayerDefn());
		feat.SetGeometry(ogeom);
		feat.SetField("id", (GIntBig) id);
		feat.SetFID(fid);
		delete ogeom;

		int err;
		{
			std::lock_guard<std::mutex> lk(__ogrmtx);
			err = layer->CreateFeature(&feat);
		}
		if(OGRERR_NONE != err)
			g_runerr("Failed to add geometry.");
	}

	~Poly() {
		//for(Geometry* g : geoms)
		//	delete g;
	}
};


// 1 - divide grid up in to "tiles"
// 2 - in parallel, flood fill each tile, extracting filled regions as polygons
// 3 - for each polygon whose bounds correspond to a tile edge, add to btree
// 4 - search btree for touching polygons with the same id, merge if found
// 5 - output
void Grid::polygonize(const std::string &filename, const std::string &layerName,
                const std::string &driver, int srid, int band, bool removeHoles, bool removeDangles,
				const std::string& mask, int maskBand, int threads,
				bool& cancel, geo::util::Status& status) {

	const GridProps& props = this->props();
	int gcols = props.cols();
	int grows = props.rows();
	double resx = props.resolutionX();
	double resy = props.resolutionY();

	int tcols = gcols;
	int trows = grows;

	if(gcols > 200 && grows > 200) {
		tcols = 100;
		trows = 100;
	}

	GridProps tprops(props);
	tprops.setSize(tcols, trows);

	GEOSContextHandle_t gctx = OGRGeometry::createGEOSContext();
	GeometryFactory::unique_ptr fact = GeometryFactory::create(new PrecisionModel());

	// Get the vector driver.
	GDALDriver *drv = GetGDALDriverManager()->GetDriverByName(driver.c_str());
	if(!drv)
		g_runerr("Failed to find driver for " << driver << ".");

	// Create an output dataset for the polygons.
	char **dopts = NULL;
	if(Util::lower(driver) == "sqlite")
		dopts = CSLSetNameValue(dopts, "SPATIALITE", "YES");
	GDALDataset* ds = drv->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, dopts);
	CPLFree(dopts);
	if(!ds)
		g_runerr("Failed to create dataset " << filename << ".");

	// Create the layer.
	OGRSpatialReference sr;
	if(srid > 0)
		sr.importFromEPSG(srid);

	char **lopts = NULL;
	if(Util::lower(driver) == "sqlite")
		lopts = CSLSetNameValue(lopts, "FORMAT", "SPATIALITE");
	OGRLayer* layer = ds->CreateLayer(layerName.c_str(), srid > 0 ? &sr : nullptr, wkbMultiPolygon, lopts);
	CPLFree(lopts);
	if(!layer)
		g_runerr("Failed to create layer " << layerName << ".");

	// There's only one field -- an ID.
	OGRFieldDefn field( "id", OFTInteger);
	layer->CreateField(&field);

	std::vector<Poly> fpolys;

	for(int tr = 0; tr < grows / trows; tr += trows) {
		for(int tc = 0; tc < gcols / tcols; tc += tcols) {

			int cols = tc + tcols >= gcols ? tc + tcols - gcols : tc + tcols;
			int rows = tr + trows >= grows ? tr + trows - grows : tr + trows;
			int minc, minr, maxc, maxr, area;

			MemRaster mr(tprops);
			mr.fillInt(0, band);

			std::unordered_map<int, Poly> polys;

			this->writeTo(mr, tcols, trows, tc, tr, 0, 0, band, 1);

			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c) {

					int v = mr.getInt(c, r, 1);

					if(v != -1 && v != 0) { // TODO Configure to allow keeping/skipping zero areas.

						TargetFillOperator<int, int> op(&mr, 1, 1, v, -2);
						mr.floodFill(c, r, op, false, &minc, &minr, &maxc, &maxr, &area);

						if(area > 0) {

							for(int rr = minr; rr <= maxr; ++rr) {

								double x0 = 0, y0 = 0, x1 = 0, y1 = 0;
								int cells = 0;

								for(int cc = minc; cc <= maxc; ++cc) {

									int vv = mr.getInt(cc, rr, 1);

									if(vv == -2) {
										if(!cells) {
											x0 = props.toX(cc + tr) + (resx > 0 ? 0.0001 : -0.0001); // TODO: Perturbation should be configurable.
											y0 = props.toY(rr + tr) + (resy > 0 ? 0.0001 : -0.0001);
										}
										++cells;
									}

									if(cells && (cc == maxc || vv != -2)) {
										x1 = props.toX(cc + tr) + (resx > 0 ? 0.0001 : -0.0001); // TODO: Perturbation should be configurable.
										y1 = props.toY(rr + tr) + resy + (resy > 0 ? 0.0001 : -0.0001);
										// Build the geometry.
										CoordinateSequence* seq = fact->getCoordinateSequenceFactory()->create(5, 2);
										seq->setAt(Coordinate(x0, y0), 0);
										seq->setAt(Coordinate(x0, y1), 1);
										seq->setAt(Coordinate(x1, y1), 2);
										seq->setAt(Coordinate(x1, y0), 3);
										seq->setAt(Coordinate(x0, y0), 4);
										LinearRing* ring = fact->createLinearRing(seq);
										Polygon* geom = fact->createPolygon(ring, NULL);
										polys[v].update(v, minc, minr, maxc, maxr, geom);
										cells = 0;
									}
								}
							}

						}

					}

				}
			}

			if(OGRERR_NONE != layer->StartTransaction())
				g_runerr("Failed to start transaction.");

			for(auto& it : polys) {
				Poly& p = it.second;
				if(p.minx <= 0 || p.maxx >= cols || p.miny <= 0 || p.maxy >= rows) {
					fpolys.push_back(std::move(p));
				} else {
					p.generate(layer, fact, gctx);
				}
			}

			if(OGRERR_NONE != layer->CommitTransaction())
				g_runerr("Failed to commit transaction.");
		}
	}

	if(OGRERR_NONE != layer->StartTransaction())
		g_runerr("Failed to start transaction.");

	for(Poly& p : fpolys) {
		p.generate(layer, fact, gctx);
	}

	if(OGRERR_NONE != layer->CommitTransaction())
		g_runerr("Failed to commit transaction.");


}


