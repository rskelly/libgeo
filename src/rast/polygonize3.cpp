/*
 * polygonize3.cpp
 *
 *  Created on: Nov 21, 2018
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

std::atomic<long> __fid(0);
std::mutex __gmtx;  // For the geoms map.
std::mutex __idmtx; // For the final list.
std::mutex __omtx;  // For OGR writes.
std::condition_variable __cv;

std::string _finishStatus(int c, int final) {
	std::stringstream ss;
	ss << "Finishing " << c << " of " << final << " polygons";
	return ss.str();
}

void _writeToFile(std::unordered_map<int, std::vector<Polygon*> >* geoms, std::set<int>* finalIds, OGRLayer* layer,
		GeometryFactory::unique_ptr* fact, GEOSContextHandle_t* gctx, bool* running, Status* status, int* finalCount) {

	while(*running || !finalIds->empty()) {

		if(!*running)
			status->update(0.5 + (float) finalIds->size() / *finalCount, _finishStatus(finalIds->size(), *finalCount));

		std::vector<Polygon*> polys;
		Geometry* geom;
		int id;

		// Get an ID and the list of polys from the queue.
		{
			std::unique_lock<std::mutex> lk(__idmtx);
			// Wait for a notification if the queue is empty.
			while(*running && finalIds->empty())
				__cv.wait(lk);
			if(finalIds->empty())
				continue;
			// Get the ID.
			id = *(finalIds->begin());
			finalIds->erase(id);
		}
		{
			std::lock_guard<std::mutex> lk(__gmtx);
			// If the geoms list is still here, grab it.
			if(geoms->find(id) == geoms->end())
				continue;
			polys = geoms->at(id);
			geoms->erase(id);
		}

		if(polys.empty())
			continue;
		// Union the polys. If the result is not a multi, make it one.
		geom = geos::operation::geounion::CascadedPolygonUnion::CascadedPolygonUnion::Union(&polys);
		polys.clear();
		if(geom->getGeometryTypeId() != GEOS_MULTIPOLYGON) {
			std::vector<Geometry*>* gs = new std::vector<Geometry*>();
			gs->push_back(geom);
			geom = (*fact)->createMultiPolygon(gs);
		}

		// Create and write the OGR geometry.
		OGRGeometry* ogeom = OGRGeometryFactory::createFromGEOS(*gctx, (GEOSGeom) geom);
		if(!geom)
			g_runerr("Null geometry.");
		OGRFeature feat(layer->GetLayerDefn());
		feat.SetGeometry(ogeom);
		feat.SetField("id", (GIntBig) id);
		feat.SetFID(++__fid);

		int err;
		{
			std::lock_guard<std::mutex> lk(__omtx);
			err = layer->CreateFeature(&feat);
			//std::cout << "poly " << std::this_thread::get_id() << "\n";
		}

		// Delete the polys.
		delete ogeom;
		for(Polygon* p : polys)
			delete p;

		if(OGRERR_NONE != err)
			g_runerr("Failed to add geometry.");
	}

	std::cout << "Finished " << std::this_thread::get_id() << "\n";
}

Polygon* _makeGeom(double x0, double y0, double x1, double y1, GeometryFactory::unique_ptr* fact) {
	// Build the geometry.
	CoordinateSequence* seq = (*fact)->getCoordinateSequenceFactory()->create(5, 2);
	seq->setAt(Coordinate(x0, y0), 0);
	seq->setAt(Coordinate(x0, y1), 1);
	seq->setAt(Coordinate(x1, y1), 2);
	seq->setAt(Coordinate(x1, y0), 3);
	seq->setAt(Coordinate(x0, y0), 4);
	LinearRing* ring = (*fact)->createLinearRing(seq);
	return (*fact)->createPolygon(ring, NULL);
}

std::tuple<GDALDataset*, OGRLayer*> _makeDataset(const std::string& filename, const std::string& driver, const std::string& layerName, int srid) {
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
	if(!layer) {
		g_runerr("Failed to create layer " << layerName << ".");
	}

	// There's only one field -- an ID.
	OGRFieldDefn field( "id", OFTInteger);
	layer->CreateField(&field);

	return std::make_tuple(ds, layer);
}

std::string _rowStatus(int r, int rows) {
	std::stringstream ss;
	ss << "Row " << (r + 1) << " of " << rows;
	return ss.str();
}

void Grid::polygonize(const std::string& filename, const std::string& layerName,
		const std::string& driver, int srid, int band, bool removeHoles, bool removeDangles,
		const std::string& mask, int maskBand, int threads,
		bool& cancel, Status& status) {

	if(!props().isInt())
		g_runerr("Only int rasters can be polygonized.");

	std::unique_ptr<MemRaster> src;
	std::unique_ptr<MemRaster> msk;

	{
		GridProps mprops(props());
		mprops.setBands(1);
		src.reset(new MemRaster(mprops, true));
		writeTo(*src, mprops.cols(), mprops.rows(), 0, 0, 0, 0, band, 1);
	}

	const GridProps& props = src->props();
	int cols = props.cols();
	int rows = props.rows();
	double resX = props.resolutionX();
	double resY = props.resolutionY();
	const Bounds& bounds = props.bounds();
	GridProps rowProps(props);
	rowProps.setSize(cols, 1);
	MemRaster buf(rowProps, false);

	double startX = resX > 0 ? bounds.minx() : bounds.maxx();
	double startY = resY > 0 ? bounds.miny() : bounds.maxy();
	double epsX = resX > 0 ? std::numeric_limits<double>::min() : -std::numeric_limits<double>::min();
	double epsY = resY > 0 ? std::numeric_limits<double>::min() : -std::numeric_limits<double>::min();

	GEOSContextHandle_t gctx = OGRGeometry::createGEOSContext();
	GeometryFactory::unique_ptr fact = GeometryFactory::create(new PrecisionModel());

	std::unique_ptr<GDALDataset> ds;
	std::unique_ptr<OGRLayer> layer;

	{
		std::tuple<GDALDataset*, OGRLayer*> t = _makeDataset(filename, driver, layerName, srid);
		ds.reset(std::get<0>(t));
		layer.reset(std::get<1>(t));
	}

	std::unordered_map<int, std::vector<Polygon*> > geoms;
	std::set<int> activeIds;
	std::set<int> finalIds;
	bool running = true;
	int finalCount = 0;
	std::list<std::thread> ths;

	for(int i = 0; i < 4; ++i)
		ths.emplace_back(&_writeToFile, &geoms, &finalIds, layer.get(), &fact, &gctx, &running, &status, &finalCount);

	for(int r = 0; r < rows; ++r) {

		status.update((float) r / rows / 2, _rowStatus(r, rows));

		src->writeTo(buf, cols, 1, 0, r, 0, 0, 1, 1);

		double x0 = startX;
		double y0 = startY + r * resY;
		double x1 = x0;
		double y1 = y0 + resY + epsY; // Perturbation for intersection to work.

		int v0 = buf.getInt(0, 0, 1);
		int v1 = -1;

		activeIds.clear();

		for(int c = 1; c < cols; ++c) {
			if((v1 = buf.getInt(c, 0, 1)) != v0) {

				x1 = startX + c * resX + epsX;

				if(v0 > 0) {
					geoms[v0].push_back(_makeGeom(x0, y0, x1, y1, &fact));
					activeIds.insert(v0);
				}

				v0 = v1;
				x0 = x1;
			}
		}

		{
			std::lock_guard<std::mutex> lk(__gmtx);
			for(auto& it : geoms) {
				if(activeIds.find(it.first) == activeIds.end())
					finalIds.insert(it.first);
			}
		}
		__cv.notify_all();

	}

	{
		std::lock_guard<std::mutex> lk(__gmtx);
		for(auto& it : geoms)
			finalIds.insert(it.first);
	}

	running = false;
	finalCount = finalIds.size();

	__cv.notify_all();

	for(std::thread& th : ths) {
		if(th.joinable())
			th.join();
	}

	status.update(1.0f, "Finished");

}

