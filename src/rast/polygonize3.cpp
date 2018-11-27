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

std::atomic<long> __fid(0);		// A feature ID for geometries.
std::mutex __gmtx;  			// For the geoms map.
std::mutex __fmtx; 			// For the final list.
std::mutex __omtx;  			// For OGR writes.
std::condition_variable __cv;	// For waiting on the queue.

// Write to the file from the map of geometry lists.
// On each loop, extracts a single finalized poly ID and loads those polys for unioning.
void _writeToFile(std::unordered_map<int, std::vector<Polygon*> >* geoms, std::set<int>* finalIds,
		OGRLayer* layer, GeometryFactory::unique_ptr* fact, GEOSContextHandle_t* gctx,
		bool removeHoles, bool removeDangles, bool* running, bool* cancel) {

	std::vector<Polygon*> polys;
	Geometry* geom;
	int id;

	while(!*cancel && (*running || !finalIds->empty())) {

		// Get an ID and the list of polys from the queue.
		{
			std::unique_lock<std::mutex> lk(__fmtx);
			// Wait for a notification if the queue is empty.
			while(!*cancel && *running && finalIds->empty())
				__cv.wait(lk);
			// If the wakeup is spurious, skip.
			if(finalIds->empty())
				continue;
			// Get the ID.
			id = *(finalIds->begin());
			finalIds->erase(id);
		}
		// Get the list of polys and remove the list from the map.
		{
			std::lock_guard<std::mutex> lk(__gmtx);
			// If the geoms list is still here, grab it.
			if(geoms->find(id) == geoms->end())
				continue;
			polys = geoms->at(id);
			geoms->erase(id);
		}

		if(polys.empty() || *cancel)
			continue;

		// Union the polys.
		geom = geos::operation::geounion::CascadedPolygonUnion::CascadedPolygonUnion::Union(&polys);
		for(Polygon* p : polys)
			delete p;
		polys.clear();

		// If we're removing dangles, throw away all but the
		// largest single polygon. If it was originally a polygon, there are no dangles.
		if(removeDangles && geom->getNumGeometries() > 1) {
			size_t idx = 0;
			double area = 0;
			for(size_t i = 0; i < geom->getNumGeometries(); ++i) {
				const Geometry* p = geom->getGeometryN(i);
				double a = p->getArea();
				if(a > area) {
					area = a;
					idx = i;
				}
			}
			Geometry *g = geom->getGeometryN(idx)->clone(); // Force copy.
			delete geom;
			geom = g;
		}

		// If we're removing holes, extract the exterior rings of all constituent polygons.
		if(removeHoles) {
			std::vector<Geometry*>* geoms0 = new std::vector<Geometry*>();
			for(size_t i = 0; i < geom->getNumGeometries(); ++i) {
				const Polygon* p = dynamic_cast<const Polygon*>(geom->getGeometryN(i));
				const LineString* l = p->getExteriorRing();
				LinearRing* r = (*fact)->createLinearRing(l->getCoordinates());
				geoms0->push_back((*fact)->createPolygon(r, nullptr));
			}
			Geometry* g = (*fact)->createMultiPolygon(geoms0); // Do not copy -- take ownership.
			delete geom;
			geom = g;
		}

		// If the result is not a multi, make it one.
		if(geom->getGeometryTypeId() != GEOS_MULTIPOLYGON) {
			std::vector<Geometry*>* gs = new std::vector<Geometry*>();
			gs->push_back(geom);
			geom = (*fact)->createMultiPolygon(gs);
		}

		if(!geom)
			g_runerr("Null geometry.");

		// Create and write the OGR geometry.
		OGRGeometry* ogeom = OGRGeometryFactory::createFromGEOS(*gctx, (GEOSGeom) geom);
		OGRFeature feat(layer->GetLayerDefn());
		feat.SetGeometry(ogeom);
		feat.SetField("id", (GIntBig) id);
		feat.SetFID(++__fid);

		// Write to the output file.
		int err;
		{
			std::lock_guard<std::mutex> lk(__omtx);
			err = layer->CreateFeature(&feat);
		}

		// Delete the polys.
		delete ogeom;
		delete geom;

		if(OGRERR_NONE != err)
			g_runerr("Failed to add geometry.");
	}
}

// Produce a rectangular polygon from the four corners.
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

// Make an OGR database to write the polygons to.
std::tuple<GDALDataset*, OGRLayer*> _makeDataset(const std::string& filename, const std::string& driver, const std::string& layerName, OGRSpatialReference* sr) {
	// Get the vector driver.
	GDALDriver *drv = GetGDALDriverManager()->GetDriverByName(driver.c_str());
	if(!drv)
		g_runerr("Failed to find driver for " << driver << ".");

	// Create an output dataset for the polygons.
	char** dopts = NULL;
	if(Util::lower(driver) == "sqlite")
		dopts = CSLSetNameValue(dopts, "SPATIALITE", "YES");
	GDALDataset* ds = drv->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, dopts);
	CPLFree(dopts);
	if(!ds)
		g_runerr("Failed to create dataset " << filename << ".");

	// Create the layer.
	char** lopts = NULL;
	if(Util::lower(driver) == "sqlite")
		lopts = CSLSetNameValue(lopts, "FORMAT", "SPATIALITE");

	OGRLayer* layer = ds->CreateLayer(layerName.c_str(), sr, wkbMultiPolygon, lopts);
	CPLFree(lopts);

	if(!layer)
		g_runerr("Failed to create layer " << layerName << ".");

	// There's only one field -- an ID.
	OGRFieldDefn field( "id", OFTInteger);
	layer->CreateField(&field);

	return std::make_tuple(ds, layer);
}

// Produce a status message from the current and total row counds.
std::string _rowStatus(int r, int rows) {
	std::stringstream ss;
	ss << "Polygonizing row " << (r + 1) << " of " << rows;
	return ss.str();
}

void Grid::polygonize(const std::string& filename, const std::string& layerName,
		const std::string& driver, const std::string& projection, int band, bool removeHoles, bool removeDangles,
		const std::string& mask, int maskBand, int threads,
		bool& cancel, Status& status) {

	if(!props().isInt())
		g_runerr("Only int rasters can be polygonized.");

	if(threads < 1)
		threads = 1;

	std::unique_ptr<MemRaster> src;
	std::unique_ptr<MemRaster> msk;

	// Create a memory-mapped version of the input raster.
	{
		GridProps mprops(props());
		mprops.setBands(1);
		src.reset(new MemRaster(mprops, true));
		writeTo(*src, mprops.cols(), mprops.rows(), 0, 0, 0, 0, band, 1);
	}

	// Extract some grid properties.
	const GridProps& props = src->props();
	int cols = props.cols();
	int rows = props.rows();
	double resX = props.resolutionX();
	double resY = props.resolutionY();
	const Bounds& bounds = props.bounds();

	// Create a buffer for the row.
	GridProps rowProps(props);
	rowProps.setSize(cols, 1);
	MemRaster buf(rowProps, false);

	// The starting corner coordinates.
	double startX = resX > 0 ? bounds.minx() : bounds.maxx();
	double startY = resY > 0 ? bounds.miny() : bounds.maxy();

	// Perturbation values for unioning geometries.
	double epsX = resX > 0 ? std::numeric_limits<double>::min() : -std::numeric_limits<double>::min();
	double epsY = resY > 0 ? std::numeric_limits<double>::min() : -std::numeric_limits<double>::min();

	// Create the output dataset
	GEOSContextHandle_t gctx = OGRGeometry::createGEOSContext();
	GeometryFactory::unique_ptr fact = GeometryFactory::create(new PrecisionModel());
	OGRSpatialReference* sr = nullptr;
	if(!projection.empty())
		sr = new OGRSpatialReference(projection.c_str());

	std::unique_ptr<GDALDataset> ds;
	std::unique_ptr<OGRLayer> layer;
	{
		std::tuple<GDALDataset*, OGRLayer*> t = _makeDataset(filename, driver, layerName, sr);
		ds.reset(std::get<0>(t));
		layer.reset(std::get<1>(t));
	}

	// Data containers.
	std::unordered_map<int, std::vector<Polygon*> > geoms;
	std::set<int> activeIds;
	std::set<int> finalIds;

	// Thread control features.
	bool running = true;
	std::list<std::thread> ths;

	// Start output threads.
	for(int i = 0; i < threads; ++i)
		ths.emplace_back(&_writeToFile, &geoms, &finalIds, layer.get(), &fact, &gctx, removeHoles, removeDangles, &running, &cancel);

	// Process raster.
	for(int r = 0; r < rows; ++r) {

		status.update((float) r / rows, _rowStatus(r, rows));

		// Load the row buffer.
		src->writeTo(buf, cols, 1, 0, r, 0, 0, 1, 1);

		// Initialize the corner coordinates.
		double x0 = startX;
		double y0 = startY + r * resY;
		double x1 = x0;
		double y1 = y0 + resY + epsY; // Perturbation for intersection to work.

		// For tracking cell values.
		int v0 = buf.getInt(0, 0, 1);
		int v1 = -1;

		// Reset the list of IDs extant in the current row.
		activeIds.clear();

		for(int c = 1; c < cols; ++c) {

			// If the current cell value differs from the previous one...
			if((v1 = buf.getInt(c, 0, 1)) != v0) {
				// Update the right x coordinate.
				x1 = startX + c * resX + epsX;
				// If the value is a valid ID, create and the geometry and save it for writing.
				if(v0 > 0) {
					geoms[v0].push_back(_makeGeom(x0, y0, x1, y1, &fact));
					activeIds.insert(v0);
				}
				// Update values for next loop.
				v0 = v1;
				x0 = x1;
			}
		}

		// IDs that are in the geoms array and not in the current row are ready to be finalized.
		{
			std::lock_guard<std::mutex> lk0(__gmtx);
			std::lock_guard<std::mutex> lk1(__fmtx);
			for(const auto& it : geoms) {
				if(activeIds.find(it.first) == activeIds.end())
					finalIds.insert(it.first);
			}
		}
		__cv.notify_all();

		while(!finalIds.empty()) {
			//std::this_thread::sleep_for(std::chrono::duration<double, std::milli>(1));
			std::this_thread::yield();
		}
	}

	// Finalize all remaining geometries.
	{
		std::lock_guard<std::mutex> lk0(__gmtx);
		std::lock_guard<std::mutex> lk1(__fmtx);
		for(auto& it : geoms)
			finalIds.insert(it.first);
	}

	// Let the threads shut down when they run out of geometries.
	running = false;
	__cv.notify_all();

	for(std::thread& th : ths) {
		if(th.joinable())
			th.join();
	}

	// Release the layer -- GDAL will take care of it. But close the dataset so that can happen.
	layer.release();
	GDALClose(ds.release());

	status.update(1.0f, "Finished polygonization");

}

