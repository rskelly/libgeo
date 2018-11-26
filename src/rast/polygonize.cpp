/*
 * rast_polygonize.cpp
 *
 *  Created on: Feb 23, 2018
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

/**
 * Class used in polygonization
 * implementation.
 */
class Poly {
private:


	Geometry* m_geom;				///<! The final geometry.
	std::vector<Geometry*> m_geoms;	///<! The list of consitiuent geometries.
	size_t m_id;					///<! Unique geometry ID.
	long m_minRow;					///<! The minimum row index from the source raster.
	long m_maxRow;					///<! The maximum row index from the source raster.
	bool m_finalized;				///<! True if this poly has already been finalized.
	bool m_collapsed;				///<! True if the geometries have been collapsed after the most recent update.

	/**
	 * Add the list of polygons to this instance with the rows covered by them.
	 *
	 * @param u The vector of Polygons.
	 * @param minRow The lowest row index occupied by this geometry.
	 * @param maxRow The highest row index occupied by this geometry.
	 */
	void update(std::vector<Polygon*>& u, int minRow, int maxRow) {
		for (Polygon* p : u)
			m_geoms.push_back(static_cast<Geometry*>(p));
		update(minRow, maxRow);
		m_collapsed = false;
	}


	/**
	 * Updates the row indices covered by this polygon.
	 *
	 * @param minRow The lowest row index occupied by this geometry.
	 * @param maxRow The highest row index occupied by this geometry.
	 */
	void update(int minRow, int maxRow) {
		if(minRow < m_minRow) m_minRow = minRow;
		if(maxRow > m_maxRow) m_maxRow = maxRow;
	}

	/**
	 * Unions the geometries in group and adds the result to geoms.
	 * If an invalid geometry is found in the output, an exception is thrown.
	 *
	 * @param group A vector of input polygons.
	 * @param geoms A vector of output geometries.
	 */
	void _collapse(std::vector<Polygon*>* group, std::vector<Geometry*>* geoms) {
		if(!group->empty()) {
			Geometry* p = geos::operation::geounion::CascadedPolygonUnion::CascadedPolygonUnion::Union(group);
			for(size_t i = 0; i < p->getNumGeometries(); ++i) {
				const Geometry* g = p->getGeometryN(i);
				if(!g || g->getGeometryTypeId() != GEOS_POLYGON)
					g_runerr("Null or invalid polygon.");
				geoms->push_back(g->clone());
			}
		}
	}

public:

	/**
	 * Create a polygon with the given ID and initial geometry,
	 * for the given rows. The factory is shared by all instances.
	 *
	 * @param id A unique ID.
	 * @param geom A Geometry. This instance takes ownership.
	 * @param minRow The minimum row.
	 * @param maxRow The maximum row.
	 */
	Poly(size_t id, Geometry* geom, long minRow, long maxRow) :
		m_geom(nullptr),
		m_id(id),
		m_minRow(std::numeric_limits<long>::max()),
		m_maxRow(std::numeric_limits<long>::min()),
		m_finalized(false),
		m_collapsed(false) {
	}

	/**
	 * Create a polygon with the given ID and initial geometry,
	 * for the given rows. The factory is shared by all instances.
	 *
	 * @param id A unique ID.
	 */
	Poly(size_t id) : Poly(id, nullptr, 0, 0) {}

	/**
	 * Create an empty polygon.
	 */
	Poly() : Poly(0) {}

	/**
	 * Add the contents of the Poly to this instance.
	 * The Poly's geometry is cloned.
	 *
	 * @param p The Poly to add.
	 */
	void update(Poly& p) {
		p.collapse();
		update(p.geom()->clone(), p.minRow(), p.maxRow());
	}

	/**
	 * Add the polygon to this instance with the
	 * rows covered by it. This Poly now owns the geometry.
	 *
	 * @param u The Polygon to add.
	 * @param minRow The lowest row index occupied by this geometry.
	 * @param maxRow The highest row index occupied by this geometry.
	 */
	void update(Geometry* u, int minRow, int maxRow) {
		m_geoms.push_back(u);
		update(minRow, maxRow);
		m_collapsed = false;
	}

	/**
	 * Add and ID and polygon to this instance with the
	 * rows covered by it. This Poly now owns the geometry.
	 *
	 * @param u The Polygon to add.
	 * @param minRow The lowest row index occupied by this geometry.
	 * @param maxRow The highest row index occupied by this geometry.
	 */
	void update(size_t id, Geometry* u, int minRow, int maxRow) {
		m_id = id;
		m_geoms.push_back(u);
		update(minRow, maxRow);
		m_collapsed = false;
	}

	/**
	 * Returns true if the range of rows given by start and end was finalized
	 * within the given thread, or by a thread whose block is completed.
	 * The checked range includes one row above and one below,
	 * which is required to guarantee that a polygon is completed.
	 *
	 * @param rowFinished The last row read.
	 * @return True if all the rows are finished.
	 */
	bool isRangeFinalized(std::vector<bool> rowsFinished) const {
		for(long i = std::max(0L, m_minRow - 1); i < std::max((long) rowsFinished.size(), m_maxRow + 1); ++i) {
			if(!rowsFinished[i])
				return false;
		}
		return true;
	}

	/**
	 * Return the pointer to the unioned polygon.
	 *
	 * @return The pointer to the unioned polygon.
	 */
	const Geometry* geom() const {
		return m_geom;
	}

	/**
	 * Return the list of constituent geometries.
	 *
	 * @return The list of constituent geometries.
	 */
	const std::vector<Geometry*>& geoms() const {
		return m_geoms;
	}

	/**
	 * Iteratively union the constituent geometries into a
	 * smaller list of single polygons (ideally).
	 */
	void collapse() {
		if(m_collapsed || m_geoms.empty())
			return;
		size_t geomCount = m_geoms.size();
		int change = 1;
		while(change) {
			std::vector<Polygon*> polys;
			for(Geometry* geom : m_geoms) {
				for(size_t i = 0; i < geom->getNumGeometries(); ++i) {
					const Geometry* g = geom->getGeometryN(i);
					if(!g || g->getGeometryTypeId() != GEOS_POLYGON)
						g_runerr("Null or invalid polygon.");
					polys.push_back(dynamic_cast<Polygon*>(g->clone()));
				}
			}
			for(Geometry* g : m_geoms)
				delete g;
			m_geoms.clear();

			//g_debug("Collapsing " << polys.size() << " polys")

			// The max group size is 1024; min 2.
			int groupSize = g_max(2, g_min(1024, polys.size() / 16));
			// If the group size is 2, just use one group.
			if(groupSize == 2)
				groupSize = polys.size();

			// Number of threads
			// TODO: Configurable.
			size_t tc = 4;
			// The number of items in one group.
			size_t groupCount = (int) std::ceil((float) polys.size() / groupSize);
			// Container for thread instances.
			std::vector<std::thread> threads(tc);
			// The output geometry groups.
			std::vector<std::vector<Geometry*> > output(tc);
			// The input polygon groups.
			std::vector<std::vector<Polygon*> > groups(tc);

			size_t i, j;
			for(size_t g = 0; g < groupCount; g += tc) {
				for(i = 0; i < tc; ++i) {
					size_t begin = (g + i) * groupSize;
					if(begin >= polys.size())
						break;
					size_t end = std::min((g + i + 1) * groupSize, polys.size());
					groups[i].clear(); // TODO: More efficient lists.
					output[i].clear();
					groups[i].assign(polys.begin() + begin, polys.begin() + end);
					threads[i] = std::thread(&Poly::_collapse, this, &groups[i], &output[i]);
				}
				for(j = 0; j < i; ++j) {
					threads[j].join();
					m_geoms.insert(m_geoms.end(), output[j].begin(), output[j].end());
				}
			}

			change = geomCount - m_geoms.size();
			geomCount = m_geoms.size();

			//g_debug("Collapse " << change << ", " << geomCount)
		}

		m_collapsed = true;
	}

	/**
	 * Generate the unioned polygon from its parts. Remove dangles and holes if required.
	 *
	 * @param fact The GeometryFactory.
	 * @param removeHoles True, if it is desired that holes in polygons be removed.
	 * @param removeDangles True, if it is desired that single-pixel artifacts, attached at the corners, be removed.
	 */
	void finalize(const GeometryFactory::unique_ptr& fact, bool removeHoles, bool removeDangles) {

		if(m_finalized)
			return;
		m_finalized = true;

		collapse();

		//g_debug("Finalize")

		Geometry* geom = fact->createMultiPolygon(m_geoms);

		// If we're removing dangles, throw away all but the
		// largest single polygon. If it was originally a polygon,
		// there are no dangles.
		if(removeDangles) {
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

		// If we're removing holes, extract the exterior rings
		// of all constituent polygons.
		if(removeHoles) {
			std::vector<Geometry*>* geoms0 = new std::vector<Geometry*>();
			for(size_t i = 0; i < geom->getNumGeometries(); ++i) {
				const Polygon* p = dynamic_cast<const Polygon*>(geom->getGeometryN(i));
				const LineString* l = p->getExteriorRing();
				LinearRing* r = fact->createLinearRing(l->getCoordinates());
				geoms0->push_back(fact->createPolygon(r, nullptr));
			}
			Geometry* g = fact->createMultiPolygon(geoms0); // Do not copy -- take ownership.
			delete geom;
			geom = g;
		}

		m_geom = geom;
	}

	/**
	 * Return the unique ID.
	 *
	 * @return The unique ID.
	 */
	size_t id() const {
		return m_id;
	}

	/**
	 * Return the minimum row index covered by the geometry.
	 *
	 * @return The minimum row index covered by the geometry.
	 */
	int minRow() const {
		return m_minRow;
	}

	/**
	 * Return the maximum row index covered by the geometry.
	 *
	 * @return The maximum row index covered by the geometry.
	 */
	int maxRow() const {
		return m_maxRow;
	}

	~Poly() {
		for(Geometry* geom : m_geoms)
			delete geom;
		if(m_geom)
			delete m_geom;
	}
};

typedef std::unordered_map<size_t, Poly> PolyMap; 	///<! A map for maintaining the collection of Polys.
typedef std::queue<Poly> PolyQueue;					///<! The Poly processing queue.

/**
 * Maintains the set of objects and values required for parallel
 * polygonization of rasters and processing of geometries.
 */
class PolyContext {
private:
	Grid* m_raster;
	int m_band;
	Grid* m_maskRaster;
	int m_maskBand;
	Status* m_status;
	int* m_block;
	bool* m_cancel;
	GDALDataset* m_ds;
	OGRLayer *m_layer;
	GeometryFactory::unique_ptr m_geomFactory;
	GEOSContextHandle_t m_gctx;
	OGRSpatialReference m_sr;
	PolyMap m_polyMap;
	PolyQueue m_polyQueue;
	std::condition_variable m_polyMapCond;
	std::condition_variable m_polyQueueCond;
	std::mutex m_polyMapMtx;    // For the PolyMap.
	std::mutex m_polyQueueMtx;  // For the PolyQueue.
	std::mutex m_finMtx;        // For the finished list.
	std::mutex m_blockMtx;      // For the block counter.
	std::mutex m_fidMtx;        // For the featureId.
	std::mutex m_ogrMtx;        // For the OGRLayer.
	size_t m_featureId;
	int m_bufSize;
	bool m_removeHoles;
	bool m_removeDangles;
	bool m_queueUpdate;
	bool m_readFinish;
	bool m_transferFinish;
	std::vector<bool> m_rowsFinished;

public:
	PolyContext(Grid* raster, Status* status, int* block, bool* cancel,
			int bufSize, int band, bool removeHoles, bool removeDangles, Grid* maskRaster, int maskBand) :
		m_raster(raster),
		m_band(band),
		m_maskRaster(maskRaster),
		m_maskBand(maskBand),
		m_status(status),
		m_block(block),
		m_cancel(cancel),
		m_ds(nullptr),
		m_layer(nullptr),
		m_geomFactory(nullptr),
		m_gctx(nullptr),
		m_featureId(0),
		m_bufSize(bufSize),
		m_removeHoles(removeHoles),
		m_removeDangles(removeDangles),
		m_queueUpdate(false),
		m_readFinish(false),
		m_transferFinish(false) {

		// Size the finished array w/ or w/o mask.
		const GridProps& props = m_raster->props();
		m_rowsFinished.resize(getRows(props, m_maskRaster));
		std::fill(m_rowsFinished.begin(), m_rowsFinished.end(), false);
	}

	~PolyContext() {
		if(m_ds)
			GDALClose(m_ds);
	}

	const GeometryFactory::unique_ptr* geomFactory() const {
		return &m_geomFactory;
	}

	void readFinish() {
		m_readFinish = true;
	}

	void transferFinish() {
		m_transferFinish = true;
	}

	void initOutput(const std::string& driver, const std::string& filename, const std::string& layerName, const std::string& projection) {
		GDALAllRegister();

		// Create the GEOS context and factory.
		m_gctx = OGRGeometry::createGEOSContext();
		m_geomFactory = GeometryFactory::create(new PrecisionModel());

		// Get the vector driver.
		GDALDriver *drv = GetGDALDriverManager()->GetDriverByName(driver.c_str());
		if(!drv)
			g_runerr("Failed to find driver for " << driver << ".");

		// Create an output dataset for the polygons.
		char **dopts = NULL;
		if(Util::lower(driver) == "sqlite")
			dopts = CSLSetNameValue(dopts, "SPATIALITE", "YES");
		m_ds = drv->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, dopts);
		CPLFree(dopts);
		if(!m_ds)
			g_runerr("Failed to create dataset " << filename << ".");

		// Create the layer.
		if(!projection.empty()) {
			char* wkt = projection.c_str();
			m_sr.importFromWkt(&wkt);
		}
		char **lopts = NULL;
		if(Util::lower(driver) == "sqlite")
			lopts = CSLSetNameValue(lopts, "FORMAT", "SPATIALITE");
		m_layer = m_ds->CreateLayer(layerName.c_str(), projection.empty() ? nullptr : &m_sr, wkbMultiPolygon, lopts);
		CPLFree(lopts);
		if(!m_layer) {
			g_runerr("Failed to create layer " << layerName << ".");
		}

		// There's only one field -- an ID.
		OGRFieldDefn field( "id", OFTInteger);
		m_layer->CreateField(&field);

		if(OGRERR_NONE != m_layer->StartTransaction())
			g_runerr("Failed to start transaction.");

	}

	void commitOutput() {
		if(OGRERR_NONE != m_layer->CommitTransaction())
			g_runerr("Failed to commit transation.");
	}

	OGRLayer* layer() {
		return m_layer;
	}

	size_t nextFeatureId() {
		std::lock_guard<std::mutex> lk(m_fidMtx);
		return ++m_featureId;
	}

	void notifyTransfer() {
		m_polyMapCond.notify_all();
	}

	void notifyWrite() {
		m_polyQueueCond.notify_all();
	}

	int getCols(const GridProps& props, Grid* mask) {
		if(mask) {
			const GridProps& mprops = mask->props();
			const Bounds& mbounds = mprops.bounds();
			return std::abs(props.toCol(mbounds.maxx()) - props.toCol(mbounds.minx())) + 1;
		} else {
			return props.cols();
		}
	}

	int getRows(const GridProps& props, Grid* mask) {
		if(mask) {
			const GridProps& mprops = mask->props();
			const Bounds& mbounds = mprops.bounds();
			return std::abs(props.toRow(mbounds.maxy()) - props.toRow(mbounds.miny())) + 1;
		} else {
			return props.rows();
		}
	}

	int getCol(const GridProps& props, Grid* mask) {
		if(mask) {
			const GridProps& mprops = mask->props();
			const Bounds& mbounds = mprops.bounds();
			return props.toCol(mprops.resolutionX() > 0 ? mbounds.minx() : mbounds.maxx());
		} else {
			return 0;
		}
	}

	int getRow(const GridProps& props, Grid* mask) {
		if(mask) {
			const GridProps& mprops = mask->props();
			const Bounds& mbounds = mprops.bounds();
			return props.toRow(mprops.resolutionY() > 0 ? mbounds.miny() : mbounds.maxy());
		} else {
			return 0;
		}
	}

	void polyReadBlocks() {

		// Get cols, rows an nodata.
		const GridProps& props = m_raster->props();
		int col = getCol(props, m_maskRaster);
		int row = getRow(props, m_maskRaster);
		int cols = getCols(props, m_maskRaster);
		int rows = getRows(props, m_maskRaster);
		size_t nd = (size_t) props.nodata();

		// The number of blocks int he raster.
		int blocks = (int) rows / m_bufSize + 1;
		// The current block index.
		int b;

		// Buffer for reading raster.
		GridProps gp(props);
		gp.setSize(cols, m_bufSize);
		MemRaster blockBuf(gp, false);

		// Buffer to hold individual row data.
		std::vector<int> rowBuf(cols);

		while(true) {
			{
				std::lock_guard<std::mutex> lk(m_blockMtx);
				b = *m_block;
				++*m_block;
			}

			// If cancelled or no more blocks, quit.
			if(*m_cancel || b >= blocks)
				break;

			// Work out the height of the current buffer. If it's too small quit.
			int bufHeight = g_min(m_bufSize, rows - b * m_bufSize);
			if(bufHeight < 1)
				break;

			m_status->update((float) b / blocks);

			// Write into the buffer.
			m_raster->writeTo(blockBuf, cols, bufHeight, col, row + b * m_bufSize, 0, 0, m_band);

			// Read over the rows in the buffer.
			for(int rr = 0; rr < bufHeight; ++rr) {

				if(*m_cancel)
					break;

				// The current overall row index.
				int r = b * m_bufSize + rr;

				// Read into the row buffer from the block.
				blockBuf.getIntRow(rr, m_band, rowBuf.data());

				std::list<std::tuple<size_t, Geometry*> > toWrite;

				// Read over the columns in the current row.
				for(int c = 0; c < cols; ++c) {

					// Get the current ID, skip if nodata.
					size_t id0 = rowBuf[c];
					if(id0 == nd || id0 == 0)
						continue;

					// Get the coord of one corner of the polygon.
					double x0 = gp.toX(c + col);
					double y0 = gp.toY(r + row);

					// Scan right...
					while(++c <= cols) {

						// Get the next ID or zero if beyond the edge.
						size_t id1 = c < cols ? rowBuf[c] : 0;

						// If the ID changes, capture and output the polygon.
						if(id0 > 0 && id1 != id0) {

							// Coord of the other corner.
							double x1 = gp.toX(c + col) + (gp.resolutionX() > 0 ? 0.0001 : -0.0001); // TODO: Perturbation should be configurable.
							double y1 = gp.toY(r + row) + gp.resolutionY() + (gp.resolutionY() > 0 ? 0.0001 : -0.0001);

							// Build the geometry.
							CoordinateSequence* seq = m_geomFactory->getCoordinateSequenceFactory()->create(5, 2);
							seq->setAt(Coordinate(x0, y0), 0);
							seq->setAt(Coordinate(x0, y1), 1);
							seq->setAt(Coordinate(x1, y1), 2);
							seq->setAt(Coordinate(x1, y0), 3);
							seq->setAt(Coordinate(x0, y0), 4);
							LinearRing* ring = m_geomFactory->createLinearRing(seq);
							Polygon* geom = m_geomFactory->createPolygon(ring, NULL);

							toWrite.emplace_back(id0, geom);

							--c; // Back up the counter by one, to start with the new ID.
							break;
						}
						id0 = id1;
					}
				}

				// Update the polygon list with the new poly.
				{
					std::lock_guard<std::mutex> lk(m_polyMapMtx);
					for(auto& it : toWrite) {
						size_t id = std::get<0>(it);
						Geometry* geom = std::get<1>(it);
						m_polyMap[id].update(id, geom, r, r);
					}
				}

				// Update the finished rows.
				m_rowsFinished[r] = true;

				// Notify the transfer threads of an update.
				notifyTransfer();
			}

		} // while

		notifyTransfer();
	}

	void polyQueueTransfer() {

		while(!*m_cancel && !(m_readFinish && m_polyMap.empty())) {

			Poly p;
			bool found = false;
			{
				// Wait for a wake-up.
				std::unique_lock<std::mutex> lk(m_polyMapMtx);
				while(!*m_cancel && !m_readFinish && m_polyMap.empty())
					m_polyMapCond.wait(lk);

				// Find a finished Poly.
				for(const auto& it : m_polyMap) {
					Poly& p0 = m_polyMap[it.first];
					if(p0.isRangeFinalized(m_rowsFinished)) {
						p = std::move(p0);
						m_polyMap.erase(it.first);
						found = true;
					}
				}
			}

			// Finalized the finished Poly and transfer to the write queue.
			if(found) {
				p.finalize(m_geomFactory, m_removeHoles, m_removeDangles);
				std::lock_guard<std::mutex> lk(m_polyQueueMtx);
				m_polyQueue.push(std::move(p));
			}

			notifyWrite();
		}

		notifyWrite();
	}

	Poly polyFromPath(int id, const std::vector<double>& path) {
		// Build the geometry.
		CoordinateSequence* seq = m_geomFactory->getCoordinateSequenceFactory()->create((size_t) 0, 2);
		std::vector<Coordinate> coords;
		for(size_t i = 0; i < path.size(); i += 2)
			coords.push_back(Coordinate(path[i], path[i+1]));
		seq->add(&coords, false);
		// Get the ring and make a polygon.
		LinearRing* ring = m_geomFactory->createLinearRing(seq);
		Polygon* geom = m_geomFactory->createPolygon(ring, NULL);
		return Poly(id, geom, 0, 0);
	}

	void enqueuePoly(Poly& poly) {
		poly.finalize(m_geomFactory, m_removeHoles, m_removeDangles);
		{
			std::lock_guard<std::mutex> lk(m_polyQueueMtx);
			m_polyQueue.push(std::move(poly));
		}
		notifyWrite();
	}

	void polyWriteQueue() {

		// Use the first thread in the group for writing polys to the DB.
		// The loop runs as long as the queue isn't "finalized" and is not empty.
		while(!*m_cancel && !(m_transferFinish && m_polyQueue.empty())) {

			Poly p;
			bool found = false;
			{
				// Wait for a wake-up.
				std::unique_lock<std::mutex> lk(m_polyQueueMtx);
				while(!*m_cancel && !m_transferFinish && m_polyQueue.empty())
					m_polyQueueCond.wait(lk);

				// Grab the element from the queue, if there is one.
				if(!m_polyQueue.empty()) {
					p = std::move(m_polyQueue.front());
					m_polyQueue.pop();
					found = true;
				}
			}

			// If no element, loop.
			if(found) {

				// Retrieve the unioned geometry and write it.
				const Geometry* geom = p.geom();

				{
					std::vector<Polygon*> polys;
					for(size_t i = 0; i < geom->getNumGeometries(); ++i) {
						const Geometry* g = geom->getGeometryN(i);
						if(!g || g->getGeometryTypeId() != GEOS_POLYGON)
							g_runerr("Null or invalid polygon.");
						polys.push_back(dynamic_cast<Polygon*>(g->clone()));
					}
					geom = geos::operation::geounion::CascadedPolygonUnion::CascadedPolygonUnion::Union(&polys);
				}

				if(geom->getGeometryTypeId() != GEOS_MULTIPOLYGON) {
					std::vector<Geometry*>* geoms = new std::vector<Geometry*>();
					geoms->push_back(geom->clone());
					geom = m_geomFactory->createMultiPolygon(geoms);
				}

				OGRGeometry* ogeom = OGRGeometryFactory::createFromGEOS(m_gctx, (GEOSGeom) geom);
				if(!geom)
					g_runerr("Null geometry.");
				OGRFeature feat(m_layer->GetLayerDefn());
				feat.SetGeometry(ogeom);
				feat.SetField("id", (GIntBig) p.id());
				feat.SetFID(nextFeatureId());
				delete ogeom;

				int err;
				{
					std::lock_guard<std::mutex> lk(m_ogrMtx);
					err = m_layer->CreateFeature(&feat);
				}
				if(OGRERR_NONE != err)
					g_runerr("Failed to add geometry.");
			}
		}
	}

};

void Grid::polygonize(const std::string& filename, const std::string& layerName,
		const std::string& driver, const std::string& projection, int band, bool removeHoles, bool removeDangles,
		const std::string& mask, int maskBand, int threads,
		bool& cancel, Status& status) {

	if(!props().isInt())
		g_runerr("Only integer rasters can be polygonized.");

	// There need to be at least three threads for reading from the raster(1),
	// writing output (1), and unioning and transferring to the write queue (n - 2).
	if(threads < 3)
		threads = 3;

	// Counter for the current block.
	int block = 0;

	// Set the buffer size.
	int bufSize = std::min(512, props().rows());

	// Remove the original file; some can't be overwritten directly.
	// This will not take care of any auxillary files (e.g. shapefiles)
	Util::rm(filename);

	// Initialize the mask if necessary.
	Grid* maskRaster = nullptr;
	if(!mask.empty())
		maskRaster = new Raster(mask);

	// Set up the shared context.
	PolyContext ctx(this, &status, &block, &cancel, bufSize, band, removeHoles, removeDangles, maskRaster, maskBand);

	// Initialize database.
	ctx.initOutput(driver, filename, layerName, projection);

	// Start the read thread.
	std::thread readT(&PolyContext::polyReadBlocks, &ctx);

	// Start the transfer threads.
	std::thread transferT(&PolyContext::polyQueueTransfer, &ctx);

	// Start the write thread.
	std::thread writeT(&PolyContext::polyWriteQueue, &ctx);

	readT.join();

	ctx.readFinish();
	ctx.notifyTransfer();

	transferT.join();

	ctx.transferFinish();
	ctx.notifyWrite();

	writeT.join();

	status.update(0.99f, "Writing polygons...");

	ctx.commitOutput();

	status.update(1.0f, "Done.");
}
