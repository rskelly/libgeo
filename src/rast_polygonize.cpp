/*
 * rast_polygonize.cpp
 *
 *  Created on: Feb 23, 2018
 *      Author: rob
 */

#include <thread>
#include <condition_variable>

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

	// The final geometry.
	geos::geom::Geometry* m_geom;
	// The list of consitiuent geometries.
	std::vector<geos::geom::Geometry*> m_geoms;
	// Unique geometry ID.
	uint64_t m_id;
	// The minimum and maximum row
	// index from the source raster.
	long m_minRow, m_maxRow;
	// True if this poly has already been finalized.
	bool m_finalized;
	// True if the geometries have been collapsed after
	// the most recent update.
	bool m_collapsed;

	/**
	 * Add the list of polygons to this instance with the rows
	 * covered by them.
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
	 * @param minRow The lowest row index occupied by this geometry.
	 * @param maxRow The highest row index occupied by this geometry.
	 */
	void update(int minRow, int maxRow) {
		if(minRow < m_minRow) m_minRow = minRow;
		if(maxRow > m_maxRow) m_maxRow = maxRow;
	}


public:

	/**
	 * Create a polygon with the given ID and initial geometry,
	 * for the given rows. The factory is shared by all
	 * instances.
	 *
	 * @param id A unique ID.
	 */
	Poly(uint64_t id) :
		m_geom(nullptr),
		m_id(id),
		m_minRow(std::numeric_limits<long>::max()),
		m_maxRow(std::numeric_limits<long>::min()),
		m_finalized(false),
		m_collapsed(false) {
	}

	/**
	 * Add the contents of the Poly to this instance.
	 * @param p The Poly to add.
	 */
	void update(Poly& p) {
		p.collapse();
		update(p.geom(), p.minRow(), p.maxRow());
	}

	/**
	 * Add the polygon to this instance with the
	 * rows covered by it.
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
	 * Returns true if the range of rows given by start and end was finalized
	 * within the given thread, or by a thread whose block is completed.
	 * The checked range includes one row above and one below,
	 * which is required to guarantee that a polygon is completed.
	 * @param rowFinished The last row read.
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
	 */
	Geometry* geom() const {
		return m_geom;
	}

	const std::vector<Geometry*>& geoms() const {
		return m_geoms;
	}

	void _collapse(std::vector<Polygon*>* group, std::vector<Geometry*>* geoms) {
		if(!group->empty()) {
			Geometry* p = geos::operation::geounion::CascadedPolygonUnion::CascadedPolygonUnion::Union(group);
			switch(p->getGeometryTypeId()) {
			case GEOS_POLYGON:
				geoms->push_back(p);
				break;
			case GEOS_MULTIPOLYGON:
				const Geometry *g;
				for(size_t i = 0; i < p->getNumGeometries(); ++i) {
					g = p->getGeometryN(i);
					if(!g || g->getGeometryTypeId() != GEOS_POLYGON)
						g_runerr("Null or invalid polygon.");
					geoms->push_back(g->clone());
				}
				break;
			default:
				break;
			}
		}
	}

	void collapse() {
		if(m_collapsed || m_geoms.empty())
			return;
		size_t geomCount = m_geoms.size();
		int change = 1;
		Geometry* g0;
		while(change) {
			std::vector<Polygon*> polys;
			for(Geometry* g : m_geoms) {
				switch(g->getGeometryTypeId()){
				case GEOS_MULTIPOLYGON:
					for(size_t i = 0; i < g->getNumGeometries(); ++i) {
						g0 = g->getGeometryN(i)->clone();
						if(!g0 || g0->getGeometryTypeId() != GEOS_POLYGON)
							g_runerr("Null or invalid polygon.");
						polys.push_back(dynamic_cast<Polygon*>(g0));
					}
					break;
				case GEOS_POLYGON:
					g0 = g->clone();
					if(!g0 || g0->getGeometryTypeId() != GEOS_POLYGON)
						g_runerr("Null or invalid polygon.");
					polys.push_back(dynamic_cast<Polygon*>(g0));
					break;
				default:
					g_runerr("Illegal geometry type: " << g->getGeometryTypeId());
					break;
				}
				delete g;
			}
			m_geoms.clear();

			std::cerr << "collapsing " << polys.size() << " polys\n";
			// The max group size is 1024; min 2.
			int groupSize = g_max(2, g_min(1024, polys.size() / 16));
			// If the group size is 2, just use one group.
			if(groupSize == 2)
				groupSize = polys.size();

			// Number of threads
			// TODO: Configurable. Causes breaks in polys; need a final merge.
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

			std::cerr << "collapse " << change << ", " << geomCount << "\n";
		}

		m_collapsed = true;
	}

	/**
	 * Generate the unioned polygon from its parts. Remove dangles and holes if required.
	 * @param fact The GeometryFactory.
	 * @param removeHoles True, if it is desired that holes in polygons be removed.
	 * @param removeDangles True, if it is desired that single-pixel artifacts, attached at the corners, be removed.
	 */
	void finalize(const GeometryFactory::unique_ptr& fact, bool removeHoles, bool removeDangles) {

		if(m_finalized)
			return;
		m_finalized = true;

		collapse();

		std::cerr << "finalize\n";

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
		std::cerr << "geom " << m_geom << "\n";
	}

	/**
	 * Return the unique ID.
	 */
	uint64_t id() const {
		return m_id;
	}

	/**
	 * Return the minimum row index covered by the geometry.
	 */
	int minRow() const {
		return m_minRow;
	}

	/**
	 * Return the maximum row index covered by the geometry.
	 */
	int maxRow() const {
		return m_maxRow;
	}

	~Poly() {
		if(m_geom)
			delete m_geom;
	}
};

typedef std::unordered_map<uint64_t, std::unique_ptr<Poly> > PolyMap;
typedef std::queue<std::unique_ptr<Poly> > PolyQueue;

class PolyContext {
private:
	Raster* m_raster;
	Raster* m_maskRaster;
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
	uint64_t m_featureId;
	int m_bufSize;
	int m_band;
	bool m_removeHoles;
	bool m_removeDangles;
	bool m_queueUpdate;
	bool m_readFinish;
	bool m_transferFinish;
	std::vector<bool> m_rowsFinished;

public:
	PolyContext(Raster* raster, Status* status, int* block, bool* cancel,
			int bufSize, int band, bool removeHoles, bool removeDangles, Raster* maskRaster = nullptr) :
		m_raster(raster),
		m_maskRaster(maskRaster),
		m_status(status),
		m_block(block),
		m_cancel(cancel),
		m_ds(nullptr),
		m_layer(nullptr),
		m_geomFactory(nullptr),
		m_gctx(nullptr),
		m_featureId(0),
		m_bufSize(bufSize),
		m_band(band),
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

	void initOutput(const std::string& driver, const std::string& filename, const std::string& layerName, int srid) {
		GDALAllRegister();

		// Create the GEOS context and factory.
		m_gctx = OGRGeometry::createGEOSContext();
		m_geomFactory = GeometryFactory::create(new PrecisionModel(1.0));

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
		m_sr.importFromEPSG(srid);
		char **lopts = NULL;
		if(Util::lower(driver) == "sqlite")
			lopts = CSLSetNameValue(lopts, "FORMAT", "SPATIALITE");
		m_layer = m_ds->CreateLayer(layerName.c_str(), &m_sr, wkbMultiPolygon, lopts);
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

	uint64_t nextFeatureId() {
		std::lock_guard<std::mutex> lk(m_fidMtx);
		return ++m_featureId;
	}

	void notifyTransfer() {
		m_polyMapCond.notify_all();
	}

	void notifyWrite() {
		m_polyQueueCond.notify_all();
	}

	int getCols(const GridProps& props, Raster* mask) {
		if(mask) {
			const GridProps& mprops = mask->props();
			const Bounds& mbounds = mprops.bounds();
			return std::abs(props.toCol(mbounds.maxx()) - props.toCol(mbounds.minx())) + 1;
		} else {
			return props.cols();
		}
	}

	int getRows(const GridProps& props, Raster* mask) {
		if(mask) {
			const GridProps& mprops = mask->props();
			const Bounds& mbounds = mprops.bounds();
			return std::abs(props.toRow(mbounds.maxy()) - props.toRow(mbounds.miny())) + 1;
		} else {
			return props.rows();
		}
	}

	int getCol(const GridProps& props, Raster* mask) {
		if(mask) {
			const GridProps& mprops = mask->props();
			const Bounds& mbounds = mprops.bounds();
			return props.toCol(mprops.resolutionX() > 0 ? mbounds.minx() : mbounds.maxx());
		} else {
			return 0;
		}
	}

	int getRow(const GridProps& props, Raster* mask) {
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
		uint64_t nd = (uint64_t) props.nodata();

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

			if(m_status)
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

				// Read over the columns in the current row.
				for(int c = 0; c < cols; ++c) {

					// Get the current ID, skip if nodata.
					uint64_t id0 = rowBuf[c];
					if(id0 == nd || id0 == 0)
						continue;

					// Get the coord of one corner of the polygon.
					double x0 = gp.toX(c + col);
					double y0 = gp.toY(r + row);

					// Scan right...
					while(++c <= cols) {

						// Get the next ID or zero if beyond the edge.
						uint64_t id1 = c < cols ? rowBuf[c] : 0;

						// If the ID changes, capture and output the polygon.
						if(id0 > 0 && id1 != id0) {

							// Coord of the other corner.
							double x1 = gp.toX(c + col);
							double y1 = gp.toY(r + row) + gp.resolutionY();

							// Build the geometry.
							CoordinateSequence* seq = m_geomFactory->getCoordinateSequenceFactory()->create(5, 2);
							seq->setAt(Coordinate(x0, y0), 0);
							seq->setAt(Coordinate(x0, y1), 1);
							seq->setAt(Coordinate(x1, y1), 2);
							seq->setAt(Coordinate(x1, y0), 3);
							seq->setAt(Coordinate(x0, y0), 4);
							LinearRing* ring = m_geomFactory->createLinearRing(seq);
							Polygon* geom = m_geomFactory->createPolygon(ring, NULL);
							// Update the polygon list with the new poly.
							{
								std::lock_guard<std::mutex> lk(m_polyMapMtx);
								if(m_polyMap.find(id0) == m_polyMap.end())
									m_polyMap[id0].reset(new Poly(id0));
								m_polyMap[id0]->update(geom, r, r);
							}

							--c; // Back up the counter by one, to start with the new ID.
							break;
						}
						id0 = id1;
					}
				}

				// Update the finished rows.
				m_rowsFinished[r] = true;

				// Notify the transfer threads of an update.
				notifyTransfer();
			}

		} // while

		std::cerr << "read block finished\n";

		notifyTransfer();
	}

	void polyQueueTransfer() {

		while(!*m_cancel && !(m_readFinish && m_polyMap.empty())) {

			std::unique_ptr<Poly> geom;
			{
				// Wait for a wake-up.
				std::unique_lock<std::mutex> lk(m_polyMapMtx);
				while(!*m_cancel && !m_readFinish && m_polyMap.empty())
					m_polyMapCond.wait(lk);

				// Process the finished items if there are any.
				for(auto it = m_polyMap.begin(); it != m_polyMap.end(); ) {
					if(it->second->isRangeFinalized(m_rowsFinished)) {
						geom.swap(it->second);
						it = m_polyMap.erase(it);
					} else {
						++it;
					}
					break;
				}
			}

			if(geom.get()) {
				geom->finalize(m_geomFactory, m_removeHoles, m_removeDangles);
				{
					std::lock_guard<std::mutex> lk(m_polyQueueMtx);
					m_polyQueue.push(std::move(geom));
				}
			}

			notifyWrite();
		}

		std::cerr << "poly transfer queue finished\n";

		notifyWrite();
	}

	std::unique_ptr<Poly> polyFromPath(int id, const std::vector<double>& path) {
		// Build the geometry.
		CoordinateSequence* seq = m_geomFactory->getCoordinateSequenceFactory()->create((size_t) 0, 2);
		std::vector<Coordinate> coords;
		for(size_t i = 0; i < path.size(); i += 2)
			coords.push_back(Coordinate(path[i], path[i+1]));
		seq->add(&coords, false);

		// Get the ring and make a polygon.
		LinearRing* ring = m_geomFactory->createLinearRing(seq);
		Polygon* geom = m_geomFactory->createPolygon(ring, NULL);
		// Create and return the poly.
		std::unique_ptr<Poly> p(new Poly(id));
		p->update(geom, 0, 0);
		return std::move(p);
	}

	void enqueuePoly(std::unique_ptr<Poly>& poly) {
		poly->finalize(m_geomFactory, m_removeHoles, m_removeDangles);
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

			std::unique_ptr<Poly> p;
			{
				// Wait for a wake-up.
				std::unique_lock<std::mutex> lk(m_polyQueueMtx);
				while(!*m_cancel && !m_transferFinish && m_polyQueue.empty())
					m_polyQueueCond.wait(lk);

				// Grab the element from the queue, if there is one.
				if(!m_polyQueue.empty()) {
					p.swap(m_polyQueue.front());
					m_polyQueue.pop();
				}
			}

			// If no element, loop.
			if(!p.get())
				continue;

			std::cerr << "write\n";

			// Retrieve the unioned geometry and write it.
			const Geometry* g = p->geom();

			{
				std::vector<Polygon*> polys;
				Geometry* p;
				switch(g->getGeometryTypeId()){
				case GEOS_MULTIPOLYGON:
					for(size_t i = 0; i < g->getNumGeometries(); ++i) {
						p = g->getGeometryN(i)->clone();
						if(!p || p->getGeometryTypeId() != GEOS_POLYGON)
							g_runerr("Null or invalid polygon.");
						polys.push_back(dynamic_cast<Polygon*>(p));
					}
					break;
				case GEOS_POLYGON:
					p = g->clone();
					if(!p || p->getGeometryTypeId() != GEOS_POLYGON)
						g_runerr("Null or invalid polygon.");
					polys.push_back(dynamic_cast<Polygon*>(p));
					break;
				default:
					g_runerr("Illegal geometry type: " << g->getGeometryTypeId());
					break;
				}
				g = geos::operation::geounion::CascadedPolygonUnion::CascadedPolygonUnion::Union(&polys);
			}

			OGRGeometry* geom = OGRGeometryFactory::createFromGEOS(m_gctx, (GEOSGeom) g);
			if(!geom)
				g_runerr("Null geometry.");
			OGRFeature feat(m_layer->GetLayerDefn());
			feat.SetGeometry(geom);
			feat.SetField("id", (GIntBig) p->id());
			feat.SetFID(nextFeatureId());
			delete geom;        // The geom is copied by the feature.

			int err;
			{
				std::lock_guard<std::mutex> lk(m_ogrMtx);
				err = m_layer->CreateFeature(&feat);
			}
			if(OGRERR_NONE != err)
				g_runerr("Failed to add geometry.");

			std::cerr << "poly\n";
		}

		std::cerr << "poly write queue finished\n";
	}

};

bool s_cancel = false;

void Raster::polygonize(const std::string& filename, const std::string& layerName,
		const std::string& driver, uint16_t srid, uint16_t band, uint16_t threads,
		int bufSize, bool removeHoles, bool removeDangles, Status *status, bool *cancel, const std::string& mask) {

	if(!m_props.isInt())
		g_runerr("Only integer rasters can be polygonized.");

	// There need to be at least three threads for readin from the raster(1),
	// writing output (1), and unioning and transferring to the write queue (n - 2).
	if(threads < 3)
		threads = 3;

	// If no cancel pointer is given, use a dummy.
	if(cancel == nullptr)
		cancel = &s_cancel;

	// Counter for the current block.
	int block = 0;

	// If the bufsize given is invalid, use the entire raster height.
	// This can use a lot of memory!
	if(bufSize <= 0 || bufSize > props().rows()) {
		g_warn("Invalid buffer size; using raster height.");
		bufSize = props().rows();
	}

	// Remove the original file; some can't be overwritten directly.
	// This will not take care of any auxillary files (e.g. shapefiles)
	Util::rm(filename);

	Raster* maskRaster = nullptr;
	if(!mask.empty())
		maskRaster = new Raster(mask);

	// Set up the shared context.
	PolyContext ctx(this, status, &block, cancel, bufSize, band, removeHoles, removeDangles, maskRaster);

	// Initialize database.
	ctx.initOutput(driver, filename, layerName, srid);

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

	g_debug("Writing...");
	if(status)
		status->update(0.99f, "Writing polygons...");

	ctx.commitOutput();

	g_debug("Done");
	if(status)
		status->update(1.0f, "Done.");
}
