/*
 * rastermatch.cpp
 *
 *  Created on: Apr 17, 2018
 *      Author: rob
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <thread>
#include <iomanip>

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
#include <geos/linearref/LengthIndexedLine.h>


#include "raster.hpp"
#include "interp/loess.hpp"

using namespace geos::geom;
using namespace geos::linearref;

using namespace geo::raster;
using namespace geo::ds;
using namespace geo::interp;

void usage() {
	std::cerr << "Usage: rastermatch <options>\n"
			<< " -a <anchor files [anchor files [...]]>\n"
			<< " -t <target file>\n"
			<< " -o [adjusted file]\n"
			<< " -f [adjustment file]\n"
			<< " -p [points file (csv)]\n\n"
			<< "  This program samples the differences between one or more 'anchor'\n"
			<< "  rasters and a 'target' raster, then calculates an adjustment to \n"
			<< "  match them by computing a local regression surface and applying\n"
			<< "  it to the target\n";
}

class Pt {
public:
	double x, y, z;
	Pt() : Pt(0, 0, 0) {}
	Pt(double x, double y, double z) :
		x(x), y(y), z(z) {}
	double operator[](size_t idx) const {
		switch(idx % 3) {
		case 1: return y;
		case 2: return z;
		default: return x;
		}
	}
	double& operator[](size_t idx) {
		switch(idx % 3) {
		case 1: return y;
		case 2: return z;
		default: return x;
		}
	}
};

void getDiffs(std::vector<Pt>& pts, 
		const std::vector<std::string>& anchors, const std::vector<int>& abands,
		const std::string& target, int tband, double difflimit) {

	std::cerr << "get diffs\n";

	Raster traster(target, false);
	GridProps tprops(traster.props());
	tprops.setBands(1);
	tprops.setWritable(true);
	MemRaster trast(tprops);
	traster.writeTo(trast, tprops.cols(), tprops.rows(), 0, 0, 0, 0, tband, 1);

	for(size_t i = 0; i < anchors.size(); ++i) {
		
		std::cerr << "anchor: " << anchors[i] << "\n";
		Raster araster(anchors[i], false);
		int aband = abands[i];
		
		GridProps aprops(araster.props());
		aprops.setBands(1);
		aprops.setWritable(true);
		MemRaster arast(aprops);
		araster.writeTo(arast, aprops.cols(), aprops.rows(), 0, 0, 0, 0, aband, 1);

		for(int row = 0; row < aprops.rows(); ++row) {
			//std::cerr << "row " << row << "\n";
			for(int col = 0; col < aprops.cols(); ++col) {
				double x = aprops.toCentroidX(col);
				double y = aprops.toCentroidY(row);
				if(tprops.hasCell(x, y) && aprops.hasCell(x, y)) {
					int tc = tprops.toCol(x);
					int tr = tprops.toRow(y);
					double a = arast.getFloat(col, row);
					if(a == aprops.nodata()) // || ar.getInt(col, row, 1) < 100)
						continue;
					double b = trast.getFloat(tc, tr);
					if(b == tprops.nodata()) // || trast.getInt(tc, tr, 1) < 100)
						continue;
					double diff = a - b;
					//if(std::abs(diff) > difflimit)
					//	continue;
					pts.emplace_back(x, y, diff);
				}
			}
		}
	}
}

void buffer(std::vector<Pt>& pts, double buf, double seg) {
	//GEOSContextHandle_t gctx = OGRGeometry::createGEOSContext();
	GeometryFactory::unique_ptr geomFactory = GeometryFactory::create(new PrecisionModel());
	std::vector<Coordinate> coords;
	for(const Pt& p : pts)
		coords.emplace_back(p.x, p.y, p.z);
	MultiPoint* mp = geomFactory->createMultiPoint(coords);
	Geometry* hull = mp->convexHull();
	Geometry* buffered = hull->buffer(buf);
	for(size_t i = 0; i < buffered->getNumGeometries(); ++i) {
		const Polygon* poly = dynamic_cast<const Polygon*>(buffered->getGeometryN(i));
		const LineString* ring = dynamic_cast<const LinearRing*>(poly->getExteriorRing());
		LengthIndexedLine lil(ring);
		double len = ring->getLength();
		for(double p = 0; p < len; p += seg) {
			Coordinate c = lil.extractPoint(p);
			pts.emplace_back(c.x, c.y, 0);
		}
	}
	delete mp;
	delete hull;
	delete buffered;
}

void doInterp(int tileSize, std::list<std::pair<int, int> >* tiles, const std::vector<Pt>* pts, Loess<Pt>* loess,
	const GridProps* adjprops, MemRaster* adjmem, std::mutex* mtx) {

	// As long as there are tiles, keep working.
	while(true) {
		int scol, srow;
		{
			// Get a tile from the list.
			std::lock_guard<std::mutex> lk(*mtx);
			if(tiles->empty())
				return;
			srow = std::get<0>(tiles->front());
			scol = std::get<1>(tiles->front());
			tiles->pop_front();
		}

		std::cerr << "scol " << scol << ", " << srow << "\n";

		// A list to collect interpolated adjustment values.
		std::list<std::tuple<double, double, double> > out;

		// Iterate over the cells in the adj raster and compute the interp value
		// for that cell.
		for(int row = srow; row < std::min(srow + tileSize, adjprops->rows()); ++row) {
			std::cerr << "row: " << row << "\n";
			for(int col = scol; col < std::min(scol + tileSize, adjprops->cols()); ++col) {
				double x = adjprops->toCentroidX(col);
				double y = adjprops->toCentroidY(row);
				double z = loess->estimate(Pt(x, y, 0));
				out.push_back(std::make_tuple(x, y, z));
			}
		}
		{
			// Write all adjustments to the raster in one go.
			std::lock_guard<std::mutex> lk(*mtx);
			for(auto& p : out)
				adjmem->setFloat(adjprops->toCol(std::get<0>(p)), adjprops->toRow(std::get<1>(p)), std::get<2>(p));
			out.clear();
		}
	}

}

void interpolate(const std::vector<Pt>& pts,
		const std::vector<std::string>& anchors, const std::vector<int>& abands,
		const std::string& target, int tband, const std::string& adjustment) {

	// Get the props and bounds for the target.
	Raster traster(target, false);
	GridProps tprops(traster.props());

	// Create props and bounds for the adjustment raster.
	GridProps adjprops(tprops);
	Bounds adjbounds = adjprops.bounds();

	// Build the list of props for anchors and extent
	// the adjustment bounds.
	std::vector<GridProps> apropss;
	for(const std::string& a : anchors) {
		Raster arast(a, false);
		const GridProps& props = arast.props();
		adjbounds.extend(props.bounds());
		apropss.push_back(props);
	}

	// Create the adjustment raster.
	adjprops.setBounds(adjbounds);
	adjprops.setWritable(true);
	adjprops.setBands(1);
	MemRaster adjmem(adjprops, false);

	double range = 500;
	Loess<Pt> loess(pts, range);

	// Each thread will work on a tile, which is a square region of the adjustment region.
	int tileSize = 256;
	std::list<std::pair<int, int> > tiles;
	for(int r = 0; r < adjprops.rows(); r += tileSize) {
		for(int c = 1024; c < 2048 /*adjprops.cols()*/; c += tileSize)
			tiles.push_back(std::make_pair(r, c));
	}

	// Start the threads.
	g_debug("starting")
	int threadCount = 8;
	std::mutex mtx;
	std::vector<std::thread> threads;
	for(int i = 0; i < threadCount; ++i)
		threads.emplace_back(doInterp, tileSize, &tiles, &pts, &loess, &adjprops, &adjmem, &mtx);

	// Wait for completion.
	for(std::thread& t : threads)
		t.join();

	// Write the adjustment raster to the output.
	Raster adjraster(adjustment, adjprops);
	adjmem.writeTo(adjraster);
}

int main(int argc, char** argv) {

	if(argc < 6) {
		usage();
		return 1;
	}

	std::vector<std::string> anchors;
	std::vector<int> abands;
	std::string target;
	std::string ptsFile = "pts.csv";
	int tband;
	std::string adjusted;
	std::string adjustment;
	int mode = 0;
	double difflimit = 1;

	for(int i = 1; i < argc; ++i) {
		std::string arg(argv[i]);
		if(arg == "-a") {
			mode = 1;
			continue;
		} else if(arg == "-t") {
			mode = 2;
			continue;
		} else if(arg == "-o") {
			mode = 3;
			continue;
		} else if(arg == "-f") {
			mode = 4;
			continue;
		} else if(arg == "-p") {
			ptsFile = argv[++i];
			continue;
		}
		switch(mode) {
		case 1:
			std::cerr << "anchor " << arg << ", " << argv[i + 1] << "\n";
			anchors.push_back(arg);
			abands.push_back(atoi(argv[++i]));
			break;
		case 2:
			std::cerr << "target " << arg << ", " << argv[i + 1] << "\n";
			target = arg;
			tband = atoi(argv[++i]);
			break;
		case 3:
			std::cerr << "adjusted " << arg << "\n";
			adjusted = arg;
			break;
		case 4:
			std::cerr << "adjustment " << arg << "\n";
			adjustment = arg;
			break;
		}
	}

	std::vector<Pt> pts;
	{
		getDiffs(pts, anchors, abands, target, tband, difflimit);
		buffer(pts, 500, 100);
		std::cerr << pts.size() << " points\n";
		std::ofstream of(ptsFile);
		of << std::setprecision(12);
		of << "x,y,diff\n";
		for(const Pt& pt : pts)
			of << pt.x << "," << pt.y << "," << pt.z << "\n";
	}

	interpolate(pts, anchors, abands, target, tband, adjustment);

	return 0;
}
