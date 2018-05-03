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
#include "ds/kdtree.hpp"
//#include "rbf.hpp"
#include "kmeans.hpp"

#include <pcl/surface/mls.h>

using namespace geo::raster;
using namespace geo::ds;
//using namespace geo::interp;

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
/*
double smoothIDW(const Pt& pt, const std::vector<Pt>& pts) {
	double d, w, s0 = 0, s1 = 0;
	for(const Pt& p0 : pts) {
		d = dist(pt, p0);
		if(d == 0)
			return p0.z;
		w = 1.0 / std::pow(d, 2.0);
		s0 += p0.z * w;
		s1 += w;
	}
	return s0 / s1;
}

double smoothGauss(const Pt& pt, const std::vector<Pt>& pts, double sigma = 0.8, double range = 100) {
	double d, s0 = 0, s1 = 0;
	for(const Pt& p0 : pts) {
		d = dist(pt, p0) / range;
		s0 += (1.0 / (sigma * std::sqrt(2 * PI))) * std::exp(-0.5 * std::pow(d / sigma, 2)) * d;
		s1 += d;
	}
	return s0;
}
*/
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

/*
std::mutex __tree_mtx;

void surfaceWork(int tileSize, std::list<std::pair<int, int> >* tiles, const std::vector<Pt>* pts, RBF<Pt>* _rbf, KDTree<Pt>* tree,
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
				Pt p(x, y, 0);
				std::list<Pt> found;
				std::list<double> dist;
				int count = 0;
				{
					std::lock_guard<std::mutex> lk(__tree_mtx);
					count = tree->radSearch(p, 1500, 999999, std::back_inserter(found), std::back_inserter(dist));
				}
				if(count < 3) {
					out.push_back(std::make_tuple(x, y, 0));
				} else {
					// Prepare the RBF instance.
					RBF<Pt> rbf(RBF<Pt>::Type::Gaussian);
					rbf.setRange(500);
					rbf.setSigma(0.8);
					rbf.add(found.begin(), found.end());
					//rbf.setSmoothing(1);
					rbf.build();
					double z = rbf.compute(p);
					out.push_back(std::make_tuple(x, y, z));
				}
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

void buildSurface(std::vector<Pt>& pts, 
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
	*/
	/*
	// Prepare the RBF instance.
	RBF<Pt> rbf(RBF<Pt>::Type::Gaussian);
	rbf.setRange(500);
	rbf.setSigma(0.8);
	//rbf.setSmoothing(1);
	//rbf.setClusters(1000);
	//rbf.setSamples(1000);

	rbf.add(pts.begin(), pts.end());
	rbf.build();
	*/
	/*
	KDTree<Pt> tree(3);
	tree.add(pts.begin(), pts.end());
	tree.build();

	// Each thread will work on a tile, which is a square region of the adjustment region.
	int tileSize = 256;
	std::list<std::pair<int, int> > tiles;
	for(int r = 512; r < 1024 adjprops.rows(); r += tileSize) {
		for(int c = 1024; c < 2048 adjprops.cols(); c += tileSize)
			tiles.push_back(std::make_pair(r, c));
	}

	// Start the threads.
	g_debug("starting")
	int threadCount = 4;
	std::mutex mtx;
	std::vector<std::thread> threads;
	for(int i = 0; i < threadCount; ++i)
		threads.emplace_back(surfaceWork, tileSize, &tiles, &pts, nullptr &rbf, &tree, &adjprops, &adjmem, &mtx);

	// Wait for completion.
	for(std::thread& t : threads)
		t.join();

	// Write the adjustment raster to the output.
	Raster adjraster(adjustment, adjprops);
	adjmem.writeTo(adjraster);
}
*/

void buildSurface2(std::vector<Pt>& pts,
		const std::vector<std::string>& anchors, const std::vector<int>& abands,
		const std::string& target, int tband, const std::string& adjustment)  {

	pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointXYZ> mls;
	mls.setPolynomialOrder(2);
	mls.setSearchRadius(1000);
	mls.setUpsamplingMethod(mls.UpsamplingMethod::VOXEL_GRID_DILATION);
	mls.setDilationVoxelSize(5);
	mls.setDilationIterations(5);

	std::list<pcl::PointXYZ> pts0;
	for(const Pt& p : pts) {
		if(p.x == 0) {
			g_debug("zero")
			continue;
		}
		pts0.emplace_back(p.x, p.y, p.z);
	}

	pcl::PointCloud<pcl::PointXYZ>::Ptr output(new pcl::PointCloud<pcl::PointXYZ>());
	pcl::PointCloud<pcl::PointXYZ>::Ptr input(new pcl::PointCloud<pcl::PointXYZ>());
	input->insert(input->end(), pts0.begin(), pts0.end());

	mls.setInputCloud(input);
	mls.process(*output);

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

	for(const pcl::PointXYZ& p : *output) {
		if(adjprops.hasCell(p.x, p.y))
			adjmem.setFloat(adjprops.toCol(p.x), adjprops.toRow(p.y), p.z);
	}

	// Write the adjustment raster to the output.
	Raster adjraster(adjustment, adjprops);
	adjmem.writeTo(adjraster);

}

using namespace geos::geom;
using namespace geos::linearref;

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

	std::vector<Pt> means;
	{
		std::vector<Pt> pts;
		getDiffs(pts, anchors, abands, target, tband, difflimit);
		{
			std::unordered_map<size_t, std::list<Pt> > clusters;
			kmeans(pts, 1000, means, clusters, 1000);
		}
		buffer(means, 500, 100);
		std::cerr << pts.size() << " points\n";
		{
		std::ofstream of(ptsFile);
		of << std::setprecision(12);
		of << "x,y,diff\n";
		for(const Pt& pt : means)
			of << pt.x << "," << pt.y << "," << pt.z << "\n";
		}
		{
			std::ofstream of("all.csv");
			of << std::setprecision(12);
			of << "x,y,diff\n";
			for(const Pt& pt : pts)
				of << pt.x << "," << pt.y << "," << pt.z << "\n";
		}
	}

	buildSurface2(means, anchors, abands, target, tband, adjustment);

	return 0;
}
