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

#include "raster.hpp"
#include "ds/kdtree.hpp"
#include "rbf.hpp"

using namespace geo::raster;
using namespace geo::ds;

void usage() {
	std::cerr << "Usage: rastermatch -a <anchor files [anchor files [...]]> -t <target file> [-o <adjusted file>] [-f <adjustment file>]\n"
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
		case 0: return x;
		case 1: return y;
		case 2: return z;
		}
		return 0;
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

	const Bounds& tbounds = tprops.bounds();

	for(int i = 0; i < anchors.size(); ++i) {
		
		std::cerr << "anchor: " << anchors[i] << "\n";
		Raster araster(anchors[i], false);
		int aband = abands[i];
		
		GridProps aprops(araster.props());
		aprops.setBands(1);
		aprops.setWritable(true);
		MemRaster arast(aprops);
		araster.writeTo(arast, aprops.cols(), aprops.rows(), 0, 0, 0, 0, aband, 1);

		const Bounds& abounds = aprops.bounds();
		
		double x, y;
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

void surfaceWork(int span, std::list<std::pair<int, int> >* spans, const std::vector<Pt>* pts,
	const GridProps* adjprops, MemRaster* adjmem, std::mutex* mtx) {

	KDTree<Pt> tree(3);
	tree.add(pts->begin(), pts->end());
	tree.build();

	RBF<Pt> rbf(RBF<Pt>::Type::ThinPlateSpline, 1);

	while(true) {
		int scol, srow;
		{
			std::lock_guard<std::mutex> lk(*mtx);
			if(spans->empty())
				return;
			srow = std::get<0>(spans->front());
			scol = std::get<1>(spans->front());
			spans->pop_front();
		}

		std::cerr << "scol " << scol << ", " << srow << "\n";

		scol *= span;
		srow *= span;

		std::vector<Pt> pts;
		std::vector<double> dist;
		
		double x = adjprops->toCentroidX(scol + span / 2);
		double y = adjprops->toCentroidY(srow + span / 2);

		int count = tree.radSearch(Pt(x, y, 0), span * std::abs(adjprops->resolutionX()) * 2, 9999, 
			std::back_inserter(pts), std::back_inserter(dist));
		std::cerr << "found " << count << " within " << (span * std::abs(adjprops->resolutionX()) * 2) << " of " << x << ", " << y << "\n";
		if(32 > count)
			continue;

		rbf.clear();
		rbf.add(pts.begin(), pts.end());
		rbf.build();

		std::list<std::tuple<double, double, double> > out;

		// Iterate over the cells in the adj raster.
		for(int row = srow; row < std::min(srow + span, adjprops->rows()); ++row) {
			std::cerr << "row: " << row << "\n";
			for(int col = scol; col < std::min(scol + span, adjprops->cols()); ++col) {
				double x = adjprops->toCentroidX(col);
				double y = adjprops->toCentroidY(row);
				double z = rbf.compute(Pt(x, y, 0));
				//std::cerr << x << ", " << y << ", " << z << "\n";
				out.push_back(std::make_tuple(x, y, z));
			}
			{
				std::lock_guard<std::mutex> lk(*mtx);
				for(auto& p : out)
					adjmem->setFloat(adjprops->toCol(std::get<0>(p)), adjprops->toRow(std::get<1>(p)), std::get<2>(p));
				out.clear();
			}
		}
	}
}

void buildSurface(std::vector<Pt>& pts, 
		const std::vector<std::string>& anchors, const std::vector<int>& abands,
		const std::string& target, int tband, const std::string& adjustment) {

	// Get the props and bounds for the target.
	Raster traster(target, false);
	GridProps tprops(traster.props());
	const Bounds& tbounds = tprops.bounds();

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

	std::mutex mtx;
	std::vector<std::thread> threads;

	int threadCount = 8;
	int span = (int) std::ceil((double) adjprops.rows() / threadCount);
	std::list<std::pair<int, int> > spans;
	for(int r = 0; r < threadCount; ++r) {
		for(int c = 0; c < threadCount; ++c)
			spans.push_back(std::make_pair(r, c));
	}
	for(int i = 0; i < threadCount; ++i)
		threads.emplace_back(surfaceWork, span, &spans, &pts, &adjprops, &adjmem, &mtx);
	for(std::thread& t : threads)
		t.join();

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
		std::vector<Pt> tmp;
		getDiffs(tmp, anchors, abands, target, tband, difflimit);
		std::random_shuffle(tmp.begin(), tmp.end());
		size_t count = std::min(tmp.size(), (size_t) 10000);
		pts.assign(tmp.begin(), tmp.begin() + count);
		std::cerr << pts.size() << " points from " << tmp.size() << "\n";
		std::ofstream of("pts.csv");
		of << "x,y,diff\n";
		for(const Pt& pt : pts)
			of << pt.x << "," << pt.y << "," << pt.z << "\n";
	}

	buildSurface(pts, anchors, abands, target, tband, adjustment);

	return 0;
}
