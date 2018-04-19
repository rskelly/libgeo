/*
 * rastermatch.cpp
 *
 *  Created on: Apr 17, 2018
 *      Author: rob
 */

#include <iostream>
#include <vector>

#include "raster.hpp"
#include "ds/kdtree.hpp"

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

	Raster trast(target, false);
	const GridProps& tprops = trast.props();
	const Bounds& tbounds = tprops.bounds();

	std::vector<Raster> arast;
	for(const std::string& a : anchors)
		arast.emplace_back(a, false);

	for(int i = 0; i < arast.size(); ++i) {
		Raster& ar  = arast[i];
		int ab = abands[i];
		const GridProps& aprops = ar.props();
		const Bounds& abounds = aprops.bounds();

		double x, y;
		for(int row = 0; row < aprops.rows(); ++row) {
			for(int col = 0; col < aprops.cols(); ++col) {
				double x = aprops.toCentroidX(col);
				double y = aprops.toCentroidY(row);
				if(tprops.hasCell(x, y) && aprops.hasCell(x, y)) {
					int tc = tprops.toCol(x);
					int tr = tprops.toRow(y);
					double a = ar.getFloat(col, row, ab);
					if(a == aprops.nodata()) // || ar.getInt(col, row, 1) < 100)
						continue;
					double b = trast.getFloat(tc, tr, tband);
					if(b == tprops.nodata()) // || trast.getInt(tc, tr, 1) < 100)
						continue;
					double diff = a - b;
					if(std::abs(diff) > difflimit)
						continue;
					pts.emplace_back(x, y, diff);
				}
			}
		}
	}
}

void buildSurface(std::vector<Pt>& pts, KDTree<Pt>& tree,
		const std::vector<std::string>& anchors, const std::vector<int>& abands,
		const std::string& target, int tband, const std::string& adjustment) {

	// Open the target raster and get properties.
	Raster trast(target, false);
	const GridProps& tprops = trast.props();
	const Bounds& tbounds = tprops.bounds();

	// Create props and bounds for the adjustment raster.
	GridProps adjprops(tprops);
	Bounds adjbounds = adjprops.bounds();

	// Create the anchor rasters list.
	std::vector<Raster> arast;
	for(const std::string& a : anchors)
		arast.emplace_back(a, false);

	// Extend the adjustment raster's bounds.
	for(Raster& ar : arast)
		adjbounds.extend(ar.props().bounds());

	// Create the adjustment raster.
	adjprops.setBounds(adjbounds);
	adjprops.setWritable(true);
	adjprops.setBands(1);
	Raster adjraster(adjustment, adjprops);
	MemRaster adjmem(adjprops, false);
	//Raster difraster("diffs.tif", adjprops);

	// Iterate over the anchors...
	std::vector<Pt> out;
	std::vector<double> dist;
	for(int i = 0; i < arast.size(); ++i) {
		Raster& ar  = arast[i];
		int ab = abands[i];
		const GridProps& aprops = ar.props();
		const Bounds& abounds = aprops.bounds();

		// Iterate over the cells in the adj raster.
		for(int row = 0; row < adjprops.rows(); ++row) {
			for(int col = 0; col < adjprops.cols(); ++col) {
				double x = adjprops.toCentroidX(col);
				double y = adjprops.toCentroidY(row);
				// If the coord is in the target and the anchor, process the cell.
				if(tprops.hasCell(x, y) && aprops.hasCell(x, y)) {
					// Get the n nearest points to the cell centre.
					tree.knn(Pt(x, y, 0), 30, std::back_inserter(out), std::back_inserter(dist));
					double ta = 0;
					double tb = 0;
					for(size_t i = 0; i < dist.size(); ++i) {
						if(dist[i] == 0) {
							ta = out[i].z;
							tb = 1;
							break;
						}
						double w = 1 / std::pow(dist[i], 2);
						ta += out[i].z * w;
						tb += w;
					}
					adjmem.setFloat(adjprops.toCol(x), adjprops.toRow(y), ta / tb);
					out.clear();
					dist.clear();
				}
			}
		}
	}

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
			anchors.push_back(arg);
			abands.push_back(atoi(argv[++i]));
			break;
		case 2:
			target = arg;
			tband = atoi(argv[++i]);
			break;
		case 3:
			adjusted = arg;
			break;
		case 4:
			adjustment = arg;
			break;
		}
	}

	std::vector<Pt> pts;
	getDiffs(pts, anchors, abands, target, tband, difflimit);

	KDTree<Pt> tree(3);
	tree.add(pts.begin(), pts.end());
	tree.build();

	buildSurface(pts, tree, anchors, abands, target, tband, adjustment);

	return 0;
}
