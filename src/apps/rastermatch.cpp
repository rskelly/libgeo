/*
 * rastermatch.cpp
 *
 *  Created on: Apr 17, 2018
 *      Author: rob
 */

#include <iostream>
#include <vector>
#include <algorithm>

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
			std::cerr << "row " << row << "\n";
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

void buildSurface(std::vector<Pt>& pts, KDTree<Pt>& tree,
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

	std::vector<Pt> out;
	std::vector<double> dist;

	/*
	std::cerr << "set up avg tree\n";
	KDTree<Pt> avgTree(3);
	std::vector<Pt> apts;
	int i = 0;
	for(const Pt& pt : pts) {
		if(i++ % 1000 == 0)
			std::cerr << " " << ((double) i++) / pts.size() << "\n";
		tree.knn(pt, 16, std::back_inserter(out), std::back_inserter(dist));
		double s = 0;
		for(Pt& p : out)
			s += p.z;
		apts.emplace_back(pt.x, pt.y, s / out.size());
		out.clear();
		dist.clear();
	}
	avgTree.add(apts.begin(), apts.end());
	avgTree.build();
	std::cerr << "done avg tree\n";
	*/

	out.clear();
	dist.clear();

	double zerodist = 1500;

	// Iterate over the anchors...
	for(int i = 0; i < apropss.size(); ++i) {
		
		const GridProps& aprops = apropss[i];
		
		// Iterate over the cells in the adj raster.
		for(int row = 0; row < adjprops.rows(); ++row) {
			std::cerr << "row " << row << "\n";
			for(int col = 0; col < adjprops.cols(); ++col) {
				double x = adjprops.toCentroidX(col);
				double y = adjprops.toCentroidY(row);
				// Get the n nearest points to the cell centre.
				//tree.knn(Pt(x, y, 0), 16, std::back_inserter(out), std::back_inserter(dist));
				tree.radSearch(Pt(x, y, 0), 500, 1000, std::back_inserter(out), std::back_inserter(dist));

				double ta = 0;
				double tb = 0;
				/*
				bool idw = true;

				if(zerodist) {
					double maxdist = 0;
					for(double d : dist)
						maxdist = d > maxdist ? d : maxdist;
					if(maxdist > zerodist) {
						ta = 0;
						tb = 1;
						idw = false;
					} else {
						out.push_back(Pt(x + (zerodist - maxdist), y, 0));
						dist.push_back(zerodist - maxdist);
					}
				}

				if(idw) {
					for(size_t i = 0; i < dist.size(); ++i) {
						if(dist[i] == 0) {
							ta = out[i].z;
							tb = 1;
							break;
						} else {
							// IDW calculation.
							double w = 1 / std::pow(dist[i], 2);
							ta += out[i].z * w;
							tb += w;
						}
					}
				}
				*/

				double val = 0;
				for(size_t i = 0; i < dist.size(); ++i) {
					//double v = out[i].z * std::sqrt(std::pow(1.0 + dist[i], 2) + 0.000001); // multiquadratic
					double v = out[i].z * std::pow(dist[i], 2) * std::log(dist[i]); // TPS
					//std::cerr << v << "\n";
					ta += v;
				}
				tb = dist.size();

				// Set the adjustment value on the output 
				adjmem.setFloat(adjprops.toCol(x), adjprops.toRow(y), ta / tb);
				out.clear();
				dist.clear();
			}
		}
	}

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
	}

	KDTree<Pt> tree(3);
	tree.add(pts.begin(), pts.end());
	tree.build();

	buildSurface(pts, tree, anchors, abands, target, tband, adjustment);

	return 0;
}
