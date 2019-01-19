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

#include "raster.hpp"
#include "ds/kdtree.hpp"

using namespace geo::raster;
using namespace geo::ds;


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
		default: return z;
		}
	}
	bool operator==(const Pt& pt) const {
		return pt.x == x && pt.y == y && pt.z == z;
	}
};

void usage() {
	std::cerr << "Usage: rastermatch <options>\n"
			<< " -a <anchor files [anchor files [...]]>\n"
			<< " -t <target file>\n"
			<< " -o [adjusted file]\n"
			<< " -f [adjustment (difference) file]\n"
			<< " -r [kernel radius in map units. the function is gaussian. Sigma 1 is 1/4 of the radius.]\n"
			<< "  This program samples the differences between one or more 'anchor'\n"
			<< "  rasters and a 'target' raster, then calculates an adjustment to \n"
			<< "  match them\n";
}

bool getDiff(MemRaster& target, std::vector<MemRaster*>& anchors,
		std::vector<bool>& mkernel, std::vector<double>& wkernel,
		double x, double y, double& z) {

	const GridProps& tprops = target.props();
	double tn = tprops.nodata();

	double ds = 0;
	int dc = 0;
	int side = (int) std::sqrt(mkernel.size());

	for(int r = 0; r < side; ++r) {
		for(int c = 0; c < side; ++c) {
			if(mkernel[r * side + c]) {
				for(MemRaster* a : anchors) {

					const GridProps& aprops = a->props();
					double an = aprops.nodata();

					int ac = aprops.toCol(x) - side / 2 + c;
					int ar = aprops.toRow(y) - side / 2 + r;
					int tc = tprops.toCol(x) - side / 2 + c;
					int tr = tprops.toRow(y) - side / 2 + r;

					double av, tv;

					if(aprops.hasCell(ac, ar) && tprops.hasCell(tc, tr)
							&& (av = a->getFloat(ac, ar, 1)) != an
							&& (tv = target.getFloat(tc, tr, 1)) != tn) {

						ds += (av - tv) * wkernel[r * side + c];
						++dc;
					}
				}
			}
		}
	}

	if(dc) {
		z = ds;
		return true;
	}

	return false;
}
/**

 */
void doInterp(MemRaster& target, KDTree<Pt>& tree, MemRaster& adjusted, MemRaster& diffs, double radius, int count) {

	const GridProps& tprops = target.props();
	int tcols = tprops.cols();
	int trows = tprops.rows();
	double tres = std::abs(target.props().resolutionX());
	double tn = tprops.nodata();
	double tv;

	std::vector<Pt*> pts;
	std::vector<double> dist;

	for(int trow = 0; trow < trows; ++trow) {
		std::cerr << "Row: " << trow << " of " << trows << "\n";
		for(int tcol = 0; tcol < tcols; ++tcol) {

			if((tv = target.getFloat(tcol, trow, 1)) == tn)
				continue;

			double x = tprops.toCentroidX(tcol);
			double y = tprops.toCentroidY(trow);

			int res = 0;
			double rad = tres;
			double cnt = count;
			{
				Pt q(x, y, 0);

				do {
					pts.clear();
					dist.clear();
					res = tree.radSearch(q, rad, std::back_inserter(pts), std::back_inserter(dist));
					rad += tres;
					cnt = (std::log2(rad / tres) + 1) * count;
				} while(res < cnt && rad <= radius);
			}

			{
				double d, w0, w = 0, adj = 0;
				for(int i = 0; i < res; ++i) {
					if((d = dist[i]) >= 0){
						const Pt& pt = *pts[i];
						w0 = 1.0 - std::max(0.0, std::min(1.0, d / radius));
						adj += pt.z * w0;
						w += 1;
					}
				}
				if(res) {
					adj /= w;
					adjusted.setFloat(tcol, trow, tv + adj, 1);
					diffs.setFloat(tcol, trow, adj, 1);
				} else {
					adjusted.setFloat(tcol, trow, tv, 1);
					diffs.setFloat(tcol, trow, 0, 1);
				}
			}
		}
	}
}

int main(int argc, char** argv) {

	if(argc < 6) {
		usage();
		return 1;
	}

	std::vector<std::string> anchors;	// Anchor filenames.
	std::vector<int> abands;			// Anchor bands.
	std::string target;					// Target filename.
	int tband = 1;						// Target band.
	std::string adjusted;				// The adjusted raster.
	std::string adjustment;				// The adjustment (differencce).
	int count = 1;
	double radius = 100;

	int mode = 0;

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
		} else if(arg == "-c") {
			count = atoi(argv[++i]);
			mode = 0;
			continue;
		} else if(arg == "-r") {
			radius = atof(argv[++i]);
			mode = 0;
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

	// Get the target image as a mem raster.
	MemRaster tgrid;
	GridProps tprops;
	// Make an adjusted raster.
	MemRaster agrid;
	// Make an adjustment (diff) raster.
	MemRaster dgrid;
	{
		Raster traster(target);
		tprops = GridProps(traster.props());
		tprops.setWritable(true);
		tprops.setBands(1);
		tgrid.init(tprops, true);
		traster.writeTo(tgrid, tprops.cols(), tprops.rows(), 0, 0, 0, 0, tband, 1);

		agrid.init(tprops, true);
		agrid.fillFloat(tprops.nodata(), 1);

		dgrid.init(tprops, true);
		dgrid.fillFloat(tprops.nodata(), 1);
	}


	KDTree<Pt> tree(2);
	// One large raster. Each valid pixel is the average of anchor pixels
	// minus the corresponding target pixel.
	{
		// Get mem rasters of each anchor and an extended bounds object.
		std::vector<MemRaster> agrids(anchors.size());
		for(size_t i = 0; i < anchors.size(); ++i) {
			Raster anchor(anchors[i]);
			GridProps props(anchor.props());
			props.setBands(1);
			props.setWritable(true);
			agrids[i].init(props, true);
			anchor.writeTo(agrids[i], props.cols(), props.rows(), 0, 0, 0, 0, abands[i], 1);
		}

		//std::ofstream tmp("tmp.csv");
		//tmp << std::setprecision(9);
		double tn = tprops.nodata();
		for(int tr = 0; tr < tprops.rows(); ++tr) {
			for(int tc = 0; tc < tprops.cols(); ++tc) {
				double tv;
				if((tv = tgrid.getFloat(tc, tr, 1)) == tn)
					continue;
				double x = tprops.toCentroidX(tc);
				double y = tprops.toCentroidY(tr);
				double sum = 0;
				int cnt = 0;
				for(size_t i = 0; i < agrids.size(); ++i) {
					const GridProps& aprops = agrids[i].props();
					int ac = aprops.toCol(x);
					int ar = aprops.toRow(y);
					double an = aprops.nodata();
					double v;
					if(aprops.hasCell(ac, ar) && (v = agrids[i].getFloat(ac, ar, 1)) != an) {
						sum += v;
						++cnt;
					}
				}
				if(cnt) {
					tree.add(new Pt(x, y, sum / cnt - tv));
					//tmp << x << "," << y << "," << (sum / cnt - tv) << "\n";
					//std::cerr << x << ", " << y << ", " << sum << ", " << cnt << ", " << tv << ", " << (sum / cnt - tv) << "\n";
				}
			}
		}

		tree.build();
	}

	doInterp(tgrid, tree, agrid, dgrid, radius, count);

	Raster araster(adjusted, agrid.props());
	Raster draster(adjustment, dgrid.props());

	agrid.writeTo(araster, agrid.props().cols(), agrid.props().rows(), 0, 0, 0, 0, 1, 1);
	dgrid.writeTo(draster, dgrid.props().cols(), dgrid.props().rows(), 0, 0, 0, 0, 1, 1);

 	return 0;
}
