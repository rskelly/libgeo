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

using namespace geo::raster;

void usage() {
	std::cerr << "Usage: rastermatch <options>\n"
			<< " -a <anchor files [anchor files [...]]>\n"
			<< " -t <target file>\n"
			<< " -o [adjusted file]\n"
			<< " -f [adjustment (difference) file]\n"
			<< " -r [kernel radius in map units]\n"
			<< " -d [decay (>1 - 2 is best)]\n"
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
void doInterp(MemRaster& target, std::vector<MemRaster*>& anchors, MemRaster& adjusted, MemRaster& diffs,
		std::vector<bool>& mkernel, std::vector<double>& wkernel) {

	const GridProps& tprops = target.props();
	int tcols = tprops.cols();
	int trows = tprops.rows();
	double tn = tprops.nodata();

	std::list<std::pair<double, double> > pts;

	for(int trow = 0; trow < trows; ++trow) {
		std::cerr << "Row: " << trow << " of " << trows << "\n";
		for(int tcol = 0; tcol < tcols; ++tcol) {

			double v = target.getFloat(tcol, trow, 1);
			if(v != tn) {
				double x = tprops.toX(tcol);
				double y = tprops.toY(trow);
				double z;
				if(getDiff(target, anchors, mkernel, wkernel, x, y, z)) {
					diffs.setFloat(tcol, trow, z, 1);
					adjusted.setFloat(tcol, trow, v + z, 1);
				} else {
					diffs.setFloat(tcol, trow, 0, 1);
					adjusted.setFloat(tcol, trow, v, 1);
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
	double radius = 100;				// In map units; converted to cells using the abs x resolution.

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

	std::vector<MemRaster*> agrids;
	for(size_t i = 0; i < anchors.size(); ++i) {
		Raster anchor(anchors[i]);
		GridProps props(anchor.props());
		props.setBands(1);
		props.setWritable(true);
		agrids.push_back(new MemRaster(props, true));
		anchor.writeTo(*agrids[i], props.cols(), props.rows(), 0, 0, 0, 0, abands[i], 1);
	}

	Raster traster(target);
	GridProps props(traster.props());
	props.setWritable(true);
	props.setBands(1);

	MemRaster tgrid(props, true);
	MemRaster agrid(props, true);
	MemRaster dgrid(props, true);

	traster.writeTo(tgrid, props.cols(), props.rows(), 0, 0, 0, 0, tband, 1);
	agrid.fillFloat(props.nodata(), 1);
	dgrid.fillFloat(props.nodata(), 1);

	int pxRadius = (int) std::ceil(radius / std::abs(props.resolutionX()));
	int rMax = pxRadius * pxRadius;
	int side = pxRadius * 2 + 1;
	std::vector<bool> mkernel(side * side);
	std::vector<double> wkernel(side * side);

	double wt = 0;
	for(int r = 0; r < side; ++r) {
		for(int c = 0; c < side; ++c) {
			double d = std::pow(pxRadius - c, 2.0) + std::pow(pxRadius - r, 2.0);
			double w = std::min(1.0, std::max(0.0, 1.0 - d / rMax));
			mkernel[r * side + c] = d <= rMax;
			wkernel[r * side + c] = w;
			wt += w;
		}
	}
	for(int r = 0; r < side; ++r) {
		for(int c = 0; c < side; ++c)
			wkernel[r * side + c] /= wt;
	}

	doInterp(tgrid, agrids, agrid, dgrid, mkernel, wkernel);

	Raster araster(adjusted, props);
	Raster draster(adjustment, props);

	agrid.writeTo(araster, props.cols(), props.rows(), 0, 0, 0, 0, 1, 1);
	dgrid.writeTo(draster, props.cols(), props.rows(), 0, 0, 0, 0, 1, 1);

 	return 0;
}
