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

/**
 * A simple point class for storing raster diffs.
 */
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
	std::cout << "Usage: rastermatch <options>\n"
			<< " -a <<anchor file> <band> [<anchor file> <band> [...]]>\n"
			<< " -t <<target file> <band>>\n"
			<< " -k <<mask file> <band>>"
			<< " -o [adjusted file]\n"
			<< " -f [adjustment (difference) file]\n"
			<< " -n Number of nearest neighbours. Default 32.\n"
			<< "  This program samples the differences between one or more 'anchor'\n"
			<< "  rasters and a 'target' raster, then calculates an adjustment to \n"
			<< "  match them\n";
}

double makeKernel(std::vector<double>& kernel, int size, double sigma) {
	kernel.resize(size * size);
	double xmid = (size - 1) / 2.0;
	double ymid = (size - 1) / 2.0;
	double spread = 1.0 / (sigma * sigma * 2.0);
	double denom = 8.0 * std::atan(1) * sigma * sigma;
	std::vector<double> gx, gy;
	gx.reserve(size);
	for(int i = 0; i < size; ++i) {
		double x = i - xmid;
		gx.push_back(std::exp(-x * x * spread));
	}
	gy.reserve(size);
	for(int i = 0; i < size; ++i) {
		double y = i - ymid;
		gy.push_back(std::exp(-y * y * spread));
	}

	double sum = 0;
	for(int r = 0; r < size; ++r) {
		for(int c = 0; c < size; ++c) {
			sum += (kernel[r * size + c] = (gx[c] * gy[r]) / denom);
		}
	}
	std::cerr << "Sum: " << sum << "\n";
	return sum;
}

double edgeDistance(int col, int row, int cols, int rows, double res) {
	double dx = std::min(col, cols - col);
	double dy = std::min(row, rows - row);
	double d = std::sqrt(dx * dx + dy * dy) * res;
	return std::round(d * res) / res;
}

void doGaussInterp(MemRaster& errors, MemRaster& target, MemRaster& adjusted, MemRaster& diffs, double sigma) {

	const GridProps& tprops = target.props();
	int tcols = tprops.cols();
	int trows = tprops.rows();
	double tn = tprops.nodata();

	const GridProps& eprops = errors.props();
	int ecols = eprops.cols();
	int erows = eprops.rows();
	double en = eprops.nodata();

	std::vector<double> kernel;

	double tv, ev;
	double res = std::abs(tprops.resolutionX());

	int size = (int) std::ceil((sigma / res) * 3); // 3 standard deviations
	if(size % 2 == 0) ++size;
	makeKernel(kernel, size, sigma);

	diffs.fillFloat(0, 1);

	for(int erow = 0; erow < erows; ++erow) {
		if(erow % 100 == 0)
			std::cout << "A Row: " << erow << " of " << erows << "\n";
		for(int ecol = 0; ecol < ecols; ++ecol) {

			double x = eprops.toCentroidX(ecol);
			double y = eprops.toCentroidY(erow);
			int tcol = tprops.toCol(x);
			int trow = tprops.toRow(y);

			if(!tprops.hasCell(tcol, trow) || (ev = errors.getFloat(ecol, erow, 1)) == en)
				continue;
			for(int r = 0; r < size; ++r) {
				for(int c = 0; c < size; ++c) {
					int cc = tcol + c - size / 2;
					int rr = trow + r - size / 2;
					if(tprops.hasCell(cc, rr))
						diffs.setFloat(cc, rr, diffs.getFloat(cc, rr, 1) + kernel[r * size + c] * ev, 1);
				}
			}
		}
	}

	for(int trow = 0; trow < trows; ++trow) {
		if(trow % 100 == 0)
			std::cout << "B Row: " << trow << " of " << trows << "\n";
		for(int tcol = 0; tcol < tcols; ++tcol) {
			if((tv = target.getFloat(tcol, trow, 1)) != tn)
				adjusted.setFloat(tcol, trow, tv + diffs.getFloat(tcol, trow, 1), 1);
		}
	}
}


void doInterp(KDTree<Pt>& tree, MemRaster& target, MemRaster& adjusted, MemRaster& diffs, int nn) {

	const GridProps& tprops = target.props();
	int tcols = tprops.cols();
	int trows = tprops.rows();
	double tn = tprops.nodata();

	double tv;
	std::vector<Pt*> pts;
	std::vector<double> dist;

	diffs.fillFloat(0, 1);

	for(int trow = 0; trow < trows; ++trow) {
		if(trow % 100 == 0)
			std::cout << "Row: " << trow << " of " << trows << "\n";
		for(int tcol = 0; tcol < tcols; ++tcol) {

			if((tv = target.getFloat(tcol, trow, 1)) == tn)
				continue;

			double x = tprops.toCentroidX(tcol);
			double y = tprops.toCentroidY(trow);

			pts.clear();
			dist.clear();
			int count = tree.knn(Pt(x, y, 0), nn, std::back_inserter(pts), std::back_inserter(dist));

			if(count) {
				double sum = 0;
				double w = 0;
				double maxDist = std::sqrt(dist[dist.size() - 1]);
				for(int i = 0; i < count; ++i) {
					double d = 1.0 - std::sqrt(dist[i]) / maxDist;
					sum += pts[i]->z * d;
					w += d;
				}
				sum /= w;
				diffs.setFloat(tcol, trow, sum, 1);
				adjusted.setFloat(tcol, trow, tv + sum, 1);
			} else {
				diffs.setFloat(tcol, trow, 0, 1);
				adjusted.setFloat(tcol, trow, tv, 1);
			}
		}
	}

}

int match(int argc, char** argv) {

	if(argc < 6) {
		usage();
		return 1;
	}

	std::vector<std::string> anchors;	// Anchor filenames.
	std::vector<int> abands;			// Anchor bands.
	std::string target;					// Target filename.
	int tband = 1;						// Target band.
	bool hasMask = false;
	std::string maskfile;
	int maskband = 1;
	std::string adjusted;				// The adjusted raster.
	std::string adjustment;				// The adjustment (differencce).
	bool mapped = false;
	int nn = 32;
	int skip = 4;
	int sigma = 100;
	std::string points; // If outputing dif points.
	bool doPoints = false;
	bool doNN = false;
	bool doGauss = false;

	for(int i = 2; i < argc; ++i) {
		std::string arg(argv[i]);
		if(arg == "-a") {
			anchors.push_back(argv[++i]);
			abands.push_back(atoi(argv[++i]));
			continue;
		} else if(arg == "-t") {
			target = argv[++i];
			tband = atoi(argv[++i]);
			continue;
		} else if(arg == "-o") {
			adjusted = argv[++i];
			continue;
		} else if(arg == "-f") {
			adjustment = argv[++i];
			continue;
		} else if(arg == "-m") {
			mapped = true;
			continue;
		} else if(arg == "-k") {
			hasMask = true;
			maskfile = argv[++i];
			maskband = atoi(argv[++i]);
			continue;
		} else if(arg == "-n") {
			nn = atoi(argv[++i]);
			doNN = true;
			continue;
		} else if(arg == "-s") {
			skip = atoi(argv[++i]);
			continue;
		} else if(arg == "-g") {
			sigma = atof(argv[++i]);
			doGauss = true;
			continue;
		} else if(arg == "-p") {
			points = argv[++i];
			doPoints = true;
			continue;
		}
	}

	std::cout << "Skip: " << skip << "; neighbours: " << nn << "; mapped: " << mapped << "; sigma: " << sigma << "\n";

	// Get the target image as a mem raster.
	MemRaster tgrid;
	GridProps tprops;
	// Make an adjusted raster.
	MemRaster agrid;
	// Make an adjustment (diff) raster.
	MemRaster dgrid;

	{
		std::cout << "Loading target\n";
		Raster traster(target);
		tprops = GridProps(traster.props());
		tprops.setWritable(true);
		tprops.setBands(1);
		tgrid.init(tprops, mapped);
		traster.writeTo(tgrid, tprops.cols(), tprops.rows(), 0, 0, 0, 0, tband, 1);
	}

	double tn = tprops.nodata();
	KDTree<Pt> tree(2);
	// A reaster for difference
	MemRaster egrid;

	{
		// Bounds for the error raster.
		Bounds ebounds;
		// Make a mask raster.
		MemRaster mgrid;
		// Properties for the mask raster.
		GridProps mprops;
		if(hasMask) {
			std::cout << "Loading mask\n";
			Raster mask(maskfile);
			mprops = mask.props();
			mprops.setWritable(true);
			mprops.setBands(1);
			mgrid.init(mprops, true);
			mask.writeTo(mgrid, mprops.cols(), mprops.rows(), 0, 0, 0, 0, maskband, 1);
		}

		// Get mem rasters of each anchor and an extended bounds object.
		std::vector<MemRaster> agrids(anchors.size());
		for(size_t i = 0; i < anchors.size(); ++i) {
			std::cout << "Loading anchor " << anchors[i] << "\n";
			Raster anchor(anchors[i]);
			GridProps props(anchor.props());
			props.setBands(1);
			props.setWritable(true);
			agrids[i].init(props, true);
			anchor.writeTo(agrids[i], props.cols(), props.rows(), 0, 0, 0, 0, abands[i], 1);
			ebounds.extend(props.bounds());
		}

		// Configure the error raster.
		GridProps eprops(tprops);
		eprops.setBounds(ebounds);
		eprops.setNoData(0);
		egrid.init(eprops, mapped);
		egrid.fillFloat(0, 1);

		double tv;

		// Fill the error grid with differences between anchors and (shifted) target.
		std::cout << "Computing errors\n";
		for(int tr = 0; tr < tprops.rows(); tr += (1 + skip)) {
			for(int tc = 0; tc < tprops.cols(); tc += (1 + skip)) {
				if((tv = tgrid.getFloat(tc, tr, 1)) == tn)
					continue;
				double x = tprops.toCentroidX(tc);
				double y = tprops.toCentroidY(tr);
				for(size_t i = 0; i < agrids.size(); ++i) {
					const GridProps& aprops = agrids[i].props();
					int ac = aprops.toCol(x);
					int ar = aprops.toRow(y);
					int mc = mprops.toCol(x);
					int mr = mprops.toRow(y);
					int ec = eprops.toCol(x);
					int er = eprops.toRow(y);
					double an = aprops.nodata();
					double v;
					if(aprops.hasCell(ac, ar)
							&& (!hasMask || (mprops.hasCell(mc, mr) && mgrid.getInt(mc, mr, 1) == 1))
							&& ((v = agrids[i].getFloat(ac, ar, 1)) != an)
							&& (egrid.getFloat(ec, er, 1) == 0)) {
						egrid.setFloat(ec, er, v - tv, 1);
					}
				}
			}
		}
		std::ofstream ptsOutput;
		if(doPoints) {
			ptsOutput.open(points);
			std::cout << "Points file: " << points << "\n";
		}

		if(doNN || doPoints) {
			double ev;
			for(int er = 0; er < eprops.rows(); ++er) {
				for(int ec = 0; ec < eprops.cols(); ++ec) {
					double x = eprops.toCentroidX(ec);
					double y = eprops.toCentroidY(er);
					if((ev = egrid.getFloat(ec, er, 1)) != 0) {
						if(doPoints) {
							ptsOutput << x << "," << y << "," << ev << "\n";
						} else if(doNN){
							tree.add(new Pt(x, y, ev));
						}
					}
				}
			}
		}
	}

	if(!doPoints) {
		agrid.init(tprops, mapped);
		agrid.fillFloat(tprops.nodata(), 1);

		dgrid.init(tprops, mapped);
		dgrid.fillFloat(tprops.nodata(), 1);

		std::cout << "Matching with " << nn << " nearest neighbours.\n";

		if(doNN) {
			tree.build();
			doInterp(tree, tgrid, agrid, dgrid, nn);
		} else if(doGauss) {
			doGaussInterp(egrid, tgrid, agrid, dgrid, sigma);
		}

		GridProps oprops(agrid.props());
		oprops.setDataType(DataType::Float32);

		std::cout << "Creating output rasters\n";
		Raster araster(adjusted, oprops);
		Raster draster(adjustment, oprops);

		std::cout << "Writing to output rasters\n";
		agrid.writeTo(araster, oprops.cols(), oprops.rows(), 0, 0, 0, 0, 1, 1);
		dgrid.writeTo(draster, oprops.cols(), oprops.rows(), 0, 0, 0, 0, 1, 1);
	}

 	return 0;
}

double pixelWeights(Grid& grid, int col, int row, int size) {
	const GridProps& props = grid.props();
	double nd = props.nodata();
	double rad = std::pow(size / 2.0, 2.0);
	int ct = 0, t = 0;
	for(int r = 0; r < size; ++r) {
		for(int c = 0; c < size; ++c) {
			double d = std::pow(c - size / 2.0, 2.0) + std::pow(r - size / 2, 2.0);
			if(d <= rad) {
				++t;
				if(props.hasCell(c + col, r + row) && grid.getFloat(c + col, r + row, 1) != nd)
					++ct;
			}
		}
	}
	return t > 0 ? (double) ct / t : 0;
}

int merge(int argc, char** argv) {

	if(argc < 6) {
		usage();
		return 1;
	}

	std::vector<char*> args;
	std::string target1;					// Target filename.
	std::string target2;					// Target filename.
	int tband1 = 1;							// Target band.
	int tband2 = 1;							// Target band.
	std::string merged;						// Merged filename
	double radius = 50;						// Kernel radius (map units)
	bool mapped = false;					// Use mapped memory.

	for(int i = 2; i < argc; ++i) {
		std::string arg(argv[i]);
		if(arg == "-r") {
			radius = atof(argv[++i]);
			continue;
		} else if(arg == "-m") {
			mapped = true;
			continue;
		} else {
			args.push_back(argv[i]);
		}
	}

	// Get files and bands.
	target1 = args[0];
	tband1 = atoi(args[1]);
	target2 = args[2];
	tband2 = atoi(args[3]);
	merged = args[4];

	// Get the target image as a mem raster.
	MemRaster tgrid1, tgrid2;
	GridProps tprops1, tprops2;
	{
		Raster traster(target1);
		tprops1 = GridProps(traster.props());
		tprops1.setWritable(true);
		tprops1.setBands(1);
		tgrid1.init(tprops1, mapped);
		traster.writeTo(tgrid1, tprops1.cols(), tprops1.rows(), 0, 0, 0, 0, tband1, 1);
	}
	{
		Raster traster(target2);
		tprops2 = GridProps(traster.props());
		tprops2.setWritable(true);
		tprops2.setBands(1);
		tgrid2.init(tprops2, mapped);
		traster.writeTo(tgrid2, tprops2.cols(), tprops2.rows(), 0, 0, 0, 0, tband2, 1);
	}

	// Get the bounds of the merger and create props.
	GridProps mprops(tprops1);
	{
		Bounds mbounds = mprops.bounds();
		mbounds.extend(tprops2.bounds());
		mprops.setBounds(mbounds);
	}
	// Create and prepare the output raster.
	MemRaster mgrid(mprops, mapped);
	mgrid.fillFloat(mprops.nodata(), 1);

	Raster mraster(merged, mprops);

	double nd1 = tprops1.nodata();
	double nd2 = tprops2.nodata();
	int mcols = mprops.cols();
	int mrows = mprops.rows();
	double v;			// Adjusted, weighted value.
	double v1, v2;		// Cell values.
	double w1, w2; 		// Valid-pixel weights.
	double x, y;
	int tcol1, trow1;
	int tcol2, trow2;	// Col, row in grid2.
	int size = (int) std::ceil(radius / std::abs(tprops1.resolutionX()));

	// Iterate over the cells of the merged raster.
	for(int mrow = 0; mrow < mrows; ++mrow) {
		if(mrow % 100 == 0)
			std::cout << "Row " << mrow << " of " << mrows << "\n";
		for(int mcol = 0; mcol < mcols; ++mcol) {

			x = mprops.toCentroidX(mcol);
			y = mprops.toCentroidY(mrow);

			v1 = nd1;
			v2 = nd2;

			// Get the col/row and value of the first target.
			tcol1 = tprops1.toCol(x);
			trow1 = tprops1.toRow(y);
			if(tprops1.hasCell(x, y))
				v1 = tgrid1.getFloat(tcol1, trow1, 1);

			// Get the col/row and value of the second target.
			tcol2 = tprops2.toCol(x);
			trow2 = tprops2.toRow(y);
			if(tprops2.hasCell(x, y))
				v2 = tgrid2.getFloat(tcol2, trow2, 1);

			if(v1 != nd1 && v2 != nd2) {
				// If both are valid, get the weights and compute a value.
				w1 = pixelWeights(tgrid1, tcol1, trow1, size);
				w2 = pixelWeights(tgrid2, tcol2, trow2, size);
				if(w1 == 0 && w2 == 0) {
					v = nd1;
				} else if(w1 == 0) {
					v = v2;
				} else if(w1 == 0) {
					v = v1;
				} else {
					v = (w1 / (w1 + w2)) * v1 + (w2 / (w1 + w2)) * v2;
				}
			} else if(v1 != nd1) {
				// If v1 is valid, use it.
				v = v1;
			} else if(v2 != nd2) {
				// If v2 is valid use it.
				v = v2;
			} else {
				// If none are valid, use ND.
				v = nd1;
			}

			// Set the output value.
			mgrid.setFloat(mcol, mrow, v, 1);
		}
	}

	mgrid.writeTo(mraster, mprops.cols(), mprops.rows(), 0, 0, 0, 0, 1, 1);

 	return 0;
}

int main(int argc, char** argv) {

	if(argc < 2) {
		usage();
		return 1;
	}

	std::string cmd = argv[1];

	if(cmd == "merge") {
		return merge(argc, argv);
	} else if(cmd == "match") {
		return match(argc, argv);
	} else {
		usage();
		return 1;
	}
}
