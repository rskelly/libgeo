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

/**
 * Produces a 2D Gaussian kernel given the mean, std. deviation and number of std. devs to extend.
 */
class Gaussian {
private:

public:
	double sigma;
	double extent;

	Gaussian(double sigma, double extent = 3) :
		sigma(sigma), extent(extent) {}

	/**
	 * Get the value.
	 */
	double operator()(double x, double y) {
		return (1.0 / (2 * M_PI * sigma * sigma)) * std::exp(-((x * x + y * y) / (2 * sigma * sigma)));
	}

	/**
	 * Populates a square kernel by modifying the given vector.
	 * The resolution determines how many cells are used (i.e. cell-per-resolution).
	 * Returns the length of one side of the kernel.
	 */
	int kernel(std::vector<double>& data, double resolution) {

		if(resolution <= 0)
			std::runtime_error("Invalid resolution.");

		int size = (sigma * extent * 2.0) / resolution;
		if(size % 2 == 0) ++size;
		data.resize(size * size);

		double xspread = 1.0 / (sigma * sigma * 2.0);
		double yspread = 1.0 / (sigma * sigma * 2.0);
		double xmid = size / 2;
		double ymid = size / 2;
		double denom = 8.0 * std::atan(1) * sigma * sigma;

		std::vector<double> gx, gy;
		gx.reserve(size);
		gy.reserve(size);

		for(int i = 0; i < size; ++i)
			gx[i] = std::exp(-std::pow(i - xmid, 2) * xspread);

		for(int i = 0; i < size; ++i)
			gy[i] = std::exp(-std::pow(i - ymid, 2) * yspread);

		double sum = 0;
		for(int r = 0; r < size; ++r) {
			for(int c = 0; c < size; ++c)
				sum += (data[r * size + c] = (gx[c] * gy[r]) / denom);
		}
		/*
		if(false && std::abs(sum - 1.0) > 0.01) {
			double sum0 = 0;
			for(size_t i = 0; i < data.size(); ++i)
				sum0 += (data[i] /= sum);
			sum = sum0;
		}
		*/
		if(true) {
			for(int i = 0; i < size; ++i)
				std::cout << data[(size / 2) * size + i] << ",";
			std::cout << "\n";
		}
		std::cerr << "Kernel sum " << sum << "\n";
		return size;
	}
};

void usage() {
	std::cerr << "Usage: rastermatch <options>\n"
			<< " -a <<anchor file> <band> [<anchor file> <band> [...]]>\n"
			<< " -t <<target file> <band>>\n"
			<< " -k <<mask file> <band>>"
			<< " -o [adjusted file]\n"
			<< " -f [adjustment (difference) file]\n"
			<< " -s [Standard deviation in map units. The function is gaussian.]\n"
			<< " -e [The number of standard deviations to extend. Defaults to 3, or SD*3 for the kernel radius.]"
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


void doInterp(MemRaster& errors, MemRaster& target, MemRaster& adjusted, MemRaster& diffs, double sigma, double extent) {

	const GridProps& eprops = errors.props();
	double en = eprops.nodata();

	const GridProps& tprops = target.props();
	int tcols = tprops.cols();
	int trows = tprops.rows();
	double tn = tprops.nodata();

	double tv, ev;
	int size;

	std::vector<double> kernel;
	{
		Gaussian ga(sigma, extent);
		size = ga.kernel(kernel, std::abs(tprops.resolutionX()));
	}

	diffs.fillFloat(0, 1);

	for(int trow = 0; trow < trows; ++trow) {
		std::cerr << "Row: " << trow << "\n";
		for(int tcol = 0; tcol < tcols; ++tcol) {

			double x = tprops.toCentroidX(tcol);
			double y = tprops.toCentroidY(trow);
			int ecol = eprops.toCol(x);
			int erow = eprops.toRow(y);

			if((tv = target.getFloat(tcol, trow, 1)) == tn || !eprops.hasCell(ecol, erow) || (ev = errors.getFloat(ecol, erow, 1)) == en || ev == 0)
				continue;

			double w;
			for(int kr = 0; kr < size; ++kr) {
				for(int kc = 0; kc < size; ++kc) {
					if((w = kernel[kr * size + kc]) != 0) {
						int cc = tcol - size / 2 + kc;
						int rr = trow - size / 2 + kr;
						if(w > 0 && tprops.hasCell(cc, rr) && target.getFloat(cc, rr, 1) != tn)
							diffs.setFloat(cc, rr, w * ev + diffs.getFloat(cc, rr, 1), 1);
					}
				}
			}
		}
	}

	for(int trow = 0; trow < trows; ++trow) {
		for(int tcol = 0; tcol < tcols; ++tcol)
			adjusted.setFloat(tcol, trow, target.getFloat(tcol, trow, 1) + diffs.getFloat(tcol, trow, 1), 1);
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
	double sigma = 100;
	double extent = 3;
	bool mapped = false;

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
		} else if(arg == "-s") {
			sigma = atof(argv[++i]);
			continue;
		} else if(arg == "-e") {
			extent = atof(argv[++i]);
			continue;
		}
	}

	// Get the target image as a mem raster.
	MemRaster tgrid;
	GridProps tprops;
	// Make an adjusted raster.
	MemRaster agrid;
	// Make an adjustment (diff) raster.
	MemRaster dgrid;
	// A reaster for difference
	MemRaster egrid;

	{
		Raster traster(target);
		tprops = GridProps(traster.props());
		tprops.setWritable(true);
		tprops.setBands(1);
		tgrid.init(tprops, mapped);
		traster.writeTo(tgrid, tprops.cols(), tprops.rows(), 0, 0, 0, 0, tband, 1);
	}


	//std::vector<KDTree<Pt> >trees;
	//for(int i = 0; i < tcount; ++i)
	//	trees.emplace_back(2);

	//std::vector<Pt> pts;
	// One large raster. Each valid pixel is the average of anchor pixels
	// minus the corresponding target pixel.
	double avgDif;
	{
		Bounds ebounds;
		// Make a mask raster.
		MemRaster mgrid;
		GridProps mprops;
		if(hasMask) {
			std::cerr << "Loading mask\n";
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
			std::cerr << "Loading anchor " << anchors[i] << "\n";
			Raster anchor(anchors[i]);
			GridProps props(anchor.props());
			props.setBands(1);
			props.setWritable(true);
			agrids[i].init(props, true);
			anchor.writeTo(agrids[i], props.cols(), props.rows(), 0, 0, 0, 0, abands[i], 1);
			ebounds.extend(props.bounds());
		}

		GridProps eprops(tprops);
		eprops.setBounds(ebounds);
		egrid.init(eprops);
		egrid.fillFloat(0, 1);

		//std::ofstream tmp("tmp.csv");
		//tmp << std::setprecision(9);
		double tn = tprops.nodata();

		double sum = 0;
		int ct = 0;
		double tv;

		for(int tr = 0; tr < tprops.rows(); ++tr) {
			for(int tc = 0; tc < tprops.cols(); ++tc) {
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
					double an = aprops.nodata();
					double v;
					if(aprops.hasCell(ac, ar)
							&& (!hasMask || !(mprops.hasCell(mc, mr) && mgrid.getInt(mc, mr, 1) == 1))
							&& (v = agrids[i].getFloat(ac, ar, 1)) != an) {
						sum += v - tv;
						++ct;
					}
				}
			}
		}

		avgDif = ct > 0 ? sum / ct : 0;

		for(int tr = 0; tr < tprops.rows(); ++tr) {
			for(int tc = 0; tc < tprops.cols(); ++tc) {
				if((tv = tgrid.getFloat(tc, tr, 1)) != tn)
					tgrid.setFloat(tc, tr, tgrid.getFloat(tc, tr, 1) + avgDif, 1);
			}
		}

		for(int tr = 0; tr < tprops.rows(); ++tr) {
			for(int tc = 0; tc < tprops.cols(); ++tc) {
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
							&& (!hasMask || !(mprops.hasCell(mc, mr) && mgrid.getInt(mc, mr, 1) == 1))
							&& (v = agrids[i].getFloat(ac, ar, 1)) != an
							&& egrid.getFloat(ec, er, 1) == 0) {
						egrid.setFloat(ec, er, v - tv, 1);
					}
				}
			}
		}

	}


	agrid.init(tprops, mapped);
	agrid.fillFloat(tprops.nodata(), 1);

	dgrid.init(tprops, mapped);
	dgrid.fillFloat(tprops.nodata(), 1);

	doInterp(egrid, tgrid, agrid, dgrid, sigma, extent);

	GridProps oprops(agrid.props());
	oprops.setDataType(DataType::Float32);
	Raster araster(adjusted, oprops);
	Raster draster(adjustment, oprops);

	agrid.writeTo(araster, oprops.cols(), oprops.rows(), 0, 0, 0, 0, 1, 1);
	dgrid.writeTo(draster, oprops.cols(), oprops.rows(), 0, 0, 0, 0, 1, 1);

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
	int tband1 = 1;						// Target band.
	int tband2 = 1;						// Target band.
	std::string merged;
	double radius = 50;

	for(int i = 2; i < argc; ++i) {
		std::string arg(argv[i]);
		if(arg == "-r") {
			radius = atof(argv[++i]);
			continue;
		} else {
			args.push_back(argv[i]);
		}
	}

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
		tgrid1.init(tprops1, true);
		traster.writeTo(tgrid1, tprops1.cols(), tprops1.rows(), 0, 0, 0, 0, tband1, 1);
	}
	{
		Raster traster(target2);
		tprops2 = GridProps(traster.props());
		tprops2.setWritable(true);
		tprops2.setBands(1);
		tgrid2.init(tprops2, true);
		traster.writeTo(tgrid2, tprops2.cols(), tprops2.rows(), 0, 0, 0, 0, tband2, 1);
	}

	GridProps mprops(tprops1);
	{
		Bounds mbounds = mprops.bounds();
		mbounds.extend(tprops2.bounds());
		mprops.setBounds(mbounds);
	}
	MemRaster mgrid(mprops, true);
	mgrid.fillFloat(mprops.nodata(), 1);

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

	for(int mrow = 0; mrow < mrows; ++mrow) {
		if(mrow % 100 == 0)
			std::cerr << "Row " << mrow << " of " << mrows << "\n";
		for(int mcol = 0; mcol < mcols; ++mcol) {

			x = mprops.toCentroidX(mcol);
			y = mprops.toCentroidY(mrow);

			v1 = nd1;
			v2 = nd2;

			tcol1 = tprops1.toCol(x);
			trow1 = tprops1.toRow(y);

			if(tprops1.hasCell(x, y))
				v1 = tgrid1.getFloat(tcol1, trow1, 1);

			tcol2 = tprops2.toCol(x);
			trow2 = tprops2.toRow(y);

			if(tprops2.hasCell(x, y))
				v2 = tgrid2.getFloat(tcol2, trow2, 1);

			if(v1 != nd1 && v2 != nd2) {
				mcol = mprops.toCol(x);
				mrow = mprops.toRow(y);
				w1 = pixelWeights(tgrid1, tcol1, trow1, size);
				w2 = pixelWeights(tgrid2, tcol2, trow2, size);
				v = (w1 / (w1 + w2)) * v1 + (w2 / (w1 + w2)) * v2;
				mgrid.setFloat(mcol, mrow, v, 1);
			} else if(v1 != nd1) {
				mgrid.setFloat(mcol, mrow, v1, 1);
			} else if(v2 != nd2) {
				mgrid.setFloat(mcol, mrow, v2, 1);
			} else {
				mgrid.setFloat(mcol, mrow, nd1, 1);
			}
		}
	}

	Raster mraster(merged, mprops);
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
