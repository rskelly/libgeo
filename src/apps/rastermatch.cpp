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

class Gaussian {
private:
	static constexpr double s2p = std::sqrt(2 * M_PI);

public:
	double sigma;
	double mean;
	Gaussian(double sigma, double mean = 0) :
		sigma(sigma), mean(mean) {}
	double operator()(double x) {
		return (1.0 / (sigma * s2p)) * std::pow(M_E, -0.5 * std::pow((x - mean) / sigma, 2.0));
	}
	void kernel(double* data, int size, double scale) {
		double maxd = std::pow(size / 2.0, 2.0);
		double sum = 0;
		for(int r = 0; r < size; ++r) {
			for(int c = 0; c < size; ++c) {
				double d = std::pow(r - size / 2.0, 2.0) + std::pow(c - size / 2.0, 2.0);
				if(d <= maxd) {
					sum += (data[r * size + c] = (*this)(std::sqrt(d)));
				} else {
					data[r * size + c] = 0;
				}
			}
		}
		double t = 0;
		for(int i = 0; i < size * size; ++i) {
			data[i] /= sum;
			t += data[i];
		}
		std::cerr << "Kernel sum " << t << "\n";
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

void doInterp_(MemRaster& target, KDTree<Pt>& tree, MemRaster& adjusted, MemRaster& diffs, int count) {

	const GridProps& tprops = target.props();
	int tcols = tprops.cols();
	int trows = tprops.rows();
	double tn = tprops.nodata();
	double tv;

	std::vector<Pt*> pts;
	std::vector<double> dist;

	for(int trow = 0; trow < trows; ++trow) {
		if(trow % 100 == 0)
			std::cerr << "Row: " << trow << " of " << trows << "\n";
		for(int tcol = 0; tcol < tcols; ++tcol) {
			if((tv = target.getFloat(tcol, trow, 1)) == tn)
				continue;

			double x = tprops.toCentroidX(tcol);
			double y = tprops.toCentroidY(trow);
			Pt pt(x, y, 0);

			pts.clear();
			dist.clear();
			int ct = tree.knn(pt, count, std::back_inserter(pts), std::back_inserter(dist));

			const double dmax = dist[count - 1];
			double sum = 0, w = 0, w0;
			for(int i = 0; i < ct; ++i) {
				const Pt* pt = pts[i];
				w0 = 1.0 - dist[i] / dmax;
				sum += pt->z * w0;
				w += w0;
			}
			sum = ct > 0 ? sum / w : 0;
			adjusted.setFloat(tcol, trow, tv + sum, 1);
			diffs.setFloat(tcol, trow, sum, 1);
		}
	}
}

std::mutex _blkmtx;
std::mutex _tmtx;
std::mutex _dmtx;
int blockSize = 2048;

void doInterp(std::list<std::pair<int, int> >* blocks, MemRaster* target, std::vector<Pt>* pts, MemRaster* adjusted, MemRaster* diffs, double radius) {

	const GridProps& tprops = target->props();
	int tcols = tprops.cols();
	int trows = tprops.rows();
	double tn = tprops.nodata();
	double tv;

	int btrow, btcol;

	GridProps bufProps(tprops);
	bufProps.setSize(blockSize, blockSize);
	MemRaster buf(bufProps, false);
	MemRaster dif(bufProps, false);
	MemRaster adj(bufProps, false);

	while(!blocks->empty()) {
		{
			std::lock_guard<std::mutex> lk(_blkmtx);
			if(blocks->empty())
				break;
			auto pr = blocks->front();
			btrow = pr.second;
			btcol = pr.first;
			blocks->pop_front();
		}
		std::cerr << "Block " << (btcol / blockSize) << ", " << (btrow / blockSize) << "\n";

		adj.fillFloat(bufProps.nodata(), 1);
		dif.fillFloat(bufProps.nodata(), 1);
		buf.fillFloat(bufProps.nodata(), 1);
		{
			std::lock_guard<std::mutex> lk(_tmtx);
			target->writeTo(buf, blockSize, blockSize, btcol, btrow, 0, 0, 1, 1);
		}

		for(int trow = btrow; trow < std::min(btrow + blockSize, trows); ++trow) {
			for(int tcol = btcol; tcol < std::min(btcol + blockSize, tcols); ++tcol) {

				if((tv = buf.getFloat(tcol - btcol, trow - btrow, 1)) == tn)
					continue;

				double x = tprops.toCentroidX(tcol);
				double y = tprops.toCentroidY(trow);

				double sum = 0, w;
				int ct = 0;
				for(const Pt& pt : *pts) {
					w = std::sqrt(std::pow(pt.x - x, 2.0) + std::pow(pt.y - y, 2.0));
					if(w <= radius) {
						sum += pt.z * (1.0 - w / radius);
						++ct;
					}
				}

				sum = ct ? sum / ct : 0;
				adj.setFloat(tcol - btcol, trow - btrow, tv + sum, 1);
				dif.setFloat(tcol - btcol, trow - btrow, sum, 1);
			}
		}

		std::lock_guard<std::mutex> lk(_dmtx);
		adj.writeTo(*adjusted, blockSize, blockSize, 0, 0, btcol, btrow, 1, 1);
		dif.writeTo(*diffs, blockSize, blockSize, 0, 0, btcol, btrow, 1, 1);
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
	int skip = 10;
	double radius = 1000;
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
		} else if(arg == "-s") {
			skip = atoi(argv[++i]);
			continue;
		} else if(arg == "-m") {
			mapped = true;
			continue;
		} else if(arg == "-k") {
			hasMask = true;
			maskfile = argv[++i];
			maskband = atoi(argv[++i]);
			continue;
		} else if(arg == "-r") {
			radius = atof(argv[++i]);
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
	{
		Raster traster(target);
		tprops = GridProps(traster.props());
		tprops.setWritable(true);
		tprops.setBands(1);
		tgrid.init(tprops, mapped);
		traster.writeTo(tgrid, tprops.cols(), tprops.rows(), 0, 0, 0, 0, tband, 1);
	}


	//KDTree<Pt> tree(2);
	std::vector<Pt> pts;
	// One large raster. Each valid pixel is the average of anchor pixels
	// minus the corresponding target pixel.
	{
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
		}

		//std::ofstream tmp("tmp.csv");
		//tmp << std::setprecision(9);
		double tn = tprops.nodata();
		int step = std::max(1, skip);
		for(int tr = 0; tr < tprops.rows(); tr += step) {
			for(int tc = 0; tc < tprops.cols(); tc += step) {
				double tv;
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
						//tree.add(new Pt(x, y, v - tv));
						pts.emplace_back(x, y, v - tv);
					}
				}
			}
		}
	}

	//tree.build();

	agrid.init(tprops, mapped);
	agrid.fillFloat(tprops.nodata(), 1);

	dgrid.init(tprops, mapped);
	dgrid.fillFloat(tprops.nodata(), 1);

	std::list<std::pair<int, int> > blocks;
	for(int r = 0; r < tprops.rows(); r += blockSize) {
		for(int c = 0; c < tprops.cols(); c += blockSize)
			blocks.emplace_back(c, r);
	}

	std::cerr << pts.size() << " points\n";

	std::vector<std::thread> threads;
	for(int t = 0; t < 4; ++t)
		threads.emplace_back(doInterp, &blocks, &tgrid, &pts, &agrid, &dgrid, radius);

	for(int t = 0; t < 4; ++t) {
		if(threads[t].joinable())
			threads[t].join();
	}

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
