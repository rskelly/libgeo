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

#include <Eigen/Core>
#include <Eigen/Dense>

#include "raster.hpp"
#include "ds/kdtree.hpp"

using namespace geo::raster;

void usage() {
	std::cerr << "Usage: rastermatch <options>\n"
			<< " -a <anchor files [anchor files [...]]>\n"
			<< " -t <target file>\n"
			<< " -o [adjusted file]\n"
			<< " -f [adjustment file]\n"
			<< " -r [search radius for interp (default 1000)]\n"
			<< " -s [skip cells in each direction (default 0)]\n"
			<< " -c [max number of points to use when interpolating (default 999999999)]\n"
			<< "  This program samples the differences between one or more 'anchor'\n"
			<< "  rasters and a 'target' raster, then calculates an adjustment to \n"
			<< "  match them by computing a local regression surface and applying\n"
			<< "  it to the target\n";
}

class Pt {
public:
	double _x, _y, _z;
	Pt() : Pt(0, 0, 0) {}
	Pt(double x, double y, double z) :
		_x(x), _y(y), _z(z) {}
	double x() const { return _x; }
	double y() const { return _y; }
	double z() const { return _z; }
	double operator[](size_t idx) const {
		switch(idx % 3) {
		case 1: return _y;
		case 2: return _z;
		default: return _x;
		}
	}
	double& operator[](size_t idx) {
		switch(idx % 3) {
		case 1: return _y;
		case 2: return _z;
		default: return _x;
		}
	}
};


/**
 * Use the average of anchor pixels for each corresponding target pixel to compute a difference and
 * write that to the adjustment raster.
 */
void getDiffs(std::vector<std::pair<std::unique_ptr<Raster>, int> >& anchors, Raster& target, int tband, geo::ds::KDTree<Pt>& tree, int skip) {


	const GridProps& tprops = target.props();
	int tcols = tprops.cols();
	int trows = tprops.rows();
	double tnd = tprops.nodata();

	std::vector<double> avs;	// Anchor values for a given coordinate.

	for(int trow = 0; trow < trows; trow += (1 + skip)) {
		for(int tcol = 0; tcol < tcols; tcol += (1 + skip)) {

			// Get the target pixel value. If it's nodata, write nodata to the adj raster and skip.
			double tv = target.getFloat(tcol, trow, tband);
			if(tv == tnd)
				continue;

			// Collect the average of values from the anchor rasters if they're valid.
			avs.clear();
			for(auto& a : anchors) {
				const GridProps& aprops = a.first->props();
				int acol = aprops.toCol(tprops.toX(tcol));
				int arow = aprops.toRow(tprops.toY(trow));
				if(aprops.hasCell(acol, arow)) {
					double av = a.first->getFloat(acol, arow, (int) a.second);
					if(av != aprops.nodata())
						avs.push_back(av);
				}
			}

			// If there were valid anchor pixels, set the average on the adj raster.
			// Otherwise set nodata.
			if(!avs.empty()) {
				double avg = 0;
				for(double av : avs)
					avg += av;
				avg /= avs.size();
				tree.add(new Pt(tprops.toX(tcol), tprops.toY(trow), tv - avg));
			}
		}
	}
	tree.build();
}

void doInterp(Raster& target, int tband, Raster& output, Raster& adj, geo::ds::KDTree<Pt>& tree, int count) {

	const GridProps& props = target.props();
	int cols = props.cols();
	int rows = props.rows();
	double nd = props.nodata();

	double radius = std::abs(props.resolutionX());
	double step = 360.0 / (2 * radius * M_PI);

	std::vector<Pt*> pts;
	std::vector<double> dist;

	for(int row = 0; row < rows; ++row) {
		for(int col = 0; col < cols; ++col) {

			double v = target.getFloat(col, row, tband);
			if(v == nd) {
				output.setFloat(col, row, nd, 1);
				continue;
			}

			double z;
			int ct;
			double rad = radius;
			double stp = step;
			double maxRad = props.cols() * std::abs(props.resolutionX());
			do {
				z = 0 ;
				ct = 0;
				for(double a = 0.0; a < 360.0; a += stp) {
					int r = (int) std::round(std::cos(a * M_PI / 180.0) * rad);
					int c = (int) std::round(std::sin(a * M_PI / 180.0) * rad);
					if(r < 0 || r >= rows || c < 0 || c > cols)
						continue;
					double v0 = target.getFloat(c, r, tband);
					if(v0 != nd) {
						z += v0;
						++ct;
					}
				}
				rad *= 2;
				stp = 360.0 / (2 * rad * M_PI);
				if(rad > maxRad)
					break;
			} while(ct < count);

			if(ct) {
				z /= nd;
				output.setFloat(col, row, v - z, 1);
				adj.setFloat(col, row, z, 1);
			} else {
				output.setFloat(col, row, nd, 1);
				adj.setFloat(col, row, nd, 1);
			}

		}
	}
}


/**
 * Use the adjustment raster to compute the best-fit plane from which to
 * interpret the pixels of the target and produce an output.
 */
void doInterp3(Raster& target, int tband, Raster& output, Raster& adj, geo::ds::KDTree<Pt>& tree, int count) {

	const GridProps& props = target.props();
	int cols = props.cols();
	int rows = props.rows();
	double nd = props.nodata();

	std::vector<Pt*> pts;
	std::vector<double> dist;

	for(int row = 0; row < rows; ++row) {
		for(int col = 0; col < cols; ++col) {

			double v = target.getFloat(col, row, tband);
			if(v == nd) {
				output.setFloat(col, row, nd, 1);
				continue;
			}

			pts.clear();
			dist.clear();
			Pt c(props.toX(col), props.toY(row), 0);
			int found = tree.knn(c, count, std::back_inserter(pts), std::back_inserter(dist));

			if(found) {
				double z = 0;
				for(int i = 0; i < found; ++i)
					z += pts[i]->_z;
				z /= found;
				output.setFloat(col, row, v - z, 1);
				adj.setFloat(col, row, z, 1);
			} else {
				output.setFloat(col, row, nd, 1);
				adj.setFloat(col, row, nd, 1);
			}

		}
	}
}


void doInterp2(Raster& target, int tband, Raster& output, Raster& adj, geo::ds::KDTree<Pt>& tree, int minCount, int maxCount, double radius) {

	typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Mtx;

	const GridProps& props = target.props();
	int cols = props.cols();
	int rows = props.rows();
	double nd = props.nodata();

	std::vector<Pt*> pts;
	std::vector<double> dist;

	for(int row = 0; row < rows; ++row) {
		for(int col = 0; col < cols; ++col) {

			double v = target.getFloat(col, row, tband);
			if(v == nd) {
				output.setFloat(col, row, nd, 1);
				continue;
			}

			double maxRad = std::max(props.bounds().width(), props.bounds().height());
			double rad = radius;
			int count = 0;

			do {
				if(rad > maxRad) {
					std::cerr << "Exceeded max radius: " << rad << ", " << count << "\n";
					break;
				}
				pts.clear();
				dist.clear();
				Pt c(props.toX(col), props.toY(row), 0);
				count = tree.radSearch(c, rad, maxCount, std::back_inserter(pts), std::back_inserter(dist));
				rad *= 2;
			} while(count < minCount);

			if(count > 2) {
				Mtx mtx(3, pts.size());

				for(int i = 0; i < count; ++i) {
					if(dist[i] > 0) {
						mtx(0, i) = pts[i]->_x;
						mtx(1, i) = pts[i]->_y;
						mtx(2, i) = pts[i]->_z;
					}
				}

				// Recenter on zero.
				Eigen::Vector3d cent = mtx.rowwise().mean();
				mtx.colwise() -= cent;

				Eigen::JacobiSVD<Eigen::MatrixXd> svd(mtx, Eigen::ComputeThinU | Eigen::ComputeThinV);
				Eigen::Vector3d norm = svd.matrixU().col(2);
				double z = norm.dot(Eigen::Vector3d(0, 0, cent[2]));
				z = cent[2] > 0 ? std::abs(z) : -std::abs(z);

				output.setFloat(col, row, v - z, 1);
				adj.setFloat(col, row, z, 1);
			} else {
				output.setFloat(col, row, nd, 1);
				adj.setFloat(col, row, nd, 1);
			}

		}
	}
}

int main(int argc, char** argv) {

	if(argc < 6) {
		usage();
		return 1;
	}

	std::vector<std::string> anchors;
	std::vector<int> abands;
	std::string target;
	int tband = 1;
	std::string adjusted;
	std::string adjustment;
	int mode = 0;
	int count = 100;
	int maxCount = 4096;
	int skip = 0;
	double radius = 1000;

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
		} else if(arg == "-s") {
			skip = atoi(argv[++i]);
			mode = 0;
			continue;
		} else if(arg == "-r") {
			radius = atof(argv[++i]);
			mode = 0;
			continue;
		} else if(arg == "-m") {
			maxCount = atoi(argv[++i]);
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

	std::vector<std::pair<std::unique_ptr<Raster>, int> > arasters;
	for(size_t i = 0; i < anchors.size(); ++i)
		arasters.emplace_back(new Raster(anchors[i]), abands[i]);

	Raster traster(target);
	GridProps aprops(traster.props());
	aprops.setWritable(true);
	aprops.setBands(1);
	geo::ds::KDTree<Pt> tree(2);

	getDiffs(arasters, traster, tband, tree, skip);

	Raster output(adjusted, aprops);
	Raster adj(adjustment, aprops);

	doInterp(traster, tband, output, adj, tree, count);


 	return 0;
}
