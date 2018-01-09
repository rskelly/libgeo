/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include <fstream>
#include <vector>
#include <unordered_map>

#include <liblas/liblas.hpp>

#include "raster.hpp"
#include "ds/kdtree.hpp"

class LASFile {
public:
	std::string filename;
	double bounds[4];

	LASFile(const std::string& filename) :
		filename(filename) {

		init();
	}

	void init() {
		std::ifstream str(filename, std::ios::in | std::ios::binary);
		liblas::ReaderFactory f;
		liblas::Reader reader = f.CreateWithStream(str);
		const liblas::Header& header = reader.GetHeader();
		bounds[0] = header.GetMinX();
		bounds[1] = header.GetMinY();
		bounds[2] = header.GetMaxX();
		bounds[3] = header.GetMaxY();
	}

	~LASFile() {
	}
};

class LASPoint {
public:
	double x; 
	double y;
	double z;

	LASPoint(const liblas::Point& pt) :
		x(pt.GetX()), y(pt.GetY()), z(pt.GetZ()) {
	}

	LASPoint(double x, double y, double z) :
		x(x), y(y), z(z) {
	}

	double operator[](int idx) const {
		switch(idx % 2) {
		case 0:
			return x;
		case 1:
			return y;
		}
		return 0;
	}

	~LASPoint() {
	}
};

using namespace geo::raster;

class Rasterizer {
public:
	std::vector<LASFile> files;
	std::unordered_set<int> currentFiles; // Files currently in tree.
	geo::ds::KDTree<LASPoint>* tree;

	Rasterizer(const std::vector<std::string> filenames) :
		tree(nullptr) {
		for(const std::string& filename : filenames)
			files.emplace_back(filename);

	}

	void updateTree(double x, double y, double radius) {
		std::unordered_set<int> requiredFiles;
		double tbounds[4] = {x - radius, y - radius, x + radius, y + radius};
		for(size_t i = 0; i < files.size(); ++i) {
			if(!(tbounds[2] < files[i].bounds[0] || tbounds[0] > files[i].bounds[2] ||
				tbounds[3] < files[i].bounds[1] || tbounds[1] > files[i].bounds[3])) {
				requiredFiles.insert(i);
			}
		}
		bool changed = false;
		if(currentFiles.empty()) {
			currentFiles.insert(requiredFiles.begin(), requiredFiles.end());
			changed = true;
		} else {
			for(int a : currentFiles) {
				if(requiredFiles.find(a) == requiredFiles.end()) {
					changed = true;
					break;
				}
			}
			for(int a : requiredFiles) {
				if(currentFiles.find(a) == currentFiles.end()) {
					changed = true;
					break;
				}
			}
			if(changed) {
				currentFiles.clear();
				currentFiles.insert(requiredFiles.begin(), requiredFiles.end());
			}
		}
		if(changed) {
			if(tree) {
				delete tree;
				tree = nullptr;
			}
			if(!currentFiles.empty()) {
				liblas::ReaderFactory fact;
				tree = new geo::ds::KDTree<LASPoint>(2);
				for(int a : currentFiles) {
					std::ifstream str(files[a].filename);
					liblas::Reader rdr = fact.CreateWithStream(str);
					while(rdr.ReadNextPoint()) {
						const liblas::Point& pt = rdr.GetPoint();
						if(pt.GetClassification().GetClass() == 2) {
							LASPoint lpt(pt);
							tree->add(lpt);
						}
					}
				}
				tree->build();
			}
		}
	}

	int getPoints(double x, double y, double radius, int count, std::list<LASPoint>& pts, std::list<double>& dists) {
		updateTree(x, y, radius);
		LASPoint pt(x, y, 0);
		int ret = tree->radSearch(pt, radius, count, std::back_inserter(pts), std::back_inserter(dists));
		return ret;
	}

	double compute(const std::string& type, const std::list<LASPoint>& pts, const std::list<double>& dists) {
		double e = 0;
		double c = 0;
		for(const LASPoint& pt : pts) {
			e += pt.z;
			++c;
		}
		return c > 0 ? e / c : 0;
	}

	void rasterize(const std::string& filename, const std::string& type, double res, 
		double easting, double northing, double radius, int srid, int threads, double percentile = 0) {

		double bounds[4] = {9999999999, 9999999999, -9999999999, -9999999999};

		for(const LASFile& f : files) {
			if(f.bounds[0] < bounds[0]) bounds[0] = f.bounds[0];
			if(f.bounds[1] < bounds[1]) bounds[1] = f.bounds[1];
			if(f.bounds[2] > bounds[2]) bounds[2] = f.bounds[2];
			if(f.bounds[3] > bounds[3]) bounds[3] = f.bounds[3];
		}

		if(easting <= 0)
			easting = ((int) (bounds[0] * res)) / res;
		if(northing <= 0)
			northing = ((int) (bounds[3] * res)) / res;

		int cols = (int) ((bounds[2] - bounds[0]) / res) + 1;
		int rows = (int) ((bounds[3] - bounds[1]) / res) + 1;

		GridProps props;
		props.setTrans(easting, res, northing, -res);
		props.setSize(cols, rows);
		props.setNoData(-9999.0);
		props.setDataType(DataType::Float32);
		props.setSrid(srid);
		props.setWritable(true);
		Raster rast(filename, props);

		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c) {
				int count = 100;
				double x = props.toX(c) + props.resolutionX() / 2;
				double y = props.toY(r) + props.resolutionY() / 2;
				std::list<LASPoint> pts;
				std::list<double> dists;
				int ret;
				while((ret = getPoints(x, y, radius, count, pts, dists)) >= count) {
					pts.clear();
					dists.clear();
					count *= 2;
					std::cerr << "count " << count << "; ret " << ret << "\n";
				}
				if(ret) {
					rast.setFloat(c, r, compute(type, pts, dists));
				} else {
					rast.setFloat(c, r, -9999.0);
				}
			}
		}
	}

	~Rasterizer() {
		if(tree)
			delete tree;
	}
};

void usage() {
	std::cerr << "Usage: las2grid <output raster> <input las [*]>\n"
			<< " -r <resolution> The output resolution in map units.\n"
			<< " -e <easting>    The top left corner horizontal alignment\n"
			<< "                 (defaults to nearest whole number).\n"
			<< " -n <northing>   The top left corner vertical alignment\n"
			<< "                 (defaults to nearest whole number).\n"
			<< " -d <radius>     The search radius for finding points.\n"
			<< " -s <srid>       The spatial reference ID. Default 0.\n"
			<< " -p <type>       The type of raster: mean (default), median, min, max, \n"
			<< "                 rugosity, variance, std. deviation and percentile.\n"
			<< "                 For percentile, use the form, 'percenile:n', where\n"
			<< "                 n is the percentile (no % sign); 1 - 99.\n"
			<< " -t <t>          The number of threads. Default 1.\n";
}

int main(int argc, char** argv) {

	using namespace geo::raster;

	if(argc < 3) {
		usage();
		return 1;
	}

	double res = 0;
	double easting = 0;
	double northing = 0;
	double radius = 0;
	uint16_t srid = 0;
	uint16_t threads = 1;
	std::string type = "mean";
	std::vector<std::string> args;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-p") {
			type = argv[++i];
		} else if(v == "-r") {
			res = atof(argv[++i]);
		} else if(v == "-d") {
			radius = atof(argv[++i]);
		} else if(v == "-s") {
			srid = atoi(argv[++i]);
		} else if(v == "-e") {
			easting = atof(argv[++i]);
		} else if(v == "-n") {
			northing = atof(argv[++i]);
		} else if(v == "-t") {
			threads = atoi(argv[++i]);
		} else {
			args.push_back(argv[i]);
		}
	}

	if(res <= 0) {
		std::cerr << "Resolution must be >0.\n";
		usage();
		return 1;
	}

	if(args.size() < 2) {
		std::cerr << "Input and output filenames required.\n";
		usage();
		return 1;
	}

	if(threads <= 0) {
		std::cerr << "Illegal thread value. Defaulting to 1.\n";
		threads = 1;
	}

	std::vector<std::string> infiles(args.begin() + 1, args.end());
	Rasterizer r(infiles);
	r.rasterize(args[0], type, res, easting, northing, radius, srid, threads, 0);

	return 0;
}


