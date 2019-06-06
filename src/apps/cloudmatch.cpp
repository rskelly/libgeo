/*
 * cloudmatch.cpp
 *
 *  Created on: May 16, 2019
 *      Author: rob
 *
 * 1) Create a grid at resoluion r for each point cloud -- mean of ground points, etc.
 * 2) Calculate the mean elevation at each cell where cells overlap.
 * 3) Build a smooth spline of differences, vertically offset all points by the local deviation at their position.
 * 4) Half resolution, goto 1. Continue until limit.
 */

#include <vector>
#include <fstream>
#include <sstream>
#include <dirent.h>
#include <sys/mman.h>

#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/io/LasReader.hpp>
#include <pdal/io/LasWriter.hpp>
#include <pdal/io/LasHeader.hpp>
#include <pdal/Options.hpp>

#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include "ds/kdtree.hpp"

constexpr int THIN = 1000;
constexpr double NODATA = -9999.0;
constexpr double MIN = std::numeric_limits<double>::lowest();
constexpr double MAX = std::numeric_limits<double>::max();

extern "C" {

	void surfit_(int* iopt, int* m, double* x, double* y, double* z, double* w,
			double* xb, double* xe, double* yb, double* ye, int* kx, int* ky,
			double* s, int* nxest, int* nyest, int* nmax, double* eps,
			int* nx, double* tx, int* ny, double* ty, double* c, double* fp,
			double* wrk1, int* lwrk1, double* wrk2, int* lwrk2, int* iwrk, int* kwrk,
			int* ier);

	void surev_(int* idim, double* tu, int* nu, double* tv, int* nv,
			double* c, double* u, int* mu, double* v, int* mv, double* f, int* mf,
			double* wrk, int* lwrk, int* iwrk, int* kwrk, int* ier);
}

namespace {
	class Kpt {
	public:
		double x, y, z;
		size_t i;
		Kpt(double x, double y, double z, size_t i) :
			x(x), y(y), z(z), i(i) {}
		double operator[](int idx) const {
			if(idx % 2 == 0) {
				return x;
			} else {
				return y;
			}
		}
	};

	class Point {
	public:
		double x, y, z;
		int cls;
	};

	class LasFile {
	public:
		pdal::PointTable m_table;
		pdal::LasReader m_reader;
		pdal::PointViewSet m_viewSet;
		pdal::PointViewPtr m_view;
		pdal::Dimension::IdList m_dims;
		pdal::LasHeader m_header;
		pdal::PointId m_ridx;
		pdal::PointId m_widx;

		LasFile(const std::string& infile) {
			pdal::Option opt("filename", infile);
			pdal::Options opts;
			opts.add(opt);
			m_reader.setOptions(opts);
			m_reader.prepare(m_table);
			m_viewSet = m_reader.execute(m_table);
			m_view = *m_viewSet.begin();
			m_dims = m_view->dims();
			m_header = m_reader.header();
			m_ridx = 0;
			m_widx = 0;
		}

		void write(const std::string& outfile) {
			pdal::Option opt("filename", outfile);
			pdal::Options opts;
			opts.add(opt);
			pdal::LasWriter writer;
			writer.setOptions(opts);
			writer.setSpatialReference(m_reader.getSpatialReference());
			writer.prepare(m_table);
			writer.execute(m_table);
		}

		void reset() {
			m_ridx = 0;
			m_widx = 0;
		}

		bool next(Point& pt) {
			if(m_ridx >= m_view->size())
				return false;
			using namespace pdal::Dimension;
			pt.x = m_view->getFieldAs<double>(Id::X, m_ridx);
			pt.y = m_view->getFieldAs<double>(Id::Y, m_ridx);
			pt.z = m_view->getFieldAs<double>(Id::Z, m_ridx);
			pt.cls = m_view->getFieldAs<int>(Id::Classification, m_ridx);
			++m_ridx;
			return true;
		}

		bool update(Point& pt) {
			if(m_widx >= m_view->size())
				return false;
			using namespace pdal::Dimension;
			m_view->setField(Id::X, m_widx, pt.x);
			m_view->setField(Id::Y, m_widx, pt.y);
			m_view->setField(Id::Z, m_widx, pt.z);
			m_view->setField(Id::Classification, m_widx, pt.cls);
			++m_widx;
			return true;
		}

		std::string projection() {
			return m_table.spatialReference().getWKT();
		}
	};

	class Item {
	public:
		std::unique_ptr<geo::ds::KDTree<Kpt>> tree;
		std::string file;

		Item(const std::string& file) : file(file) {}

		void unload() {
			if(tree.get())
				tree->destroy();
		}

		void load() {
			tree.reset(new geo::ds::KDTree<Kpt>(2));
			LasFile las(file);
			Point pt;
			size_t i = 0;
			while(las.next(pt)) {
				if(pt.cls == 2)
					tree->add(new Kpt(pt.x, pt.y, pt.z, i));
				++i;
			}
			tree->build();
		}

		double wavg(std::vector<Kpt*> pts, std::vector<double> dist, double rad) {
			double s = 0, w = 0;
			for(size_t i  = 0; i < pts.size(); ++i) {
				s += pts[0]->z * (dist[i] / rad);
				w += dist[i] / rad;
			}
			return w > 0 ? s / w : std::nan("");
		}

		void adjustTo(Item& ref, const std::string& outfile) {
			std::vector<Kpt*> a;
			std::vector<double> da;
			std::vector<Kpt*> b;
			std::vector<double> db;
			LasFile las(file);
			Point pt;
			size_t i = 0;
			double rad = 100;
			while(las.next(pt)) {
				Kpt kpt(pt.x, pt.y, pt.z, i);
				tree->radSearch(kpt, rad, 0, std::back_inserter(a), std::back_inserter(da));
				ref.tree->radSearch(kpt, rad, 0, std::back_inserter(b), std::back_inserter(db));
				double za = wavg(a, da, rad);
				double zb = wavg(b, db, rad);
				if(std::isnan(za) || std::isnan(zb)) {
					std::cerr << "nan\n";
					pt.cls = 7;
				} else {
					pt.z += (zb - za);
				}
				las.update(pt);
			}
			las.write(outfile);
		}

		~Item() {
		}
	};

	std::string makeFile(const std::string& outdir, const std::string& tpl, double res, int i, const std::string& ext) {
		std::stringstream ss;
		ss << outdir << "/" << tpl << "_" << res << "_" << i << ext;
		return ss.str();
	}

	std::vector<std::string> getFiles(const std::string& dirname) {
		std::vector<std::string> files;
		DIR *dir;
		struct dirent *ent;
		if ((dir = opendir (dirname.c_str())) != NULL) {
		  /* print all the files and directories within directory */
		  while ((ent = readdir (dir)) != NULL) {
			  std::string file = dirname + "/" + ent->d_name;
			  if(file.substr(file.size() - 4, std::string::npos) == ".las")
				  files.push_back(file);
		  }
		  closedir (dir);
		} else {
		  /* could not open directory */
		  perror ("");
		}
		return files;
	}
}

int main(int argc, char** argv) {

	std::string outfile;
	std::string adjfile;
	std::string reffile;

	for(int i = 1; i < argc - 3; ++i) {
		std::string arg = argv[i];
		if(arg == "-a") {
			adjfile = argv[++i];
			continue;
		} else if(arg == "-r") {
			reffile = argv[++i];
			continue;
		} else if(arg == "-o") {
			outfile = argv[++i];
		}
	}

	Item ref(reffile);
	ref.load();

	Item adj(adjfile);
	adj.load();

	adj.adjustTo(ref, outfile);
}
