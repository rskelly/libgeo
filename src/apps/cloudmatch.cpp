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
		std::string file;
		std::string projection;
		std::vector<double> grid;
		std::vector<double> weight;
		double xmin, xmax, ymin, ymax, zmin, zmax;
		int cols, rows;
		double res;

		Item(const std::string& file, double res) :
			file(file), res(res) {
			reset();
		}

		void unload() {
			grid.resize(0);
			weight.resize(0);
		}

		void reset() {
			xmin = ymin = zmin = std::numeric_limits<double>::max();
			ymax = ymax = zmax = std::numeric_limits<double>::lowest();
			cols = rows = 0;
		}

		void snap() {
			xmin = std::floor(xmin / res) * res;
			xmax = std::ceil(xmax / res) * res;
			ymin = std::floor(ymin / res) * res;
			ymax = std::ceil(ymax / res) * res;
			cols = (int) std::ceil(xmax - xmin) / res;
			rows = (int) std::ceil(ymax - ymin) / res;
		}

		void load() {
			LasFile las(file);
			Point pt;
			while(las.next(pt)) {
				if(pt.x < xmin) xmin = pt.x;
				if(pt.x > xmax) xmax = pt.x;
				if(pt.y < ymin) ymin = pt.y;
				if(pt.y > ymax) ymax = pt.y;
				if(pt.z < zmin) zmin = pt.z;
				if(pt.z > zmax) zmax = pt.z;
			}

			projection = las.projection();

			snap();

			grid.resize(cols * rows);
			weight.resize(cols * rows);
			fillWeight(0);
			fillGrid(0);

			int col, row;
			double w, s, d, cx, cy;
			las.reset();
			while(las.next(pt)) {
				if(pt.cls != 2)
					continue;
				col = (int) (pt.x - xmin) / res;
				row = (int) (pt.y - ymin) / res;
				for(int r = row - 1; r < row + 2; ++r) {
					for(int c = col - 1; c < col + 2; ++c) {
						if(c < 0 || r < 0 || r >= rows || c >= cols)
							continue;
						cx = xmin + col * res + res * 0.5;
						cy = ymin + col * res + res * 0.5;
						d = std::pow(cx - pt.x, 2.0) + std::pow(cy - pt.y, 2.0);
						if(d > res * res)
							continue;
						w = d / (res * res);
						grid[r * cols + c] += pt.z * w;
						weight[r * cols + c] += w;
					}
				}
			}

			for(size_t i = 0; i < grid.size(); ++i) {
				if(weight[i] > 0) {
					grid[i] /= weight[i];
				} else {
					grid[i] = NODATA;
				}
			}
		}

		void extend(Item& ref) {
			xmin = std::min(ref.xmin, xmin);
			ymin = std::min(ref.ymin, ymin);
			xmax = std::max(ref.xmax, xmax);
			ymax = std::max(ref.ymax, ymax);
			snap();
		}

		void fillGrid(double v) {
			std::fill(grid.begin(), grid.end(), v);
		}

		void fillWeight(double v) {
			std::fill(weight.begin(), weight.end(), v);
		}

		double toX(int col) {
			return xmin + col * res + res * 0.5;
		}

		double toY(int row) {
			return ymin + row * res + res * 0.5;
		}

		double toCol(double x) {
			return (int) (x - xmin) / res;
		}

		double toRow(double y) {
			return (int) (y - ymin) / res;
		}

		double get(int c, int r) {
			return grid[r * cols + c];
		}

		double get(double x, double y) {
			return get(toCol(x), toRow(y));
		}

		void set(int c, int r, double v) {
			grid[r * cols + r] = v;
		}

		void set(double x, double y, double v) {
			set(toCol(x), toRow(y), v);
		}

		void add(Item& other, double mult = 1) {
			int cc, rr;
			double v;
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c) {
					cc = other.toCol(toX(c));
					rr = other.toRow(toY(r));
					if((v = other.get(cc, rr)) != NODATA) {
						if(get(c, r) == NODATA) {
							set(c, r, v * mult);
						} else {
							set(c, r, get(c, r) + v * mult);
						}
					}
				}
			}
		}

		std::vector<double> m_x;
		std::vector<double> m_y;
		std::vector<double> m_z;
		std::vector<double> m_w;
		std::vector<double> m_tx;
		std::vector<double> m_ty;
		std::vector<double> m_c;
		double* m_wrk1;
		int m_lwrk1;
		double* m_wrk2;
		int m_lwrk2;
		int* m_iwrk;
		int m_liwrk;
		int m_nx;
		int m_ny;

		void smoothSpline(double sm = 0) {

			if(sm == 0) {
				for(double d : grid) {
					if(d != NODATA)
						sm += 1;
				}
			}

			int iopt = 0;
			int kx = 3;
			int ky = 3;
			double eps = std::numeric_limits<double>::min();

			m_x.resize(0);
			m_y.resize(0);
			m_z.resize(0);
			m_w.resize(0);

			int m = 0;
			double gs, gv;
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c) {
					if((gv = get(c, r)) != NODATA) {
						m_x.push_back(toX(c));
						m_y.push_back(toY(r));
						m_z.push_back(gv);
						m_w.push_back(1.0);
						++m;
					}
				}
			}

			int nxest = (int) std::ceil(kx + 1.0 + std::sqrt(m / 2.0));
			int nyest = (int) std::ceil(ky + 1.0 + std::sqrt(m / 2.0));
			int nmax = std::max(std::max(m, nxest), nyest);

			m_c.resize((nxest - kx - 1) * (nyest - ky - 1));
			m_tx.resize(m);
			m_ty.resize(m);

			double fp;

			int u = nxest - kx - 1;
			int v = nyest - ky - 1;
			int km = std::max(kx, ky) + 1;
			int ne = std::max(nxest, nyest);
			int bx = kx * v + ky + 1;
			int by = ky * u + kx +1;
			int b1, b2;
			if(bx <= by) {
				b1 = bx;
				b2 = b1 + v - ky;
			} else {
				b1 = by;
				b2 = b1 + u - kx;
			}

			if(m_wrk1)
				munmap(m_wrk1, m_lwrk1);
			m_lwrk1 = u * v * (2 + b1 + b2) + 2 * (u + v + km * (m + ne) + ne - kx - ky) + b2 + 1;
			m_wrk1 = (double*) mmap(0, m_lwrk1 * sizeof(double), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

			if(m_wrk2)
				munmap(m_wrk2, m_lwrk2);
			m_lwrk2 = u * v * (b2 + 1) + b2;
			m_wrk2 = (double*) mmap(0, m_lwrk2 * sizeof(double), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

			if(m_iwrk)
				munmap(m_iwrk, m_liwrk);
			m_liwrk = m + (nxest - 2 * kx - 1) * (nyest - 2 * ky - 1);
			m_iwrk = (int*) mmap(0, m_liwrk * sizeof(int), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

			int ier, iter = 0;
			double x0 = toX(0) - res / 2, x1 = toX(cols - 1) + res / 2;
			double y0 = toY(0) - res / 2, y1 = toY(rows - 1) + res / 2;
			do {
				surfit_(&iopt, &m, m_x.data(), m_y.data(), m_z.data(), m_w.data(),
						&x0, &x1, &y0, &y1, &kx, &ky,
						&sm, &nxest, &nyest, &nmax, &eps,
						&m_nx, m_tx.data(), &m_ny, m_ty.data(), m_c.data(), &fp,
						m_wrk1, (int*) &m_lwrk1, m_wrk2, (int*) &m_lwrk2, m_iwrk, (int*) &m_liwrk,	// Dangerous converstion to int.
						&ier);
				std::cerr << "ier : " << ier << "\n";
				if(ier == 1) {
					sm *= 10;
				} else if(ier == 4) {
					sm *= 10;
				} else if(ier < 0) {
					sm *= 0.1;
				} else if(ier > 5) {
					throw std::runtime_error("Smoothing failed");
				} else if(ier == 0) {
					break;
					std::cerr << "ier : " << ier << "\n";
				}
				++iter;
			} while(iter < 100);

		}

		void smooth() {

			int idim = 1;

			m_x.resize(cols); // A single row.
			m_y.resize(rows);

			for(int c = 0; c < cols; ++c)
				m_x[c] = toX(c);
			for(int r = 0; r < rows; ++r)
				m_y[r] = toY(r);

			int mf = cols * rows * idim;
			m_z.resize(mf);

			if(m_wrk1)
				munmap(m_wrk1, m_lwrk1);
			m_lwrk1 = cols * rows * 4;
			m_wrk1 = (double*) mmap(0, m_lwrk1 * sizeof(double), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

			if(m_iwrk)
				munmap(m_iwrk, m_liwrk);
			m_liwrk = cols * rows;
			m_iwrk = (int*) mmap(0, m_liwrk * sizeof(int), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

			int ier;

			surev_(&idim, m_tx.data(), &m_nx, m_ty.data(), &m_ny,
					m_c.data(), m_x.data(), &cols, m_y.data(), &rows, m_z.data(), &mf,
					m_wrk1, (int*) &m_lwrk1, m_iwrk, (int*) &m_liwrk, &ier);

			fillGrid(NODATA);
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c)
					set(c, r, m_z[c * rows + r]); // Note: z is transposed from usual.
			}
		}

		void adjustTo(Item& ref, const std::string& outfile) {

			Item out(outfile, res);
			out.extend(*this);
			out.extend(ref);
			out.fillGrid(NODATA);
			out.add(ref);
			out.add(*this, -1);
			out.smoothSpline();
			out.smooth();
			out.save(grid, "smooth.tif");
			/*
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
			*/
		}

		void save(std::vector<double>& data, const std::string& outfile) {
			GDALAllRegister();
			GDALDriverManager* dm = GetGDALDriverManager();
			GDALDriver* drv = dm->GetDriverByName("GTiff");
			GDALDataset* ds = drv->Create(outfile.c_str(), cols, rows, 1, GDT_Float32, 0);
			const char* proj = projection.c_str();
			ds->SetProjection(proj);
			GDALRasterBand* b = ds->GetRasterBand(1);
			b->SetNoDataValue(NODATA);
			if(CE_None != b->RasterIO(GF_Write, 0, 0, cols, rows, data.data(), cols, rows, GDT_Float32, 0, 0, 0))
				throw std::runtime_error("Failed to write band 1.");
			GDALClose(ds);
		}

		~Item() {
			if(m_wrk1)
				munmap(m_wrk1, m_lwrk1);
			if(m_wrk2)
				munmap(m_wrk2, m_lwrk2);
			if(m_iwrk)
				munmap(m_iwrk, m_liwrk);
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

	Item ref(reffile, 20);
	ref.load();

	Item adj(adjfile, 20);
	adj.load();

	adj.adjustTo(ref, outfile);
}
