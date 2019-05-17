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

#include <liblas/liblas.hpp>
#include <vector>
#include <fstream>
#include <sstream>

#include <gdal_priv.h>
#include <ogr_spatialref.h>

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


class Item {
public:
	std::string file;
	double xmin, xmax, ymin, ymax;
	double res;
	int cols, rows;
	std::vector<double> grid;
	std::vector<int> counts;
	std::string projection;

	Item(const std::string& file, double xmin, double ymin, double xmax, double ymax, double res) :
		file(file),
		xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax),
		res(res) {
		calcBounds();
	}

	Item(const std::string& file) : Item(file, MAX, MAX, MIN, MIN, 0) {}

	Item() : Item("") {}

	Item(const Item& item) : Item(item.file, item.xmin, item.ymin, item.xmax, item.ymax, item.res) {
		grid.assign(item.grid.begin(), item.grid.end());
		counts.assign(item.counts.begin(), item.counts.end());
	}

	void extend(Item& other) {
		if(other.xmin < xmin) xmin = other.xmin;
		if(other.ymin < ymin) ymin = other.ymin;
		if(other.xmax > xmax) xmax = other.xmax;
		if(other.ymax > ymax) ymax = other.ymax;
	}

	void reset() {
		xmin = ymin = MAX;
		xmax = ymax = MIN;
	}

	void calcBounds() {
		xmin = std::floor(xmin / res) * res - res;
		ymin = std::floor(ymin / res) * res - res;
		xmax = std::ceil(xmax / res) * res + res;
		ymax = std::ceil(ymax / res) * res + res;
		cols = (int) ((xmax - xmin) / res);
		rows = (int) ((ymax - ymin) / res);
		grid.resize(cols * rows);
		counts.resize(cols * rows);
	}

	void load() {
		std::ifstream in(file);
		liblas::ReaderFactory rf;
		liblas::Reader rdr = rf.CreateWithStream(in);
		const liblas::Header& hdr = rdr.GetHeader();

		projection = hdr.GetSRS().GetWKT(liblas::SpatialReference::WKTModeFlag::eCompoundOK);

		reset();
		while(rdr.ReadNextPoint()) {
			const liblas::Point& pt = rdr.GetPoint();
			double x = pt.GetX();
			double y = pt.GetY();
			if(x < xmin) xmin = x;
			if(y < ymin) ymin = y;
			if(x > xmax) xmax = x;
			if(y > ymax) ymax = y;
		}

		calcBounds();

		fill(0);
		fillCounts(0);

		while(rdr.ReadNextPoint()) {
			const liblas::Point& pt = rdr.GetPoint();
			if(pt.GetClassification().GetClass() != 2)
				continue;
			double x = pt.GetX();
			double y = pt.GetY();
			double z = pt.GetZ();
			set(x, y, get(x, y) + z);
			setCount(x, y, getCount(x, y) + 1);
		}

		for(size_t i = 0; i < grid.size(); ++i) {
			if(counts[i] == 0)
				grid[i] = -9999.0;
		}

	}

	void fill(double v) {
		std::fill(grid.begin(), grid.end(), v);
	}

	void fillCounts(int v) {
		std::fill(counts.begin(), counts.end(), v);
	}

	void add(Item& other) {
		if(res != other.res)
			throw std::runtime_error("Resolutions do not match.");
		double x1 = std::max(xmin, other.xmin);
		double x2 = std::min(xmax, other.xmax);
		double y1 = std::max(ymin, other.ymin);
		double y2 = std::min(ymax, other.ymax);

		for(double y = y1; y < y2; y += res) {
			for(double x = x1; x < x2; x += res) {
				if(other.get(x, y) == -9999.0)
					continue;
				if(get(x, y) == -9999.0) {
					set(x, y, other.get(x, y));
					setCount(x, y, 1);
				} else {
					set(x, y, get(x, y) + other.get(x, y));
					setCount(x, y, getCount(x, y) + 1);
				}
			}
		}
	}

	void diff(Item& other) {
		std::vector<double> tmp(grid.size());
		std::fill(tmp.begin(), tmp.end(), -9999.0);
		double z, oz;
		for(double y = ymin; y < ymax; y += res) {
			for(double x = xmin; x < xmax; x += res) {
				if((oz = other.get(x, y)) == -9999.0 || (z = get(x, y)) == -9999.0)
					continue;
				int c = toCol(x);
				int r = toRow(y);
				tmp[r * cols + c] = oz - z;
			}
		}
		grid.swap(tmp);
	}

	void average() {
		for(size_t i = 0; i < grid.size(); ++i) {
			if(counts[i] > 0)
				grid[i] /= counts[i];
		}
	}

	int toCol(double x) {
		return (int) ((x - xmin) / res);
	}

	int toRow(double y) {
		return (int) ((y - ymin) / res);
	}

	int toX(int col) {
		return col * res + xmin + res * 0.5;
	}

	int toY(int row) {
		return row * res + ymin + res * 0.5;
	}

	void set(double x, double y, double v) {
		set(toCol(x), toRow(y), v);
	}

	void set(int col, int row, double v) {
		if(col >= 0 && col < cols && row >= 0 && row < rows)
			grid[row * cols + col] = v;
	}

	void setCount(double x, double y, int v) {
		setCount(toCol(x), toRow(y), v);
	}

	void setCount(int col, int row, int v) {
		if(col >= 0 && col < cols && row >= 0 && row < rows)
			counts[row * cols + col] = v;
	}

	double get(double x, double y) {
		return get(toCol(x), toRow(y));
	}

	double get(int col, int row) {
		if(col >= 0 && col < cols && row >= 0 && row < rows) {
			return grid[row * cols + col];
		} else {
			return -9999.0;
		}
	}

	int getCount(double x, double y) {
		return getCount(toCol(x), toRow(y));
	}

	int getCount(int col, int row) {
		if(col >= 0 && col < cols && row >= 0 && row < rows) {
			return counts[row * cols + col];
		} else {
			return 0;
		}
	}

	void splineSmooth(double smooth) {

		int iopt = 0;
		int kx = 3;
		int ky = 3;
		double eps = std::numeric_limits<double>::min();

		std::vector<double> x;
		std::vector<double> y;
		std::vector<double> z;
		std::vector<double> w;

		int m = 0;
		for(double yy = ymin; yy < ymax; yy += res) {
			for(double xx = xmin; xx < xmax; xx += res) {
				double zz = get(xx, yy);
				if(zz != -9999.0) {
					x.push_back(xx);
					y.push_back(yy);
					z.push_back(zz);
					w.push_back(1);
					++m;
				}
			}
		}

		int nmax = m;
		int nxest = (int) std::ceil(kx + 1.0 + std::sqrt(m / 2.0));
		int nyest = (int) std::ceil(ky + 1.0 + std::sqrt(m / 2.0));

		std::vector<double> c((nxest - kx - 1) * (nyest - ky - 1));
		std::vector<double> tx(m);
		std::vector<double> ty(m);

		int nx;
		int ny;
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

		int lwrk1 = u * v * (2 + b1 + b2) + 2 * (u + v + km * (m + ne) + ne - kx - ky) + b2 + 1;
		std::vector<double> wrk1(lwrk1);
		int lwrk2 = u * v * (b2 + 1) + b2;
		std::vector<double> wrk2(lwrk2);
		int kwrk = m + (nxest - 2 * kx - 1) * (nyest - 2 * ky - 1);
		std::vector<int> iwrk(kwrk);

		int ier;

		surfit_(&iopt, &m, x.data(), y.data(), z.data(), w.data(),
				&xmin, &xmax, &ymin, &ymax, &kx, &ky,
				&smooth, &nxest, &nyest, &nmax, &eps,
				&nx, tx.data(), &ny, ty.data(), c.data(), &fp,
				wrk1.data(), &lwrk1, wrk2.data(), &lwrk2, iwrk.data(), &kwrk,
				&ier);

		if(ier > 0)
			throw std::runtime_error("Smoothing failed");

		int idim = 1;

		int mf = cols * rows * idim;

		x.resize(cols); // A single row.
		y.resize(rows);
		z.resize(mf);

		lwrk1 = cols * rows * 4;
		wrk1.resize(lwrk1);
		kwrk = cols * rows;
		iwrk.resize(kwrk);

		surev_(&idim, tx.data(), &nx, ty.data(), &ny,
				c.data(), x.data(), &cols, y.data(), &rows, z.data(), &mf,
				wrk1.data(), &lwrk1, iwrk.data(), &kwrk, &ier);

		fill(-9999.0);
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c) {
				set(c, r, z[r * cols + c]);
			}
		}
	}

	void save(const std::string& outfile) {
		GDALAllRegister();
		GDALDriverManager* dm = GetGDALDriverManager();
		GDALDriver* drv = dm->GetDriverByName("GTiff");
		GDALDataset* ds = drv->Create(outfile.c_str(), cols, rows, 2, GDT_Float32, 0);
		const char* proj = projection.c_str();
		ds->SetProjection(proj);
		GDALRasterBand* b1 = ds->GetRasterBand(1);
		b1->SetNoDataValue(0);
		if(CE_None != b1->RasterIO(GF_Write, 0, 0, cols, rows, counts.data(), cols, rows, GDT_Int32, 0, 0, 0))
			throw std::runtime_error("Failed to write band 1.");
		GDALRasterBand* b2 = ds->GetRasterBand(2);
		b2->SetNoDataValue(0);
		if(CE_None != b2->RasterIO(GF_Write, 0, 0, cols, rows, grid.data(), cols, rows, GDT_Float64, 0, 0, 0))
			throw std::runtime_error("Failed to write band 2.");
	}

	void adjust(const std::string& outfile) {

		std::ofstream out(outfile);
		std::ifstream in(file);

		liblas::ReaderFactory rf;
		liblas::Reader rdr = rf.CreateWithStream(in);
		const liblas::Header& hdr = rdr.GetHeader();

		liblas::WriterFactory wf;
		liblas::Writer wtr = wf.CreateWithStream(out, hdr);

		while(rdr.ReadNextPoint()) {
			const liblas::Point& pt = rdr.GetPoint();
			liblas::Point npt(pt);
			double x = pt.GetX();
			double y = pt.GetY();
			npt.SetZ(pt.GetZ() + get(x, y));
			wtr.WritePoint(npt);
		}

	}
};

std::string makeFile(const std::string& outdir, const std::string& tpl, double res, int i, const std::string& ext) {
	std::stringstream ss;
	ss << outdir << "/" << tpl << "_" << res << "_" << i << ext;
	return ss.str();
}

int main(int argc, char** argv) {

	double res = atof(argv[argc - 3]);
	int steps = atoi(argv[argc - 2]);
	std::string outdir = argv[argc - 1];
	std::vector<Item> items;

	double res0 = res * std::pow(2, steps);

	std::vector<std::string> las;
	for(int i = 1; i < argc - 2; ++i)
		las.push_back(argv[i]);

	while(res0 >= res) {

		for(const std::string& l : las)
			items.emplace_back(l);

		Item avg;
		avg.projection = items[0].projection;

		for(size_t i = 0; i < items.size(); ++i) {
			Item& item = items[i];
			item.load();
			item.save(makeFile(outdir, "item", res, i, ".tif"));
			avg.extend(item);
			avg.add(item);
		}

		avg.average();
		avg.save(makeFile(outdir, "avg", res, 0, ".tif"));

		std::string s;
		std::vector<std::string> slas;

		for(size_t i = 0; i < items.size(); ++i) {
			Item& item = items[i];
			item.diff(avg);
			item.splineSmooth(10);
			s = makeFile(outdir, "adj", res, i, ".las");
			slas.push_back(s);
			item.adjust(s);
		}

		res /= 2;
		las.swap(slas);
	}
}
