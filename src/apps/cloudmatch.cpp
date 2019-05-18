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

constexpr double NO_DATA = -9999.0;
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

	std::vector<double> m_wrk1;
	std::vector<double> m_wrk2;
	std::vector<int> m_iwrk;

	std::vector<double> m_x;
	std::vector<double> m_y;
	std::vector<double> m_z;
	std::vector<double> m_w;

	std::vector<double> m_c;
	std::vector<double> m_tx;
	std::vector<double> m_ty;

	int m_nx;
	int m_ny;

	Item(const std::string& file, double xmin, double ymin, double xmax, double ymax, double res) :
		file(file),
		xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax),
		res(res),
		cols(0), rows(0),
		m_nx(0), m_ny(0) {
	}

	Item(const std::string& file) : Item(file, MAX, MAX, MIN, MIN, 0) {}

	Item(const std::string& file, double res) : Item(file, MAX, MAX, MIN, MIN, res) {}

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

		rdr.Reset();
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
				grid[i] = NO_DATA;
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
				if(other.get(x, y) == NO_DATA)
					continue;
				if(get(x, y) == NO_DATA) {
					set(x, y, other.get(x, y));
					setCount(x, y, other.getCount(x, y));
				} else {
					set(x, y, get(x, y) + other.get(x, y));
					setCount(x, y, getCount(x, y) + other.getCount(x, y));
				}
			}
		}
	}

	void diff(Item& other) {
		std::vector<double> tmp(grid.size());
		std::fill(tmp.begin(), tmp.end(), NO_DATA);
		double z, oz;
		for(double y = ymin; y < ymax; y += res) {
			for(double x = xmin; x < xmax; x += res) {
				if((oz = other.get(x, y)) == NO_DATA || (z = get(x, y)) == NO_DATA)
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
			if(counts[i] > 0) {
				grid[i] /= counts[i];
			} else {
				grid[i] = NO_DATA;
			}
		}
	}

	int toCol(double x) {
		return (int) ((x - xmin) / res);
	}

	int toRow(double y) {
		return (int) ((ymax - y) / res);
	}

	int toX(int col) {
		return col * res + xmin + res * 0.5;
	}

	int toY(int row) {
		return ymax - row * res - res * 0.5;
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
			return NO_DATA;
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

	void splineSmooth(double& smooth) {

		m_x.resize(0);
		m_y.resize(0);
		m_z.resize(0);
		m_w.resize(0);

		int m = 0;
		for(double yy = ymin; yy <= ymax; yy += res) {										// Include the corners to be sure the spline covers all points.
			for(double xx = xmin; xx <= xmax; xx += res) {
				double zz = get(xx == xmax ? xx - res : xx, yy == ymax ? yy - res : yy);	// The right/top corner takes the value of the cell with the previous index.
				if(zz != NO_DATA) {
					m_x.push_back(xx);
					m_y.push_back(yy);
					m_z.push_back(zz);
					m_w.push_back(1);
					++m;
				}
			}
		}

		int iopt = 0;
		int kx = 3;
		int ky = 3;
		int nmax = m;
		int nxest = (int) std::ceil(kx + 1.0 + std::sqrt(m / 2.0));
		int nyest = (int) std::ceil(ky + 1.0 + std::sqrt(m / 2.0));
		double eps = std::pow(10.0, -6.0); //std::numeric_limits<double>::min();

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

		int lwrk1 = u * v * (2 + b1 + b2) + 2 * (u + v + km * (m + ne) + ne - kx - ky) + b2 + 1;
		m_wrk1.resize(lwrk1);
		int lwrk2 = u * v * (b2 + 1) + b2;
		m_wrk2.resize(lwrk2);
		int kwrk = m + (nxest - 2 * kx - 1) * (nyest - 2 * ky - 1);
		m_iwrk.resize(kwrk);

		int ier;
		int tries = 4;

		while(--tries >= 0) {
			surfit_(&iopt, &m, m_x.data(), m_y.data(), m_z.data(), m_w.data(),
					&xmin, &xmax, &ymin, &ymax, &kx, &ky,
					&smooth, &nxest, &nyest, &nmax, &eps,
					&m_nx, m_tx.data(), &m_ny, m_ty.data(), m_c.data(), &fp,
					m_wrk1.data(), &lwrk1, m_wrk2.data(), &lwrk2, m_iwrk.data(), &kwrk,
					&ier);
			if(ier == 0)
				break;
			smooth *= 2;
		}

		if(ier > 0)
			throw std::runtime_error("Smoothing failed");
	}

	void smoothGrid() {

		int idim = 1;
		int mf = cols * rows * idim;

		m_x.resize(mf); // A single col/row.
		m_y.resize(mf);
		m_z.resize(mf);

		// Populate the cell-centre coordinates.
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c) {
				m_x[r * cols + c] = xmin + c * res + res * 0.5;
				m_y[r * cols + c] = ymin + r * res + res * 0.5;
			}
		}

		int lwrk1 = cols * rows * 4;
		m_wrk1.resize(lwrk1);
		int kwrk = cols * rows;
		m_iwrk.resize(kwrk);

		int ier;

		surev_(&idim, m_tx.data(), &m_nx, m_ty.data(), &m_ny,
				m_c.data(), m_x.data(), &cols, m_y.data(), &rows, m_z.data(), &mf,
				m_wrk1.data(), &lwrk1, m_iwrk.data(), &kwrk, &ier);

		// Apply smoothed values to raster.
		fill(NO_DATA);
		for(int i = 0; i < mf; ++i)
			set(m_x[i], m_y[i], m_z[i]);
	}

	double smoothAt(double x, double y) {

		int idim = 1;
		int mx = 1;
		int my = 1;
		int mf = mx * my * idim;

		m_x.resize(mx); // A single col/row.
		m_y.resize(my);
		m_z.resize(mf);

		m_x[0] = x;
		m_y[0] = y;

		int lwrk1 = 4 * (mx + my);
		m_wrk1.resize(lwrk1);
		int kwrk = mx + my;
		m_iwrk.resize(kwrk);

		int ier;

		surev_(&idim, m_tx.data(), &m_nx, m_ty.data(), &m_ny,
				m_c.data(), m_x.data(), &mx, m_y.data(), &my, m_z.data(), &mf,
				m_wrk1.data(), &lwrk1, m_iwrk.data(), &kwrk, &ier);

		if(ier == 10) // TODO: What about other non-zero returns?
			throw std::runtime_error("Failed to find smoothed value.");

		return m_z[0];
	}

	void save(const std::string& outfile) {
		GDALDriverManager* dm = GetGDALDriverManager();
		GDALDriver* drv = dm->GetDriverByName("GTiff");
		GDALDataset* ds = drv->Create(outfile.c_str(), cols, rows, 2, GDT_Float32, 0);
		double trans[] = {xmin, res, 0, ymax, 0, -res};
		ds->SetGeoTransform(trans);
		const char* proj = projection.c_str();
		ds->SetProjection(proj);
		GDALRasterBand* b1 = ds->GetRasterBand(1);
		b1->SetNoDataValue(NO_DATA);
		if(CE_None != b1->RasterIO(GF_Write, 0, 0, cols, rows, counts.data(), cols, rows, GDT_Int32, 0, 0, 0))
			throw std::runtime_error("Failed to write band 1.");
		GDALRasterBand* b2 = ds->GetRasterBand(2);
		b2->SetNoDataValue(NO_DATA);
		if(CE_None != b2->RasterIO(GF_Write, 0, 0, cols, rows, grid.data(), cols, rows, GDT_Float64, 0, 0, 0))
			throw std::runtime_error("Failed to write band 2.");
		GDALClose(ds);
	}

	void adjust(const std::string& outfile, Item& reference) {

		std::ofstream out(outfile, std::ios::binary|std::ios::out);
		std::ifstream in(file, std::ios::binary|std::ios::in);

		liblas::ReaderFactory rf;
		liblas::Reader rdr = rf.CreateWithStream(in);
		const liblas::Header& hdr = rdr.GetHeader();

		liblas::Header whd(hdr);
		liblas::WriterFactory wf;
		liblas::Writer wtr = liblas::Writer(out, hdr);

		double minx = MAX, miny = MAX, minz = MAX, maxx = MIN, maxy = MIN, maxz = MIN;
		while(rdr.ReadNextPoint()) {
			const liblas::Point& pt = rdr.GetPoint();
			double x = pt.GetX();
			double y = pt.GetY();
			double z = pt.GetZ() + (reference.smoothAt(x, y) - smoothAt(x, y));
			liblas::Point npt(pt);
			npt.SetCoordinates(x, y, z);
			wtr.WritePoint(npt);
			if(x < minx) minx = x;
			if(y < miny) miny = y;
			if(z < minz) minz = z;
			if(x > maxx) maxx = x;
			if(y > maxy) maxy = y;
			if(z > maxz) maxz = z;
		}

		whd.SetMin(minx, miny, minz);
		whd.SetMax(maxx, maxy, maxz);
		wtr.SetHeader(whd);
		wtr.WriteHeader();
		out.close();
		in.close();
	}
};

std::string makeFile(const std::string& outdir, const std::string& tpl, double res, int i, const std::string& ext) {
	std::stringstream ss;
	ss << outdir << "/" << tpl << "_" << res << "_" << i << ext;
	return ss.str();
}

int main(int argc, char** argv) {

	double smooth1 =1000;
	double smooth2 = 200;

	GDALAllRegister();

	double res = atof(argv[argc - 3]);
	int steps = atoi(argv[argc - 2]);
	std::string outdir = argv[argc - 1];
	std::vector<Item> items;

	double res0 = res * std::pow(2, steps);

	std::vector<std::string> las;
	for(int i = 1; i < argc - 3; ++i)
		las.push_back(argv[i]);

	while(res0 >= res) {

		for(const std::string& l : las)
			items.emplace_back(l, res);

		Item avg("", res);

		// Load each point cloud, add it to the overall average.
		// Then average the item and save it.
		for(size_t i = 0; i < items.size(); ++i) {
			Item& item = items[i];
			item.load();
			avg.extend(item);
		}
		avg.calcBounds();
		avg.fill(NO_DATA);
		for(size_t i = 0; i < items.size(); ++i) {
			Item& item = items[i];
			avg.add(item);
			item.average();
			item.save(makeFile(outdir, "item_avg", res, i, ".tif"));
		}

		// Make the average of all point clouds, then produce a smoothed raster surface.
		avg.projection = items[0].projection;
		avg.average();
		avg.save(makeFile(outdir, "avg", res, 0, ".tif"));
		avg.splineSmooth(smooth1);
		avg.smoothGrid();
		avg.save(makeFile(outdir, "smooth", res, 0, ".tif"));

		// Adjust each point cloud by adding the difference between the smoothed overall
		// average, and the smoothed individual average.
		std::string s;
		std::vector<std::string> slas;
		for(size_t i = 0; i < items.size(); ++i) {
			Item& item = items[i];
			item.splineSmooth(smooth2);
			item.smoothGrid();
			s = makeFile(outdir, "adj", res, i, ".las");
			slas.push_back(s);
			item.adjust(s, avg);
		}

		res0 /= 2;
		las.swap(slas);
	}
}
