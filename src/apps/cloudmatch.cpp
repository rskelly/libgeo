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

//#include <liblas/liblas.hpp>

#include <gdal_priv.h>
#include <ogr_spatialref.h>

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
	std::vector<std::string> files;
	double xmin, xmax, ymin, ymax;
	double res;
	int cols, rows;
	int validCount;
	std::vector<double> grid;
	std::vector<int> counts;
	std::vector<double> std;
	std::vector<double> diffs;
	std::string projection;

	//std::vector<double> m_wrk1;
	//std::vector<double> m_wrk2;
	//std::vector<int> m_iwrk;

	double* m_wrk1;
	size_t m_lwrk1;
	double* m_wrk2;
	size_t m_lwrk2;
	int* m_iwrk;
	size_t m_liwrk;

	std::vector<double> m_x;
	std::vector<double> m_y;
	std::vector<double> m_z;
	std::vector<double> m_w;

	std::vector<double> m_c;
	std::vector<double> m_tx;
	std::vector<double> m_ty;

	int m_nx;
	int m_ny;

	Item(const std::vector<std::string>& files, double xmin, double ymin, double xmax, double ymax, double res) :
		files(files),
		xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax),
		res(res),
		cols(0), rows(0),
		validCount(0),
		m_wrk1(nullptr), m_lwrk1(0),
		m_wrk2(nullptr), m_lwrk2(0),
		m_iwrk(nullptr), m_liwrk(0),
		m_nx(0), m_ny(0) {
	}

	Item(const std::vector<std::string>& files) : Item(files, MAX, MAX, MIN, MIN, 0) {}

	Item(const std::string& file) : Item(std::vector<std::string>({file}), MAX, MAX, MIN, MIN, 0) {}

	Item() : Item(std::vector<std::string>({})) {}

	Item(const Item& item) : Item(item.files, item.xmin, item.ymin, item.xmax, item.ymax, item.res) {
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

		reset();

		{
			LasFile las(files[0]);
			projection = las.projection();
		}

		int i = 0;
		for(const std::string& file : files) {
			std::cout << "Loading " << ++i << " of " << files.size() << "\n";
			LasFile las(file);
			Point pt;
			while(las.next(pt)) {
				if(pt.x < xmin) xmin = pt.x;
				if(pt.y < ymin) ymin = pt.y;
				if(pt.x > xmax) xmax = pt.x;
				if(pt.y > ymax) ymax = pt.y;
			}
		}

		calcBounds();

		fill(0);
		fillCounts(0);

		for(const std::string& file : files) {
			LasFile las(file);
			Point pt;
			while(las.next(pt)) {
				if(pt.cls != 2)
					continue;
				double x = pt.x;
				double y = pt.y;
				double z = pt.z;
				set(x, y, get(x, y) + z);
				setCount(x, y, getCount(x, y) + 1);
			}
		}

		for(size_t i = 0; i < grid.size(); ++i) {
			if(counts[i] == 0) {
				grid[i] = NODATA;
			} else {
				grid[i] /= counts[i];
			}
		}
	}

	void fill(double v) {
		std::fill(grid.begin(), grid.end(), v);
	}

	void fillCounts(int v) {
		std::fill(counts.begin(), counts.end(), v);
	}

	void fillStd(double v) {
		std::fill(std.begin(), std.end(), v);
	}

	bool intersects(const Item& other) {
		return !(other.xmax < xmin || other.xmin > xmax || other.ymax < ymin || other.ymin > ymax);
	}

	void diff(Item& other) {
		diffs.resize(cols * rows);
		std::fill(diffs.begin(), diffs.end(), NODATA);
		double z, oz;
		for(double y = ymin; y < ymax; y += res) {
			for(double x = xmin; x < xmax; x += res) {
				if((oz = other.get(x, y)) == NODATA || (z = get(x, y)) == NODATA)
					continue;
				int c = toCol(x);
				int r = toRow(y);
				diffs[r * cols + c] = oz - z;
			}
		}
	}

	void stats() {
		std.resize(cols * rows);

		fillStd(0);
		validCount = 0;

		for(const std::string& file : files) {
			LasFile las(file);
			Point pt;

			while(las.next(pt)) {
				if(pt.cls != 2)
					continue;
				double x = pt.x;
				double y = pt.y;
				double z = pt.z;
				setStd(x, y, getStd(x, y) + std::pow(z - get(x, y), 2.0));
			}
		}

		for(size_t i = 0; i < grid.size(); ++i) {
			if(counts[i] == 0) {
				std[i] = NODATA;
			} else {
				std[i] = std::sqrt(std[i]);
				++validCount;
			}
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

	void setStd(double x, double y, double v) {
		setStd(toCol(x), toRow(y), v);
	}

	void setStd(int col, int row, double v) {
		if(col >= 0 && col < cols && row >= 0 && row < rows)
			std[row * cols + col] = v;
	}

	double get(double x, double y) {
		return get(toCol(x), toRow(y));
	}

	double get(int col, int row) {
		if(col >= 0 && col < cols && row >= 0 && row < rows) {
			return grid[row * cols + col];
		} else {
			return NODATA;
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

	double getStd(double x, double y) {
		return getStd(toCol(x), toRow(y));
	}

	double getStd(int col, int row) {
		if(col >= 0 && col < cols && row >= 0 && row < rows) {
			return std[row * cols + col];
		} else {
			return 0;
		}
	}

	void splineSmooth(double smooth, std::vector<double>& target) {

		stats();

		if(smooth == 0)
			smooth = validCount * 2;

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
				if((gv = target[r * cols + c]) != NODATA && (gs = getStd(c, r)) > 0) {
					m_x.push_back(toX(c));
					m_y.push_back(toY(r));
					m_z.push_back(gv);
					m_w.push_back(1.0 / gs);
					++m;
				}
			}
		}

		int nxest = (int) std::ceil(kx + 1.0 + std::sqrt(m / 2.0)) * 2;
		int nyest = (int) std::ceil(ky + 1.0 + std::sqrt(m / 2.0)) * 2;
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
		//m_wrk1.resize(lwrk1);
		m_wrk1 = (double*) mmap(0, m_lwrk1 * sizeof(double), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

		if(m_wrk2)
			munmap(m_wrk2, m_lwrk2);
		m_lwrk2 = u * v * (b2 + 1) + b2;
		//m_wrk2.resize(lwrk2);
		m_wrk2 = (double*) mmap(0, m_lwrk2 * sizeof(double), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

		if(m_iwrk)
			munmap(m_iwrk, m_liwrk);
		m_liwrk = m + (nxest - 2 * kx - 1) * (nyest - 2 * ky - 1);
		//m_iwrk.resize(kwrk);
		m_iwrk = (int*) mmap(0, m_liwrk * sizeof(int), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

		int ier, iter = 0;
		double x0 = toX(0) - res / 2, x1 = toX(cols - 1) + res / 2;
		double y0 = toY(0) - res / 2, y1 = toY(rows - 1) + res / 2;
		do {
			surfit_(&iopt, &m, m_x.data(), m_y.data(), m_z.data(), m_w.data(),
					&x0, &x1, &y0, &y1, &kx, &ky,
					&smooth, &nxest, &nyest, &nmax, &eps,
					&m_nx, m_tx.data(), &m_ny, m_ty.data(), m_c.data(), &fp,
					m_wrk1, (int*) &m_lwrk1, m_wrk2, (int*) &m_lwrk2, m_iwrk, (int*) &m_liwrk,	// Dangerous converstion to int.
					&ier);
			std::cerr << "ier : " << ier << "\n";
			if(ier == 1) {
				smooth *= 10;
			} else if(ier == 4) {
				smooth *= 10;
			} else if(ier < 0) {
				smooth *= 0.1;
			} else if(ier > 5) {
				throw std::runtime_error("Smoothing failed");
			} else if(ier == 0) {
				break;
				std::cerr << "ier : " << ier << "\n";
			}
			++iter;
		} while(iter < 100);

	}

	void smooth(std::vector<double>& target) {

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
		//m_wrk1.resize(lwrk1);
		m_wrk1 = (double*) mmap(0, m_lwrk1 * sizeof(double), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

		if(m_iwrk)
			munmap(m_iwrk, m_liwrk);
		m_liwrk = cols * rows;
		//m_iwrk.resize(kwrk);
		m_iwrk = (int*) mmap(0, m_liwrk * sizeof(int), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);

		int ier;

		surev_(&idim, m_tx.data(), &m_nx, m_ty.data(), &m_ny,
				m_c.data(), m_x.data(), &cols, m_y.data(), &rows, m_z.data(), &mf,
				m_wrk1, (int*) &m_lwrk1, m_iwrk, (int*) &m_liwrk, &ier);

		fill(NODATA);
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c)
				target[r * cols + c] = m_z[c * rows + r]; // Note: z is transposed from usual.
		}
	}

	void save(const std::string& outfile) {
		GDALAllRegister();
		GDALDriverManager* dm = GetGDALDriverManager();
		GDALDriver* drv = dm->GetDriverByName("GTiff");
		int bands = 0;
		if(!grid.empty()) bands++;
		if(!std.empty()) bands++;
		if(!counts.empty()) bands++;
		if(!diffs.empty()) bands++;
		GDALDataset* ds = drv->Create(outfile.c_str(), cols, rows, bands, GDT_Float32, 0);
		const char* proj = projection.c_str();
		ds->SetProjection(proj);
		int b = 1;
		GDALRasterBand* b1 = ds->GetRasterBand(b++);
		b1->SetNoDataValue(NODATA);
		if(CE_None != b1->RasterIO(GF_Write, 0, 0, cols, rows, counts.data(), cols, rows, GDT_Int32, 0, 0, 0))
			throw std::runtime_error("Failed to write band 1.");
		if(!std.empty()) {
			GDALRasterBand* b2 = ds->GetRasterBand(b++);
			b2->SetNoDataValue(NODATA);
			if(CE_None != b2->RasterIO(GF_Write, 0, 0, cols, rows, std.data(), cols, rows, GDT_Float64, 0, 0, 0))
				throw std::runtime_error("Failed to write band 2.");
		}
		if(!diffs.empty()) {
			GDALRasterBand* b2 = ds->GetRasterBand(b++);
			b2->SetNoDataValue(NODATA);
			if(CE_None != b2->RasterIO(GF_Write, 0, 0, cols, rows, diffs.data(), cols, rows, GDT_Float64, 0, 0, 0))
				throw std::runtime_error("Failed to write band 2.");
		}
		GDALRasterBand* b3 = ds->GetRasterBand(b++);
		b3->SetNoDataValue(NODATA);
		if(CE_None != b3->RasterIO(GF_Write, 0, 0, cols, rows, grid.data(), cols, rows, GDT_Float64, 0, 0, 0))
			throw std::runtime_error("Failed to write band 2.");
		GDALClose(ds);
	}

	void adjust(Item& ref, const std::string& outfile) {
		adjust(ref, std::vector<std::string>({outfile}));
	}

	void adjust(Item& ref, const std::vector<std::string>& outfiles) {
		for(size_t i = 0; i < files.size(); ++i) {
			LasFile las(files[i]);
			Point pt;
			while(las.next(pt)) {
				pt.z += (get(pt.x, pt.y) - ref.get(pt.x, pt.y));
				las.update(pt);
			}
			las.write(outfiles[i]);
		}
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

int main(int argc, char** argv) {

	double res = 50;
	std::string outdir;
	std::vector<std::string> las;

	for(int i = 1; i < argc - 3; ++i) {
		std::string arg = argv[i];
		if(arg == "-f") {
			std::vector<std::string> files = getFiles(std::string(argv[++i]));
			for(const std::string& file : files)
				las.emplace_back(file);
			continue;
		} else if(arg == "-r") {
			res = atof(argv[++i]);
			continue;
		} else if(arg == "-o") {
			outdir = argv[++i];
		}
	}

	Item ref(las);
	ref.res = res;
	ref.load();
	ref.splineSmooth(0, ref.grid);
	ref.smooth(ref.grid);
	ref.save(makeFile(outdir, "ref_smooth", res, 0, ".tif"));

	int i = 0;
	for(const std::string& f : las) {
		Item item(f);
		item.res = res;
		item.load();
		item.splineSmooth(0, item.grid);
		item.smooth(item.grid);
		item.save(makeFile(outdir, "item_smooth", res, ++i, ".tif"));
		item.adjust(ref, makeFile(outdir, "las", res, ++i, ".las"));
	}

}
