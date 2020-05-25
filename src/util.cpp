/*
 * contrem_util.cpp
 *
 *  Created on: Jun 4, 2019
 *      Author: rob
 */

#include "geo.hpp"

#ifdef _WIN32
#include <Windows.h>
#include <io.h>
constexpr char pathsep = '\\';
#else
#include <sys/time.h>
constexpr char pathsep = '/';
#endif

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <sstream>
#include <regex>
#include <fstream>
#include <random>

#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include "fitpack.hpp"
#include "util.hpp"


using namespace geo::util;

#if __cplusplus == 201704L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else 
#include <filesystem>
namespace fs = std::filesystem;
#endif

namespace {

	const std::string defaultChars = "abcdefghijklmnaoqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890";

	//https://codereview.stackexchange.com/a/118957
	std::string randomString(size_t len = 15) {
		std::mt19937_64 gen { std::random_device()() };
	    std::uniform_int_distribution<size_t> dist { 0, defaultChars.length()-1 };
	    std::string ret;
	    std::generate_n(std::back_inserter(ret), len, [&] {
	    	return defaultChars[dist(gen)];
	    });
	    return ret;
	}

} // anon


FileType geo::util::getFileType(const std::string& filename) {
	std::string ext;
	{
		size_t p = filename.find('.');
		if(p < std::string::npos) {
			std::string ext0 = filename.substr(p);
			std::transform(ext0.begin(), ext0.end(), std::back_inserter(ext), ::tolower);
		}
	}
	if(ext == ".csv" || ext == ".txt") {
		return FileType::CSV;
	} else if(ext == ".roi") {
		return FileType::ROI;
	} else {
		GDALAllRegister();
		GDALDataset* ds = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_READONLY, 0, 0, 0));
		if(ds) {
			std::string drv(ds->GetDriverName());
			FileType type = FileType::Unknown;
			if(drv == "GTiff") {
				type = FileType::GTiff;
			} else if(drv == "ENVI") {
				type = FileType::ENVI;
			} else if(drv == "ESRI Shapefile") {
				type = FileType::SHP;
			} else if(drv == "SQLite") {
				type = FileType::SQLITE;
			}
			GDALClose(ds);
			return type;
		}
	}
	return FileType::Unknown;
}

std::string geo::util::fileTypeAsString(FileType type) {
	switch(type) {
	case FileType::GTiff: return "GTiff";
	case FileType::ENVI: return "ENVI";
	case FileType::ROI: return "ENVI ROI";
	case FileType::SHP: return "Shapefile";
	case FileType::CSV: return "CSV";
	default: return "";
	}
}

FileType geo::util::fileTypeFromString(const std::string& type) {
	if(type == "GTiff") {
		return FileType::GTiff;
	} else if(type == "ENVI") {
		return FileType::ENVI;
	} else if(type == "ENVI ROI") {
		return FileType::ROI;
	} else if(type == "Shapefile" || type == "ESRI Shapefile") {
		return FileType::SHP;
	} else if(type == "CSV") {
		return FileType::CSV;
	} else {
		return FileType::Unknown;
	}
}

std::string geo::util::normMethodAsString(NormMethod method) {
	switch(method) {
	case NormMethod::ConvexHull:
		return "Convex Hull";
	case NormMethod::ConvexHullLongestSeg:
		return "Convex Hull, Longest Segment";
	case NormMethod::Line:
		return "Line";
	case NormMethod::Unknown:
	default:
		return "Unknown";
	}
}

NormMethod geo::util::normMethodFromString(const std::string& method) {
	if(method == normMethodAsString(NormMethod::ConvexHull)) {
		return NormMethod::ConvexHull;
	} else if(method == normMethodAsString(NormMethod::ConvexHullLongestSeg)) {
		return NormMethod::ConvexHullLongestSeg;
	} else if(method == normMethodAsString(NormMethod::Line)) {
		return NormMethod::Line;
	} else {
		return NormMethod::Unknown;
	}
}

bool geo::util::isnonzero(const double& v) {
	return v != 0;
}

bool geo::util::isdir(const std::string& path) {
	return fs::is_directory(path);
}

bool geo::util::isfile(const std::string& path) {
	return fs::is_regular_file(path);
}

bool geo::util::rem(const std::string& dir) {
	try {
		fs::path p(dir);
		fs::remove_all(p);
	}
	catch (const std::exception& ex) {
		g_warn(ex.what());
		return false;
	}
	return true;
}

std::string geo::util::parent(const std::string& path) {
	std::string _p = path;

	while(_p.size() > 0 && _p.back() == pathsep)
		_p = _p.substr(0, _p.size() - 1);

	if(_p.empty())
		return _p;

	size_t a = _p.find_last_of(pathsep);
	if(a == std::string::npos)
		return "";

	return _p.substr(0, a);

}

bool geo::util::rename(const std::string& from, const std::string& to) {
	if (isdir(to))
		g_runerr(to << " is a directory.")
	return fs::copy_file(from, to);
}

std::string geo::util::join(const std::string& a, const std::string& b) {
	if(b.empty()) {
		return a.empty() ? "" : a;
	} else if(a.empty()) {
		return b.empty() ? "" : b;
	}
	std::string _a, _b;
	for(size_t i = a.size() - 1; i < std::string::npos; --i) {
		if(a[i] != pathsep) {
			_a = a.substr(0, i + 1);
			break;
		}
	}
	for(size_t i = 0; i < b.size(); ++i) {
		if(b[i] != pathsep) {
			_b = b.substr(i);
			break;
		}
	}
	return _a + pathsep + _b;
}

std::string geo::util::basename(const std::string& path) {
	std::string _p = path;

	while(_p.size() > 0 && _p.back() == pathsep)
		_p = _p.substr(0, _p.size() - 1);

	if(_p.empty())
		return _p;

	size_t a = _p.find_last_of(pathsep);
	size_t b = _p.find_last_of('.');
	return _p.substr(a + 1, b - a - 1);
}

std::string geo::util::extension(const std::string& path) {
	size_t pos = path.find_last_of('.');
	if(pos < std::string::npos)
		return path.substr(pos, std::string::npos);
	return path;
}

std::string geo::util::gettmpdir() {
	fs::path p = fs::temp_directory_path();
	return p.string();
}

int geo::util::pid() {
#ifdef _WIN32
	return GetCurrentProcessId();
#else
	return getpid();
#endif
}

std::string geo::util::tmpdir(const std::string& prefix, const std::string& dir) {
	// Assemble the target directory, check and attempt to create if needed.
	std::string tdir = join(gettmpdir(), dir);
	if(!isdir(tdir)) {
		if(!makedir(tdir))
			g_runerr("Failed to make target dir: " << tdir);
	}
	int tries = 16;
	while(--tries) {
		std::stringstream ss;
		ss << join(tdir, prefix) << '_' << pid() << '_' << randomString();
		std::string path = ss.str();
		if(makedir(path))
			return path;
	}
	g_runerr("Failed to make temporary directory.");
}

std::string geo::util::tmpfile(const std::string& prefix, const std::string& dir) {
	// Assemble the target directory, check and attempt to create if needed.
	std::string tdir = join(gettmpdir(), dir);
	if(!isdir(tdir)) {
		if(!makedir(tdir))
			g_runerr("Failed to make target directory: " << tdir);
	}
	int tries = 16;
	while(--tries) {
		// Assemble the file path.
		std::stringstream ss;
		ss << join(tdir, prefix) << '_' << pid() << '_' << randomString();
		std::string path = ss.str();
		if(!isfile(path))
			return path;
	}
	g_runerr("Failed to create non-extant filename.");
}

bool geo::util::makedir(const std::string& path) {
	bool res = fs::create_directories(path);
	if (res)
		fs::permissions(path, fs::perms::owner_all | fs::perms::group_all);
	return res;
}

std::string geo::util::sanitize(const std::string& str) {
	std::regex repl("([^0-9A-Za-z]+)");
	std::stringstream ss;
	std::regex_replace(std::ostreambuf_iterator<char>(ss), str.begin(), str.end(), repl, "_");
	return ss.str();
}

std::string geo::util::lowercase(const std::string& str) {
	std::string out;
	std::transform(str.begin(), str.end(), std::back_inserter(out), ::tolower);
	return out;
}

int geo::util::gdalTypeSize(GDALDataType type) {
	switch(type) {
	case GDT_Float32:
	case GDT_Int32:
	case GDT_UInt32: 	return 4;
	case GDT_Int16:
	case GDT_UInt16: 	return 2;
	case GDT_Float64: 	return 8;
	case GDT_Byte: 		return 1;
	default:
		throw std::runtime_error("Unknown GDAL data type: " + std::to_string((int) type));
	}
}

std::string geo::util::projectionFromSRID(int srid) {
	std::string out;
	OGRSpatialReference osr;
	if((OGRERR_NONE == osr.importFromEPSG(srid))) {
		char* wkt;
		osr.exportToWkt(&wkt);
		out = wkt;
		CPLFree(wkt);
	}
	return out;
}

TmpFile::TmpFile(size_t size) :
	fd(0), size(0) {
	filename = geo::util::tmpfile("geo_util");
	fd = ::open(filename.c_str(), O_CREAT|O_RDWR, 0777);
	if(fd <= 0)
		g_runerr("Failed to open temp file: " << strerror(errno));
	resize(size);
}

void TmpFile::resize(size_t newSize) {
	if(newSize > size) {
		if(fd <= 0)
			throw std::runtime_error(std::string("File is not open: ") + strerror(errno));

		if((size_t) lseek(fd, newSize - 1, SEEK_SET) != (newSize - 1))
			throw std::runtime_error(std::string("Failed to create temporary file for mapping.") + strerror(errno));

		if(write(fd, "", 1) < 1)
			throw std::runtime_error(std::string("Failed to create temporary file for mapping.") + strerror(errno));

		size = newSize;
	}
}

void TmpFile::close() {
	::close(fd);
}

TmpFile::~TmpFile() {
	::close(fd);
	unlink(filename.c_str());
}

constexpr uint32_t MAX_UINT32 = 32768;

uint64_t geo::util::morton(uint32_t x, uint32_t y) {
	if(x >= MAX_UINT32)
		g_runerr("x coordinate is too large for bit shuffling: " << x)
	if(y >= MAX_UINT32)
		g_runerr("y coordinate is too large for bit shuffling: " << y)
	uint32_t xa = bitsplit2(x);
	uint32_t ya = bitsplit2(y);
	return xa | (ya << 1);
}

double geo::util::random(double min, double max) {
	return min + ((double) rand() / RAND_MAX) * (max - min);
}

uint64_t geo::util::microtime() {
#ifdef _WIN32
	FILETIME ft;
	GetSystemTimePreciseAsFileTime(&ft);
	return ((uint64_t) ft.dwHighDateTime << 32) | ft.dwLowDateTime;
#else
	struct timeval t;
	gettimeofday(&t, NULL);
	return (uint64_t) t.tv_sec * 1000000 + t.tv_usec;
#endif
}


double BivariateSpline::stddev(const std::vector<double>& v) const {
	double mean = 0, var = 0;
	for(const double& vv : v)
		mean += vv;
	mean /= v.size();
	for(const double& vv : v)
		var += std::pow(vv - mean, 2.0);
	return std::sqrt(var / v.size());
}

int BivariateSpline::init(double& smooth, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
		std::vector<double>& weights,
		double x0, double y0, double x1, double y1) {

	if(smooth <= 0) {
		smooth = x.size();
		std::cout << "Setting smooth value to " << smooth << "\n";
	}

	if(weights.empty()) {
		std::cout << "Using std. dev. for weights.\n";
		double s = 1.0 / stddev(z);
		weights.resize(x.size());
		for(size_t i = 0; i < x.size(); ++i)
			weights[i] = s;
	}

	int iopt = 0;
	int kx = 3;
	int ky = 3;
	double eps = std::pow(10.0, -15.0);

	int m = x.size();
	if(m < (kx + 1) * (ky + 1))
		throw std::runtime_error("x array size must be greater than or equal to (kx + 1) * (ky + 1).");

	int nxest = (int) std::ceil(kx + 1.0 + std::sqrt(m / 2.0));
	if(nxest < 2 * kx + 2)
		throw std::runtime_error("nxest too small.");
	int nyest = (int) std::ceil(ky + 1.0 + std::sqrt(m / 2.0));
	if(nyest < 2 * ky + 2)
		throw std::runtime_error("nyest too small.");
	int nmax = geo::max(geo::max(m, nxest), nyest);

	m_c.resize((nxest - kx - 1) * (nyest - ky - 1));
	m_tx.resize(m);
	m_ty.resize(m);

	double fp;

	int u = nxest - kx - 1;
	int v = nyest - ky - 1;
	int km = geo::max(kx, ky) + 1;
	int ne = geo::max(nxest, nyest);
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

	int ier, iter = 0;
	std::cout << "X: " << x.size() << "; Y: " << y.size() << "\n";
	do {
		std::cout << "Smoothing with " << smooth << "\n";
		surfit_(&iopt, &m, x.data(), y.data(), z.data(), weights.data(),
				&x0, &x1, &y0, &y1, &kx, &ky,
				&smooth, &nxest, &nyest, &nmax, &eps,
				&m_nx, m_tx.data(), &m_ny, m_ty.data(), m_c.data(), &fp,
				wrk1.data(), (int*) &lwrk1, wrk2.data(), (int*) &lwrk2, iwrk.data(), (int*) &kwrk,	// Dangerous converstion to int.
				&ier);

		if(ier == 1) {
			std::cerr << "Smoothing parameter too small. Increasing: " << smooth << "->" << smooth * 2 << "\n";
			smooth *= 2;
		} else if(ier == 2) {
			std::cerr << "Smoothing parameter probably too small. Increasing: " << smooth << "->" << smooth * 2 << "\n";
			smooth *= 2;
		} else if(ier == 4) {
			std::cerr << "Too many knots. Increased smoothing parameter: " << smooth << "->" << smooth * 2 << "\n";
			smooth *= 2;
		} else if(ier == 5) {
			std::cerr << "Can't add more knots. Increased smoothing parameter: " << smooth << "->" << smooth * 2 << "\n";
			smooth *= 2;
		} else if(ier == -2) {
			std::cerr << "Output is the least squares fit. Smoothing should be no larger than " << fp << " (" << smooth << ")\n";
			smooth *= 0.5;
		} else if(ier < 0) {
			std::cerr << "The coefficients are the minimal norm least-squares solution of a rank deficient system (" << ier << ")\n";
			break;
		} else if(ier > 5) {
			throw std::runtime_error("Smoothing failed");
		} else if(ier == 0) {
			break;
		}
		++iter;
	} while(iter < 100);

	return ier;

}

int BivariateSpline::evaluate(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& z) {
	return evaluate(x.data(), x.size(), y.data(), y.size(), z.data(), z.size());
}

int BivariateSpline::evaluate(const double* x, int nx, const double* y, int ny, double* z, int nz) {

	int idim = 2;
	int mf = nx* ny * idim;

	if(nz < mf)
		g_runerr("Z array is too small at " << nz << "; " << mf << " required.")

	int lwrk1 = (nx + ny) * 4;
	std::vector<double> wrk1(lwrk1);

	int liwrk = nx + ny;
	std::vector<int> iwrk(liwrk);

	int ier;

	surev_(&idim, m_tx.data(), &m_nx, m_ty.data(), &m_ny,
			m_c.data(), x, &nx, y, &ny, z, &mf,
			wrk1.data(), (int*) &lwrk1, iwrk.data(), (int*) &liwrk, &ier);

	return ier;
}


double SmoothingSpline::stddev(const std::vector<double>& v) const {
	double mean = 0, var = 0;
	for(const double& vv : v)
		mean += vv;
	mean /= v.size();
	for(const double& vv : v)
		var += std::pow(vv - mean, 2.0);
	return std::sqrt(var / v.size());
}

int SmoothingSpline::init(double& smooth, const std::vector<double>& x, const std::vector<double>& y,
		std::vector<double>& weights,
		double x0, double x1) {

	if(smooth <= 0) {
		smooth = x.size();
		std::cout << "Setting smooth value to " << smooth << "\n";
	}

	if(weights.empty()) {
		std::cout << "Using std. dev. for weights.\n";
		double s = 1.0 / stddev(y);
		weights.resize(x.size());
		for(size_t i = 0; i < x.size(); ++i)
			weights[i] = s;
	}

	int iopt = 0;
	int k = 3;
	//double eps = std::pow(10.0, -15.0);

	int m = x.size();
	if(m < k + 1)
		throw std::runtime_error("x array size must be greater than or equal to (kx + 1).");

	int nest = m + k + 1;
	//int nmax = geo::max(m, nest);

	m_c.resize(nest);
	m_tx.resize(nest);

	double fp;

	int lwrk = (m * (k + 1) + nest * (7 + 3 * k));
	std::vector<double> wrk(lwrk);

	int kwrk = nest;
	std::vector<int> iwrk(kwrk);

	int n;
	int ier, iter = 0;
	do {
		std::cout << "Smoothing with " << smooth << "\n";
		curfit_(&iopt, &m, x.data(), y.data(), weights.data(),
				&x0, &x1, &k, &smooth, &nest, &n, m_tx.data(), m_c.data(),
				&fp, wrk.data(), &lwrk, iwrk.data(), &ier);

		m_tx.resize(n);
		m_c.resize(n - k - 1);

		if(ier == 1) {
			std::cerr << "Smoothing parameter too small. Increasing: " << smooth << "->" << smooth * 2 << "\n";
			smooth *= 2;
		} else if(ier == 2) {
			std::cerr << "Smoothing parameter probably too small. Increasing: " << smooth << "->" << smooth * 2 << "\n";
			smooth *= 2;
		} else if(ier == 4) {
			std::cerr << "Too many knots. Increased smoothing parameter: " << smooth << "->" << smooth * 2 << "\n";
			smooth *= 2;
		} else if(ier == 5) {
			std::cerr << "Can't add more knots. Increased smoothing parameter: " << smooth << "->" << smooth * 2 << "\n";
			smooth *= 2;
		} else if(ier == -2) {
			std::cerr << "Output is the least squares fit. Smoothing should be no larger than " << fp << " (" << smooth << ")\n";
			smooth *= 0.5;
		} else if(ier < 0) {
			std::cerr << "The coefficients are the minimal norm least-squares solution of a rank deficient system (" << ier << ")\n";
			break;
		} else if(ier > 5) {
			throw std::runtime_error("Smoothing failed");
		} else if(ier == 0) {
			break;
		}
		++iter;
	} while(iter < 100);

	return ier;

}

int SmoothingSpline::evaluate(const std::vector<double>& x, std::vector<double>& y) {

	int idim = 1;
	int nx = x.size();
	int n = m_tx.size();
	int nc = m_c.size();
	int k = 3;
	int mx = nx * idim;

	int ier;

	y.resize(x.size());

	curev_(&idim, m_tx.data(), &n, m_c.data(), &nc, &k,
			x.data(), &nx, y.data(), &mx, &ier);


	return ier;
}

const std::vector<double>& SmoothingSpline::knots() const {
	return m_tx;
}

const std::vector<double>& SmoothingSpline::coefficients() const {
	return m_c;
}

using namespace geo::util::csv;

int CSVValue::asInt() const {
	if(t == Int) {
		return i;
	} else {
		return (int) d;
	}
}

double CSVValue::asDouble() const {
	if(t == Double) {
		return d;
	} else {
		return (double) i;
	}
}

const std::string& CSVValue::asString() const {
	return s;
}


bool CSV::isdouble(const std::string& s) {
	if(s == "inf" || s == "-inf" || s == "NaN")
		return true;
	for(size_t i = 0; i < s.size(); ++i) {
		if(!std::isdigit(s[i], std::locale()) && s[i] != '.' && s[i] != '+' && s[i] != '-' && s[i] != 'e')
			return false;
	}
	return true;
}

bool CSV::isint(const std::string& s) {
	for(size_t i = 0; i < s.size(); ++i) {
		if(!std::isdigit(s[i], std::locale()))
			return false;
	}
	return true;
}

CSV::CSV(const std::string& file, bool header) {
	if(!file.empty())
		load(file, header);
}

void CSV::load(const std::string& file, bool header) {
	std::ifstream in(file);
	std::string line;
	std::string cell;
	std::vector<std::string> names;
	std::vector<CSVType> types;
	std::vector<std::vector<std::string>> values;
	int colCount = 0;
	bool doNames = true;
	while(std::getline(in, line)) {
		std::stringstream ss(line);
		if(header) {
			while(std::getline(ss, cell, ','))
				names.push_back(cell);
			colCount = names.size();
			values.resize(colCount);
			header = false;
			doNames = false;
		} else {
			int idx = 0;
			while(std::getline(ss, cell, ',')) {
				if(idx == colCount)
					break;
				if(doNames) {
					names.push_back("col_" + std::to_string(++colCount));
					values.resize(colCount);
				}
				values[idx++].push_back(cell);
			}
			doNames = false;
		}
	}
	types.resize(names.size());
	for(size_t i = 0; i < names.size(); ++i) {
		int t = 2;
		for(size_t j = 0; j < values[i].size(); ++j) {
			if(t == 2 && !isint(values[i][j])) {
				--t;
			} else if(t == 1 && !isdouble(values[i][j])) {
				--t;
				break;
			}
		}
		switch(t) {
		case 2: types[i] = Int; break;
		case 1: types[i] = Double; break;
		default: types[i] = String; break;
		}
	}
	m_values.resize(names.size());
	for(size_t i = 0; i < names.size(); ++i) {
		m_values[i].name = names[i];
		m_values[i].type = types[i];
		m_values[i].values.resize(values[i].size());
		for(size_t j = 0; j < values[i].size(); ++j) {
			m_values[i].values[j].t = types[i];
			switch(types[i]) {
			case Double:
				m_values[i].values[j].d = atof(values[i][j].c_str());
				break;
			case Int:
				m_values[i].values[j].i = atoi(values[i][j].c_str());
				break;
			default:
				m_values[i].values[j].s = values[i][j];
				break;
			}
		}
		values[i].clear();
	}
}

std::vector<std::string> CSV::columnNames() const {
	std::vector<std::string> names;
	for(const CSVColumn& c : m_values)
		names.push_back(c.name);
	return names;
}

std::vector<CSVValue> CSV::row(size_t idx) const {
	if(idx >= m_values[0].values.size())
		throw std::runtime_error("Index is too large.");
	std::vector<CSVValue> row;
	for(size_t i = 0; i < m_values.size(); ++i)
		row.push_back(m_values[i].values[idx]);
	return row;
}

std::vector<CSVValue> CSV::column(const std::string& name) const {
	for(size_t i = 0; i < m_values.size(); ++i) {
		if(m_values[i].name == name)
			return m_values[i].values;
	}
	throw std::runtime_error("No column named " + name);
}

CSVType CSV::columnType(const std::string& name) const {
	for(size_t i = 0; i < m_values.size(); ++i) {
		if(m_values[i].name == name)
			return m_values[i].type;
	}
	throw std::runtime_error("No column named " + name);
}

std::vector<CSVValue> CSV::column(size_t i) const {
	if(i < m_values.size())
		return m_values[i].values;
	throw std::runtime_error("No column with index " + i);
}

CSVType CSV::columnType(size_t i) const {
	if(i < m_values.size())
		return m_values[i].type;
	throw std::runtime_error("No column with index " + i);
}

void Stopwatch::start() {
	m_begin = std::chrono::steady_clock::now();
}

void Stopwatch::reset() {
	start();
}

std::string Stopwatch::time() {
	int t = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now()- m_begin).count();
	std::stringstream ss;
	ss << std::setfill('0');
	ss << std::setw(2) << (t / 3600) << ':';
	ss << std::setw(2) << (t / 60) << ':';
	ss << std::setw(2) << (t % 60);
	return ss.str();
}


void geo::util::saveGrid(const std::string& file, const std::vector<double> grid,
		int cols, int rows, double minx, double miny, double xres, double yres,
		const std::string& proj) {
	GDALAllRegister();
	GDALDriverManager* dm = GetGDALDriverManager();
	GDALDriver* drv = dm->GetDriverByName("GTiff");
	GDALDataset* ds = drv->Create(file.c_str(), cols, rows, 1, GDT_Float32, 0);
	double trans[] = {minx, xres, 0, miny, 0, yres};
	ds->SetGeoTransform(trans);
	ds->SetProjection(proj.c_str());
	if(CE_None != ds->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, cols, rows, (void*) grid.data(), cols, rows, GDT_Float64, 0, 0, 0))
		std::cerr << "Failed to write to raster.\n";
	ds->GetRasterBand(1)->SetNoDataValue(-9999.0);
	GDALClose(ds);
}
