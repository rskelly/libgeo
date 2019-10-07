/*
 * contrem_util.cpp
 *
 *  Created on: Jun 4, 2019
 *      Author: rob
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <ftw.h>

#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <sstream>
#include <regex>

#include <gdal_priv.h>

#include "util.hpp"

using namespace geo::util;

namespace {

	// Used in rem.
	// https://stackoverflow.com/questions/2256945/removing-a-non-empty-directory-programmatically-in-c-or-c
	int rmFiles(const char *pathname, const struct stat*, int, struct FTW*) {
		if(remove(pathname) < 0) {
			perror("ERROR: remove");
			return -1;
		}
		return 0;
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
	struct stat st;
	if(!stat(path.c_str(), &st))
		return S_ISDIR(st.st_mode);
	return false;
}

bool geo::util::isfile(const std::string& path) {
	struct stat st;
	if(!stat(path.c_str(), &st))
		return S_ISREG(st.st_mode);
	return false;
}

bool geo::util::rem(const std::string& dir) {
	if(isfile(dir)) {
		return !::unlink(dir.c_str());
	} else if (isdir(dir) && nftw(dir.c_str(), rmFiles,10, FTW_DEPTH|FTW_MOUNT|FTW_PHYS) < 0) {
		perror("ERROR: ntfw");
		return false;
	}
	return true;
}

#ifdef _WIN32
	const char pathsep = '\\';
#else
	const char pathsep = '/';
#endif

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

std::string geo::util::join(const std::string& a, const std::string& b) {
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

std::string geo::util::tmpdir(const std::string& tpl) {
	char tname[PATH_MAX];
	std::strncpy(tname, tpl.c_str(), tpl.size());
	tname[tpl.size()] = '\0';
	if(!mkdtemp(tname))
		g_runerr("Failed to make temporary dir " << tname << ": " << strerror(errno));
	return std::string(tname);
}

bool geo::util::makedir(const std::string& filename) {
	std::stringstream path(filename);
	std::stringstream inter;
	std::string part, current;
	if(filename[0] == '/')
		inter << '/';
	while(std::getline(path, part, '/')) {
		inter << part;
		current = inter.str();
		if(!isdir(current) && !isfile(current)) {
			if(mkdir(current.c_str(), 0755))
				return false;
		}
		if(!part.empty())
			inter << '/';
	}
	return true;
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


TmpFile::TmpFile(size_t size) :
	fd(0), size(0) {
	char tpl[] = {"geo_util_XXXXXX"};
	fd = mkstemp(tpl);
	filename = tpl;
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
	struct timeval t;
	gettimeofday(&t, NULL);
	return (uint64_t) t.tv_sec * 1000000 + t.tv_usec;
}

Bounds::Bounds() :
	m_minx(G_DBL_MAX_POS), m_miny(G_DBL_MAX_POS), m_maxx(G_DBL_MAX_NEG), m_maxy(G_DBL_MAX_NEG), m_minz(G_DBL_MAX_POS), m_maxz(G_DBL_MAX_NEG) {
}

Bounds::Bounds(double minx, double miny, double maxx, double maxy) :
		m_minx(std::min(minx, maxx)), m_miny(std::min(miny, maxy)),
		m_maxx(std::max(minx, maxx)), m_maxy(std::max(miny, maxy)),
		m_minz(G_DBL_MAX_POS), m_maxz(G_DBL_MAX_NEG) {
}

Bounds::Bounds(double minx, double miny, double maxx, double maxy, double minz, double maxz) :
		m_minx(std::min(minx, maxx)), m_miny(std::min(miny, maxy)),
		m_maxx(std::max(minx, maxx)), m_maxy(std::max(miny, maxy)),
		m_minz(std::min(minz, maxz)), m_maxz(std::min(minz, maxz)) {
}

void Bounds::assign(const Bounds& bounds) {
	set(bounds.minx(), bounds.miny() , bounds.maxx(), bounds.maxy(), bounds.minz(), bounds.maxz());
}

void Bounds::set(double minx, double miny, double maxx, double maxy, double minz, double maxz) {
	m_minx = std::min(minx, maxx);
	m_miny = std::min(miny, maxy);
	m_minz = std::min(minz, maxz);
	m_maxx = std::max(minx, maxx);
	m_maxy = std::max(miny, maxy);
	m_maxz = std::max(minz, maxz);
}

bool Bounds::contains(double x, double y) const {
	return x >= m_minx && x < m_maxx && y >= m_miny && y < m_maxy;
}

bool Bounds::contains(double x, double y, double z) const {
	return contains(x, y) && z >= m_minz && z < m_maxz;
}

bool Bounds::contains(const Bounds &b, int dims) const {
	if (dims == 3) {
		return contains(b.minx(), b.miny(), b.minz())
				&& contains(b.maxx(), b.maxy(), b.maxz());
	} else {
		return contains(b.minx(), b.miny()) && contains(b.maxx(), b.maxy());
	}
}

bool Bounds::intersects(const Bounds &b, int dims) const {
	if (dims == 3) {
		return !(b.maxx() < minx() || b.maxy() < miny() || b.minx() > maxx()
				|| b.miny() > maxy() || b.minz() > maxz() || b.maxz() < minz());
	} else {
		return !(b.maxx() < minx() || b.maxy() < miny() || b.minx() > maxx()
				|| b.miny() > maxy());
	}
}

Bounds Bounds::intersection(const Bounds &other) const {
	return Bounds(g_max(minx(), other.minx()), g_max(miny(), other.miny()),
			g_min(maxx(), other.maxx()), g_min(maxy(), other.maxy()));
}

void Bounds::cube() {
	double max = g_max(width(), height());
	if(m_maxz != G_DBL_MAX_POS) {
		max = g_max(max, depth()) / 2.0;
		set(midx() - max, midy() - max, midx() + max, midy() + max, midz() - max, midz() + max);
	} else {
		max /= 2.0;
		set(midx() - max, midy() - max, midx() + max, midy() + max);
	}
}

double Bounds::midx() const {
	return m_minx + (m_maxx - m_minx) / 2.0;
}

double Bounds::midy() const {
	return m_miny + (m_maxy - m_miny) / 2.0;
}

double Bounds::midz() const {
	return m_minz + (m_maxz - m_minz) / 2.0;
}


double Bounds::minx() const {
	return m_minx;
}

void Bounds::minx(double minx) {
	m_minx = minx;
}

double Bounds::miny() const {
	return m_miny;
}

void Bounds::miny(double miny) {
	m_miny = miny;
}

double Bounds::minz() const {
	return m_minz;
}

void Bounds::minz(double minz) {
	m_minz = minz;
}

double Bounds::maxx() const {
	return m_maxx;
}

void Bounds::maxx(double maxx) {
	m_maxx = maxx;
}

double Bounds::maxy() const {
	return m_maxy;
}

void Bounds::maxy(double maxy) {
	m_maxy = maxy;
}

double Bounds::maxz() const {
	return m_maxz;
}

void Bounds::maxz(double maxz) {
	m_maxz = maxz;
}

double Bounds::width() const {
	return maxx() - minx();
}

double Bounds::height() const {
	return maxy() - miny();
}

double Bounds::depth() const {
	return maxz() - minz();
}

double Bounds::volume() const {
	return width() * height() * depth();
}

int Bounds::maxCol(double resolution) const {
	return (int) g_abs(width() / resolution);
}

int Bounds::maxRow(double resolution) const {
	return (int) g_abs(height() / resolution);
}

int Bounds::toCol(double x, double resolution) const {
    if(width() == 0.0)
        return 0;
    if(resolution > 0) {
        return (int) ((x - m_minx) / width() * (width() / resolution));
    } else {
        return (int) ((x - m_maxx) / width() * (width() / resolution));
    }
}

int Bounds::toRow(double y, double resolution) const {
    if(height() == 0.0)
        return 0;
    if(resolution > 0) {
        return (int) ((y - m_miny) / height() * (height() / resolution));
    } else {
        return (int) ((y - m_maxy) / height() * (height() / resolution));
    }
}

double Bounds::toX(int col, double resolution) const {
	if(resolution > 0) {
		return m_minx + resolution * col;
	} else {
		return m_maxx + resolution * col;
	}
}

double Bounds::toY(int row, double resolution) const {
	if(resolution > 0) {
		return m_miny + resolution * row;
	} else {
		return m_maxy + resolution * row;
	}
}

void Bounds::extend(const Bounds &b) {
	m_minx = g_min(b.minx(), m_minx);
	m_maxx = g_max(b.maxx(), m_maxx);
	m_miny = g_min(b.miny(), m_miny);
	m_maxy = g_max(b.maxy(), m_maxy);
	m_minz = g_min(b.minz(), m_minz);
	m_maxz = g_max(b.maxz(), m_maxz);
}

void Bounds::extendX(double x) {
	m_minx = g_min(x, m_minx);
	m_maxx = g_max(x, m_maxx);
}

void Bounds::extendY(double y) {
	m_miny = g_min(y, m_miny);
	m_maxy = g_max(y, m_maxy);
}

void Bounds::extendZ(double z) {
	m_minz = g_min(z, m_minz);
	m_maxz = g_max(z, m_maxz);
}

void Bounds::extend(double x, double y) {
	extendX(x);
	extendY(y);
}

void Bounds::extend(double x, double y, double z) {
	extend(x, y);
	extendZ(z);
}

double Bounds::operator[](size_t pos) const {
	switch (pos) {
	case 0:
		return m_minx;
	case 1:
		return m_miny;
	case 2:
		return m_maxx;
	case 3:
		return m_maxy;
	case 4:
		return m_minz;
	case 5:
		return m_maxz;
	default:
		g_argerr("Illegal position: " << pos);
	}
}

void Bounds::snap(double resolution) {
	minx(std::floor(minx() / resolution) * resolution);
	miny(std::floor(miny() / resolution) * resolution);
	maxx(std::floor(maxx() / resolution) * resolution + resolution);
	maxy(std::floor(maxy() / resolution) * resolution + resolution);
}

void Bounds::collapse(int dims) {
	minx (G_DBL_MAX_POS);
	miny(G_DBL_MAX_POS);
	maxx (G_DBL_MAX_NEG);
	maxy(G_DBL_MAX_NEG);
	if (dims == 3) {
		minz(G_DBL_MAX_POS);
		maxz(G_DBL_MAX_NEG);
	}
}

std::string Bounds::print() const {
	std::stringstream s;
	print(s);
	return s.str();
}

void Bounds::print(std::ostream &str) const {
	str << "[Bounds: " << minx() << ", " << miny() << ", " << minz() << "; "
			<< maxx() << ", " << maxy() << ", " << maxz() << "]";
}

void Bounds::fromString(const std::string &str) {
	std::vector<std::string> parts;
	split(std::back_inserter(parts), str);
	if (parts.size() < 4)
		g_runerr("Bounds string must be 4 or 6 comma-separated doubles.");
	m_minx = atof(parts[0].c_str());
	m_miny = atof(parts[1].c_str());
	m_maxx = atof(parts[2].c_str());
	m_maxy = atof(parts[3].c_str());
	if (parts.size() >= 6) {
		m_maxz = atof(parts[4].c_str());
		m_maxz = atof(parts[5].c_str());
	}
}

std::string Bounds::toString() const {
	std::stringstream ss;
	ss << minx() << "," << miny() << "," << maxx() << "," << maxy() << ","
			<< minz() << "," << maxz();
	return ss.str();
}

void Bounds::align(double x, double y, double xres, double yres) {
	xres = g_abs(xres);
	yres = g_abs(yres);
	while (x < m_minx)
		x += xres;
	while (x > m_minx)
		x -= xres;
	m_minx = x;
	while (x < m_maxx)
		x += xres;
	m_maxx = x;
	while (y < m_miny)
		y += yres;
	while (y > m_miny)
		y -= yres;
	m_miny = y;
	while (y < m_maxy)
		y += yres;
	m_maxy = y;
}

