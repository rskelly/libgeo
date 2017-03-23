#include "util.hpp"
#include "geo.hpp"

#include <ogr_spatialref.h>

#include <set>
#include <list>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <map>
#include <cmath>
#include <string>
#include <tuple>

using namespace geo::util;

void Callbacks::stepCallback(float status) const {
	g_debug("Step: " << (int) (status * 100.0f) << "%");
}

void Callbacks::overallCallback(float status) const {
	g_debug("Overall: " << (int) (status * 100.0f) << "%");
}

void Callbacks::statusCallback(const std::string &msg) const {
	g_debug("Status: " << msg);
}

Callbacks::~Callbacks() {
}

Status::Status(Callbacks *callbacks, float start, float end) :
	callbacks(callbacks), start(start), end(end) {
}
void Status::update(float s) {
	callbacks->stepCallback(start + (end - start) * s);
}

Point::Point(double x, double y, double z) :
		x(x), y(y), z(z) {
}

Point::Point(double x, double y, double z,
		const std::map<std::string, std::string> &fields) :
		x(x), y(y), z(z) {
	for (auto it : fields)
		this->fields[it.first] = it.second;
}

Point::Point() :
		x(0), y(0), z(0) {
}

Bounds::Bounds() :
		m_minx(G_DBL_MAX_POS), m_miny(G_DBL_MAX_POS), m_minz(G_DBL_MAX_POS), m_maxx(
				G_DBL_MAX_NEG), m_maxy(G_DBL_MAX_NEG), m_maxz(G_DBL_MAX_NEG) {
}

Bounds::Bounds(double minx, double miny, double maxx, double maxy) :
		m_minx(minx), m_miny(miny), m_minz(G_DBL_MAX_NEG), m_maxx(maxx), m_maxy(
				maxy), m_maxz(G_DBL_MAX_POS) {
}

Bounds::Bounds(double minx, double miny, double maxx, double maxy, double minz,
		double maxz) :
		m_minx(minx), m_miny(miny), m_minz(minz), m_maxx(maxx), m_maxy(maxy), m_maxz(
				maxz) {
}

void Bounds::assign(const Bounds& bounds) {
	set(bounds.minx(), bounds.miny() , bounds.maxx(), bounds.maxy(), bounds.minz(), bounds.maxz());
}

void Bounds::set(double _minx, double _miny, double _maxx, double _maxy, double _minz, double _maxz) {
	m_minx = _minx;
	m_miny = _miny;
	m_minz = _minz;
	m_maxx = _maxx;
	m_maxy = _maxy;
	m_maxz = _maxz;
}

bool Bounds::contains(double x, double y) const {
	return x > m_minx && x < m_maxx && y > m_miny && y < m_maxy;
}

bool Bounds::contains(double x, double y, double z) const {
	return contains(x, y) && z > m_minz && z < m_maxz;
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
	Util::splitString(std::back_inserter(parts), str);
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

double Util::computeArea(double x1, double y1, double z1, double x2, double y2,
		double z2, double x3, double y3, double z3) {
	double side0 = std::sqrt(
			std::pow(x1 - x2, 2.0) + std::pow(y1 - y2, 2.0)
					+ std::pow(z1 - z2, 2.0));
	double side1 = std::sqrt(
			std::pow(x2 - x3, 2.0) + std::pow(y2 - y3, 2.0)
					+ std::pow(z2 - z3, 2.0));
	double side2 = std::sqrt(
			std::pow(x3 - x1, 2.0) + std::pow(y3 - y1, 2.0)
					+ std::pow(z3 - z1, 2.0));
	double s = (side0 + side1 + side2) / 2.0;
	return std::sqrt(s * (s - side0) * (s - side1) * (s - side2));
}

void Util::copyfile(std::string &srcfile, std::string &dstfile) {
	std::ifstream src(srcfile.c_str(), std::ios::binary);
	std::ofstream dst(dstfile.c_str(), std::ios::binary);
	dst << src.rdbuf();
}

bool Util::exists(const std::string &name) {
	boost::filesystem::path p(name);
	return boost::filesystem::exists(p);
}

std::string Util::pathJoin(const std::string& a, const std::string& b) {
	boost::filesystem::path pa(a);
	boost::filesystem::path pb(b);
	return (pa / pb).string();
}

bool Util::pathExists(const std::string &name) {
	boost::filesystem::path p(name);
	return boost::filesystem::exists(p.remove_filename());
}

bool Util::rm(const std::string &name) {
	using namespace boost::filesystem;
	path p(name);
	return remove_all(p) > 0;
}

bool Util::mkdir(const std::string &dir) {
	using namespace boost::filesystem;
	path bdir(dir);
	if (!boost::filesystem::exists(bdir))
		return create_directories(bdir);
	return true;
}

std::string Util::parent(const std::string& file) {
	using namespace boost::filesystem;
	path p(file);
	return p.parent_path().string();
}

std::string Util::extension(const std::string &filename) {
	using namespace boost::filesystem;
	path p(filename);
	std::string ext = p.extension().string();
	if(ext.size())
		ext = ext.substr(1);
	//lower(ext);
	return ext;
}

std::string& Util::lower(std::string &str) {
	std::transform(str.begin(), str.end(), str.begin(), ::tolower);
	return str;
}

std::string& Util::upper(std::string &str) {
	std::transform(str.begin(), str.end(), str.begin(), ::toupper);
	return str;
}

std::string Util::lower(const std::string &str) {
	std::string n(str);
	std::transform(n.begin(), n.end(), n.begin(), ::tolower);
	return n;
}

std::string Util::upper(const std::string &str) {
	std::string n(str);
	std::transform(n.begin(), n.end(), n.begin(), ::toupper);
	return n;
}

std::string Util::tmpDir() {
	return boost::filesystem::temp_directory_path().string();
}

uint64_t Util::diskSpace(const std::string& path) {
	using namespace boost::filesystem;
	space_info si = space(path);
	return si.available;
}

std::string Util::tmpFile(const std::string &root) {
	using namespace boost::filesystem;
	path p = unique_path();
	if (!root.empty()) {
		path r(root);
		return (r / p).string();
	}
	p = temp_directory_path() / p;
	return p.string(); // Windows can have wide string paths.
}

using namespace boost::interprocess;

MappedFile::MappedFile(const std::string &filename, uint64_t size, bool remove) :
		m_filename(filename),
		m_size(size),
		m_remove(remove),
		m_mapping(nullptr),
		m_region(nullptr) {

	using namespace boost::interprocess;
	using namespace boost::filesystem;

	{
		std::filebuf fbuf;
		if (size > 0) {
			std::string dirname = Util::parent(filename);
			Util::mkdir(dirname);
			if(Util::diskSpace(dirname) < size)
				g_runerr("There is not enough free space for the requested memory-mapped file.");
			if(!fbuf.open(filename,
					std::ios_base::in | std::ios_base::out
							| std::ios_base::trunc | std::ios_base::binary))
				g_runerr("Failed to open file for memory-mapping.");
			long res = fbuf.pubseekoff(size - 1, std::ios_base::beg);
			if(res < 0 || (uint64_t) res < size - 1)
				g_runerr("Failed to reserve space for memory-mapped file.");
			if(EOF == fbuf.sputc(0))
				g_runerr("Failed to reserve space for memory-mapped file.");
		}
	}

	m_mapping = new file_mapping(filename.c_str(), read_write);
	m_region = new mapped_region(*m_mapping, read_write);
}

void* MappedFile::data() {
	return m_region->get_address();
}

uint64_t MappedFile::size() {
	return m_size;
}

size_t MappedFile::pageSize() {
	return m_region->get_page_size();
}

MappedFile::~MappedFile() {
	m_mapping->remove(m_filename.c_str());
	delete m_region;
	delete m_mapping;
	if (m_remove) // TODO: Check if shared? Also, remove_file_on_destroy
		Util::rm(m_filename);
}

std::unique_ptr<MappedFile> Util::mapFile(const std::string &filename,
		uint64_t size, bool remove) {
	std::unique_ptr<MappedFile> mf(new MappedFile(filename, size, remove));
	return std::move(mf);
}

std::string CRS::epsg2Proj4(int crs) const {
	OGRSpatialReference ref;
	char *wkt;
	ref.importFromEPSG(crs);
	ref.exportToProj4(&wkt);
	return std::string(wkt);
}

std::string CRS::epsg2WKT(int crs) const {
	OGRSpatialReference ref;
	char *wkt;
	ref.importFromEPSG(crs);
	ref.exportToWkt(&wkt);
	return std::string(wkt);
}

