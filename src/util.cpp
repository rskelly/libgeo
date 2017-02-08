#include <set>
#include <list>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <map>
#include <cmath>
#include <string>
#include <tuple>

#include <ogr_spatialref.h>

#include "geotools.hpp"
#include "util.hpp"

using namespace geotools::util;

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
	Util::splitString(str, parts);
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

void Util::parseRanges(std::set<double> &values, const char *str, double step) {
	std::stringstream ss;
	double first = 0, second = 0;
	bool range = false;
	int i = 0;
	char c = str[i++];
	while (true) {
		if (c == '-') {
			range = true;
			first = atof(ss.str().c_str());
			ss.str(std::string());
		} else if (c == ',' || c == '\0') {
			if (!range) {
				values.insert(atof(ss.str().c_str()));
				ss.str(std::string());
			} else {
				second = atof(ss.str().c_str());
				ss.str(std::string());
				step = g_abs(step);
				if (first > second) {
					double tmp = second;
					second = first;
					first = tmp;
				}
				for (double i = first; i <= second; i += step)
					values.insert(i);
				range = false;
			}
			if (c == '\0')
				break;
		} else {
			ss << c;
		}
		c = str[i++];
	}
}

void Util::parseRanges(std::set<int> &values, const char *str) {
	std::stringstream ss;
	int first = 0, second = 0;
	bool range = false;
	int i = 0;
	char c = str[i++];
	while (true) {
		if (c == '-') {
			range = true;
			first = atoi(ss.str().c_str());
			ss.str(std::string());
		} else if (c == ',' || c == '\0') {
			if (!range) {
				values.insert(atoi(ss.str().c_str()));
				ss.str(std::string());
			} else {
				second = atoi(ss.str().c_str());
				ss.str(std::string());
				if (first > second) {
					int tmp = second;
					second = first;
					first = tmp;
				}
				for (int i = first; i <= second; ++i)
					values.insert(i);
				range = false;
			}
			if (c == '\0')
				break;
		} else {
			ss << c;
		}
		c = str[i++];
	}
}

void Util::splitString(const std::string &str, std::list<std::string> &lst) {
	std::stringstream ss(str);
	std::string item;
	while (std::getline(ss, item, ','))
		lst.push_back(item);
}

void Util::splitString(const std::string &str, std::vector<std::string> &lst) {
	std::stringstream ss(str);
	std::string item;
	while (std::getline(ss, item, ','))
		lst.push_back(item);
}

std::string Util::join(const std::vector<std::string> &lst, const std::string &delim) {
	return boost::algorithm::join(lst, delim);
}

/**
 // Split a comma-delimited string into a set of unique integers.
 */
void Util::intSplit(std::set<int> &values, const char *str) {
	std::stringstream ss(str);
	std::string item;
	while (std::getline(ss, item, ','))
		values.insert(atoi(item.c_str()));
}

// TODO: Template

void Util::intSplit(std::set<uint8_t> &values, const char *str) {
	std::stringstream ss(str);
	std::string item;
	while (std::getline(ss, item, ','))
		values.insert((uint8_t) atoi(item.c_str()));
}

/**
 // Split a comma-delimited string into a set of unique integers.
 */
void Util::intSplit(std::list<int> &values, const char *val) {
	std::stringstream ss(val);
	std::string item;
	while (std::getline(ss, item, ','))
		values.push_back(atoi(item.c_str()));
}

/**
 // Split a comma-delimited string into a set of unique integers.
 */
void Util::intSplit(std::vector<int> &values, const char *str) {
	std::stringstream ss(str);
	std::string item;
	while (std::getline(ss, item, ','))
		values.push_back(atoi(item.c_str()));
}

/**
 // Return true if the integer is in the set, or the set is empty.
 */
bool Util::inList(std::set<int> &values, int value) {
	return values.size() == 0 || values.find(value) != values.end();
}

bool Util::inList(std::vector<int> &values, int value) {
	return std::find(values.begin(), values.end(), value) != values.end();
}

void Util::copyfile(std::string &srcfile, std::string &dstfile) {
	std::ifstream src(srcfile.c_str(), std::ios::binary);
	std::ofstream dst(dstfile.c_str(), std::ios::binary);
	dst << src.rdbuf();
}

// Load the samples from a csv file. The file must have x, y and z headers.

void Util::loadXYZSamples(std::string &datafile,
		std::vector<std::tuple<double, double, double> > &samples) {
	g_runerr("Not implemented");
	/*
	 io::CSVReader<3> in(datafile.c_str());
	 in.read_header(io::ignore_extra_column, "x", "y", "z");
	 double x, y, z;
	 while (in.read_row(x, y, z))
	 samples.push_back(std::make_tuple(x, y, z));
	 */
}

void Util::loadIDXYZSamples(std::string &datafile,
		std::vector<std::tuple<std::string, double, double, double> > &samples) {
	g_runerr("Not implemented");
	/*
	 io::CSVReader<4> in(datafile.c_str());
	 in.read_header(io::ignore_extra_column, "id", "x", "y", "z");
	 std::string id;
	 double x, y, z;
	 while (in.read_row(id, x, y, z))
	 samples.push_back(std::make_tuple(id, x, y, z));
	 */
}

void Util::status(int step, int of, const std::string &message, bool end) {
#pragma omp critical(__status)
	{
		if (step < 0)
			step = 0;
		if (of <= 0)
			of = 1;
		if (step > of)
			of = step;
		float status = (float) (step * 100) / of;
		std::stringstream out;
		out << "Status: " << std::fixed << std::setprecision(2) << status
				<< "% " << message << std::right << std::setw(100)
				<< std::setfill(' ');
		if (end)
			out << std::endl;
		else
			out << '\r';
		std::cerr << out.str();
		std::cerr.flush();
	}
}

bool Util::exists(const std::string &name) {
	boost::filesystem::path p(name);
	return boost::filesystem::exists(p);
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
		return create_directory(bdir);
	return true;
}

std::string Util::extension(const std::string &filename) {
	using namespace boost::filesystem;
	path p(filename);
	std::string ext = p.extension().string();
	if(ext.size())
		ext = ext.substr(1);
	lower(ext);
	return ext;

}

size_t Util::dirlist(const std::string &dir, std::vector<std::string> &files,
		const std::string &ext) {
	using namespace boost::filesystem;
	using namespace boost::algorithm;
	if (is_regular_file(dir)) {
		files.push_back(dir);
	} else {
		directory_iterator end;
		directory_iterator di(dir);
		for (; di != end; ++di) {
			if (!ext.empty()) {
				std::string p(di->path().string());
				to_lower(p);
				if (ends_with(p, ext))
					files.push_back(p);
			} else {
				files.push_back(di->path().string());
			}
		}
	}
	return files.size();
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

const std::string Util::tmpFile() {
	return Util::tmpFile("");
}

const std::string Util::tmpFile(const std::string &root) {
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
		m_filename(filename), m_size(size), m_remove(remove) {

	using namespace boost::interprocess;
	using namespace boost::filesystem;

	{
		std::filebuf fbuf;
		if (size > 0) {
			fbuf.open(filename,
					std::ios_base::in | std::ios_base::out
							| std::ios_base::trunc | std::ios_base::binary);
			fbuf.pubseekoff(size - 1, std::ios_base::beg);
			fbuf.sputc(0);
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

MappedFile::~MappedFile() {
	if (m_remove)
		file_mapping::remove(m_filename.c_str());
	Util::rm(m_filename);
	delete m_region;
	delete m_mapping;
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


