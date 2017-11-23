#include <iterator>
#include <fstream>
#include <algorithm>

#include <ogr_spatialref.h>

#include <boost/filesystem.hpp>

#include "crypto/md5.hpp"
#include "crypto/uuid.hpp"
#include "util.hpp"

using namespace geo::util;

Stopwatch::Stopwatch() : m_reset(false), m_running(false) {
	m_stop = m_start = std::chrono::system_clock::now();
	reset();
}

void Stopwatch::reset() {
	m_reset = true;
}

void Stopwatch::start() {
	if(m_reset)
		m_start = std::chrono::system_clock::now();
	m_running = true;
}

void Stopwatch::stop() {
	m_stop = std::chrono::system_clock::now();
	m_running = false;
}
std::string Stopwatch::time() {
	uint64_t sec = millis() / 1000;
	std::stringstream ss;
	uint64_t h = sec / 3600;
	uint64_t m = (sec / 60) % 60;
	uint64_t s = sec % 60;
	ss << h << ":" << (m > 9 ? "" : "0") << m << ":" << (s > 9 ? "" : "0") << s;
	return ss.str();
}

uint64_t Stopwatch::millis() {
	std::chrono::time_point<std::chrono::system_clock> now = m_running ? std::chrono::system_clock::now() : m_stop;
	uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - m_start).count();
	return ms;
}

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
	m_callbacks(callbacks), m_start(start), m_end(end) {
}

Callbacks* Status::callbacks() const {
	return m_callbacks;
}

float Status::start() const {
	return m_start;
}

float Status::end()  const {
	return m_end;
}

void Status::update(float s, const std::string& msg) {
	if(!m_callbacks) {
		g_debug("Status: " << s << ", " << msg);
	} else {
		m_callbacks->stepCallback(m_start + (m_end - m_start) * s);
		if(!msg.empty())
			m_callbacks->statusCallback(msg);
	}
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

std::vector<std::pair<int, int> > Util::circularKernel(int outerRadius, int innerRadius, bool includeCenter) {
	double outr = g_sq((double) outerRadius) + 1.0;
	double inr = g_sq((double) innerRadius);
	std::vector<std::pair<int, int> > offsets;
	for(int r = -outerRadius; r < outerRadius + 1; ++r) {
		for(int c = -outerRadius; c < outerRadius + 1; ++c) {
			double d0 = g_sq((double) c) + g_sq((double) r);
			if((d0 <= outr && d0 >= inr) || (includeCenter && r == 0 && c == 0)) {
				offsets.push_back(std::make_pair(c, r));
				std::cerr << "x ";
			} else {
				std::cerr << "  ";
			}
		}
		std::cerr << "\n";
	}
	return offsets;
}

std::vector<std::pair<int, int> > Util::squareKernel(int size, bool includeCenter) {
	std::vector<std::pair<int, int> > offsets;
	for(int r = -size / 2; r < size / 2 + 1; ++r) {
		for(int c = -size / 2; c < size / 2 + 1; ++c) {
			if(includeCenter || !(r == 0 && c == 0))
				offsets.push_back(std::make_pair(c, r));
		}
	}
	return offsets;
}

void Util::copyfile(const std::string &srcfile, const std::string &dstfile) {
	std::ifstream src(srcfile.c_str(), std::ios::binary);
	std::ofstream dst(dstfile.c_str(), std::ios::binary);
	dst << src.rdbuf();
}

bool Util::exists(const std::string &name) {
	boost::filesystem::path p(name);
	return boost::filesystem::exists(p);
}

std::string Util::basename(const std::string &filename) {
	boost::filesystem::path p(filename);
	return boost::filesystem::basename(filename);
}

std::string Util::filename(const std::string &filename) {
	boost::filesystem::path p(filename);
	return p.filename().string();
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

uint64_t Util::filesize(const std::string& name) {
    struct stat stat_buf;
    int rc = stat(name.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}

bool Util::mkdir(const std::string &dir) {
	using namespace boost::filesystem;
	path bdir(dir);
	if (!boost::filesystem::exists(bdir))
		return create_directories(bdir);
	return true;
}

bool Util::isDir(const std::string& path) {
	boost::filesystem::path p(path);
	return boost::filesystem::is_directory(p);

}

bool Util::isFile(const std::string& path) {
	boost::filesystem::path p(path);
	return boost::filesystem::is_regular(p);

}

std::string Util::parent(const std::string& file) {
	using namespace boost::filesystem;
	path p(file);
	return p.parent_path().string();
}

std::string Util::extension(const std::string &filename) {
	using namespace boost::filesystem;
	path p(filename);
	return p.extension().string();
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

std::string Util::md5(const std::string& input) {
	geo::crypto::MD5 m(input);
	m.finalize();
	return m.hexdigest();
}

#include "openssl/sha.h"

std::string Util::sha256(const std::string& input) {
	const char* str = input.c_str();
	char outputBuffer[65];
	unsigned char hash[SHA256_DIGEST_LENGTH];
	SHA256_CTX sha256;
	SHA256_Init(&sha256);
	SHA256_Update(&sha256, str, strlen(str));
	SHA256_Final(hash, &sha256);
	for(int i = 0; i < SHA256_DIGEST_LENGTH; i++)
		sprintf(outputBuffer + (i * 2), "%02x", hash[i]);
	outputBuffer[64] = 0;
	return std::string(outputBuffer);
}

std::string Util::sha256File(const std::string& file) {
	const char *path = file.c_str();
	char outputBuffer[65];
	std::FILE *f = std::fopen(path, "rb");
	if(!f)
		throw std::runtime_error("Failed to open file for hashing.");
	unsigned char hash[SHA256_DIGEST_LENGTH];
	SHA256_CTX sha256;
	SHA256_Init(&sha256);
	const int bufSize = 32768;
	char *buffer = (char*) malloc(bufSize);
	if(!buffer)
		throw std::runtime_error("Failed to allocate buffer for hashing.");
	int bytesRead = 0;
	while((bytesRead = std::fread(buffer, 1, bufSize, f)))
		SHA256_Update(&sha256, buffer, bytesRead);
	std::fclose(f);
	SHA256_Final(hash, &sha256);
	for(int i = 0; i < SHA256_DIGEST_LENGTH; i++)
		sprintf(outputBuffer + (i * 2), "%02x", hash[i]);
	free(buffer);
	return std::string(outputBuffer);
}

MappedFile::MappedFile(const std::string& name, uint64_t size) :
	m_size(0),
	m_name(name),
	m_region(nullptr),
	m_shm(nullptr) {

	if (size > 0)
		reset(size);
}

MappedFile::MappedFile(uint64_t size) :
	m_size(0),
	m_name(geo::crypto::UUID::uuid()),
	m_region(nullptr),
	m_shm(nullptr) {

	if (size > 0)
		reset(size);
}

const std::string& MappedFile::name() const {
	return m_name;
}

size_t MappedFile::fixSize(size_t size) {
	size_t page = MappedFile::pageSize();
	if(size % page == 0) {
		return size;
	} else {
		return (size / page + 1) * page;
	}
}

bool MappedFile::write(void* input, uint64_t position, uint64_t length) {
    if(size() < position + length)
        reset(position + length);
    std::memcpy((char*) data() + position, input, length);
    return true;
}

bool MappedFile::read(void* output, uint64_t position, uint64_t length) {
    if(position + length > size())
        return false;
    std::memcpy(output, (char*) data() + position, length);
    return true;
}

void MappedFile::reset(uint64_t size) {
	using namespace boost::interprocess;
	#pragma omp critical(__mapped_file_reset__)
	{
		size = fixSize(size);
		if(size != 0 && size != m_size) {
			m_size = size;
			if (m_region) {
				m_region->flush();
				delete m_region;
			}
			if (m_shm)
				delete m_shm;
			m_shm = new shared_memory_object(open_or_create, m_name.c_str(), read_write);
			m_shm->truncate(m_size);
			m_region = new mapped_region(*m_shm, read_write);
		}
	}
}

void* MappedFile::data() {
	return m_region ? m_region->get_address() : nullptr;
}


uint64_t MappedFile::size() const {
	return m_size;
}

size_t MappedFile::pageSize() {
	return boost::interprocess::mapped_region::get_page_size();
}

MappedFile::~MappedFile() {
	if(m_region) {
		m_region->flush();
		delete m_region;
		delete m_shm;
	}
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
