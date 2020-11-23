
#define _USE_MATH_DEFINES

#include <fstream>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <memory>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdint>
#include <list>
#include <condition_variable>
#include <thread>

#include <errno.h>
#include <math.h>
#include <string.h>

#include <pdal/StageFactory.hpp>
#include <pdal/io/BufferReader.hpp>

#include "util.hpp"
#include "grid.hpp"
#include "pointcloud.hpp"
#include <ds/mqtree.hpp>
#include "pc_computer.hpp"

using namespace geo::grid;
using namespace geo::util;
using namespace geo::pc;

bool geo::pc::pointSort(const geo::pc::Point& a, const geo::pc::Point& b) {
	return a.value() < b.value();
}

PCFile::PCFile(const std::string& filename, double x, double y, double size, double buffer) :
	m_x(x), m_y(y),
	m_fileBounds{
		geo::maxvalue<double>(), geo::maxvalue<double>(),
		geo::minvalue<double>(), geo::minvalue<double>(),
		geo::maxvalue<double>(), geo::minvalue<double>()
	},
	m_bounds{x, y, x + size, y + size},
	m_bufferedBounds{x - buffer, y - buffer, x + size + buffer, y + size + buffer},
	m_pointCount(0),
	m_inited(false),
	m_index(0),
	m_source(nullptr) {

	m_filenames.push_back(filename);
}

PCFile::PCFile(const std::vector<std::string>& filenames, double x, double y, double size, double buffer) :
	m_x(x), m_y(y),
	m_fileBounds{
		geo::maxvalue<double>(), geo::maxvalue<double>(),
		geo::minvalue<double>(), geo::minvalue<double>(),
		geo::maxvalue<double>(), geo::minvalue<double>()
	},
	m_bounds{x, y, x + size, y + size},
	m_bufferedBounds{x - buffer, y - buffer, x + size + buffer, y + size + buffer},
	m_pointCount(0),
	m_inited(false),
	m_index(0),
	m_source(nullptr),
	m_filenames(filenames) {
}

void PCFile::resize(double x, double y, double size, double buffer) {
	m_x = x;
	m_y = y;
	m_bounds[0] = x;
	m_bounds[1] = y;
	m_bounds[2] = x + size;
	m_bounds[3] = y + size;
	m_bufferedBounds[0] = x - buffer;
	m_bufferedBounds[1] = y - buffer;
	m_bufferedBounds[2] = x + size + buffer;
	m_bufferedBounds[3] = y + size + buffer;
}

double PCFile::x() const {
	return m_x;
}

double PCFile::y() const {
	return m_y;
}

size_t PCFile::pointCount() const {
	return m_pointCount;
}

void PCFile::fileBounds(double* bounds) const {
	for(int i = 0; i < 6; ++i)
		bounds[i] = m_fileBounds[i];
}

void PCFile::bounds(double* bounds) const {
	for(int i = 0; i < 4; ++i)
		bounds[i] = m_bounds[i];
}

void PCFile::bufferedBounds(double* bounds) const {
	for(int i = 0; i < 4; ++i)
		bounds[i] = m_bufferedBounds[i];
}

const std::vector<std::string>& PCFile::filenames() const {
	return m_filenames;
}

void PCFile::init(bool useHeader) {
	if(m_inited)
		return;
	m_pointCount = 0;

	for(const std::string& filename : m_filenames) {
		g_debug("Initializing: " << filename);
		PDALSource src(filename);
		if(useHeader) {
			double minx = src.hdr.minX();
			double miny = src.hdr.minY();
			double minz = src.hdr.minZ();
			double maxx = src.hdr.maxX();
			double maxy = src.hdr.maxY();
			double maxz = src.hdr.maxZ();
			if(minx < m_fileBounds[0]) m_fileBounds[0] = minx;
			if(miny < m_fileBounds[1]) m_fileBounds[1] = miny;
			if(maxx > m_fileBounds[2]) m_fileBounds[2] = maxx;
			if(maxy > m_fileBounds[3]) m_fileBounds[3] = maxy;
			if(minz < m_fileBounds[4]) m_fileBounds[4] = minz;
			if(maxz > m_fileBounds[5]) m_fileBounds[5] = maxz;
			m_pointCount += src.hdr.pointCount();
		} else {
			using namespace pdal::Dimension;
			for (pdal::PointId idx = 0; idx < src.view->size(); ++idx) {
				double x = src.view->getFieldAs<double>(Id::X, idx);
				double y = src.view->getFieldAs<double>(Id::Y, idx);
				double z = src.view->getFieldAs<double>(Id::Z, idx);
				if(x < m_fileBounds[0]) m_fileBounds[0] = x;
				if(y < m_fileBounds[1]) m_fileBounds[1] = y;
				if(x > m_fileBounds[2]) m_fileBounds[2] = x;
				if(y > m_fileBounds[3]) m_fileBounds[3] = y;
				if(z < m_fileBounds[4]) m_fileBounds[4] = z;
				if(z > m_fileBounds[5]) m_fileBounds[5] = z;
				++m_pointCount;
			}
		}
	}
	m_inited = true;
}

bool PCFile::contains(double x, double y) const {
	return x >= m_bounds[0] && x < m_bounds[2]  && y >= m_bounds[1] && y < m_bounds[3];
}

bool PCFile::intersects(double* b) const {
	const double* a = m_fileBounds;
	double aa[4] {std::min(a[0], a[2]), std::min(a[1], a[3]), std::max(a[0], a[2]), std::max(a[1], a[3])};
	double bb[4] {std::min(b[0], b[2]), std::min(b[1], b[3]), std::max(b[0], b[2]), std::max(b[1], b[3])};
	return !(aa[2] <= bb[0] || aa[0] >= bb[2] || aa[3] <= bb[1] || aa[1] >= bb[3]);
}

bool PCFile::next(geo::pc::Point& pt) {
	if((isReaderOpen() || openReader()) && m_source->nextPoint(pt)) {
		return true;
	}
	return false;
}

bool PCFile::openReader() {
	closeReader();
	if(m_index < m_filenames.size()) {
		m_source = new PDALSource(m_filenames[m_index]);
		++m_index;
		return true;
	}
	return false;
}

void PCFile::closeReader() {
	if(m_source) {
		delete m_source;
		m_source = nullptr;
	}
}

bool PCFile::isReaderOpen() const {
	return m_source != nullptr;
}

bool PCFile::containsBuffered(double x, double y) const {
	return x >= m_bufferedBounds[0] && x < m_bufferedBounds[2] && y >= m_bufferedBounds[1] && y < m_bufferedBounds[3];
}

PCFile::~PCFile() {
	closeReader();
}



const long maxPoints = 20000000;

PCWriter::PCWriter(const std::string& /*filename*/, const pdal::LasHeader& /*hdr*/, double /*x*/, double /*y*/, double /*size*/, double /*buffer*/) :
		m_writer(nullptr), m_header(nullptr),
		m_idx(0) {
	/*
	m_fileIdx(0),
	m_returns(0), m_retNum{0,0,0,0,0},
	m_totalReturns(0),
	m_bounds{x, y, x + size, y + size},
	m_bufferedBounds{x - buffer, y - buffer, x + size + buffer, y + size + buffer},
	m_outBounds{
		geo::maxvalue<double>(), geo::maxvalue<double>(),
		geo::minvalue<double>(), geo::minvalue<double>(),
		geo::maxvalue<double>(), geo::minvalue<double>()
	},
	m_x(x), m_y(y),
	m_filename(filename),
	m_writer(nullptr), m_header(nullptr),
	m_dod(true),
	m_buffer(buffer),
	m_size(size) {

	m_header = new pdal::LasHeader(hdr);

	open();
	*/
}

PCWriter::PCWriter(const std::string& filename, const std::string& tpl) :
	m_writer(nullptr), m_header(nullptr),
	m_idx(0) {
	/*
	m_fileIdx(0),
	m_returns(0), m_retNum{0,0,0,0,0},
	m_totalReturns(0),
	m_bounds{x, y, x + size, y + size},
	m_outBounds{
		geo::maxvalue<double>(), geo::maxvalue<double>(),
		geo::minvalue<double>(), geo::minvalue<double>(),
		geo::maxvalue<double>(), geo::minvalue<double>()
	},
	m_x(x), m_y(y),
	m_filename(filename),
	m_writer(nullptr), m_header(nullptr),
	m_dod(true),
	m_buffer(0),
	m_size(size) {
	*/
	m_tpl = tpl;
	m_filename = filename;

	open();
}

//PCWriter::PCWriter(PCWriter&& other) = default;

double PCWriter::x() const {
	return 0; //m_x;
}

double PCWriter::y() const {
	return 0; //m_y;
}

void PCWriter::outBounds(double* /*bounds*/) const {
	/*
	for(int i = 0; i < 6; ++i)
		bounds[i] = m_outBounds[i];
		*/
}

void PCWriter::bufferedBounds(double* /*bounds*/) const {
	/*
	for(int i = 0; i < 4; ++i)
		bounds[i] = m_bufferedBounds[i];
		*/
}

void PCWriter::bounds(double* /*bounds*/) const {
	/*
	for(int i = 0; i < 4; ++i)
		bounds[i] = m_bufferedBounds[i];
		*/
}

/*
const std::vector<std::string>& PCWriter::filenames() const {
	return m_filenames;
}
*/

/*
std::string PCWriter::nextFile() {
	std::stringstream ss;
	ss << m_filename << "_" << ++m_fileIdx << ".las";
	std::string filename = ss.str();
	m_filenames.push_back(filename);
	return filename;
}
*/

void PCWriter::open() {
	close();
	PDALSource src(m_tpl);
	pdal::StageFactory fact;
	pdal::Options opts;
	opts.add("filename", m_filename);
	m_writer = fact.createStage("writers.las");
	m_writer->setOptions(opts);
	m_writer->setSpatialReference(src.reader.getSpatialReference());
	m_table.layout()->registerDims(src.table.layout()->dims());
	m_view.reset(new pdal::PointView(m_table));
}

void PCWriter::deleteOnDestruct(bool /*dod*/) {
	//m_dod = dod;
}

void PCWriter::close() {
	try {
		if(m_writer) {
			pdal::BufferReader rdr;
			rdr.addView(m_view);
			m_writer->setInput(rdr);
			m_writer->prepare(m_table);
			m_writer->execute(m_table);
			delete m_writer;
			m_writer = nullptr;
		}
	} catch(const std::exception& ex) {
		std::cerr << "Fail in close: " << ex.what() << "\n";
	}
}

void PCWriter::addPoint(const geo::pc::Point& pt) {
	if(!m_writer)
		open();

	pdal::PointRef ptr(m_table, m_idx);
	pt.getPoint(ptr);
	++m_idx;

}

size_t PCWriter::count() const {
	return m_idx;
}

double PCWriter::size() const {
	return 0; //m_size;
}

bool PCWriter::contains(double /*x*/, double /*y*/) const {
	return false;//x >= m_bounds[0] && x < m_bounds[2]  && y >= m_bounds[1] && y < m_bounds[3];
}

bool PCWriter::containsBuffered(double /*x*/, double /*y*/) const {
	return false; //x >= m_bufferedBounds[0] && x < m_bufferedBounds[2] && y >= m_bufferedBounds[1] && y < m_bufferedBounds[3];
}

PCWriter::~PCWriter() {
	close();
}

int even(int num) {
	if(num % 2 == 1)
		++num;
	return num;
}

bool intersects(double* a, double* b) {
	double aa[4] {std::min(a[0], a[2]), std::min(a[1], a[3]), std::max(a[0], a[2]), std::max(a[1], a[3])};
	double bb[4] {std::min(b[0], b[2]), std::min(b[1], b[3]), std::max(b[0], b[2]), std::max(b[1], b[3])};
	return !(aa[2] <= bb[0] || aa[0] >= bb[2] || aa[3] <= bb[1] || aa[1] >= bb[3]);
}


geo::pc::Point::Point(double x, double y, double z, double intensity, double angle,
		int cls, int returnNum, int numReturns, bool isEdge) :
	m_x(x),
	m_y(y),
	m_z(z),
	m_intensity(intensity),
	m_angle(angle),
	m_cls(cls),
	m_edge(isEdge),
	m_numReturns(numReturns),
	m_returnNum(returnNum) {
}

geo::pc::Point::Point(const geo::pc::Point& pt) :
		geo::pc::Point(pt.x(), pt.y(), pt.z(), pt.intensity(), pt.scanAngle(),
				pt.classId(), pt.returnNum(), pt.numReturns(), pt.isEdge()) {}

geo::pc::Point::Point(const pdal::PointRef& pt) {
	setPoint(pt);
}

void geo::pc::Point::setPoint(const pdal::PointRef& pt) {
	using namespace pdal::Dimension;
	m_x = pt.getFieldAs<double>(Id::X);
	m_y = pt.getFieldAs<double>(Id::Y);
	m_z = pt.getFieldAs<double>(Id::Z);
	m_intensity = pt.getFieldAs<short>(Id::Z);
	m_angle = pt.getFieldAs<char>(Id::ScanAngleRank);
	m_cls = pt.getFieldAs<char>(Id::Classification);
	m_edge = pt.getFieldAs<bool>(Id::EdgeOfFlightLine);
	m_numReturns = pt.getFieldAs<char>(Id::NumberOfReturns);
	m_returnNum = pt.getFieldAs<char>(Id::ReturnNumber);
}

void geo::pc::Point::getPoint(pdal::PointRef& pt) const {
	using namespace pdal::Dimension;
	pt.setField(Id::X, m_x);
	pt.setField(Id::Y, m_y);
	pt.setField(Id::Z, m_z);
	pt.setField(Id::Z, m_intensity);
	pt.setField(Id::ScanAngleRank, m_angle);
	pt.setField(Id::Classification, m_cls);
	pt.setField(Id::EdgeOfFlightLine, m_edge);
	pt.setField(Id::NumberOfReturns, m_numReturns);
	pt.setField(Id::ReturnNumber, m_returnNum);
}

geo::pc::Point::Point(double x, double y, double z) :
	geo::pc::Point() {
	m_x = x;
	m_y = y;
	m_z = z;
}

geo::pc::Point::Point() :
	m_x(0),
	m_y(0),
	m_z(0),
	m_intensity(0),
	m_angle(0),
	m_cls(0),
	m_edge(0),
	m_numReturns(0),
	m_returnNum(0) {
}

int geo::pc::Point::classId() const {
	return m_cls;
}

void geo::pc::Point::classId(int cls) {
	m_cls = cls;
}

double geo::pc::Point::operator[](int idx) const {
	if(idx % 2 == 0) {
		return m_x;
	} else {
		return m_y;
	}
}

double geo::pc::Point::x() const {
	return m_x;
}

double geo::pc::Point::y() const {
	return m_y;
}

double geo::pc::Point::z() const {
	return m_z;
}

void geo::pc::Point::x(double x) {
	m_x = x;
}

void geo::pc::Point::y(double y) {
	m_y = y;
}

void geo::pc::Point::z(double z) {
	m_z = z;
}

double geo::pc::Point::value() const {
	return z();
}

double geo::pc::Point::intensity() const {
	return m_intensity;
}

void geo::pc::Point::intensity(double intensity) {
	m_intensity = intensity;
}

double geo::pc::Point::scanAngle() const {
	return m_angle;
}

void geo::pc::Point::scanAngle(double angle) {
	m_angle = angle;
}

bool geo::pc::Point::isEdge() const {
	return m_edge == 1;
}

void geo::pc::Point::isEdge(bool isEdge) {
	m_edge = isEdge;
}

bool geo::pc::Point::isLast() const {
	return returnNum() == numReturns();
}

bool geo::pc::Point::isFirst() const {
	return returnNum() == 1;
}

int geo::pc::Point::returnNum() const {
	return m_returnNum;
}

void geo::pc::Point::returnNum(int returnNum) {
	m_returnNum = returnNum;
}

int geo::pc::Point::numReturns() const {
	return m_numReturns;
}

void geo::pc::Point::numReturns(int numReturns) {
	m_numReturns = numReturns;
}

bool geo::pc::Point::operator<(const Point& other) const {
	return x() == other.x() ? y() < other.y() : x() < other.x();
}

geo::pc::Point::~Point() {
}



bool pctBoundsSort(const PCFile& a, const PCFile& b) {
	double ba[6];
	double bb[6];
	a.bounds(ba);
	b.bounds(bb);
	return ba[1] < bb[1] && ba[0] < bb[0];
}

PCTreeIterator::PCTreeIterator(const std::vector<std::string>& files, double size, double buffer) :
	m_idx(-1), m_cols(-1), m_rows(-1),
	m_minX(0), m_minY(0),
	m_size(size),
	m_buffer(buffer) {
	for(const std::string& filename : files)
		m_files.emplace_back(filename);
}

void PCTreeIterator::init() {
	std::sort(m_files.begin(), m_files.end(), pctBoundsSort);
	m_minX = geo::maxvalue<double>();
	m_minY = geo::maxvalue<double>();
	double maxX = -geo::maxvalue<double>();
	double maxY = -geo::maxvalue<double>();
	double bounds[6];
	for(PCFile& file : m_files) {
		file.init();
		file.fileBounds(bounds);
		if(bounds[0] < m_minX) m_minX = bounds[0];
		if(bounds[1] < m_minY) m_minY = bounds[1];
		if(bounds[2] > maxX) maxX = bounds[2];
		if(bounds[3] > maxY) maxY = bounds[3];
	}
	m_cols = (int) std::ceil((maxX - m_minX) / m_size);
	m_rows = (int) std::ceil((maxY - m_minY) / m_size);
}

void PCTreeIterator::reset() {
	m_idx = -1;
}

bool PCTreeIterator::next(geo::ds::mqtree<geo::pc::Point>& tree) {
	tree.clear();
	if(++m_idx >= m_cols * m_rows)
		return false;
	int col = m_idx % m_cols;
	int row = m_idx / m_cols;
	int minX = m_minX + col * m_size;
	int minY = m_minY + row * m_size;
	double tileBounds[4] { minX - m_buffer, minY - m_buffer, minX + m_size + m_buffer, minY + m_size + m_buffer };
	geo::pc::Point pt;
	for(PCFile& file : m_files) {
		if(file.intersects(tileBounds)) {
			while(file.next(pt)) {
				double x = pt.x();
				double y = pt.y();
				if(x >= tileBounds[0] && x <= tileBounds[2] && y >= tileBounds[1] && y <= tileBounds[2])
					tree.add(pt);
			}
		}
	}
	return true;
}

