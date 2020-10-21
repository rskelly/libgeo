
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
	m_fileBounds{geo::maxvalue<double>(), geo::maxvalue<double>(), geo::minvalue<double>(), geo::minvalue<double>(), geo::maxvalue<double>(), geo::minvalue<double>()},
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
	m_fileBounds{geo::maxvalue<double>(), geo::maxvalue<double>(), geo::minvalue<double>(), geo::minvalue<double>(), geo::maxvalue<double>(), geo::minvalue<double>()},
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

	liblas::ReaderFactory f;
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

PCWriter::PCWriter(const std::string& filename, const liblas::Header& hdr, double x, double y, double size, double buffer) :
	m_fileIdx(0),
	m_returns(0), m_retNum{0,0,0,0,0},
	m_totalReturns(0),
	m_bounds{x, y, x + size, y + size},
	m_bufferedBounds{x - buffer, y - buffer, x + size + buffer, y + size + buffer},
	m_outBounds{geo::maxvalue<double>(), geo::maxvalue<double>(), geo::minvalue<double>(), geo::minvalue<double>(), geo::maxvalue<double>(), geo::minvalue<double>()},
	m_x(x), m_y(y),
	m_filename(filename),
	m_writer(nullptr), m_header(nullptr),
	m_dod(true),
	m_buffer(buffer),
	m_size(size) {

	m_header = new liblas::Header(hdr);

	open();
}

PCWriter::PCWriter(PCWriter&& other) = default;

double PCWriter::x() const {
	return m_x;
}

double PCWriter::y() const {
	return m_y;
}

void PCWriter::outBounds(double* bounds) const {
	for(int i = 0; i < 6; ++i)
		bounds[i] = m_outBounds[i];
}

void PCWriter::bufferedBounds(double* bounds) const {
	for(int i = 0; i < 4; ++i)
		bounds[i] = m_bufferedBounds[i];
}

void PCWriter::bounds(double* bounds) const {
	for(int i = 0; i < 4; ++i)
		bounds[i] = m_bufferedBounds[i];
}

const std::vector<std::string>& PCWriter::filenames() const {
	return m_filenames;
}

std::string PCWriter::nextFile() {
	std::stringstream ss;
	ss << m_filename << "_" << ++m_fileIdx << ".las";
	std::string filename = ss.str();
	m_filenames.push_back(filename);
	return filename;
}

void PCWriter::open() {
	close();
	std::string filename = nextFile();
	m_str.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
	m_writer = new liblas::Writer(m_str, *m_header);
}

void PCWriter::deleteOnDestruct(bool dod) {
	m_dod = dod;
}

void PCWriter::close() {
	try {
		if(m_writer) {
			m_header->SetMin(m_outBounds[0], m_outBounds[1], m_outBounds[4]);
			m_header->SetMax(m_outBounds[2], m_outBounds[3], m_outBounds[5]);
			for(int i = 0; i < 5; ++i)
				m_header->SetPointRecordsByReturnCount(i, m_retNum[i]);
			m_header->SetPointRecordsCount(m_returns);
			m_writer->SetHeader(*m_header);
			m_writer->WriteHeader();
			m_str.close();
			delete m_writer;
			m_writer = nullptr;

			m_returns = 0;
			for(int i = 0; i < 5; ++i)
				m_retNum[i] = 0;
		}
	} catch(const std::exception& ex) {
		std::cerr << "Fail in close: " << ex.what() << "\n";
	}
}

void PCWriter::addPoint(const liblas::Point& pt) {
	if(m_returns >= maxPoints || !m_writer)
		open();

	double x = pt.GetX();
	double y = pt.GetY();
	double z = pt.GetZ();

	++m_returns;
	++m_totalReturns;
	m_retNum[pt.GetReturnNumber() - 1]++;

	if(x < m_outBounds[0]) m_outBounds[0] = x;
	if(x > m_outBounds[2]) m_outBounds[2] = x;
	if(y < m_outBounds[1]) m_outBounds[1] = y;
	if(y > m_outBounds[3]) m_outBounds[3] = y;
	if(z < m_outBounds[4]) m_outBounds[4] = z;
	if(z > m_outBounds[5]) m_outBounds[5] = z;

	try {
		m_writer->WritePoint(pt);
	} catch(const std::exception& ex) {
		std::cerr << "Fail in addPoint: " << ex.what() << "\n";
	}
}

size_t PCWriter::count() const {
	return m_totalReturns;
}

double PCWriter::size() const {
	return m_size;
}

bool PCWriter::contains(double x, double y) const {
	return x >= m_bounds[0] && x < m_bounds[2]  && y >= m_bounds[1] && y < m_bounds[3];
}

bool PCWriter::containsBuffered(double x, double y) const {
	return x >= m_bufferedBounds[0] && x < m_bufferedBounds[2] && y >= m_bufferedBounds[1] && y < m_bufferedBounds[3];
}

PCWriter::~PCWriter() {
	close();
	delete m_header;
	if(!m_totalReturns || m_dod) {
		for(const std::string& f : m_filenames)
			rem(f);
	}
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


Tiler::Tiler(const std::vector<std::string> filenames) {
	for(const std::string& filename : filenames)
		files.emplace_back(filename);
}

void Tiler::tile(const std::string& outdir, double size, double buffer, int,
	double easting, double northing, int maxFileHandles) {

	if(buffer < 0)
		g_runerr("Negative buffer is not allowed. Use easting, northing and tile size to crop tiles.");

	// Calculate the overall bounds of the file set.
	double allBounds[6] = {geo::maxvalue<double>(), geo::maxvalue<double>(), geo::minvalue<double>(), geo::minvalue<double>(), geo::maxvalue<double>(), geo::minvalue<double>()};
	double fBounds[6];
	for(PCFile& f : files) {
		f.init();
		f.fileBounds(fBounds);
		if(fBounds[0] < allBounds[0]) allBounds[0] = fBounds[0];
		if(fBounds[1] < allBounds[1]) allBounds[1] = fBounds[1];
		if(fBounds[2] > allBounds[2]) allBounds[2] = fBounds[2];
		if(fBounds[3] > allBounds[3]) allBounds[3] = fBounds[3];
		if(fBounds[4] < allBounds[4]) allBounds[4] = fBounds[4];
		if(fBounds[5] > allBounds[5]) allBounds[5] = fBounds[5];
	}

	// If the easting and northing aren't given, calculate as a
	// multiple of size.
	if(std::isnan(easting) || std::isnan(northing)) {
		easting = ((int) (allBounds[0] / size)) * size - size;
		northing = ((int) (allBounds[1] / size)) * size - size;

		// Reset the bounds using easting, northing and multiples of size.
		allBounds[0] = easting;
		allBounds[1] = northing;
		allBounds[2] = std::ceil(allBounds[2] / size) * size + size;
		allBounds[3] = std::ceil(allBounds[3] / size) * size + size;
	} else {
		allBounds[0] = easting;
		allBounds[1] = northing;
		allBounds[2] = easting + size;
		allBounds[3] = northing + size;
	}

	// Compute the number of columns and rows of tiles.
	int cols = (int) (allBounds[2] - allBounds[0]) / size;
	int rows = (int) (allBounds[3] - allBounds[1]) / size;

	// Double the tile size until few enough writers are used.
	double size0 = size;
	int cols0 = cols;
	int rows0 = rows;
	while(cols0 * rows0 > maxFileHandles && cols0 > 1 && rows0 > 1) {
		size0 *= 2.0;
		cols0 = even(std::ceil(cols0 * 0.5));
		rows0 = even(std::ceil(rows0 * 0.5));
	}

	liblas::ReaderFactory rfact;
	std::vector<std::unique_ptr<PCWriter> > writers;

	do {

		if(!writers.empty()) {
			// This is run n>0, so we now read from the intermediate tiles.

			// To prevent the previous set of tiles from getting deleted before
			// we're done reading them.
			std::vector<std::unique_ptr<PCWriter> > tmpWriters;

			// Clear the readers list and rebuild with files from the previous
			// level of tiles.
			files.clear();
			for(std::unique_ptr<PCWriter>& writer : writers) {
				writer->close();
				//if(!writer->count()) {
				//	for(const std::string& filename : writer->filenames())
					//	Util::rm(filename);
				//} else {
					// Create and emplace a LASReader.
					files.emplace_back(writer->filenames(), writer->x(), writer->y(), writer->size(), buffer);
					// Close the PCWriter.
					tmpWriters.emplace_back(writer.release());
				//}
			}

			// Clear the tiles list to start rebuilding it with the next level.
			writers.clear();

			for(PCFile& file : files) {
				file.init();
				{
					// Get a liblas::Header from the first file to use as a template
					// for the PCWriters.
					std::ifstream istr(file.filenames()[0], std::ios::in | std::ios::binary);
					liblas::Reader irdr = rfact.CreateWithStream(istr);
					const liblas::Header& ihdr = irdr.GetHeader();

					// Create the writers for the intermediate tiles. There are four tiles
					// associated with each LASReader so only four handles at a time can be open.
					for(int r = 0; r < 2; ++r) {
						for(int c = 0; c < 2; ++c) {
							std::stringstream ss;
							double x = file.x() + c * size0;
							double y = file.y() + r * size0;
							ss << "tile_" << (int) x << "_" << (int) y << "_" << size0;
							std::string outfile = join(outdir, ss.str());
							writers.emplace_back(new PCWriter(outfile, ihdr, x, y, size0, buffer));
						}
					}
				}

				// Now run over the files and reat their points into the new tiles.
				for(const std::string& filename : file.filenames()) {
					std::ifstream istr(filename, std::ios::in | std::ios::binary);
					liblas::Reader irdr = rfact.CreateWithStream(istr);
					while(irdr.ReadNextPoint()) {
						const liblas::Point& pt = irdr.GetPoint();
						for(std::unique_ptr<PCWriter>& writer : writers) {
							double x = pt.GetX();
							double y = pt.GetY();
							if(file.contains(x, y) && writer->containsBuffered(x, y))
								writer->addPoint(pt);
						}
					}
				}

			}

			// Destroy the previous level of tiles.
			tmpWriters.clear();

		} else {
			// This is the first run, so we read all the input files into the first set of
			// intermediate files.

			{
				// Get the header from the first file to use as a template.
				std::ifstream istr(files[0].filenames()[0], std::ios::in | std::ios::binary);
				liblas::Reader irdr = rfact.CreateWithStream(istr);
				const liblas::Header& ihdr = irdr.GetHeader();

				// Create the writers for the intermediate tiles.
				for(int r = 0; r < rows0; ++r) {
					for(int c = 0; c < cols0; ++c) {
						std::stringstream ss;
						double x = allBounds[0] + c * size0;
						double y = allBounds[1] + r * size0;
						ss << "tile_" << (int) x << "_" << (int) y << "_" << size0;
						std::string outfile = join(outdir, ss.str());
						writers.emplace_back(new PCWriter(outfile, ihdr, x, y, size0, buffer));
					}
				}
			}

			// Iterate over the files, filling the tiles.
			for(PCFile& file : files) {

				double fbounds[4];
				file.fileBounds(fbounds);

				double wbounds[4];
				std::vector<PCWriter*> wtrs;
				for(std::unique_ptr<PCWriter>& writer : writers) {
					writer->bufferedBounds(wbounds);
					if(intersects(fbounds, wbounds))
						wtrs.push_back(writer.get());
				}

				for(const std::string& filename : file.filenames()) {
					std::ifstream istr(filename, std::ios::in | std::ios::binary);
					liblas::Reader irdr = rfact.CreateWithStream(istr);
					while(irdr.ReadNextPoint()) {
						const liblas::Point& pt = irdr.GetPoint();
						for(PCWriter* writer : wtrs) {
							double x = pt.GetX();
							double y = pt.GetY();
							if(writer->containsBuffered(x, y))
								writer->addPoint(pt);
						}
					}
				}
			}
		}

		// At the end of each loop, we halve the tile size and
		// double the col/row size to move down the pyramid.
		size0 *= 0.5;
		cols0 *= 2;
		rows0 *= 2;

	} while(size0 >= size);

	// We don't want to delete the last batch of tiles, so this is
	// set to false.
	for(std::unique_ptr<PCWriter>& writer : writers)
		writer->deleteOnDestruct(false);

	// Destroy the last batch of writers.
	writers.clear();
}

Tiler::~Tiler() {}

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

geo::pc::Point::Point(const liblas::Point& pt) :
	geo::pc::Point(pt.GetX(), pt.GetY(), pt.GetZ(),
			pt.GetIntensity(),
			pt.GetScanAngleRank(),
			pt.GetClassification().GetClass(),
			pt.GetReturnNumber(),
			pt.GetNumberOfReturns(),
			pt.GetFlightLineEdge()) {
}

void geo::pc::Point::setPoint(const liblas::Point& pt) {
	m_x = pt.GetX();
	m_y = pt.GetY();
	m_z = pt.GetZ();
	m_intensity = pt.GetIntensity();
	m_angle = pt.GetScanAngleRank();
	m_cls = pt.GetClassification().GetClass();
	m_edge = pt.GetFlightLineEdge();
	m_numReturns = pt.GetNumberOfReturns();
	m_returnNum = pt.GetReturnNumber();
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
	m_minX = DBL_MAX;
	m_minY = DBL_MAX;
	double maxX = -DBL_MAX;
	double maxY = -DBL_MIN;
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

