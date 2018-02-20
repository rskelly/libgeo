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

#include <errno.h>
#include <math.h>
#include <sys/mman.h>
#include <string.h>

#include <liblas/liblas.hpp>

#include "util.hpp"
#include "raster.hpp"
#include "pointcloud.hpp"
#include "ds/kdtree.hpp"
#include "pc_computer.hpp"

#define MIN_DBL std::numeric_limits<double>::lowest()
#define MAX_DBL std::numeric_limits<double>::max()
#define NODATA -9999.0

using namespace geo::raster;
using namespace geo::util;
using namespace geo::pc;

PCFile::PCFile(const std::string& filename, double x, double y, double size, double buffer) :
	m_x(x), m_y(y),
	m_fileBounds{MAX_DBL, MAX_DBL, MIN_DBL, MIN_DBL, MAX_DBL, MIN_DBL},
	m_bounds{x, y, x + size, y + size},
	m_bufferedBounds{x - buffer, y - buffer, x + size + buffer, y + size + buffer},
	m_pointCount(0),
	m_inited(false) {

	m_filenames.push_back(filename);
}

PCFile::PCFile(const std::vector<std::string>& filenames, double x, double y, double size, double buffer) :
	m_x(x), m_y(y),
	m_fileBounds{MAX_DBL, MAX_DBL, MIN_DBL, MIN_DBL, MAX_DBL, MIN_DBL},
	m_bounds{x, y, x + size, y + size},
	m_bufferedBounds{x - buffer, y - buffer, x + size + buffer, y + size + buffer},
	m_pointCount(0),
	m_inited(false),
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
		std::ifstream str(filename, std::ios::in | std::ios::binary);
		liblas::Reader reader = f.CreateWithStream(str);
		if(useHeader) {
			const liblas::Header& hdr = reader.GetHeader();
			double minx = hdr.GetMinX();
			double miny = hdr.GetMinY();
			double minz = hdr.GetMinZ();
			double maxx = hdr.GetMaxX();
			double maxy = hdr.GetMaxY();
			double maxz = hdr.GetMaxZ();
			if(minx < m_fileBounds[0]) m_fileBounds[0] = minx;
			if(miny < m_fileBounds[1]) m_fileBounds[1] = miny;
			if(maxx > m_fileBounds[2]) m_fileBounds[2] = maxx;
			if(maxy > m_fileBounds[3]) m_fileBounds[3] = maxy;
			if(minz < m_fileBounds[4]) m_fileBounds[4] = minz;
			if(maxz > m_fileBounds[5]) m_fileBounds[5] = maxz;
			m_pointCount += hdr.GetPointRecordsCount();
		} else {
			while(reader.ReadNextPoint()) {
				const liblas::Point& pt = reader.GetPoint();
				double x = pt.GetX();
				double y = pt.GetY();
				double z = pt.GetZ();
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

bool PCFile::containsBuffered(double x, double y) const {
	return x >= m_bufferedBounds[0] && x < m_bufferedBounds[2] && y >= m_bufferedBounds[1] && y < m_bufferedBounds[3];
}

PCFile::~PCFile() {}

PCPointFilter::PCPointFilter() :
	minScanAngle(-90),
	maxScanAngle(90),
	keepEdges(false),
	minZ(DBL_MIN),
	maxZ(DBL_MAX),
	minIntensity(DBL_MIN),
	maxIntensity(DBL_MAX),
	lastOnly(false),
	firstOnly(false) {
}

void PCPointFilter::printHelp(std::ostream& str) {
	str << " Point filtering parameters:\n"
		<< " -p:c <class(es)>     Comma-delimited list of classes to keep.\n"
		<< " -p:minz <z>          The minimum height threshold.\n"
		<< " -p:maxz <z>          The maximum height threshold.\n"
		<< " -p:mini <intensity>  The minimum intensity.\n"
		<< " -p:maxi <intensity>  The maximum intensity.\n"
		<< " -p:mina <angle>      The minimum scan angle.\n"
		<< " -p:maxa <angle>      The maximum scan angle.\n"
		<< " -p:f                 First returns only\n"
		<< " -p:l                 Last returns only\n";
}

bool PCPointFilter::parseArgs(int& idx, char** argv) {
	std::string v = argv[idx];
	bool found = false;
	if(v == "-p:c") {
		std::vector<std::string> tmp;
		std::string cls = argv[++idx];
		Util::splitString(std::back_inserter(tmp), cls);
		for(const std::string& t : tmp)
			classes.push_back(atoi(t.c_str()));
		found = true;
	} else if(v == "-p:l") {
		lastOnly = true;
		found = true;
	} else if(v == "-p:f") {
		firstOnly = true;
		found = true;
	} else if(v == "-p:minz") {
		minZ = atof(argv[++idx]);
		found = true;
	} else if(v == "-p:maxz") {
		maxZ = atof(argv[++idx]);
		found = true;
	} else if(v == "-p:mini") {
		minIntensity = atof(argv[++idx]);
		found = true;
	} else if(v == "-p:maxi") {
		maxIntensity = atof(argv[++idx]);
		found = true;
	} else if(v == "-p:mina") {
		minScanAngle = atof(argv[++idx]);
		found = true;
	} else if(v == "-p:maxa") {
		maxScanAngle = atof(argv[++idx]);
		found = true;
	} else if(v == "-p:e") {
		keepEdges = true;
		found = true;
	}
	return found;
}

bool PCPointFilter::keep(const geo::pc::Point& pt) const {
	if(lastOnly && !pt.isLast())
		return false;
	if(firstOnly && !pt.isFirst())
		return false;
	double z = pt.z();
	if(z < minZ || z > maxZ ||
			pt.intensity() < minIntensity || pt.intensity() > maxIntensity ||
			pt.scanAngle() < minScanAngle || pt.scanAngle() > maxScanAngle ||
			(pt.isEdge() && !keepEdges))
		return false;
	if(!classes.empty()) {
		int cls = pt.classId();
		for(size_t i = 0; i < classes.size(); ++i) {
			if(cls == classes[i])
				return true;
		}
		return false;
	}
	return true;
}


const long maxPoints = 20000000;

PCWriter::PCWriter(const std::string& filename, const liblas::Header& hdr, double x, double y, double size, double buffer) :
	m_fileIdx(0),
	m_returns(0), m_retNum{0,0,0,0,0},
	m_totalReturns(0),
	m_bounds{x, y, x + size, y + size},
	m_bufferedBounds{x - buffer, y - buffer, x + size + buffer, y + size + buffer},
	m_outBounds{MAX_DBL, MAX_DBL, MIN_DBL, MIN_DBL, MAX_DBL, MIN_DBL},
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
	std::cerr << filename << "\n";
	m_str.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
	m_writer = new liblas::Writer(m_str, *m_header);
}

void PCWriter::deleteOnDestruct(bool dod) {
	m_dod = dod;
}

void PCWriter::close() {
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

	m_writer->WritePoint(pt);
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
			Util::rm(f);
	}
}

int even(int num) {
	if(num % 2 == 1)
		++num;
	return num;
}

bool intersects(double* a, double* b) {
	double aa[4] {g_min(a[0], a[2]), g_min(a[1], a[3]), g_max(a[0], a[2]), g_max(a[1], a[3])};
	double bb[4] {g_min(b[0], b[2]), g_min(b[1], b[3]), g_max(b[0], b[2]), g_max(b[1], b[3])};
	return !(aa[2] <= bb[0] || aa[0] >= bb[2] || aa[3] <= bb[1] || aa[1] >= bb[3]);
}


Tiler::Tiler(const std::vector<std::string> filenames) {
	for(const std::string& filename : filenames)
		files.emplace_back(filename);
}

void Tiler::tile(const std::string& outdir, double size, double buffer, int srid,
	double easting, double northing, int maxFileHandles) {

	std::cerr << std::setprecision(12);

	if(buffer < 0)
		g_runerr("Negative buffer is not allowed. Use easting, northing and tile size to crop tiles.");

	// Calculate the overall bounds of the file set.
	double allBounds[6] = {MAX_DBL, MAX_DBL, MIN_DBL, MIN_DBL, MAX_DBL, MIN_DBL};
	double fBounds[6];
	for(PCFile& f : files) {
		f.init();
		f.fileBounds(fBounds);
		std::cerr << "file bounds " << fBounds[0] << ", " << fBounds[1] << "; " << fBounds[2] << ", " << fBounds[3] << "\n";
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

		std::cerr << "cols " << cols0 << "; rows " << rows0 << "; size " << size0 << "\n";
		std::cerr << "bounds " << allBounds[0] << ", " << allBounds[1] << "; " << allBounds[2] << ", " << allBounds[3] << "\n";

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
							std::string outfile = Util::pathJoin(outdir, ss.str());
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
						std::string outfile = Util::pathJoin(outdir, ss.str());
						std::cerr << outfile << "\n";
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


geo::pc::Point::Point(double x, double y, double z, double intensity, double angle, int cls, int returnNum, int numReturns, bool isEdge) :
	m_x(x),
	m_y(y),
	m_z(z),
	m_intensity(intensity),
	m_angle(angle),
	m_cls(cls),
	m_returnNum(returnNum),
	m_numReturns(numReturns),
	m_isEdge(isEdge),
	m_point(nullptr) {
}

geo::pc::Point::Point(const liblas::Point& pt) :
	geo::pc::Point(pt.GetX(), pt.GetY(), pt.GetZ(),
			pt.GetIntensity(),
			pt.GetScanAngleRank(),
			pt.GetClassification().GetClass(),
			pt.GetReturnNumber(),
			pt.GetNumberOfReturns(),
			pt.GetFlightLineEdge()) {

	m_point = new liblas::Point(pt);
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
	m_returnNum(0),
	m_numReturns(0),
	m_isEdge(false),
	m_point(nullptr) {
}

int geo::pc::Point::classId() const {
	return m_cls;
}

double geo::pc::Point::operator[](int idx) const {
	switch(idx % 2) {
	case 0:
		return m_x;
	case 1:
		return m_y;
	}
	return 0;
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

double geo::pc::Point::value() const {
	return m_z;
}

double geo::pc::Point::intensity() const {
	return m_intensity;
}

double geo::pc::Point::scanAngle() const {
	return m_angle;
}

bool geo::pc::Point::isEdge() const {
	return m_isEdge;
}

bool geo::pc::Point::isLast() const {
	return m_returnNum == m_numReturns;
}

bool geo::pc::Point::isFirst() const {
	return m_returnNum == 1;
}

int geo::pc::Point::returnNum() const {
	return m_returnNum;
}

int geo::pc::Point::numReturns() const {
	return m_numReturns;
}

geo::pc::Point::~Point() {
	if(m_point)
		delete m_point;
}

/**
 * Struct used to store the line headers in the point cloud cache file.
 */
typedef struct {
	size_t idx;			// The index of the cell: row * cols + col
	size_t nextLine;	// The address in the mapped memory of the next row containing data for this cell.
	size_t count;		// The number of points in the cell.
} MappedLine;

/**
 * Struct used to store point in the point cache file.
 */
typedef struct {
	float x;
	float y;
	float z;
	float intensity;
	float angle;
	int cls;
	int retNum;
	int numRets;
	int edge;
} MappedPoint;

class MemGrid {
private:
	geo::util::MappedFile m_mapped;
	size_t m_lineCount;
	size_t m_cellCount;
	size_t m_lineLength;
	size_t m_totalLength;
	std::string m_mapFile;

	size_t m_currentLine;						// The next available line index.
	std::list<size_t> m_finalized; 				// The line offsets of finalized cells available for re-use.
	std::unordered_map<size_t, size_t> m_map;   // Maps cell indices to line indices.

	void resize(size_t size) {
		std::cerr << "resize " << size << "\n";
		m_mapped.reset(size);
		m_totalLength = size;
	}

	// row structure: idx (size_t -> 8 bytes) | next line (size_t -> 8 bytes) | count (size_t -> 8 bytes) | points (3 * double * n)
	// row length = 8 + 8 + 8 + (3 * 8) * n; n = m_lineCount
	// if the count >= line count, the next row must be set.

	void writePoint(size_t idx, float x, float y, float z, float intensity, float angle, int cls, int retNum, int numRets, int edge) {

		char* data = (char*) m_mapped.data();
		size_t offset;
		bool newLine;

		if(m_map.find(idx) != m_map.end()) {
			offset = m_map[idx];
			newLine = false;
		} else {
			if(!m_finalized.empty()) {
				offset = m_finalized.front();
				m_finalized.pop_front();
			} else {
				offset = m_currentLine++;
			}
			m_map[idx] = offset;
			newLine = true;
		}
		offset *= m_lineLength;

		if(offset + m_lineLength >= m_totalLength) {
			resize(m_totalLength + std::max(m_lineLength, (size_t) (m_totalLength * 0.25)));
			data = (char*) m_mapped.data();
		}

		MappedPoint mp = {x, y, z, intensity, angle, cls, retNum, numRets, edge};
		MappedLine ml = {0, 0, 0};

		if(newLine) {

			// There is no record in the points file; start one.
			ml.idx = idx;
			ml.count = 1;
			ml.nextLine = 0;
			std::memcpy(data + offset, &ml, sizeof(MappedLine));
			std::memcpy(data + offset + sizeof(MappedLine), &mp, sizeof(MappedPoint));

		} else {

			// Retrieve the line header at the offset.
			std::memcpy(&ml, data + offset, sizeof(MappedLine));

			// Find the last line for this cell.
			while(ml.nextLine) {
				offset = ml.nextLine * m_lineLength;
				std::memcpy(&ml, data + offset, sizeof(MappedLine));
			}

			// If the line is full, update the nextLine and grab a finalized one or start a new one.
			if(ml.count == m_lineCount) {

				if(!m_finalized.empty()) {
					ml.nextLine = m_finalized.front();
					m_finalized.pop_front();
				} else {
					ml.nextLine = m_currentLine++;
				}

				// If the new offset is too far, resize the memory.
				if(ml.nextLine * m_lineLength + m_lineLength >= m_totalLength) {
					resize(m_totalLength + std::max(m_lineLength, (size_t) (m_totalLength * 0.25)));
					data = (char*) m_mapped.data();
				}

				// Copy the line header to the previous line
				std::memcpy(data + offset, &ml, sizeof(MappedLine));

				// Create the new line.
				offset = ml.nextLine * m_lineLength;
				ml.nextLine = 0;
				ml.count = 1;
				std::memcpy(data + offset, &ml, sizeof(MappedLine));
				std::memcpy(data + offset + sizeof(MappedLine), &mp, sizeof(MappedPoint));

			} else {

				// Add the point toe the line.
				size_t count = ml.count++;
				std::memcpy(data + offset, &ml, sizeof(MappedLine));
				offset += sizeof(MappedLine) + count * sizeof(MappedPoint);
				std::memcpy(data + offset, &mp, sizeof(MappedPoint));

			}
		}
	}

	void finalize(size_t idx) {
		size_t offset = m_map[idx];
		m_finalized.push_back(offset);
		char* data = (char*) m_mapped.data();
		MappedLine ml = {0, 0, 0};
		std::memcpy(&ml, data + offset * m_lineLength, sizeof(MappedLine));
		while(ml.nextLine) {
			m_finalized.push_back(ml.nextLine);
			std::memcpy(&ml, data + ml.nextLine * m_lineLength, sizeof(MappedLine));
		}
		m_map.erase(idx);
	}

	size_t readPoints(size_t idx, std::vector<geo::pc::Point>& pts) {
		if(!hasUnread(idx))
			return 0; //g_runerr("No unread pixel for that index.");
		char* data = (char*) m_mapped.data();
		size_t offset = m_map[idx] * m_lineLength;
		size_t count = 0;
		MappedLine ml;
		MappedPoint mp;
		do {
			char* buf = data + offset;
			std::memcpy(&ml, buf, sizeof(MappedLine));
			if(ml.idx != idx)
				break;
			buf += sizeof(MappedLine);
			for(size_t i = 0; i < ml.count; ++i) {
				std::memcpy(&mp, buf, sizeof(MappedPoint));
				pts.emplace_back(mp.x, mp.y, mp.z, mp.intensity, mp.angle, mp.cls, mp.retNum, mp.numRets, mp.edge);
				buf += sizeof(MappedPoint);
				++count;
			}
			offset = ml.nextLine * m_lineLength;
		} while(ml.nextLine);
		finalize(idx);
		return count;
	}

public:

	MemGrid() :
		m_lineCount(0),
		m_cellCount(0),
		m_lineLength(0),
		m_totalLength(0),
		m_currentLine(0) {
	}

	MemGrid(size_t cellCount, size_t lineCount) :
		MemGrid() {
		init(cellCount, lineCount);
	}

	/**
	 * If the line count is zero, will calculate the line count to fit in a page.
	 */
	void init(const std::string& mapFile, size_t cellCount, size_t lineCount) {
		m_lineCount = lineCount == 0 ? (sysconf(_SC_PAGESIZE) - sizeof(MappedLine)) / sizeof(MappedPoint) : lineCount;
		m_cellCount = cellCount;
		m_lineLength = sizeof(MappedLine) + sizeof(MappedPoint) * m_lineCount;
		m_totalLength = cellCount * m_lineLength;
		m_mapFile = mapFile;

		std::cerr << "MemGrid: cells: " << cellCount << "; line count: " << lineCount << "; total size: " << m_totalLength << "; line length: " << m_lineLength << "\n";

		if(mapFile.empty()) {
			m_mapped.init(m_totalLength, true);
		} else {
			m_mapped.init(mapFile, m_totalLength, true);
		}
	}

	void init(size_t cellCount, size_t lineCount) {
		init("", cellCount, lineCount);
	}

	bool hasUnread(size_t idx) {
		return m_map.find(idx) != m_map.end();
	}

	const std::unordered_map<size_t, size_t>& indexMap() const {
		return m_map;
	}

	void add(size_t idx, double x, double y, double z, double intensity, double angle, int cls, int retNum, int numRets, int edge) {
		writePoint(idx, x, y, z, intensity, angle, cls, retNum, numRets, edge);
	}

	size_t get(size_t idx, std::vector<geo::pc::Point>& out) {
		return readPoints(idx, out);
	}

	void flush() {
		m_mapped.flush();
	}

	~MemGrid() {
		if(m_mapFile.empty())
			Util::rm(m_mapped.name());
	}
};

const std::unordered_map<std::string, std::string> computerNames = {
		{"min", "The minimum value"},
		{"max", "The maximum value"},
		{"percentile-5", "The 5th percentile"},
		{"decile-1", "The 1st decile"},
		{"decile-2", "The 2nd decile"},
		{"quartile-1", "The 1st quartile"},
		{"decile-3", "The 3rd decile"},
		{"decile-4", "The 4th decile"},
		{"decile-5", "The 4th decile"},
		{"quartile-2", "The 2nd quartile"},
		{"median", "The median value"},
		{"decile-6", "The 6th decile"},
		{"decile-7", "The 7th decile"},
		{"quantile-3", "The 3rd quantile"},
		{"decile-8", "The 8th decile"},
		{"decile-9", "The 9th decile"},
		{"percentile-95", "The 95th percentile"},
		{"mean", "The mean value"},
		{"variance", "The variance with n-1"},
		{"std-dev", "The standard deviation with n-1"},
		{"rugosity-acr", "The arc-chord rugosity (DuPreez, 2004)"},
		{"idw-2", "Inverse distance weighting; coefficient 2"},
		{"hlrg-bio", "HLRG biometrics set"}
};

Computer* getComputer(const std::string& name) {
	if(name == "min") { 						return new MinComputer();
	} else if(name == "min") { 					return new MaxComputer();
	} else if(name == "percentile-5") { 		return new PercentileComputer(0.05);
	} else if(name == "decile-1") { 			return new PercentileComputer(0.1);
	} else if(name == "decile-2") { 			return new PercentileComputer(0.2);
	} else if(name == "quartile-1") { 			return new PercentileComputer(0.25);
	} else if(name == "decile-3") { 			return new PercentileComputer(0.3);
	} else if(name == "decile-4") { 			return new PercentileComputer(0.4);
	} else if(name == "decile-5") { 			return new PercentileComputer(0.5);
	} else if(name == "quartile-2") { 			return new PercentileComputer(0.5);
	} else if(name == "median") { 				return new PercentileComputer(0.5);
	} else if(name == "decile-6") { 			return new PercentileComputer(0.6);
	} else if(name == "decile-7") { 			return new PercentileComputer(0.7);
	} else if(name == "quartile-3") { 			return new PercentileComputer(0.75);
	} else if(name == "decile-8") { 			return new PercentileComputer(0.8);
	} else if(name == "decile-9") { 			return new PercentileComputer(0.9);
	} else if(name == "percentile-95") { 		return new PercentileComputer(0.95);
	} else if(name == "mean") { 				return new MeanComputer();
	} else if(name == "variance") { 			return new VarianceComputer();
	} else if(name == "std-dev") { 				return new StdDevComputer();
	} else if(name == "rugosity-acr") { 		return new RugosityComputer();
	} else if(name == "idw-2") {				return new IDWComputer();
	} else if(name == "hlrg-bio") {				return new HLRGBiometricsComputer(20, 75, 2.0);
	}
	g_runerr("Unknown computer name (" << name << ")");
}

Rasterizer::Rasterizer(const std::vector<std::string> filenames) :
	m_filter(nullptr) {
	for(const std::string& filename : filenames)
		m_files.emplace_back(filename);
}

const std::unordered_map<std::string, std::string>& Rasterizer::availableComputers() {
	return computerNames;
}

double Rasterizer::density(double resolution, double radius) {
	size_t count = 0;
	double w, h, cells, sum = 0;
	double fBounds[6];
	for(PCFile& f: m_files) {
		f.init();
		f.fileBounds(fBounds);
		w = fBounds[2] - fBounds[0];
		h = fBounds[3] - fBounds[1];
		cells = (w * h) / (resolution * resolution);
		if(cells > 0) {
			sum += f.pointCount() / cells;
			++count;
		}
	}
	std::cerr << sum << ", " << count << "\n";
	double cell = radius > 0 ? (M_PI * radius * radius) / (resolution * resolution) : 1;
	return (sum / count) * 1.5 * cell;
}

void Rasterizer::rasterize(const std::string& filename, const std::vector<std::string>& _types,
		double resX, double resY, double easting, double northing, double radius, int srid, int density, double ext) {

	if(std::isnan(resX) || std::isnan(resY))
		g_runerr("Resolution not valid.");

	std::vector<std::string> types(_types);
	if(types.empty()) {
		g_warn("No methods given; defaulting to mean.");
		types.push_back("mean");
	}

	std::vector<std::unique_ptr<Computer> > computers;
	for(const std::string& name : types)
		computers.emplace_back(getComputer(name));

	double bounds[4] = {MAX_DBL, MAX_DBL, MIN_DBL, MIN_DBL};
	{
		double fBounds[6];
		for(PCFile& f: m_files) {
			f.init();
			f.fileBounds(fBounds);
			if(fBounds[0] < bounds[0]) bounds[0] = fBounds[0];
			if(fBounds[1] < bounds[1]) bounds[1] = fBounds[1];
			if(fBounds[2] > bounds[2]) bounds[2] = fBounds[2];
			if(fBounds[3] > bounds[3]) bounds[3] = fBounds[3];
		}
	}

	if(std::isnan(easting))
		easting = ((int) ((resX > 0 ? bounds[0] : bounds[2]) / resX)) * resX;
	if(std::isnan(northing))
		northing = ((int) ((resY > 0 ? bounds[1] : bounds[3]) / resY)) * resY;

	bounds[resX > 0 ? 0 : 2] = easting;
	bounds[resY > 0 ? 1 : 3] = northing;

	if(radius < 0)
		radius = std::sqrt(std::pow(resX / 2, 2) * 2);

	int cols = (int) ((bounds[2] - bounds[0]) / std::abs(resX)) + 1;
	int rows = (int) ((bounds[3] - bounds[1]) / std::abs(resY)) + 1;

	std::cerr << "cols " << cols << "; rows " << rows << "\n";
	std::cerr << "bounds " << bounds[0] << ", " << bounds[1] << ", " << bounds[2] << ", " << bounds[3] << "\n";

	int bandCount = 0;
	for(const std::unique_ptr<Computer>& comp : computers)
		bandCount += comp->bandCount();

	GridProps props;
	props.setTrans(easting, resX, northing, resY);
	props.setSize(cols, rows);
	props.setNoData(NODATA);
	props.setDataType(DataType::Float32);
	props.setSrid(srid);
	props.setWritable(true);
	props.setBands(bandCount + 1);
	Raster rast(filename, props);
	rast.fillFloat(0, 1);
	rast.fillFloat(NODATA, 2);

	liblas::ReaderFactory fact;
	MemGrid grid;

	// The squared radius for comparison.
	double rad0 = radius * radius;
	// The radius of the "box" of pixels to check for a claim on the current point.
	int radpx = (int) std::ceil(radius / std::abs(props.resolutionX()));
	int i = 0;

	// Initialize the grid with some starting slots.
	grid.init(100000, 0);

	// As we go through the m_files list, we'll remove pointers from
	// this list and use it to calculate bounds for finalizing cells.
	std::unordered_set<PCFile*> files;
	for(PCFile& file : m_files)
		files.insert(&file);

	std::vector<geo::pc::Point> values;
	std::vector<double> out;

	for(PCFile& file : m_files) {
		std::cerr << "Reading file " << i++ << " of " << m_files.size() << ".\n";
		for(const std::string& filename : file.filenames()) {
			std::ifstream str(filename);
			liblas::Reader rdr = fact.CreateWithStream(str);
			while(rdr.ReadNextPoint()) {
				const geo::pc::Point pt(rdr.GetPoint());
				double x = pt.x();
				double y = pt.y();
				double z = pt.z();
				int col = props.toCol(x);
				int row = props.toRow(y);
				if(radius == 0) {
					grid.add(row * cols + col, x, y, z, pt.intensity(), pt.scanAngle(), pt.classId(), pt.returnNum(), pt.numReturns(), pt.isEdge());
				} else {
					for(int r = g_max(0, row - radpx); r < g_min(rows, row + radpx + 1); ++r) {
						for(int c = g_max(0, col - radpx); c < g_min(cols, col + radpx + 1); ++c) {
							double cx = props.toCentroidX(c);
							double cy = props.toCentroidY(r);
							double dist = std::pow(cx - x, 2.0) + std::pow(cy - y, 2.0);
							if(dist <= rad0)
								grid.add(r * cols + c, x, y, z, pt.intensity(), pt.scanAngle(), pt.classId(), pt.returnNum(), pt.numReturns(), pt.isEdge());
						}
					}
				}
			}
		}

		// Get the bounds of the finished file and then remove it
		// from the set.
		double fbounds[4];
		file.fileBounds(fbounds);
		files.erase(&file);

		int r0 = props.toRow(fbounds[1]);
		int r1 = props.toRow(fbounds[3]);
		int c0 = props.toCol(fbounds[0]);
		int c1 = props.toCol(fbounds[2]);
		double tbounds[6], cbounds[4];
		double rX = radius > 0 ? radius : std::abs(resX) * 0.5;
		double rY = radius > 0 ? radius : std::abs(resY) * 0.5;
		for(int r = g_max(0, g_min(r0, r1) - 1); r < g_min(rows, g_max(r0, r1) + 2); ++r) {
			for(int c = g_max(0, g_min(c0, c1) - 1); c < g_min(cols, g_max(c0, c1) + 2); ++c) {
				bool final = true;
				if(!files.empty()) {
					double x = props.toCentroidX(c);
					double y = props.toCentroidY(r);
					cbounds[0] = x - rX;
					cbounds[1] = y - rY;
					cbounds[2] = x + rX;
					cbounds[3] = y + rY;
					for(PCFile* f : files) {
						f->fileBounds(tbounds);
						if(intersects(tbounds, cbounds)) {
							final = false;
							break;
						}
					}
				}
				if(final) {
					size_t count = grid.get(r * cols + c, values);
					rast.setFloat((int) c, (int) r, count, 1);
					int band = 2;
					if(count) {
						double x = props.toCentroidX(c);
						double y = props.toCentroidY(r);
						for(size_t i = 0; i < computers.size(); ++i) {
							computers[i]->compute(x, y, values, radius, out, m_filter);
							for(double val : out)
								rast.setFloat((int) c, (int) r, std::isnan(val) ? NODATA : val, band++);
							out.clear();
						}
					} else {
						for(size_t i = 0; i < computers.size(); ++i) {
							for(int j = 0; j < computers[i]->bandCount(); ++j)
								rast.setFloat((int) c, (int) r, NODATA, band++);
						}
					}
					values.clear();
				}
			}
		}
	}

	std::unordered_map<size_t, size_t> mp(grid.indexMap()); // Copy to avoid invalidating the iterator.
	for(const auto& item : mp) {
		int r = item.first / cols;
		int c = item.first % cols;
		size_t count = grid.get(r * cols + c, values);
		rast.setFloat(c, r, count, 1);
		int band = 2;
		if(count) {
			double x = props.toCentroidX(c);
			double y = props.toCentroidY(r);
			for(size_t i = 0; i < computers.size(); ++i) {
				computers[i]->compute(x, y, values, radius, out, m_filter);
				for(double val : out)
					rast.setFloat(c, r, std::isnan(val) ? NODATA : val, band++);
				out.clear();
			}
		} else {
			for(size_t i = 0; i < computers.size(); ++i) {
				for(int j = 0; j < computers[i]->bandCount(); ++j)
					rast.setFloat(c, r, NODATA, band++);
			}
		}
		values.clear();
	}

}

void Rasterizer::setFilter(const PCPointFilter& filter) {
	if(m_filter)
		delete m_filter;
	m_filter = new PCPointFilter(filter);
}

Rasterizer::~Rasterizer() {
	if(m_filter)
		delete m_filter;
}


Normalizer::Normalizer(const std::vector<std::string> filenames) :
		m_filenames(filenames),
		m_filter(nullptr) {
}

void Normalizer::setFilter(const PCPointFilter& filter) {
	m_filter = new PCPointFilter(filter);
}

Normalizer::~Normalizer() {
	if(m_filter)
		delete m_filter;
}

void Normalizer::normalize(const std::string& dtmpath, const std::string& outdir, int band) {

	Raster dtm(dtmpath);
	const GridProps& props = dtm.props();
	double nodata = props.nodata();

	liblas::ReaderFactory fact;

	for(const std::string& filename : m_filenames) {

		std::string outfile = Util::pathJoin(outdir, Util::basename(filename) + ".las");

		std::ifstream str(filename);
		liblas::Reader rdr = fact.CreateWithStream(str);

		std::cerr << outfile << "\n";

		std::ofstream ostr(outfile, std::ios::binary | std::ios::trunc | std::ios::out);
		liblas::Header hdr(rdr.GetHeader());
		liblas::Writer wtr(ostr, hdr);

		double minZ = DBL_MAX;
		double maxZ = DBL_MIN;
		double minX = DBL_MAX;
		double maxX = DBL_MIN;
		double minY = DBL_MAX;
		double maxY = DBL_MIN;

		while(rdr.ReadNextPoint()) {

			const liblas::Point& pt = rdr.GetPoint();

			if(m_filter && !m_filter->keep(pt))
				continue;

			double x = pt.GetX();
			double y = pt.GetY();

			int col = props.toCol(x);
			int row = props.toRow(y);

			if(col < 0 || col >= props.cols() || row < 0 || row >= props.rows()) {
				std::cerr << x << ", " << y << "; " << props.toCol(x) << ", " << props.toRow(y) << "; " << props.cols() << "; " << props.rows() << "\n";
				continue;
			}

			double t = dtm.getFloat(props.toCol(x), props.toRow(y), band);

			if(t == nodata)
				continue;

			double z = pt.GetZ();
			double z0 = z - t;

			if(z0 < 0) z0 = 0;

			if(z0 < minZ) minZ = z0;
			if(z0 > maxZ) maxZ = z0;
			if(x < minX) minX = x;
			if(x > maxX) maxX = x;
			if(y < minY) minY = y;
			if(y > maxY) maxY = y;

			liblas::Point npt(pt);
			npt.SetZ(z0);
			wtr.WritePoint(npt);
		}

		hdr.SetMin(minX, minY, minZ);
		hdr.SetMax(maxX, maxY, maxZ);
		wtr.SetHeader(hdr);

	}

}
