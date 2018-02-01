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
	m_bufferedBounds{x - buffer, y - buffer, x + size + buffer, y + size + buffer} {

	m_filenames.push_back(filename);
}

PCFile::PCFile(const std::vector<std::string>& filenames, double x, double y, double size, double buffer) :
	m_x(x), m_y(y),
	m_fileBounds{MAX_DBL, MAX_DBL, MIN_DBL, MIN_DBL, MAX_DBL, MIN_DBL},
	m_bounds{x, y, x + size, y + size},
	m_bufferedBounds{x - buffer, y - buffer, x + size + buffer, y + size + buffer},
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
			}
		}
	}
}

bool PCFile::contains(double x, double y) const {
	return x >= m_bounds[0] && x < m_bounds[2]  && y >= m_bounds[1] && y < m_bounds[3];
}

bool PCFile::containsBuffered(double x, double y) const {
	return x >= m_bufferedBounds[0] && x < m_bufferedBounds[2] && y >= m_bufferedBounds[1] && y < m_bufferedBounds[3];
}

PCFile::~PCFile() {}


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
	return !(a[2] <= b[0] || a[0] >= b[2] || a[3] <= b[1] || a[1] >= b[3]);
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


geo::pc::Point::Point(const liblas::Point& pt) :
	m_x(pt.GetX()),
	m_y(pt.GetY()),
	m_z(pt.GetZ()),
	m_point(new liblas::Point(pt)) {
}

geo::pc::Point::Point(double x, double y, double z) :
	m_x(x),
	m_y(y),
	m_z(z),
	m_point(nullptr) {
}

geo::pc::Point::Point() :
	m_x(0),
	m_y(0),
	m_z(0),
	m_point(nullptr) {
}

int geo::pc::Point::classId() const {
	if(m_point)
		return m_point->GetClassification().GetClass();
	return 0;
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
	return m_point ? m_point->GetIntensity() : 0;
}

double geo::pc::Point::scanAngle() const {
	return m_point ? m_point->GetScanAngleRank() : 0;
}

bool geo::pc::Point::isEdge() const {
	return m_point ? m_point->GetFlightLineEdge() : false;
}

bool geo::pc::Point::isLast() const {
	return m_point ? m_point->GetReturnNumber() == m_point->GetNumberOfReturns() : false;
}

bool geo::pc::Point::isFirst() const {
	return m_point ? m_point->GetReturnNumber() == 1 : false;
}

int geo::pc::Point::returnNum() const {
	return m_point ? m_point->GetReturnNumber() : 0;
}


geo::pc::Point::~Point() {
	if(m_point)
		delete m_point;
}

typedef struct {
	size_t idx;
	size_t nextLine;
	size_t count;
} MappedLine;

typedef struct {
	double x;
	double y;
	double z;
} MappedPoint;

class MemGrid {
private:
	geo::util::MappedFile m_mapped;
	size_t m_lineCount;
	size_t m_cellCount;
	size_t m_currentLine;
	size_t m_lineLength;
	size_t m_totalLength;
	std::string m_mapFile;

	void resize(size_t size) {
		std::cerr << "resize " << size << "\n";
		m_mapped.reset(size);
		m_totalLength = size;
	}

	// row structure: idx (size_t -> 8 bytes) | next line (size_t -> 8 bytes) | count (size_t -> 8 bytes) | points (3 * double * n)
	// row length = 8 + 8 + 8 + (3 * 8) * n; n = m_lineCount
	// if the count >= line count, the next row must be set.

	void writePoint(size_t idx, double x, double y, double z) {
		char* data = (char*) m_mapped.data();
		size_t offset = idx * m_lineLength;
		MappedPoint mp = {x, y, z};
		MappedLine ml;
		std::memcpy(&ml, data + offset, sizeof(MappedLine));
		if(ml.idx != idx) {
			// There is no record in the points file; start one.
			ml.idx = idx;
			ml.count = 1;
			ml.nextLine = 0;
			std::memcpy(data + offset, &ml, sizeof(MappedLine));
			std::memcpy(data + offset + sizeof(MappedLine), &mp, sizeof(MappedPoint));
		} else {
			// Get the current line start and try to add to it.
			while(ml.nextLine) {
				offset = ml.nextLine * m_lineLength;
				std::memcpy(&ml, data + offset, sizeof(MappedLine));
			}
			if(ml.count == m_lineCount) {
				// If the line is full, update the nextLine and start a new one.
				ml.nextLine = m_currentLine++;
				std::memcpy(data + offset, &ml, sizeof(MappedLine));
				if(ml.nextLine * m_lineLength >= m_totalLength) {
					resize(m_totalLength * 2);
					data = (char*) m_mapped.data();
				}
				offset = ml.nextLine * m_lineLength;
				ml.nextLine = 0;
				ml.count = 1;
				std::memcpy(data + offset, &ml, sizeof(MappedLine));
				std::memcpy(data + offset + sizeof(MappedLine), &mp, sizeof(MappedPoint));
			} else {
				size_t count = ml.count++;
				std::memcpy(data + offset, &ml, sizeof(MappedLine));
				offset += sizeof(MappedLine) + count * sizeof(MappedPoint);
				std::memcpy(data + offset, &mp, sizeof(MappedPoint));
			}
		}
	}

	size_t readPoints(size_t idx, std::vector<geo::pc::Point>& pts) {
		char* data = (char*) m_mapped.data();
		size_t offset = idx * m_lineLength;
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
				pts.emplace_back(mp.x, mp.y, mp.z);
				buf += sizeof(MappedPoint);
				++count;
			}
			offset = ml.nextLine * m_lineLength;
		} while(ml.nextLine);
		return count;
	}

public:

	MemGrid() :
		m_lineCount(0),
		m_cellCount(0),
		m_currentLine(0),
		m_lineLength(0),
		m_totalLength(0) {
	}

	MemGrid(size_t cellCount, size_t lineCount = 128) {
		init(cellCount, lineCount);
	}

	MemGrid(const std::string& mapFile, size_t cellCount, size_t lineCount = 128) {
		init(mapFile, cellCount, lineCount);
	}

	void init(const std::string& mapFile, size_t cellCount, size_t lineCount = 128) {
		m_lineCount = lineCount;
		m_cellCount = cellCount;
		m_currentLine = cellCount;
		m_lineLength = sizeof(size_t) * 3 + sizeof(double) * 3 * m_lineCount;
		m_totalLength = cellCount * m_lineLength;
		m_mapFile = mapFile;

		if(mapFile.empty()) {
			m_mapped.init(m_totalLength, true);
		} else {
			m_mapped.init(mapFile, m_totalLength, true);
		}
	}

	void init(size_t cellCount, size_t lineCount = 128) {
		m_cellCount = cellCount;
		m_lineCount = lineCount;
		m_currentLine = cellCount;
		m_lineLength = sizeof(size_t) * 3 + sizeof(double) * 3 * m_lineCount;
		m_totalLength = cellCount * m_lineLength;

		m_mapped.init(m_totalLength, true);
	}

	void add(size_t idx, double x, double y, double z) {
		writePoint(idx, x, y, z);
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
		{"min", "The maximum value"},
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

bool Rasterizer::filter(const geo::pc::Point& pt) const {
	if(m_filter) {
		return m_filter->keep(pt);
	} else {
		return true;
	}
}

void Rasterizer::rasterize(const std::string& filename, const std::vector<std::string>& types,
		double res,	double easting, double northing, double radius, int srid, int density, double ext,
		const std::string& mapFile) {

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

	if(easting <= 0)
		easting = ((int) (bounds[0] / res)) * res;
	if(northing <= 0)
		northing = ((int) (bounds[3] / res)) * res;

	int cols = (int) ((bounds[2] - bounds[0]) / res) + 1;
	int rows = (int) ((bounds[3] - bounds[1]) / res) + 1;

	std::cerr << "cols " << cols << "; rows " << rows << "\n";
	std::cerr << "bounds " << bounds[0] << ", " << bounds[1] << ", " << bounds[2] << ", " << bounds[3] << "\n";

	int bandCount = 0;
	for(const std::unique_ptr<Computer>& comp : computers)
		bandCount += comp->bandCount();

	GridProps props;
	props.setTrans(easting, res, northing, -res);
	props.setSize(cols, rows);
	props.setNoData(NODATA);
	props.setDataType(DataType::Float32);
	props.setSrid(srid);
	props.setWritable(true);
	props.setBands(bandCount + 1);
	Raster rast(filename, props);

	liblas::ReaderFactory fact;
	MemGrid grid;

	if(mapFile.empty()) {
		grid.init(cols * rows, density);

		// The squared radius for comparison.
		double rad0 = radius * radius;
		// The radius of the "box" of pixels to check for a claim on the current point.
		int radpx = (int) std::ceil(radius / std::abs(props.resolutionX()));
		double xOffset = props.resolutionX() * 0.5;
		double yOffset = props.resolutionY() * 0.5;
		int i = 0;
		for(PCFile& file : m_files) {
			std::cerr << "Reading file " << i++ << " of " << m_files.size() << ".\n";
			for(const std::string& filename : file.filenames()) {
				std::ifstream str(filename);
				liblas::Reader rdr = fact.CreateWithStream(str);
				while(rdr.ReadNextPoint()) {
					const geo::pc::Point pt(rdr.GetPoint());
					if(filter(pt)) {
						double x = pt.x();
						double y = pt.y();
						double z = pt.z();
						int col = props.toCol(x);
						int row = props.toRow(y);
						if(radius == 0) {
							grid.add(row * cols + col, x, y, z);
						} else {
							for(int r = row - radpx; r < row + radpx + 1; ++r) {
								if(r < 0 || r >= rows) continue;
								for(int c = col - radpx; c < col + radpx + 1; ++c) {
									if(c < 0 || c >= cols) continue;
									double cx = props.toX(c) + xOffset;
									double cy = props.toY(r) + yOffset;
									double dist = std::pow(cx - x, 2.0) + std::pow(cy - y, 2.0);
									if(dist <= rad0)
										grid.add(r * cols + c, x, y, z);
								}
							}
						}
					}
				}
			}
		}
	} else {
		grid.init(mapFile, cols * rows, density);
	}

	std::vector<geo::pc::Point> values;
	std::vector<double> out;
	for(size_t r = 0; r < (size_t) rows; ++r) {
		if(r % 100 == 0)
			std::cerr << "Row " << r << " of " << rows << "\n";
		for(size_t c = 0; c < (size_t) cols; ++c) {
			size_t count = grid.get(r * cols + c, values);
			rast.setFloat((int) c, (int) r, count, 1);
			if(count) {
				double x = props.toX(c) + props.resolutionX();
				double y = props.toY(r) + props.resolutionY();
				int band = 2;
				for(size_t i = 0; i < computers.size(); ++i) {
					computers[i]->compute(x, y, values, radius, out);
					for(double val : out)
						rast.setFloat((int) c, (int) r, val == NODATA ? NODATA : val, band++);
					out.clear();
				}
				values.clear();
			} else {
				//for(size_t i = 0; i < computers.size(); ++i)
				//	rast.setFloat((int) c, (int) r, NODATA, i + 2);
			}
		}
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

void Normalizer::normalize(const std::string& dtmpath, const std::string& outdir) {

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

			double t = dtm.getFloat(props.toCol(x), props.toRow(y));

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
