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

using namespace geo::raster;
using namespace geo::util;
using namespace geo::pc;

PCFile::PCFile(const std::string& filename, double x, double y) :
	m_x(x), m_y(y),
	m_bounds{99999999., 99999999., -999999999., -99999999., 99999999., -99999999.} {

	m_filenames.push_back(filename);
}

PCFile::PCFile(const std::vector<std::string>& filenames, double x, double y) :
	m_x(x), m_y(y),
	m_bounds{99999999., 99999999., -999999999., -99999999., 99999999., -99999999.},
	m_filenames(filenames) {
}

double PCFile::x() const {
	return m_x;
}

double PCFile::y() const {
	return m_y;
}

void PCFile::bounds(double* bounds) const {
	for(int i = 0; i < 6; ++i)
		bounds[i] = m_bounds[i];
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
			if(minx < m_bounds[0]) m_bounds[0] = minx;
			if(miny < m_bounds[1]) m_bounds[1] = miny;
			if(maxx > m_bounds[2]) m_bounds[2] = maxx;
			if(maxy > m_bounds[3]) m_bounds[3] = maxy;
			if(minz < m_bounds[4]) m_bounds[4] = minz;
			if(maxz > m_bounds[5]) m_bounds[5] = maxz;
		} else {
			while(reader.ReadNextPoint()) {
				const liblas::Point& pt = reader.GetPoint();
				double x = pt.GetX();
				double y = pt.GetY();
				double z = pt.GetZ();
				if(x < m_bounds[0]) m_bounds[0] = x;
				if(y < m_bounds[1]) m_bounds[1] = y;
				if(x > m_bounds[2]) m_bounds[2] = x;
				if(y > m_bounds[3]) m_bounds[3] = y;
				if(z < m_bounds[4]) m_bounds[4] = z;
				if(z > m_bounds[5]) m_bounds[5] = z;
			}
		}
	}
}

PCFile::~PCFile() {}


const long maxPoints = 20000000;

PCWriter::PCWriter(const std::string& filename, const liblas::Header& hdr, double x, double y) :
	m_fileIdx(0),
	m_returns(0), m_retNum{0,0,0,0,0},
	m_totalReturns(0),
	m_outBounds{99999999.0,99999999.0,-99999999.0,-99999999.0,99999999.0,-99999999.0},
	m_x(x), m_y(y),
	m_filename(filename),
	m_writer(nullptr), m_header(nullptr),
	m_dod(true) {

	m_header = new liblas::Header(hdr);

	open();
}

double PCWriter::x() const {
	return m_x;
}

double PCWriter::y() const {
	return m_y;
}

void PCWriter::bounds(double* bounds) const {
	for(int i = 0; i < 6; ++i)
		bounds[i] = m_outBounds[i];
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
	m_str.open(filename, std::ios::out | std::ios::binary);
	m_writer = new liblas::Writer(m_str, *m_header);
}

void PCWriter::deleteOnDestruct(bool dod) {
	m_dod = dod;
}

void PCWriter::close() {
	if(m_writer) {
		m_header->SetMin(m_outBounds[0], m_outBounds[2], m_outBounds[4]);
		m_header->SetMax(m_outBounds[1], m_outBounds[3], m_outBounds[5]);
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

geo::pc::Tile::Tile(double minx, double miny, double maxx, double maxy, double buffer) {
	m_bounds[0] = minx;
	m_bounds[1] = miny;
	m_bounds[2] = maxx;
	m_bounds[3] = maxy;
	m_bufferedBounds[0] = m_bounds[0] - buffer;
	m_bufferedBounds[1] = m_bounds[1] - buffer;
	m_bufferedBounds[2] = m_bounds[2] + buffer;
	m_bufferedBounds[3] = m_bounds[3] + buffer;
}

PCWriter* geo::pc::Tile::writer(bool release) {
	if(release) {
		return m_writer.release();
	} else {
		return m_writer.get();
	}
}

void geo::pc::Tile::writer(PCWriter* wtr) {
	m_writer.reset(wtr);
}

bool geo::pc::Tile::contains(double x, double y) {
	return x >= m_bounds[0] && x < m_bounds[2] 
		&& y >= m_bounds[1] && y < m_bounds[3];
}

bool geo::pc::Tile::containsBuffered(double x, double y) {
	return x >= m_bufferedBounds[0] && x < m_bufferedBounds[2] 
		&& y >= m_bufferedBounds[1] && y < m_bufferedBounds[3];
}


Tiler::Tiler(const std::vector<std::string> filenames) {
	for(const std::string& filename : filenames)
		files.emplace_back(filename);
}

void Tiler::tile(const std::string& outdir, double size, double buffer, int srid,
	double easting, double northing, int maxFileHandles) {

	std::cerr << std::setprecision(12);

	// Calculate the overall bounds of the file set.
	double allBounds[6] = {9999999999., 9999999999., -9999999999., -9999999999., 999999999., -999999999.};
	double fBounds[6];
	for(PCFile& f : files) {
		f.init();
		f.bounds(fBounds);
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
	if(std::isnan(easting))
		easting = ((int) (allBounds[0] / size)) * size - size;
	if(std::isnan(northing))
		northing = ((int) (allBounds[1] / size)) * size - size;

	// Reset the bounds using easting, northing and multiples of size.
	allBounds[0] = easting;
	allBounds[1] = northing;
	allBounds[2] = std::ceil(allBounds[2] / size) * size + size;
	allBounds[3] = std::ceil(allBounds[3] / size) * size + size;

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
	std::vector<std::unique_ptr<Tile> > tiles;

	do {

		std::cerr << "cols " << cols0 << "; rows " << rows0 << "; size " << size0 << "\n";
		std::cerr << "bounds " << allBounds[0] << ", " << allBounds[1] << "; " << allBounds[2] << ", " << allBounds[3] << "\n";

		if(!tiles.empty()) {
			// This is run n>0, so we now read from the intermediate tiles.

			// To prevent the previous set of tiles from getting deleted before
			// we're done reading them.
			std::vector<std::unique_ptr<PCWriter> > tmpWriters;

			// Clear the readers list and rebuild with files from the previous
			// level of tiles.
			files.clear();
			for(std::unique_ptr<Tile>& tile : tiles) {
				PCWriter* wtr = tile->writer(true);
				// Create and emplace a LASReader.
				files.emplace_back(wtr->filenames(), wtr->x(), wtr->y());
				// Close the PCWriter.
				wtr->close();
				// Add to tmpWriters to prevent immediate destruction.
				tmpWriters.emplace_back(wtr);
			}

			// Clear the tiles list to start rebuilding it with the next level.
			tiles.clear();

			for(PCFile& file : files) {
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
							std::unique_ptr<Tile> tile(new Tile(x, y, x + size0, y + size0, buffer));
							tile->writer(new PCWriter(outfile, ihdr, x, y));
							tiles.push_back(std::move(tile));
						}
					}
				}

				// Now run over the files and reat their points into the new tiles.
				for(const std::string& filename : file.filenames()) {
					std::ifstream istr(filename, std::ios::in | std::ios::binary);
					liblas::Reader irdr = rfact.CreateWithStream(istr);
					while(irdr.ReadNextPoint()) {
						const liblas::Point& pt = irdr.GetPoint();
						for(std::unique_ptr<Tile>& tile : tiles) {
							if(tile->containsBuffered(pt.GetX(), pt.GetY()))
								tile->writer()->addPoint(pt);
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
						std::unique_ptr<Tile> tile(new Tile(x, y, x + size0, y + size0, buffer));
						tile->writer(new PCWriter(outfile, ihdr, x, y));
						tiles.push_back(std::move(tile));
					}
				}
			}

			// Iterate over the files, filling the tiles.
			for(PCFile& file : files) {
				for(const std::string& filename : file.filenames()) {
					std::ifstream istr(filename, std::ios::in | std::ios::binary);
					liblas::Reader irdr = rfact.CreateWithStream(istr);
					while(irdr.ReadNextPoint()) {
						const liblas::Point& pt = irdr.GetPoint();
						for(std::unique_ptr<Tile>& tile : tiles) {
							if(tile->containsBuffered(pt.GetX(), pt.GetY()))
								tile->writer()->addPoint(pt);
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
	for(std::unique_ptr<Tile>& tile : tiles)
		tile->writer()->deleteOnDestruct(false);

	// Destroy the last batch of tiles.
	tiles.clear();
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

geo::pc::Point::~Point() {
}

class MemGrid {
private:
	geo::util::MappedFile m_mapped;
	size_t m_lineCount;
	size_t m_cellCount;
	size_t m_currentLine;
	size_t m_lineLength;
	size_t m_totalLength;
	std::unordered_map<size_t, std::vector<size_t> > m_lines;	///< The current position in the data for idx. May be spread across segments.
	std::unordered_map<size_t, size_t> m_positions;				///< The start indices of the lines for an index. Modulus to get the line index.

	void resize(size_t size) {
		m_mapped.reset(size);
		m_totalLength = size;
	}

	void writePoint(size_t idx, double x, double y, double z) {
		size_t line;
		size_t pos;
		if(m_positions.find(idx) == m_positions.end()) {
			// If there's no entry for this pixel, start one.
			pos = 0;
			line = m_currentLine++;
			m_lines[idx].push_back(line);
			m_positions[idx] = pos;
		} else {
			pos = m_positions[idx];
			size_t i = pos / m_lineCount;
			std::vector<size_t>& lines = m_lines[idx];
			if(i >= lines.size()) {
				// If the current position is beyond the end of the
				// last line, start a new one.
				line = m_currentLine++;
				lines.push_back(line);
			} else {
				// Otherwise use the line corresponding to the position.
				line = lines[i];
			}
		}
		// Skip to the line, then advance by the position.
		size_t offset = line * m_lineCount + pos % m_lineCount;
		// If the offset is beyond the end of the file, extend the buffer.
		if(offset * sizeof(double) >= m_totalLength)
			resize(m_totalLength * 2);
		// Get the address at offset.
		double* buf = ((double*) m_mapped.data()) + offset;
		*buf = x;
		*(buf + 1) = y;
		*(buf + 2) = z;
		m_positions[idx] = pos + 3;
	}

	size_t readPoints(size_t idx, std::vector<geo::pc::Point>& pts) {
		if(m_lines.find(idx) == m_lines.end())
			return 0;
		double* data = (double*) m_mapped.data();
		std::vector<size_t>& lines = m_lines[idx];
		size_t count = 0;
		size_t maxPos = m_positions[idx];
		size_t pos = 0;
		while(pos < maxPos) {
			size_t offset = lines[pos / m_lineCount] * m_lineCount;
			double* buf = data + offset;
			for(size_t p = 0;p < m_lineCount && (pos + p) < maxPos; p += 3) {
				pts.emplace_back(*buf, *(buf + 1), *(buf + 2));
				buf += 3;
				++count;
			}
			pos += m_lineCount;
		}
		return count;
	}

public:

	MemGrid(size_t cellCount, size_t lineCount = 512) :
		m_lineCount(lineCount * 3),
		m_cellCount(cellCount),
		m_currentLine(0),
		m_lineLength(m_lineCount * sizeof(double)),
		m_totalLength(cellCount * m_lineLength) {

		m_mapped.init(m_totalLength, true);
	}

	void add(int col, int row, double x, double y, double z) {
		size_t idx = ((size_t) row << 32) | col;
		writePoint(idx, x, y, z);
	}

	size_t get(size_t col, size_t row, std::vector<geo::pc::Point>& out) {
		size_t idx = (row << 32) | col;
		return readPoints(idx, out);
	}

	~MemGrid() {
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
		{"rugosity-acr", "The arc-chord rugosity (DuPreez, 2004)"}
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
	}
	g_runerr("Unknown computer name: " << name);
}

Rasterizer::Rasterizer(const std::vector<std::string> filenames) {
	for(const std::string& filename : filenames)
		m_files.emplace_back(filename);
}

const std::unordered_map<std::string, std::string>& Rasterizer::availableComputers() {
	return computerNames;
}

bool Rasterizer::filter(const geo::pc::Point& pt) const {
	return pt.classId() == 2;
}

void Rasterizer::rasterize(const std::string& filename, const std::vector<std::string>& types,
		double res,	double easting, double northing, double radius, int srid, int density, double ext) {

	std::vector<std::unique_ptr<Computer> > computers;
	for(const std::string& name : types)
		computers.emplace_back(getComputer(name));

	double bounds[4] = {9999999999, 9999999999, -9999999999, -9999999999};
	{
		double fBounds[6];
		for(PCFile& f: m_files) {
			f.init();
			f.bounds(fBounds);
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

	GridProps props;
	props.setTrans(easting, res, northing, -res);
	props.setSize(cols, rows);
	props.setNoData(-9999.0);
	props.setDataType(DataType::Float32);
	props.setSrid(srid);
	props.setWritable(true);
	props.setBands(computers.size() + 1);
	Raster rast(filename, props);

	liblas::ReaderFactory fact;
	MemGrid grid(cols * rows, density);

	// The squared radius for comparison.
	double rad0 = radius * radius;
	// The radius of the "box" of pixels to check for a claim on the current point.
	int radpx = (int) std::ceil(radius / std::abs(props.resolutionX()));
	int i = 0;
	for(PCFile& file : m_files) {
		std::cerr << "Reading file " << i++ << " of " << m_files.size() << ".\n";
		for(const std::string& filename : file.filenames()) {
			std::ifstream str(filename);
			liblas::Reader rdr = fact.CreateWithStream(str);
			while(rdr.ReadNextPoint()) {
				const geo::pc::Point pt(rdr.GetPoint());
				if(filter(pt)) { // TODO: Configurable.
					double x = pt.x();
					double y = pt.y();
					double z = pt.z();
					int col = props.toCol(x);
					int row = props.toRow(y);
					for(int r = row - radpx; r < row + radpx + 1; ++r) {
						for(int c = col - radpx; c < col + radpx + 1; ++c) {
							double cx = props.toX(c);
							double cy = props.toY(r);
							double dist = std::pow(cx - x, 2.0) + std::pow(cy - y, 2.0);
							if(dist <= rad0)
								grid.add(c, r, x, y, z);
						}
					}
				}
			}
		}
	}

	std::vector<geo::pc::Point> values;
	for(size_t r = 0; r < (size_t) rows; ++r) {
		if(r % 100 == 0)
			std::cerr << "Row " << r << " of " << rows << "\n";
		for(size_t c = 0; c < (size_t) cols; ++c) {
			size_t count = grid.get(c, r, values);
			rast.setFloat((int) c, (int) r, count, 1);
			if(count) {
				double x = props.toX(c) + props.resolutionX();
				double y = props.toY(r) + props.resolutionY();
				for(size_t i = 0; i < computers.size(); ++i) {
					double val = computers[i]->compute(x, y, values, radius);
					rast.setFloat((int) c, (int) r, val, i + 2);
				}
				values.clear();
			} else {
				for(size_t i = 0; i < computers.size(); ++i)
					rast.setFloat((int) c, (int) r, -9999.0, i + 2);
			}
		}
	}
}

Rasterizer::~Rasterizer() {
}
