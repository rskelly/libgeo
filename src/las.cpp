#include <fstream>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <memory>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdint>

#include <liblas/liblas.hpp>

#include "util.hpp"
#include "las.hpp"

using namespace geo::util;
using namespace geo::las;

LASFile::LASFile(const std::string& filename, double x, double y) :
	m_x(x), m_y(y),
	m_bounds{99999999., 99999999., -999999999., -99999999., 99999999., -99999999.} {

	m_filenames.push_back(filename);
}

LASFile::LASFile(const std::vector<std::string>& filenames, double x, double y) :
	m_x(x), m_y(y),
	m_bounds{99999999., 99999999., -999999999., -99999999., 99999999., -99999999.},
	m_filenames(filenames) {
}

double LASFile::x() const {
	return m_x;
}

double LASFile::y() const {
	return m_y;
}

void LASFile::bounds(double* bounds) const {
	for(int i = 0; i < 6; ++i)
		bounds[i] = m_bounds[i];
}

const std::vector<std::string>& LASFile::filenames() const {
	return m_filenames;
}

void LASFile::init(bool useHeader) {
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

LASFile::~LASFile() {}


const long maxPoints = 20000000;

LASWriter::LASWriter(const std::string& filename, const liblas::Header& hdr, double x, double y) :
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

double LASWriter::x() const {
	return m_x;
}

double LASWriter::y() const {
	return m_y;
}

void LASWriter::bounds(double* bounds) const {
	for(int i = 0; i < 6; ++i)
		bounds[i] = m_outBounds[i];
}

const std::vector<std::string>& LASWriter::filenames() const {
	return m_filenames;
}

std::string LASWriter::nextFile() {
	std::stringstream ss;
	ss << m_filename << "_" << ++m_fileIdx << ".las";
	std::string filename = ss.str();
	m_filenames.push_back(filename);
	return filename;
}

void LASWriter::open() {
	close();
	std::string filename = nextFile();
	m_str.open(filename, std::ios::out | std::ios::binary);
	m_writer = new liblas::Writer(m_str, *m_header);
}

void LASWriter::deleteOnDestruct(bool dod) {
	m_dod = dod;
}

void LASWriter::close() {
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

void LASWriter::addPoint(const liblas::Point& pt) {
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

LASWriter::~LASWriter() {
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

Tile::Tile(double minx, double miny, double maxx, double maxy, double buffer) {
	m_bounds[0] = minx;
	m_bounds[1] = miny;
	m_bounds[2] = maxx;
	m_bounds[3] = maxy;
	m_bufferedBounds[0] = m_bounds[0] - buffer;
	m_bufferedBounds[1] = m_bounds[1] - buffer;
	m_bufferedBounds[2] = m_bounds[2] + buffer;
	m_bufferedBounds[3] = m_bounds[3] + buffer;
}

LASWriter* Tile::writer(bool release) {
	if(release) {
		return m_writer.release();
	} else {
		return m_writer.get();
	}
}

void Tile::writer(LASWriter* wtr) {
	m_writer.reset(wtr);
}

bool Tile::contains(double x, double y) {
	return x >= m_bounds[0] && x < m_bounds[2] 
		&& y >= m_bounds[1] && y < m_bounds[3];
}

bool Tile::containsBuffered(double x, double y) {
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
	for(LASFile& f : files) {
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
			std::vector<std::unique_ptr<LASWriter> > tmpWriters;

			// Clear the readers list and rebuild with files from the previous
			// level of tiles.
			files.clear();
			for(std::unique_ptr<Tile>& tile : tiles) {
				LASWriter* wtr = tile->writer(true);
				// Create and emplace a LASReader.
				files.emplace_back(wtr->filenames(), wtr->x(), wtr->y());
				// Close the LASWriter.
				wtr->close();
				// Add to tmpWriters to prevent immediate destruction.
				tmpWriters.emplace_back(wtr);
			}

			// Clear the tiles list to start rebuilding it with the next level.
			tiles.clear();

			for(LASFile& file : files) {
				{
					// Get a liblas::Header from the first file to use as a template
					// for the LASWriters.
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
							tile->writer(new LASWriter(outfile, ihdr, x, y));
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
						tile->writer(new LASWriter(outfile, ihdr, x, y));
						tiles.push_back(std::move(tile));
					}
				}
			}

			// Iterate over the files, filling the tiles.
			for(LASFile& file : files) {
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