/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

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

using namespace geo::util;

const bool useHeader = true;

class LASFile {
public:
	std::vector<std::string> filenames;
	double bounds[6];
	double x; // Min corner coords.
	double y;

	LASFile(const std::string& filename, double x = 0, double y = 0) :
		bounds{99999999., 99999999., -999999999., -99999999., 99999999., -99999999.},
		x(x), y(y) {
		filenames.push_back(filename);
	}

	LASFile(const std::vector<std::string>& filenames, double x = 0, double y = 0) :
		filenames(filenames),
		bounds{99999999., 99999999., -999999999., -99999999., 99999999., -99999999.},
		x(x), y(y) {
	}

	void init() {
		liblas::ReaderFactory f;
		for(const std::string& filename : filenames) {

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

				if(minx < bounds[0]) bounds[0] = minx;
				if(miny < bounds[1]) bounds[1] = miny;
				if(maxx > bounds[2]) bounds[2] = maxx;
				if(maxy > bounds[3]) bounds[3] = maxy;
				if(minz < bounds[4]) bounds[4] = minz;
				if(maxz > bounds[5]) bounds[5] = maxz;
			} else {
				while(reader.ReadNextPoint()) {
					const liblas::Point& pt = reader.GetPoint();
					double x = pt.GetX();
					double y = pt.GetY();
					double z = pt.GetZ();
					if(y < 6500000.) {
						std::cerr << "oob " << x << ", " << y << ", " << z << "\n";
					}
					if(x < bounds[0]) bounds[0] = x;
					if(y < bounds[1]) bounds[1] = y;
					if(x > bounds[2]) bounds[2] = x;
					if(y > bounds[3]) bounds[3] = y;
					if(z < bounds[4]) bounds[4] = z;
					if(z > bounds[5]) bounds[5] = z;
				}
			}
		}
	}

	~LASFile() {
	}
};

const long maxPoints = 20000000;

class LASWriter {
public:
	int fileIdx;
	int returns;
	int retNum[5];
	long totalReturns;
	double outBounds[6];
	double x;
	double y;
	std::string filename;
	std::vector<std::string> filenames;
	liblas::Writer* writer;
	liblas::Header* header;
	std::ofstream str;
	bool dod;

	LASWriter(const std::string& filename, const liblas::Header& hdr, double x = 0, double y = 0) :
		fileIdx(0),
		returns(0), retNum{0,0,0,0,0},
		totalReturns(0),
		outBounds{99999999.0,99999999.0,-99999999.0,-99999999.0,99999999.0,-99999999.0},
		x(x), y(y),
		filename(filename),
		writer(nullptr), header(nullptr),
		dod(true) {

		header = new liblas::Header(hdr);

		open();
	}

	std::string nextFile() {
		std::stringstream ss;
		ss << filename << "_" << ++fileIdx << ".las";
		std::string filename = ss.str();
		filenames.push_back(filename);
		return filename;
	}

	void open() {
		close();
		std::string filename = nextFile();
		str.open(filename, std::ios::out | std::ios::binary);
		writer = new liblas::Writer(str, *header);
	}

	void deleteOnDestruct(bool dod) {
		this->dod = dod;
	}

	void close() {
		if(writer) {
			header->SetMin(outBounds[0], outBounds[2], outBounds[4]);
			header->SetMax(outBounds[1], outBounds[3], outBounds[5]);
			for(int i = 0; i < 5; ++i)
				header->SetPointRecordsByReturnCount(i, retNum[i]);
			header->SetPointRecordsCount(returns);
			writer->SetHeader(*header);
			writer->WriteHeader();
			str.close();
			delete writer;
			writer = nullptr;

			returns = 0;
			for(int i = 0; i < 5; ++i)
				retNum[i] = 0;
		}
	}

	void addPoint(const liblas::Point& pt) {
		if(returns >= maxPoints || !writer)
			open();

		double x = pt.GetX();
		double y = pt.GetY();
		double z = pt.GetZ();

		++returns;
		++totalReturns;
		retNum[pt.GetReturnNumber() - 1]++;

		if(x < outBounds[0]) outBounds[0] = x;
		if(x > outBounds[2]) outBounds[2] = x;
		if(y < outBounds[1]) outBounds[1] = y;
		if(y > outBounds[3]) outBounds[3] = y;
		if(z < outBounds[4]) outBounds[4] = z;
		if(z > outBounds[5]) outBounds[5] = z;

		writer->WritePoint(pt);
	}

	~LASWriter() {
		close();
		delete header;
		if(!totalReturns || dod) {
			for(const std::string& f : filenames)
				Util::rm(f);
		}
	}
};

class Tile {
public:
	double bounds[4];
	double bufferedBounds[4];
	std::unique_ptr<LASWriter> writer;

	Tile(double minx, double miny, double maxx, double maxy, double buffer = 0) {
		bounds[0] = minx;
		bounds[1] = miny;
		bounds[2] = maxx;
		bounds[3] = maxy;
		bufferedBounds[0] = bounds[0] - buffer;
		bufferedBounds[1] = bounds[1] - buffer;
		bufferedBounds[2] = bounds[2] + buffer;
		bufferedBounds[3] = bounds[3] + buffer;
	}

	bool contains(double x, double y) {
		return x >= bounds[0] && x < bounds[2] && y >= bounds[1] && y < bounds[3];
	}

	bool containsBuffered(double x, double y) {
		return x >= bufferedBounds[0] && x < bufferedBounds[2] && y >= bufferedBounds[1] && y < bufferedBounds[3];
	}

};

int even(int num) {
	if(num % 2 == 1)
		++num;
	return num;
}

class Tiler {
public:
	std::vector<LASFile> files;

	Tiler(const std::vector<std::string> filenames) {
		for(const std::string& filename : filenames)
			files.emplace_back(filename);
	}

	void tile(const std::string& outdir, double size,
		double easting, double northing, double buffer, int srid, int maxFiles = 32) {

		std::cerr << std::setprecision(12);

		double allBounds[6] = {9999999999., 9999999999., -9999999999., -9999999999., 999999999., -999999999.};

		for(LASFile& f : files) {
			f.init();
			std::cerr << "file bounds " << f.bounds[0] << ", " << f.bounds[1] << "; " << f.bounds[2] << ", " << f.bounds[3] << "\n";
			if(f.bounds[0] < allBounds[0]) allBounds[0] = f.bounds[0];
			if(f.bounds[1] < allBounds[1]) allBounds[1] = f.bounds[1];
			if(f.bounds[2] > allBounds[2]) allBounds[2] = f.bounds[2];
			if(f.bounds[3] > allBounds[3]) allBounds[3] = f.bounds[3];
			if(f.bounds[4] < allBounds[4]) allBounds[4] = f.bounds[4];
			if(f.bounds[5] > allBounds[5]) allBounds[5] = f.bounds[5];
		}

		if(easting <= 0)
			easting = ((int) (allBounds[0] / size)) * size - size;
		if(northing <= 0)
			northing = ((int) (allBounds[1] / size)) * size - size;

		std::cerr << "bounds a " << allBounds[0] << ", " << allBounds[1] << "; " << allBounds[2] << ", " << allBounds[3] << "\n";
		allBounds[0] = easting;
		allBounds[1] = northing;
		allBounds[2] = std::ceil(allBounds[2] / size) * size + size;
		allBounds[3] = std::ceil(allBounds[3] / size) * size + size;
		std::cerr << "bounds b " << allBounds[0] << ", " << allBounds[1] << "; " << allBounds[2] << ", " << allBounds[3] << "\n";

		int cols = (int) (allBounds[2] - allBounds[0]) / size;
		int rows = (int) (allBounds[3] - allBounds[1]) / size;

		// Double the tile size until few enough writers are used.
		double size0 = size;
		int cols0 = cols;
		int rows0 = rows;
		while(cols0 * rows0 > maxFiles && cols0 > 1 && rows0 > 1) {
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
				// This is run n>0, so we delete the original las files and now read from the
				// intermediate tiles.

				std::vector<std::unique_ptr<LASWriter> > tmpWriters;

				files.clear();
				for(std::unique_ptr<Tile>& tile : tiles) {
					files.emplace_back(tile->writer->filenames, tile->writer->x, tile->writer->y);
					tile->writer->close();
					tmpWriters.push_back(std::move(tile->writer));
				}

				tiles.clear();

				for(LASFile& file : files) {
					{
						std::ifstream istr(file.filenames[0], std::ios::in | std::ios::binary);
						liblas::Reader irdr = rfact.CreateWithStream(istr);
						const liblas::Header& ihdr = irdr.GetHeader();

						// Create the writers for the intermediate tiles.
						for(int r = 0; r < 2; ++r) {
							for(int c = 0; c < 2; ++c) {
								std::stringstream ss;
								double x = file.x + c * size0;
								double y = file.y + r * size0;
								ss << "tile_" << (int) x << "_" << (int) y << "_" << size0;
								std::string outfile = Util::pathJoin(outdir, ss.str());
								std::unique_ptr<Tile> tile(new Tile(x, y, x + size0, y + size0, buffer));
								tile->writer.reset(new LASWriter(outfile, ihdr, x, y));
								tiles.push_back(std::move(tile));
							}
						}
					}

					for(const std::string& filename : file.filenames) {

						std::ifstream istr(filename, std::ios::in | std::ios::binary);
						liblas::Reader irdr = rfact.CreateWithStream(istr);

						while(irdr.ReadNextPoint()) {
							const liblas::Point& pt = irdr.GetPoint();
							for(std::unique_ptr<Tile>& tile : tiles) {
								if(tile->containsBuffered(pt.GetX(), pt.GetY()))
									tile->writer->addPoint(pt);
							}
						}
					}

				}

				tmpWriters.clear();

			} else {
				// This is the first run, so we read all the input files into the first set of
				// intermediate files.

				{
					// Get the header from the first file to use as a template.
					std::ifstream istr(files[0].filenames[0], std::ios::in | std::ios::binary);
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
							tile->writer.reset(new LASWriter(outfile, ihdr, x, y));
							tiles.push_back(std::move(tile));
						}
					}
				}

				// Iterate over the files, filling the tiles.
				for(LASFile& file : files) {

					for(const std::string& filename : file.filenames) {

						std::ifstream istr(filename, std::ios::in | std::ios::binary);
						liblas::Reader irdr = rfact.CreateWithStream(istr);

						while(irdr.ReadNextPoint()) {
							const liblas::Point& pt = irdr.GetPoint();
							for(std::unique_ptr<Tile>& tile : tiles) {
								if(tile->containsBuffered(pt.GetX(), pt.GetY()))
									tile->writer->addPoint(pt);
							}
						}
					}
				}
			}

			size0 *= 0.5;
			cols0 *= 2;
			rows0 *= 2;

		} while(size0 >= size);

		for(std::unique_ptr<Tile>& tile : tiles)
			tile->writer->deleteOnDestruct(false);

		tiles.clear();

	}

	~Tiler() {
	}
};

void usage() {
	std::cerr << "Usage: las2grid <output dir> <input las [*]>\n"
			<< " -t <size>       The tile size in map units (default 1000; square).\n"
			<< " -e <easting>    The top left corner horizontal alignment\n"
			<< "                 (defaults to nearest whole multiple of size).\n"
			<< " -n <northing>   The top left corner vertical alignment\n"
			<< "                 (defaults to nearest whole multiple of size).\n"
			<< " -s <srid>       The spatial reference ID. Default 0.\n"
			<< " -b <buffer>     Buffer for each tile. Default 0.\n";
}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	double size = 1000;
	double easting = 0;
	double northing = 0;
	double buffer = 0;
	uint16_t srid = 0;
	std::vector<std::string> args;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-t") {
			size = atof(argv[++i]);
		} else if(v == "-s") {
			srid = atoi(argv[++i]);
		} else if(v == "-b") {
			buffer = atof(argv[++i]);
		} else if(v == "-e") {
			easting = atof(argv[++i]);
		} else if(v == "-n") {
			northing = atof(argv[++i]);
		} else {
			args.push_back(argv[i]);
		}
	}

	if(size <= 0) {
		std::cerr << "Size must be >0.\n";
		usage();
		return 1;
	}

	if(args.size() < 2) {
		std::cerr << "Input and output filenames required.\n";
		usage();
		return 1;
	}

	std::vector<std::string> infiles(args.begin() + 1, args.end());
	Tiler r(infiles);
	r.tile(args[0], size, easting, northing, buffer, srid);

	return 0;
}


