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

#include <liblas/liblas.hpp>

#include "util.hpp"

using namespace geo::util;

class LASFile {
public:
	std::string filename;
	double bounds[6];
	double gridBounds[6];
	double x; // Min corner coords.
	double y;

	LASFile(const std::string& filename, double x = 0, double y = 0) :
		filename(filename),
		x(x), y(y) {

		init();
	}

	void init() {
		std::ifstream str(filename, std::ios::in | std::ios::binary);
		liblas::ReaderFactory f;
		liblas::Reader reader = f.CreateWithStream(str);
		const liblas::Header& header = reader.GetHeader();
		bounds[0] = header.GetMinX();
		bounds[1] = header.GetMinY();
		bounds[2] = header.GetMaxX();
		bounds[3] = header.GetMaxY();
		bounds[4] = header.GetMinZ();
		bounds[5] = header.GetMaxZ();
		str.close();
	}

	~LASFile() {
	}
};


class LASWriter {
public:
	std::string filename;
	std::vector<std::string> filenames;
	int fileIdx;
	liblas::Writer* writer;
	liblas::Header* header;
	std::ofstream str;
	int returns;
	int retNum[5];
	double outBounds[6];
	long totalReturns;
	double x;
	double y;

	LASWriter(const std::string& filename, const liblas::Header& hdr, double x = 0, double y = 0) :
		filename(filename),
		fileIdx(0),
		writer(nullptr), header(nullptr),
		returns(0), retNum{0,0,0,0,0},
		outBounds{99999999.0,99999999.0,-99999999.0,-99999999.0,99999999.0,-99999999.0},
		totalReturns(0),
		x(x), y(y) {

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
		if(writer)
			close();
		std::string filename = nextFile();
		str.open(filename, std::ios::out | std::ios::binary);
		writer = new liblas::Writer(str, *header);
	}

	void close() {
		header->SetMin(outBounds[0], outBounds[2], outBounds[4]);
		header->SetMax(outBounds[1], outBounds[3], outBounds[5]);
		for(int i = 0; i < 5; ++i)
			header->SetPointRecordsByReturnCount(i, retNum[i]);
		header->SetPointRecordsCount(returns);
		writer->SetHeader(*header);
		str.close();
		delete writer;
		writer = nullptr;

		returns = 0;
		for(int i = 0; i < 5; ++i)
			retNum[i] = 0;
	}

	void addPoint(const liblas::Point& pt) {
		if(returns >= 20000000)
			close();
		if(!writer)
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
		double easting, double northing, int srid, int maxFiles = 32) {

		double allBounds[6] = {9999999999, 9999999999, -9999999999, -9999999999, 999999999, -999999999};

		for(const LASFile& f : files) {
			if(f.bounds[0] < allBounds[0]) allBounds[0] = f.bounds[0];
			if(f.bounds[1] < allBounds[1]) allBounds[1] = f.bounds[1];
			if(f.bounds[2] > allBounds[2]) allBounds[2] = f.bounds[2];
			if(f.bounds[3] > allBounds[3]) allBounds[3] = f.bounds[3];
			if(f.bounds[4] < allBounds[4]) allBounds[4] = f.bounds[4];
			if(f.bounds[5] > allBounds[5]) allBounds[5] = f.bounds[5];
		}

		if(easting <= 0)
			easting = ((int) (allBounds[0] / size)) * size;
		if(northing <= 0)
			northing = ((int) (allBounds[1] / size)) * size;

		allBounds[0] = easting;
		allBounds[1] = northing;
		allBounds[2] = easting + (((int) ((allBounds[2] - easting) / size)) + 1) * size;
		allBounds[3] = northing + (((int) ((allBounds[3] - northing) / size)) + 1) * size;

		int cols = even((int) (allBounds[2] - allBounds[0]) / size);
		int rows = even((int) (allBounds[3] - allBounds[1]) / size);


		// Double the tile size until few enough are used.
		double size0 = size;
		int cols0 = cols, rows0 = rows;
		while(cols0 * rows0 > maxFiles && cols0 > 1 && rows0 > 1) {
			size0 *= 2.0;
			cols0 = even(((int) (allBounds[2] - allBounds[0]) / size0) + 1);
			rows0 = even(((int) (allBounds[3] - allBounds[1]) / size0) + 1);
		}

		liblas::ReaderFactory rfact;
		std::unordered_map<uint64_t, std::unique_ptr<LASWriter> > writers;

		do {

			std::cerr << "cols " << cols0 << "; rows " << rows0 << "\n";

			if(!writers.empty()) {
				// This is run n>0, so we delete the original las files and now read from the
				// intermediate tiles.

				files.clear();
				for(auto& w: writers) {
					for(const std::string& filename : w.second->filenames)
						files.emplace_back(filename, w.second->x, w.second->y);
				}

				writers.clear();

				for(LASFile& file : files) {

					std::ifstream istr(file.filename, std::ios::in | std::ios::binary);
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
							std::cout << outfile << "\n";
							uint64_t idx = ((uint64_t) r << 32) | c;
							writers[idx].reset(new LASWriter(outfile, ihdr));
						}
					}

					while(irdr.ReadNextPoint()) {
						const liblas::Point& pt = irdr.GetPoint();
						int c = (int) ((pt.GetX() - file.x) / size0);
						int r = (int) ((pt.GetY() - file.y) / size0);
						uint64_t idx = ((uint64_t) r << 32) | c;
						if(c < 0 || c > 1 || r < 0 || r > 1) {
							std::cerr << "warning: point out of bounds: " << pt.GetX() << ", " << pt.GetY() << "; " << file.x << ", " << file.y << ", " << file.filename << "; " << c << ", " << r << "\n";
							if(c < 0) c = 0;
							if(c > 1) c = 1;
							if(r < 0) r = 0;
							if(r > 1) r = 1;
							idx = ((uint64_t) r << 32) | c;
						}
						writers[idx]->addPoint(pt);
					}

					writers.clear();

					istr.close();
					Util::rm(file.filename);

				}

			} else {
				// This is the first run, so we read all the input files into the first set of
				// intermediate files.

				// Get the header from the first file to use as a template.
				std::ifstream istr(files[0].filename, std::ios::in | std::ios::binary);
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
						std::cout << outfile << "\n";
						uint64_t idx = ((uint64_t) r << 32) | c;
						writers[idx].reset(new LASWriter(outfile, ihdr, x, y));
					}
				}

				// Iterate over the files, filling the tiles.
				for(LASFile& file : files) {
					std::ifstream istr(file.filename, std::ios::in | std::ios::binary);
					liblas::Reader irdr = rfact.CreateWithStream(istr);
					while(irdr.ReadNextPoint()) {
						const liblas::Point& pt = irdr.GetPoint();
						double x = pt.GetX();
						double y = pt.GetY();
						int c = (int) ((x - allBounds[0]) / size0);
						int r = (int) ((y - allBounds[1]) / size0);
						uint64_t idx = ((uint64_t) r << 32) | c;
						writers[idx]->addPoint(pt);
					}
				}
			}

			size0 *= 0.5;
			cols0 = even(((int) (allBounds[2] - allBounds[0]) / size0) + 1);
			rows0 = even(((int) (allBounds[3] - allBounds[1]) / size0) + 1);

		} while(size0 >= size);

		writers.clear();


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
			<< " -s <srid>       The spatial reference ID. Default 0.\n";
}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	double size = 1000;
	double easting = 0;
	double northing = 0;
	uint16_t srid = 0;
	std::vector<std::string> args;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-t") {
			size = atof(argv[++i]);
		} else if(v == "-s") {
			srid = atoi(argv[++i]);
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
	r.tile(args[0], size, easting, northing, srid);

	return 0;
}


