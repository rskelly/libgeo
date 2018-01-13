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

	LASFile(const std::string& filename) :
		filename(filename) {

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

class Tiler {
public:
	std::vector<LASFile> files;

	Tiler(const std::vector<std::string> filenames) {
		for(const std::string& filename : filenames)
			files.emplace_back(filename);
	}

	void tile(const std::string& outdir, double size,
		double easting, double northing, int srid) {

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
			northing = ((int) (allBounds[3] / size)) * size;

		int cols = (int) ((allBounds[2] - allBounds[0]) / size) + 1;
		int rows = (int) ((allBounds[3] - allBounds[1]) / size) + 1;

		std::cerr << "cols " << cols << "; rows " << rows << "\n";

		liblas::ReaderFactory rfact;

		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c) {
	
				std::stringstream ss;
				ss << "tile_" << (int) (c * size + allBounds[0]) << "_" << (int) (r * size + allBounds[1] + size) << ".las";
				std::string outfile = Util::pathJoin(outdir, ss.str());

				std::cout << outfile << "\n";

				std::ofstream ostr(outfile, std::ios::out | std::ios::binary);
				std::unique_ptr<liblas::Writer> owtr;

				double bnds[6] = {9999999999, 9999999999, -9999999999, -9999999999, 999999999, -999999999};
				int returns = 0;
				int retNum[5] = {0, 0, 0, 0, 0};

				double minx = (c * size) + easting;
				double maxx = minx + size;
				double maxy = northing - (r * size);
				double miny = maxy - size;

				for(LASFile& file : files) {

					if(file.bounds[0] > maxx || file.bounds[2] < minx || file.bounds[1] > maxy || file.bounds[3] < miny)
						continue;

					// Get the header from the first file to use as a template.
					std::ifstream istr(file.filename, std::ios::in | std::ios::binary);
					liblas::Reader irdr = rfact.CreateWithStream(istr);
					const liblas::Header& ihdr = irdr.GetHeader();

					if(!owtr.get()) {
						liblas::Header ohdr(ihdr);
						owtr.reset(new liblas::Writer(ostr, ohdr));
					}

					const liblas::Header& ohdr = owtr->GetHeader();

					while(irdr.ReadNextPoint()) {
						
						liblas::Point pt(irdr.GetPoint());
						
						double x = pt.GetX();
						double y = pt.GetY();

						if(x < minx || x >= maxx || y < miny || y >= maxy)
							continue;

						double z = pt.GetZ();
						
						pt.SetHeader(&ohdr);
						owtr->WritePoint(pt);
						
						++returns;
						retNum[pt.GetReturnNumber() - 1]++;

						if(x < bnds[0]) bnds[0] = x;
						if(x > bnds[2]) bnds[2] = x;
						if(y < bnds[1]) bnds[1] = y;
						if(y > bnds[3]) bnds[3] = y;
						if(z < bnds[4]) bnds[4] = z;
						if(z > bnds[5]) bnds[5] = z;
					}
				}

				if(returns == 0) {
					owtr.reset();
					ostr.close();
					Util::rm(outfile);
				} else {
					liblas::Header ohdr(owtr->GetHeader());
					ohdr.SetMin(bnds[0], bnds[2], bnds[4]);
					ohdr.SetMax(bnds[1], bnds[3], bnds[5]);
					for(int i = 0; i < 5; ++i)
						ohdr.SetPointRecordsByReturnCount(i, retNum[i]);
					ohdr.SetPointRecordsCount(returns);
					owtr->SetHeader(ohdr);
				}
			}
		}

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


