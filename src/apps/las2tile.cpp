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

		std::unordered_map<int, double*> bounds;
		std::unordered_map<int, std::ofstream> streams;
		std::unordered_map<int, std::unique_ptr<liblas::Writer> > writers;
		std::unordered_map<int, std::string> outfiles;
		std::unordered_map<int, int> pos;
		int currentOpen = -1;

		liblas::ReaderFactory rfact;

		{
			// Get the header from the first file to use as a template.
			std::ifstream str(files[0].filename, std::ios::in | std::ios::binary);
			liblas::Reader rdr = rfact.CreateWithStream(str);
			const liblas::Header& hdr = rdr.GetHeader();

			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c) {
					
					int idx = (r << 16) | c;
					
					std::stringstream ss;
					ss << "tile_" << (int) (c * size + allBounds[0]) << "_" << (int) (r * size + allBounds[3]) << ".las";
					
					double bnds[6] = {9999999999, 9999999999, -9999999999, -9999999999, 999999999, -999999999};

					liblas::Header whdr(hdr);
					std::ofstream of;

					// Set the bounds.
					bounds.insert(std::make_pair(idx, bnds));
					// Set the filename.
					outfiles[idx] = Util::pathJoin(outdir, ss.str());
					// Open the stream for Writer creation.
					streams[idx].open(outfiles[idx]);
					// Create the Writer.
					writers[idx].reset(new liblas::Writer(streams[idx], whdr));
					// Get the stream position and close it.
					pos[idx] = streams[idx].tellp();
					streams[idx].close();
				}
			}
		}

		for(const LASFile& f : files) {

			std::ifstream str(f.filename, std::ios::in | std::ios::binary);
			liblas::Reader rdr = rfact.CreateWithStream(str);

			while(rdr.ReadNextPoint()) {
				
				liblas::Point pt(rdr.GetPoint());
				
				double x = pt.GetX();
				double y = pt.GetY();
				double z = pt.GetZ();
				
				int col = (int) ((x - easting) / size);
				int row = (int) ((northing - y) / size);
				int idx = (row << 16) | col;
				
				if(idx != currentOpen) {
					if(currentOpen > -1 && streams[currentOpen].is_open()) {
						// If the currentOpen index is valid, get the 
						// open stream's position and close it.
						std::ofstream& of = streams[currentOpen]; 
						pos[currentOpen] = of.tellp();
						of.close();
					}
					// Open and seek the new stream, record active
					// index.
					std::ofstream& of = streams[idx];
					of.open(outfiles[idx]);
					of.seekp(pos[idx]);
					currentOpen = idx;
				}

				liblas::Writer* wtr = writers[idx].get();
				if(!wtr)
					continue;
				const liblas::Header& hdr = wtr->GetHeader();
				
				pt.SetHeader(&hdr);
				wtr->WritePoint(pt);
				
				double* bnds = bounds[idx];
				if(x < bnds[0]) bnds[0] = x;
				if(x > bnds[2]) bnds[2] = x;
				if(y < bnds[1]) bnds[1] = y;
				if(y > bnds[3]) bnds[3] = y;
				if(z < bnds[4]) bnds[4] = z;
				if(z > bnds[5]) bnds[5] = z;
			}
		}

		for(auto& p : writers) {
			double* bnds = bounds[p.first];
			liblas::Writer* wtr = writers[p.first].get();
			if(!wtr)
				continue;
			liblas::Header nhdr(wtr->GetHeader());
			nhdr.SetMax(bnds[2], bnds[3], bnds[5]);
			nhdr.SetMin(bnds[0], bnds[1], bnds[4]);
			wtr->SetHeader(nhdr);
			streams[p.first].close();
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


