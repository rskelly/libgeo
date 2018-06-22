/**
 * las2tile re-tiles a set of LAS point cloud files into 
 * new tiles representing a regular grid. A user-configurable
 * buffer around each tile is optional. A user might want to
 * include a buffer for, for example, performing a ground classification,
 * then cropping off the buffer.
 *
 * las2tile works by iteratively creating "pyramids" of tiles. So,
 * if the desired tile size is 1000m, the program might start by tiling to
 * 8000m, then 4000m, then 2000m, then finally 1000m. This is done to
 * reduce the number of open file handles at any given time, as well
 * as saving memory. The program automatically removes the previous
 * level's tiles, and automatically handles intermediate tiles that 
 * exceed the maximum size of a LAS file.
 *
 *  Created on: Apr 13, 2017
 *  Author: rob
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
#include "pointcloud.hpp"

/**
 * Prints the useage message.
 */
void usage() {
	std::cerr << "Usage: las2grid <output dir> <input las [*]>\n"
			<< " -t <size>       The tile size in map units (default 1000; square).\n"
			<< " -e <easting>    The top left corner horizontal alignment\n"
			<< "                 (defaults to nearest whole multiple of size).\n"
			<< " -n <northing>   The top left corner vertical alignment\n"
			<< "                 (defaults to nearest whole multiple of size).\n"
			<< " -s <srid>       The spatial reference ID. Default 0.\n"
			<< " -b <buffer>     Buffer for each tile. Default 0.\n"
			<< " -h <handles>    Maximum number of open file handles. Default 64.\n";
}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	double size = 1000;
	double easting = std::nan("");
	double northing = std::nan("");
	double buffer = 0;
	uint16_t srid = 0;
	int handles = 64;

	std::vector<std::string> args;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-t") {
			size = atof(argv[++i]);
		} else if(v == "-s") {
			srid = atoi(argv[++i]);
		} else if(v == "-h") {
			handles = atoi(argv[++i]);
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
	geo::pc::Tiler r(infiles);
	r.tile(args[0], size, buffer, srid, easting, northing, handles);

	return 0;
}


