/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include "raster.hpp"

void usage() {
	std::cerr << "Usage: voidfill [options] <input raster> <output raster>\n"
			<< " -b <band>       The band. Default 1.\n"
			<< " -t <method>     One of idw, nearest.\n"
			<< " -m <mask>       Filename of mask layer.\n";
}

int main(int argc, char** argv) {

	using namespace geo::raster;

	if(argc < 3) {
		usage();
		return 1;
	}

	uint16_t band = 1;
	std::string method;
	std::string mask;
	std::vector<std::string> args;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-b") {
			band = atoi(argv[++i]);
		} else if(v == "-t") {
			method = argv[++i];
		} else if(v == "-m") {
			mask = argv[++i];
		} else {
			args.push_back(argv[i]);
		}
	}

	if(args.size() < 2) {
		std::cerr << "Input and output filenames required.\n";
		usage();
		return 1;
	}

	if(args.size() < 3)
		args.push_back("dn");

	if(threads <= 0) {
		std::cerr << "Illegal thread value. Defaulting to 1.\n";
		threads = 1;
	}

	Raster test(args[0]);
	test.voidfill(args[1], band, method, mask);

	return 0;
}


