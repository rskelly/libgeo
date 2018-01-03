/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include "raster.hpp"

void usage() {
	std::cerr << "Usage: polygonize <input raster> <output vector> [field name (dn)]\n"
			<< " -f <format>     The output format. Any of the available\n"
			<< "                 OGR drivers. Default SQLite.\n"
			<< " -s <srid>       The spatial reference ID. Default 0.\n"
			<< " -b <band>       The band. Default 1.\n"
			<< " -d              Remove dangles.\n"
			<< " -h              Remove holes.\n"
			<< " -t <t>          The number of threads. Default 1.\n";
}

int main(int argc, char** argv) {

	using namespace geo::raster;

	if(argc < 3) {
		usage();
		return 1;
	}

	std::string driver = "SQLite";
	uint16_t srid = 0;
	uint16_t band = 1;
	uint16_t threads = 1;
	bool holes = false;
	bool dangles = false;
	std::vector<std::string> args;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-f") {
			driver = argv[++i];
		} else if(v == "-s") {
			srid = atoi(argv[++i]);
		} else if(v == "-b") {
			band = atoi(argv[++i]);
		} else if(v == "-h") {
			holes = true;
		} else if(v == "-d") {
			dangles = true;
		} else if(v == "-t") {
			threads = atoi(argv[++i]);
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
	test.polygonize(args[1], args[2], driver, srid, band, threads, holes, dangles, nullptr, nullptr);

	return 0;
}


