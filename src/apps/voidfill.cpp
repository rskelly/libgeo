/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include "raster.hpp"

void usage() {
	std::cerr << "Usage: voidfill [options] <input raster> <output raster>\n"
			<< " -b  <band>       The band. Default 1.\n"
			<< " -m  <mask>       Filename of mask layer.\n"
			<< " -mb <mask band>  Band for mask.\n"
			<< " -r  <radius>     Radius for search.\n"
			<< " -e  <exponent>   IWD exponent.\n"
			<< " -c  <count>      Min number of pixels required to calculate a value for a cell.\n";
}

int main(int argc, char** argv) {

	using namespace geo::raster;

	if(argc < 3) {
		usage();
		return 1;
	}

	double radius = 10;
	int count = 4;
	double exp = 2.0;
	int band = 1;
	std::string mask;
	int maskBand = 1;
	std::vector<std::string> args;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-r") {
			radius = atof(argv[++i]);
		} else if(v == "-b") {
			band = atoi(argv[++i]);
		} else if(v == "-mb") {
			maskBand = atoi(argv[++i]);
		} else if(v == "-e") {
			exp = atof(argv[++i]);
		} else if(v == "-c") {
			count = atoi(argv[++i]);
		} else if(v == "-m") {
			mask = argv[++i];
		} else {
			args.push_back(argv[i]);
		}
	}

	if(count < 3)
		std::cerr << "Count should probably be more than " << count << ". Continuing anyway.\n";

	if(band < 1) {
		std::cerr << "Illegal band number: " << band << "\n";
		usage();
		return 1;
	}

	if(maskBand < 1) {
		std::cerr << "Illegal mask band number: " << band << "\n";
		usage();
		return 1;
	}

	if(exp <= 0) {
		std::cerr << "Illegal exponent: " << exp << "\n";
		usage();
		return 1;
	}

	if(args.size() < 2) {
		std::cerr << "Input and output filenames required.\n";
		usage();
		return 1;
	}

	Raster test(args[0]);
	test.voidFillIDW(args[1], radius, count, exp, band, mask, maskBand);

	return 0;
}


