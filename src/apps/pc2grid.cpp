/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include <pointcloud.hpp>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <iostream>

#include "pointcloud.hpp"

void usage() {
	std::cerr << "Usage: pc2grid [options] <output raster> <input las [*]>\n"
			<< " -r <resolution> The output resolution in map units.\n"
			<< " -e <easting>    The top left corner horizontal alignment\n"
			<< "                 (defaults to nearest multiple of resolution).\n"
			<< " -n <northing>   The top left corner vertical alignment\n"
			<< "                 (defaults to nearest whole multiple of resolution).\n"
			<< " -d <radius>     The search radius for finding points.\n"
			<< " -s <srid>       The spatial reference ID. Default 0.\n"
			<< " -p <type>       The type of raster: mean (default), median, min, max, \n"
			<< "                 rugosity, variance, std. deviation and percentile.\n"
			<< "                 For percentile, use the form, 'percenile:n', where\n"
			<< "                 n is the percentile (no % sign); 1 - 99.\n"
			<< " -t <t>          The number of threads. Default 1.\n";
}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	double res = 0;
	double easting = 0;
	double northing = 0;
	double radius = 0;
	uint16_t srid = 0;
	uint16_t threads = 1;
	std::string type = "mean";
	std::vector<std::string> args;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-p") {
			type = argv[++i];
		} else if(v == "-r") {
			res = atof(argv[++i]);
		} else if(v == "-d") {
			radius = atof(argv[++i]);
		} else if(v == "-s") {
			srid = atoi(argv[++i]);
		} else if(v == "-e") {
			easting = atof(argv[++i]);
		} else if(v == "-n") {
			northing = atof(argv[++i]);
		} else if(v == "-t") {
			threads = atoi(argv[++i]);
		} else {
			args.push_back(argv[i]);
		}
	}

	if(res <= 0) {
		std::cerr << "Resolution must be >0.\n";
		usage();
		return 1;
	}

	if(args.size() < 2) {
		std::cerr << "Input and output filenames required.\n";
		usage();
		return 1;
	}

	if(threads <= 0) {
		std::cerr << "Illegal thread value. Defaulting to 1.\n";
		threads = 1;
	}

	std::vector<std::string> infiles(args.begin() + 1, args.end());
	geo::pc::Rasterizer r(infiles);
	r.rasterize(args[0], type, res, easting, northing, radius, srid, threads, 0);

	return 0;
}


