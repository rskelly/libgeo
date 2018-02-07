/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include <fstream>
#include <vector>
#include <unordered_set>
#include <iostream>

#include <pointcloud.hpp>

#include "pointcloud.hpp"

using namespace geo::pc;

void usage() {
	std::cerr << "Usage: pc2grid [options] <output raster> <input las [*]>\n"
			<< " -r <resolution> The output resolution in map units.\n"
			<< " -e <easting>    The top left corner horizontal alignment\n"
			<< "                 (defaults to nearest multiple of resolution).\n"
			<< " -n <northing>   The top left corner vertical alignment\n"
			<< "                 (defaults to nearest whole multiple of resolution).\n"
			<< " -d <radius>     The search radius for finding points. Set to zero\n"
			<< "                 to use the rectangular cell bounds. Defaults to root\n"
			<< "                 of 1/2 cell diagonal (a circle that touches the corners.\n"
			<< " -s <srid>       The spatial reference ID. Default 0.\n"
			<< " -m <methods>    Comma-separated list of statistics to compute; one on each \n"
			<< "                 layer (see below.) Default mean.\n"
			<< "                 rugosity, variance, std. deviation and percentile.\n"
			<< "                 For percentile, use the form, 'percenile:n', where\n"
			<< "                 n is the percentile (no % sign); 1 - 99.\n"
			<< " -i <density>    The estimated number of points per cell underestimating this\n"
			<< "                 saves disk space at the cost of efficiency.\n"
			<< " -t              Estimate the per-cell point density. Do nothing else\n\n";

	PCPointFilter::printHelp(std::cerr);

	std::cerr << "\n Available computers: \n";
	for(auto& item : geo::pc::Rasterizer::availableComputers())
		std::cerr << " - " << item.first << ": " << item.second << ".\n";

}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	double res = 0;
	double easting = 0;
	double northing = 0;
	double radius = -1;
	int density = 64;
	bool estimate = false;
	uint16_t srid = 0;
	std::string mapFile;
	std::vector<std::string> types;
	std::vector<std::string> args;
	PCPointFilter filter;

	for(int i = 1; i < argc; ++i) {
		if(filter.parseArgs(i, argv))
			continue;
		std::string v = argv[i];
		if(v == "-m") {
			std::string type = argv[++i];
			Util::splitString(std::back_inserter(types), Util::lower(type), ",");
		} else if(v == "-r") {
			res = atof(argv[++i]);
		} else if(v == "-i") {
			density = atoi(argv[++i]);
		} else if(v == "-d") {
			radius = atof(argv[++i]);
		} else if(v == "-s") {
			srid = atoi(argv[++i]);
		} else if(v == "-e") {
			easting = atof(argv[++i]);
		} else if(v == "-n") {
			northing = atof(argv[++i]);
		} else if(v == "-t") {
			estimate = true;
		} else {
			args.push_back(argv[i]);
		}
	}

	if(types.empty()) {
		std::cerr << "No methods given; defaulting to mean.\n";
		types.push_back("mean");
	}
	if(res <= 0) {
		std::cerr << "Resolution must be >0.\n";
		usage();
		return 1;
	}

	if(radius < 0)
		radius = std::sqrt(std::pow(res / 2, 2) * 2);

	std::cerr << "Radius: " << radius << "\n";

	if(args.size() < 2) {
		std::cerr << "Input and output filenames required.\n";
		usage();
		return 1;
	}

	std::vector<std::string> infiles(args.begin() + 1, args.end());

	try {
		Rasterizer r(infiles);
		r.setFilter(filter);
		if(estimate) {
			density = (int) std::ceil(r.density(res, radius));
			std::cerr << "Estimated point density: " << density << "\n";
		} else {
			std::cerr << "Point density: " << density << "\n";
		}
		r.rasterize(args[0], types, res, easting, northing, radius, srid, density, 0, mapFile);
	} catch(const std::exception& ex) {
		std::cerr << ex.what() << "\n";
		usage();
	}
	return 0;
}


