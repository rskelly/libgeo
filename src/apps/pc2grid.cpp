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

#include "pointcloud.hpp"

using namespace geo::pc;

void usage() {
	std::cerr << "Usage: pc2grid [options] <output raster> <input las [*]>\n"
			<< " -rx <resolution> The output x resolution in map units.\n"
			<< " -ry <resolution> The output y resolution in map units.\n"
			<< " -e <easting>     The top left corner horizontal alignment\n"
			<< "                  (defaults to nearest multiple of resolution).\n"
			<< " -n <northing>    The top left corner vertical alignment\n"
			<< "                  (defaults to nearest whole multiple of resolution).\n"
			<< " -d <radius>      The search radius for finding points. Set to zero\n"
			<< "                  to use the rectangular cell bounds. Defaults to root\n"
			<< "                  of 1/2 cell diagonal (a circle that touches the corners.\n"
			<< " -s <srid>        The spatial reference ID. Default 0.\n"
			<< " -m <methods>     Comma-separated list of statistics to compute; one on each \n"
			<< "                  layer (see below.) Default mean.\n"
			<< "                  rugosity, variance, std. deviation and percentile.\n"
			<< "                  For percentile, use the form, 'percenile:n', where\n"
			<< "                  n is the percentile (no % sign); 1 - 99.\n"
			<< " -v               Verbose. Enable debug and warning messages.\n"
			<< " -h               If given *do not* trust the LAS file headers to contain\n"
			<< "                  good bounds (etc.) info\n"
			<< " -t <min>         Normalize the point density in each cell. Any cell with more than this\n"
			<< "                  number of points will be randomly thinned. Any cell with less will be\n"
			<< "                  reduced to zero. This occurs after filtering.\n"
			<< " -d <nodata>      A nodata value. The default is -9999.0\n"
			<< " -o               Fill voids. Does this by doubling the search radius iteratively.\n"
			<< "                  the point count in this cells remains at zero.\n"
			<< " -l <bytes>       The limit of memory devoted to point data that will trigger the use\n"
			<< "                  of file-backed memory. This will be slow but less likely to crash.\n"
			<< " -c <scale>       A scale value to scale coordinates for storage in the quadtree.\n"
			<< "                  must be non-zero, usually a multiple or fraction of 10. Default 1.\n";

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

	double resX = std::nan("");
	double resY = std::nan("");
	double easting = std::nan("");
	double northing = std::nan("");
	double radius = std::nan("");
	uint16_t srid = 0;
	bool useHeader = true;
	std::vector<std::string> types;
	std::vector<std::string> args;
	PCPointFilter filter;
	int thin = 0;
	double nodata = -9999;
	bool voids = false;
	size_t limit = 0;
	double scale = 1;

	for(int i = 1; i < argc; ++i) {
		if(filter.parseArgs(i, argv))
			continue;
		std::string v = argv[i];
		if(v == "-m") {
			std::string type = argv[++i];
			split(std::back_inserter(types), lowercase(type), ",");
		} else if(v == "-h") {
			useHeader = false;
		} else if(v == "-v") {
			g_loglevel(G_LOG_TRACE);
		} else if(v == "-rx") {
			resX = atof(argv[++i]);
		} else if(v == "-ry") {
			resY = atof(argv[++i]);
		} else if(v == "-d") {
			radius = atof(argv[++i]);
		} else if(v == "-s") {
			srid = atoi(argv[++i]);
		} else if(v == "-e") {
			easting = atof(argv[++i]);
		} else if(v == "-n") {
			northing = atof(argv[++i]);
		} else if(v == "-t") {
			thin = atoi(argv[++i]);
		} else if(v == "-l") {
			limit = atoi(argv[++i]);
		} else if(v == "-d") {
			nodata = atof(argv[++i]);
		} else if(v == "-c") {
			scale = atof(argv[++i]);
		} else if(v == "-o") {
			voids = true;
		} else {
			args.push_back(argv[i]);
		}
	}

	filter.print();
	
	if(args.size() < 2) {
		std::cerr << "Input and output filenames required.\n";
		usage();
		return 1;
	}

	std::vector<std::string> infiles(args.begin() + 1, args.end());

	try {
		Rasterizer r(infiles);
		r.setFilter(&filter);
		r.setThin(thin);
		r.setNoData(nodata);
		r.setMemLimit(limit);
		r.setScale(scale);
		r.rasterize(args[0], types, resX, resY, easting, northing, radius, srid, useHeader, voids);
	} catch(const std::exception& ex) {
		std::cerr << ex.what() << "\n";
		usage();
	}
	return 0;
}


