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
	std::cerr << "Usage: pcnorm [options] <output dir> <terrain model> <input las [*]>\n"
			<< " Point filtering parameters:\n"
			<< " -c <class(es)>  Comma-delimited list of classes to keep.\n"
			<< " -t <threshold>  The minimum height threshold.\n";
}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	std::vector<std::string> args;
	std::vector<int> classes;
	geo::pc::PCPointFilter filter;
	int band = 1;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-c") {
			std::vector<std::string> tmp;
			std::string cls = argv[++i];
			Util::splitString(std::back_inserter(tmp), cls);
			for(const std::string& t : tmp)
				classes.push_back(atoi(t.c_str()));
			filter.addClassFilter(classes);
		} else if(v == "-t") {
			filter.addZRangeFilter(atof(argv[++i]), DBL_MAX);
		} else if(v == "-b") {
			band = atoi(argv[++i]);
		} else {
			args.push_back(argv[i]);
		}
	}

	if(args.size() < 3) {
		std::cerr << "Output dir, terrain model and input files are required.\n";
		usage();
		return 1;
	}

	std::vector<std::string> infiles(args.begin() + 2, args.end());

	try {
		geo::pc::Normalizer n(infiles);
		n.setFilter(filter);
		n.normalize(args[1], args[0], band);
	} catch(const std::exception& ex) {
		std::cerr << ex.what() << "\n";
		usage();
	}
	return 0;
}


