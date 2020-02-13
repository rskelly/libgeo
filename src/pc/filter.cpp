/*
 * pc_filter.cpp
 *
 *  Created on: Mar 17, 2018
 *      Author: rob
 */

#include <vector>

#include "pointcloud.hpp"
#include "pc_filter.hpp"

using namespace geo::pc;
using namespace geo::pc::filter;


PCPointFilter::PCPointFilter() {
}

PCPointFilter::~PCPointFilter() {
	for(PointFilter* f : m_filters)
		delete f;
	for(auto& it : m_computerFilters) {
		for(PointFilter* f : it.second)
			delete f;
	}
}

void PCPointFilter::printHelp(std::ostream& str) {
	str << " Point filtering parameters:\n"
		<< " -p:c <class(es)>     Comma-delimited list of classes to keep.\n"
		<< " -p:minz <z>          The minimum height threshold.\n"
		<< " -p:maxz <z>          The maximum height threshold.\n"
		<< " -p:z    <min,max>    The minimum and maximum height threshold.\n"
		<< " -p:mini <intensity>  The minimum intensity.\n"
		<< " -p:maxi <intensity>  The maximum intensity.\n"
		<< " -p:i    <min,max>    The minimum and maximum intensity threshold.\n"
		<< " -p:mina <angle>      The minimum scan angle.\n"
		<< " -p:maxa <angle>      The maximum scan angle.\n"
		<< " -p:a    <min,max>    The minimum and maximum scan angle threshold.\n"
		<< " -p:f                 First returns only.\n"
		<< " -p:l                 Last returns only.\n"
		<< " -p:e                 Drop edges.\n"
		<< " Computer-specific filtering:\n"
		<< " -f:[method]:[filter] Use the filter part from the point filtering parameters;\n"
		<< "                      The filter will be applied inside the computer on points\n"
		<< "                      that have been pre-filtered by the general filter (or not\n"
		<< "                      filtered if no general filter is given.\n";
}

bool PCPointFilter::parseArgs(int& idx, char** argv) {
	std::string v = argv[idx];
	bool found = false;
	if(v.substr(0, 3) == "-p:") {
		addFilter(getFilter(v.substr(3), argv, idx));
		found = true;
	} else if(v.substr(0, 3) == "-f:") {
		addComputerFilter(v.substr(3), argv, idx);
		found = true;
	}
	return found;
}

PointFilter* PCPointFilter::getFilter(const std::string& v, char** argv, int& idx) {
	if(v == "c") {
		std::vector<std::string> tmp;
		std::vector<int> classes;
		std::string cls = argv[++idx];
		split(std::back_inserter(tmp), cls);
		for(const std::string& t : tmp)
			classes.push_back(atoi(t.c_str()));
		return new PointClassFilter(classes);
	} else if(v == "l") {
		return new PointLastOnlyFilter(true);
	} else if(v == "f") {
		return new PointFirstOnlyFilter(true);
	} else if(v == "minz") {
		double minZ = atof(argv[++idx]);
		return new PointZRangeFilter(minZ, G_DBL_MAX_POS);
	} else if(v == "maxz") {
		double maxZ = atof(argv[++idx]);
		return new PointZRangeFilter(-G_DBL_MAX_NEG, maxZ);
	} else if(v == "z") {
		std::string arg = argv[++idx];
		std::vector<std::string> parts;
		split(std::back_inserter(parts), arg);
		return new PointZRangeFilter(atof(parts[0].c_str()), atof(parts[1].c_str()));
	} else if(v == "mini") {
		double minIntensity = atof(argv[++idx]);
		return new PointIntensityFilter(minIntensity, DBL_MAX);
	} else if(v == "maxi") {
		double maxIntensity = atof(argv[++idx]);
		return new PointIntensityFilter(-DBL_MAX, maxIntensity);
	} else if(v == "i") {
		std::string arg = argv[++idx];
		std::vector<std::string> parts;
		split(std::back_inserter(parts), arg);
		return new PointIntensityFilter(atof(parts[0].c_str()), atof(parts[1].c_str()));
	} else if(v == "mina") {
		double minScanAngle = atof(argv[++idx]);
		return new PointScanAngleFilter(minScanAngle, DBL_MAX);
	} else if(v == "maxa") {
		double maxScanAngle = atof(argv[++idx]);
		return new PointScanAngleFilter(-DBL_MAX, maxScanAngle);
	} else if(v == "a") {
		std::string arg = argv[++idx];
		std::vector<std::string> parts;
		split(std::back_inserter(parts), arg);
		return new PointScanAngleFilter(atof(parts[0].c_str()), atof(parts[1].c_str()));
	} else if(v == "e") {
		return new PointEdgeFilter();
	} else {
		g_argerr("Unknown filter: " << v);
	}
}

void PCPointFilter::addComputerFilter(const std::string& arg, char** argv, int& idx) {
	std::string method = arg.substr(0, arg.find(':'));
	std::string filter = arg.substr(arg.find(':') + 1);
	m_computerFilters[method].push_back(getFilter(filter, argv, idx));
}

const std::vector<PointFilter*> PCPointFilter::computerFilters(const std::string& name) const {
	if(m_computerFilters.find(name) == m_computerFilters.end()) {
		return {};
	} else{
		return m_computerFilters.at(name);
	}
}

void PCPointFilter::print() const {
	for(PointFilter* f : m_filters)
		f->print();
}

bool PCPointFilter::keep(const geo::pc::Point& pt) const {
	for(PointFilter* f : m_filters) {
		if(!f->keep(pt))
			return false;
	}
	return true;
}

void PCPointFilter::addFilter(PointFilter* filter) {
	m_filters.push_back(filter);
}

double PCPointFilter::minZRange() const {
	PointZRangeFilter* zf;
	double z = std::numeric_limits<double>::max();
	bool found = false;
	for(PointFilter* f : m_filters) {
		if((zf = dynamic_cast<PointZRangeFilter*>(f)) != nullptr) {
			if(z > zf->minZ) {
				z = zf->minZ;
				found = true;
			}
		}
	}
	return found ? z : G_DBL_MAX_NEG; // NOTE: Compromise -- if there's no filter, it's the minimum. This ensures that when no filter is set all >0 points are kept (with maxZRange).
}

double PCPointFilter::maxZRange() const {
	PointZRangeFilter* zf;
	double z = std::numeric_limits<double>::lowest();
	bool found = false;
	for(PointFilter* f : m_filters) {
		if((zf = dynamic_cast<PointZRangeFilter*>(f)) != nullptr) {
			if(z < zf->maxZ) {
				z = zf->maxZ;
				found = true;
			}
		}
	}
	return found ? z : G_DBL_MAX_POS; // NOTE: Compromise -- if there's no filter, it's zero. This ensures that when no filter is set all >0 points are kept (with maxZRange).
}



