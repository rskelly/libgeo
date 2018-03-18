/*
 * pc_filter.cpp
 *
 *  Created on: Mar 17, 2018
 *      Author: rob
 */

#include <vector>

#include "pointcloud.hpp"

using namespace geo::pc;

class PointClassFilter : public PointFilter {
public:
	std::vector<int> classes;
	PointClassFilter(int cls) {
		classes.push_back(cls);
	}
	PointClassFilter(const std::vector<int>& classes) :
		classes(classes) {
	}
	bool keep(const geo::pc::Point& pt) const {
		int cls = pt.classId();
		for(size_t i = 0; i < classes.size(); ++i) {
			if(classes[i] == cls)
				return true;
		}
		return false;
	}
	void print() const {
		std::stringstream ss;
		for(int c : classes)
			ss << c << ", ";
		g_debug("Class filter: " << ss.str());
	}
};

class PointScanAngleFilter : public PointFilter  {
public:
	double minAngle;
	double maxAngle;
	PointScanAngleFilter(double minAngle = -DBL_MAX, double maxAngle = DBL_MAX) :
		minAngle(minAngle), maxAngle(maxAngle) {
	}
	bool keep(const geo::pc::Point& pt) const {
		double a = pt.scanAngle();
		if(a < minAngle || a > maxAngle)
			return false;
		return true;
	}
	void print() const {
		g_debug("Scan angle filter: " << minAngle << ", " << maxAngle);
	}
};

class PointEdgeFilter : public PointFilter  {
public:
	bool keepEdge;
	PointEdgeFilter(bool keepEdge = false) :
		keepEdge(keepEdge) {
	}
	bool keep(const geo::pc::Point& pt) const {
		return keepEdge || !pt.isEdge();
	}
	void print() const {
		g_debug("Keep edge filter: " << keepEdge);
	}
};

class PointZRangeFilter : public PointFilter  {
public:
	double minZ;
	double maxZ;
	PointZRangeFilter(double minZ = -DBL_MAX, double maxZ = DBL_MAX) :
		minZ(minZ), maxZ(maxZ) {
	}
	bool keep(const geo::pc::Point& pt) const {
		double z = pt.z();
		return z >= minZ && z <= maxZ;
	}
	void print() const {
		g_debug("Z range filter: " << minZ << ", " << maxZ);
	}
};

class PointIntensityFilter : public PointFilter  {
public:
	double minIntensity;
	double maxIntensity;
	PointIntensityFilter(double minIntensity = -DBL_MAX, double maxIntensity = DBL_MAX) :
		minIntensity(minIntensity), maxIntensity(maxIntensity) {
	}
	bool keep(const geo::pc::Point& pt) const {
		double i = pt.intensity();
		return i >= minIntensity && i <= maxIntensity;
	}
	void print() const {
		g_debug("Intensity filter: " << minIntensity << ", " << maxIntensity);
	}
};

class PointFirstOnlyFilter : public PointFilter  {
public:
	bool keepFirst;
	PointFirstOnlyFilter(bool keepFirst = false) :
		keepFirst(keepFirst) {
	}
	bool keep(const geo::pc::Point& pt) const {
		return !keepFirst || pt.isFirst();
	}
	void print() const {
		g_debug("First only filter: " << keepFirst);
	}
};

class PointLastOnlyFilter : public PointFilter  {
public:
	bool keepLast;
	PointLastOnlyFilter(bool keepLast = false) :
		keepLast(keepLast) {
	}
	bool keep(const geo::pc::Point& pt) const {
		return !keepLast || pt.isLast();
	}
	void print() const {
		g_debug("Last only filter: " << keepLast);
	}
};

PCPointFilter::PCPointFilter() {
}

PCPointFilter::~PCPointFilter() {
	for(PointFilter* f : filters)
		delete f;
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
		<< " -p:f                 First returns only\n"
		<< " -p:l                 Last returns only\n";
}

bool PCPointFilter::parseArgs(int& idx, char** argv) {
	std::string v = argv[idx];
	bool found = false;
	if(v == "-p:c") {
		std::vector<std::string> tmp;
		std::vector<int> classes;
		std::string cls = argv[++idx];
		Util::splitString(std::back_inserter(tmp), cls);
		for(const std::string& t : tmp)
			classes.push_back(atoi(t.c_str()));
		addClassFilter(classes);
		found = true;
	} else if(v == "-p:l") {
		addKeepLastFilter();
		found = true;
	} else if(v == "-p:f") {
		addKeepFirstFilter();
		found = true;
	} else if(v == "-p:minz") {
		double minZ = atof(argv[++idx]);
		addZRangeFilter(minZ, DBL_MAX);
		found = true;
	} else if(v == "-p:maxz") {
		double maxZ = atof(argv[++idx]);
		addZRangeFilter(-DBL_MAX, maxZ);
		found = true;
	} else if(v == "-p:a") {
		std::string arg = argv[++idx];
		std::vector<std::string> parts;
		Util::splitString(std::back_inserter(parts), arg);
		addZRangeFilter(atof(parts[0].c_str()), atof(parts[1].c_str()));
		found = true;
	} else if(v == "-p:mini") {
		double minIntensity = atof(argv[++idx]);
		addIntensityFilter(minIntensity, DBL_MAX);
		found = true;
	} else if(v == "-p:maxi") {
		double maxIntensity = atof(argv[++idx]);
		addIntensityFilter(-DBL_MAX, maxIntensity);
		found = true;
	} else if(v == "-p:i") {
		std::string arg = argv[++idx];
		std::vector<std::string> parts;
		Util::splitString(std::back_inserter(parts), arg);
		addIntensityFilter(atof(parts[0].c_str()), atof(parts[1].c_str()));
		found = true;
	} else if(v == "-p:mina") {
		double minScanAngle = atof(argv[++idx]);
		addScanAngleFilter(minScanAngle, DBL_MAX);
		found = true;
	} else if(v == "-p:maxa") {
		double maxScanAngle = atof(argv[++idx]);
		addScanAngleFilter(-DBL_MAX, maxScanAngle);
		found = true;
	} else if(v == "-p:a") {
		std::string arg = argv[++idx];
		std::vector<std::string> parts;
		Util::splitString(std::back_inserter(parts), arg);
		addScanAngleFilter(atof(parts[0].c_str()), atof(parts[1].c_str()));
		found = true;
	} else if(v == "-p:e") {
		addRejectEdgeFilter();
		found = true;
	}
	return found;
}

void PCPointFilter::print() const {
	for(PointFilter* f : filters)
		f->print();
}

bool PCPointFilter::keep(const geo::pc::Point& pt) const {
	for(PointFilter* f : filters) {
		if(!f->keep(pt))
			return false;
	}
	return true;
}

void PCPointFilter::addClassFilter(int cls) {
	filters.push_back(new PointClassFilter(cls));
}

void PCPointFilter::addClassFilter(const std::vector<int>& cls) {
	filters.push_back(new PointClassFilter(cls));
}

void PCPointFilter::addIntensityFilter(double min, double max) {
	filters.push_back(new PointIntensityFilter(min, max));
}

void PCPointFilter::addZRangeFilter(double min, double max) {
	filters.push_back(new PointZRangeFilter(min, max));
}

void PCPointFilter::addScanAngleFilter(double min, double max) {
	filters.push_back(new PointIntensityFilter(min, max));
}

void PCPointFilter::addKeepLastFilter() {
	filters.push_back(new PointLastOnlyFilter(true));
}

void PCPointFilter::addKeepFirstFilter() {
	filters.push_back(new PointFirstOnlyFilter(true));
}

void PCPointFilter::addRejectEdgeFilter() {
	filters.push_back(new PointEdgeFilter(false));
}



