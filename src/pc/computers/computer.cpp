#include <vector>

#include "pointcloud.hpp"
#include "pc_computer.hpp"

using namespace geo::pc;
using namespace geo::pc::compute;

bool geo::pc::compute::pointSort(const geo::pc::Point& a, const geo::pc::Point& b) {
	return a.value() < b.value();
}

int Computer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int DensityComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int DensityComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	if(!filtered.empty()) {
		double area = radius * radius * M_PI;
		out.push_back(pts.size() / area);
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int DensityComputer::bandCount() const {
	return 1;
}



