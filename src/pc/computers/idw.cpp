/*
 * idw.cpp
 *
 *  Created on: Apr 3, 2018
 *      Author: rob
 */

#include <vector>

#include "pointcloud.hpp"
#include "pc_computer.hpp"

using namespace geo::pc::compute;

IDWComputer::IDWComputer(double exponent) : exponent(exponent) {}

int IDWComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int IDWComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	if(!filtered.empty()) {
		double result = 0;
		double div = 0;
		for(const geo::pc::Point& pt : filtered) {
			double d = std::pow(x - pt.x(), 2.0) + std::pow(y - pt.y(), 2.0);
			if(d == 0) {
				result = pt.z();
				div = 1;
				break;
			} else {
				double w = std::pow(1 / d, exponent);
				result += w * pt.z();
				div += w;
			}
		}
		out.push_back(result / div);
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int IDWComputer::bandCount() const {
	return 1;
}
