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

int IDWComputer::compute(double x, double y, const std::vector<geo::pc::Point>&, const std::vector<geo::pc::Point>& filtered, double, std::vector<double>& out) {
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

std::vector<std::pair<std::string, std::string>> IDWComputer::bandMeta() const {
	return {{"idx_" + std::to_string(exponent), "idw - exponent: " + std::to_string(exponent)}};
}

