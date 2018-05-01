/*
 * rbf.cpp
 *
 *  Created on: Apr 21, 2018
 *      Author: rob
 */

/*
#include "pointcloud.hpp"
#include "rbf.hpp"
#include "pc_computer.hpp"

using namespace geo::pc;
using namespace geo::pc::compute;

class RBFComputer : public geo::pc::Computer {
public:
	RBFComputer(geo::interp::RBF<geo::pc::Point>::Type type, double smoothing);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};


int RBFComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int RBFComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	// TODO: Implement me!
	if(!filtered.empty()) {
		double area = radius * radius * M_PI;
		out.push_back(pts.size() / area);
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int RBFComputer::bandCount() const {
	return 1;
}
*/
