/*
 * percentile.cpp
 *
 *  Created on: Apr 3, 2018
 *      Author: rob
 */

#include <vector>

#include "pointcloud.hpp"
#include "pc_computer.hpp"

using namespace geo::pc::compute;

PercentileComputer::PercentileComputer(double percentile) :
	m_percentile(percentile) {
	if(m_percentile <= 0 || m_percentile >= 1)
		g_runerr("Percentile must be between 0 and 1.");
}

void PercentileComputer::setPercentile(double percentile) {
	m_percentile = percentile;
}

int PercentileComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int PercentileComputer::compute(double, double, const std::vector<geo::pc::Point>&, const std::vector<geo::pc::Point>& filtered, double, std::vector<double>& out) {
	if(!filtered.empty()) {
		std::vector<geo::pc::Point> _pts(filtered);
		std::sort(_pts.begin(), _pts.end(), pointSort);
		int count = _pts.size();
		if(count % 2 == 0) {
			size_t idx = (size_t) (count * m_percentile) - 1;
			out.push_back((_pts[idx + 1].value() + _pts[idx].value()) / 2.0);
		} else {
			out.push_back(_pts[(size_t) (count * m_percentile)].value());
		}
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int PercentileComputer::bandCount() const {
	return 1;
}

