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

int PercentileComputer::compute(double, double, const std::vector<geo::pc::Point>&, const std::vector<geo::pc::Point>& filtered, double, std::vector<double>& out) {
	if(!filtered.empty()) {
		static std::vector<geo::pc::Point> _pts;
		_pts.assign(filtered.begin(), filtered.end());
		std::sort(_pts.begin(), _pts.end(), pointSort);
		int idx = (_pts.size() - 1) * m_percentile;
		out.push_back(_pts[idx].z());
		_pts.clear();
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int PercentileComputer::bandCount() const {
	return 1;
}

std::vector<std::string> PercentileComputer::bandMeta() const {
	return {"percentile: " + std::to_string(m_percentile)};
}

