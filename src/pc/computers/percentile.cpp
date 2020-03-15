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
		std::vector<geo::pc::Point> _pts(filtered);
		std::sort(_pts.begin(), _pts.end(), pointSort);
		int count = _pts.size();
		if(count % 2 == 0) {
			int idx1 = (int) (count * m_percentile) - 1;
			int idx2 = idx1 + 1;
			if(idx1 < 0) ++idx1;
			if(idx2 >= (int) _pts.size()) --idx2;
			out.push_back((_pts[idx1].value() + _pts[idx2].value()) / 2.0);
		} else {
			int idx = (int) (count * m_percentile) - 1;
			if(idx < 0) ++idx;
			if(idx >= (int) _pts.size()) --idx;
			out.push_back(_pts[idx].value());
		}
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

