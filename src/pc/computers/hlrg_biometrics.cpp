/*
 * hlrg_biometrics.cpp
 *
 * Bands:
 * 1) count
 * 2) std. dev
 * 3) gap fraction
 * 4) 85th percentile
 * 5-9) L-Moments
 * 10-29) LHQ
 * 30-50) CCF
 *
 *  Created on: Apr 3, 2018
 *      Author: rob
 */

#include <vector>

#include "pointcloud.hpp"
#include "pc_computer.hpp"
#include "pc_filter.hpp"

using namespace geo::pc::compute;
using namespace geo::pc::filter;

namespace {

	int comb(int n, int k) {
		if(k > n || n < 0 || k < 0)
			return 0;
		int v = 1;
		for(int j = 0; j < std::min(k, n - k); ++j)
			v = (v * (n - j)) / (j + 1);
		return v;
	}

	int hlrgLHQ(const std::vector<geo::pc::Point>& pts, int bands, std::vector<double>& out) {

		if(pts.empty()) {
			for(int band = 0; band <= bands; ++band)
				out.push_back(std::nan(""));
		} else {

			out.push_back(pts[0].value());

			for(int band = 1; band < bands; ++band) {
				double rank = ((double) band / bands) * (pts.size() + 1);
				int base = (int) std::floor(rank);
				double diff = rank - base;
				int rank1 = base - 1; if(rank1 < 0) rank1 = 0;
				int rank2 = base;
				out.push_back(((1.0 - diff) * pts[rank1].value()) + (diff * pts[rank2].value()));
			}

			out.push_back(pts[pts.size() - 1].value());
		}

		return bands + 1;

	}

	int getCCFCount(const std::vector<geo::pc::Point>& pts, double thresh) {
		int count = 0;
		for(const geo::pc::Point& pt : pts) {
			if(pt.value() > thresh)
				++count;
		}
		return count;
	}

	int hlrgCCF(const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filteredPts,
			int bands, double threshold, std::vector<double>& out) {

		if(pts.empty() || filteredPts.empty()) {
			for(int band = 0; band <= bands; ++band)
				out.push_back(std::nan(""));
		} else {
			double htIncrement = (filteredPts.back().value() - threshold) / bands;
			double curHeight = threshold;

			for(int band = 0; band <= bands; ++band) {
				int ccfCount = getCCFCount(filteredPts, curHeight);
				double ccfPercent = (double) ccfCount / pts.size();
				out.push_back(ccfPercent);
				curHeight += htIncrement;
			}
		}

		return bands + 1;
	}

	int hlrgLMoments(const std::vector<geo::pc::Point>& pts, std::vector<double>& out) {

		// Stolen from here: https://pypi.org/project/lmoments/0.1.0/#files

		if(pts.empty()) {
			for(int i = 0; i < 4; ++i)
				out.push_back(std::nan(""));
		} else {

			size_t n = pts.size();

			std::vector<double> com1; com1.reserve(n);
			std::vector<double> com2; com2.reserve(n);
			std::vector<double> com3; com3.reserve(n);
			std::vector<double> com4; com4.reserve(n);
			std::vector<double> com5; com5.reserve(n);
			std::vector<double> com6; com6.reserve(n);

			for(size_t i = 1; i < n + 1; ++i) {
				com1.push_back(comb(i - 1, 1));
				com2.push_back(comb(n - i, 1));
				com3.push_back(comb(i - 1, 2));
				com4.push_back(comb(n - i, 2));
				com5.push_back(comb(i - 1, 3));
				com6.push_back(comb(n - i, 3));
			}

			double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;

			for(const geo::pc::Point& p : pts)
				sum1 += p.z();

			double co1 = 1.0 / comb(n, 1);
			double co2 = 0.5 * 1.0 / comb(n, 2);
			double co3 = 1.0 / 3.0 * 1.0 / comb(n, 3);
			double co4 = 1.0 / 4.0 * 1.0 / comb(n, 5);
			double v;

			for(size_t i = 0; i < n; ++i) {
				v = pts[i].value();
				double tmp2 = com1[i] - com2[i];
				double tmp3 = com3[i] - 2 * com1[i] * com2[i] + com4[i];
				double tmp4 = com5[i] - 3 * com3[i] * com2[i] + 3 * com1[i] * com4[i] - com6[i];
				sum2 += tmp2 * v;
				sum3 += tmp3 * v;
				sum4 += tmp4 * v;
			}

			double l1 = co1 * sum1;
			double l2 = co2 * sum2;
			double l3 = co3 * sum3 / l2;
			double l4 = co4 * sum4 / l2;

			out.push_back(l1);
			out.push_back(l2);
			out.push_back(l3);
			out.push_back(l4);
		}

		return 4;
	}

}

HLRGBiometricsComputer::HLRGBiometricsComputer(int bands, int minCount) :
		m_bands(bands),
		m_minCount(minCount) {
	m_stdDev.setBias(-1);
	m_perc.setPercentile(.85);
}


int HLRGBiometricsComputer::compute(double x, double y, const std::vector<geo::pc::Point>& /*pts*/, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {

	// TODO: Canopy metrics are performed on the filtered set, split by a threshold
	// value. Thus the comparison is between filtered points>threshold vs. all filtered points.
	// I do not agree with this strategy but it is what it is.
	// (Inspired by a problem with a dataset where the canopy density varied differently from
	// the ground density, even though the cover seemed to be completely uniform due to a
	// postprocessing error.

	double thresh = std::nan("");
	for(const PointFilter* filter : filters()) {
		const PointZRangeFilter* f;
		if((f = dynamic_cast<const PointZRangeFilter*>(filter)) != nullptr)
			thresh = f->minZ;
	}
	if(std::isnan(thresh))
		g_runerr("A threshold must be given for HRLG biometrics, even if it's zero. Use -f:hlrg-bio:minz <z>");

	std::vector<geo::pc::Point> _allpts(filtered);
	int count = _allpts.size();
	std::sort(_allpts.begin(), _allpts.end(), pointSort);
	std::vector<geo::pc::Point> _threshpts;
	filter(_allpts, _threshpts);

	// "Rugosity"
	m_stdDev.compute(x, y, _allpts, _threshpts, radius, out);

	// "Gap fraction."
	if(count) {
		out.push_back(1.0 - ((double) _threshpts.size() / _allpts.size()));
	} else {
		out.push_back(1.0);
	}

	// "85th percentile"
	m_perc.compute(x, y, _allpts, _threshpts, radius, out);

	// "L-Moments"
	if(count) {
		hlrgLMoments(_threshpts, out);
	} else {
		for(int i = 0; i < 4; ++i)
			out.push_back(std::nan(""));
	}

	// "LHQ"
	if(count) {
		hlrgLHQ(_threshpts, m_bands, out);
	} else {
		for(int i = 0; i <= m_bands; ++i)
			out.push_back(std::nan(""));
	}

	// "Canopy closure fraction."
	if(count) {
		hlrgCCF(_allpts, _threshpts, m_bands, thresh, out);
	} else {
		for(int i = 0; i <= m_bands; ++i)
			out.push_back(0);
	}

	return bandCount();
}

std::vector<std::pair<std::string, std::string>> HLRGBiometricsComputer::bandMeta() const {
	std::vector<std::pair<std::string, std::string>> ret;
	ret.emplace_back("hlrg_rugosity", "Rugosity (Std. Dev.)");
	ret.emplace_back("hlrg_gap", "Gap Fraction");
	ret.emplace_back("hlrg_p85", "85th Percentile");
	for(int i = 0; i < 4; ++i) {
		std::string l = std::to_string(i + 1);
		ret.emplace_back("hlrg_lmoment_" + l, "L-Moment " + l);
	}
	float slice = 100.0 / m_bands;
	for(int i = 0; i <= m_bands; ++i) {
		std::string s = std::to_string((int) (i * slice));
		ret.emplace_back("hlrg_lhq_" + s, "LHQ " + s + "%");
	}
	for(int i = 0; i <= m_bands; ++i) {
		std::string s = std::to_string((int) (i * slice));
		ret.emplace_back("hlrg_ccf_" + std::to_string((int) (i * slice)), "CCF " + s + "%");
	}
	return ret;
}

int HLRGBiometricsComputer::bandCount() const {
	return 7 + (m_bands + 1) * 2;
}

