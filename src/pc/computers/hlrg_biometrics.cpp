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

		double htIncrement = (filteredPts[filteredPts.size() - 1].value() - threshold) / bands;
		double curHeight = threshold;

		for(int band = 0; band <= bands; ++band) {
			int ccfCount = getCCFCount(filteredPts, curHeight);
			double ccfPercent = (double) ccfCount / pts.size();
			out.push_back(ccfPercent);
			curHeight += htIncrement;
		}

		return bands + 1;
	}

	int hlrgLMoments(const std::vector<geo::pc::Point>& pts, std::vector<double>& out) {

		// Stolen from here: https://pypi.org/project/lmoments/0.1.0/#files

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

	double thresh = G_DBL_MAX_NEG;
	for(const PointFilter* filter : filters()) {
		const PointZRangeFilter* f;
		if((f = dynamic_cast<const PointZRangeFilter*>(filter)) != nullptr)
			thresh = f->minZ;
	}

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

std::vector<std::string> HLRGBiometricsComputer::bandMeta() const {
	std::vector<std::string> ret;
	ret.push_back("Rugosity (Std. Dev.)");
	ret.push_back("Gap Fraction");
	ret.push_back("85th Percentile");
	for(int i = 0; i < 4; ++i)
		ret.push_back("L-Moment " + std::to_string(i + 1));
	float slice = 100.0 / m_bands;
	for(int i = 0; i <= m_bands; ++i)
		ret.push_back("LHQ " + std::to_string((int) (i * slice)) + "%");
	for(int i = 0; i <= m_bands; ++i)
		ret.push_back("CCF " + std::to_string((int) (i * slice)) + "%");
	return ret;
}

int HLRGBiometricsComputer::bandCount() const {
	return 7 + (m_bands + 1) * 2;
}



/*	if(zValues.empty() || zValues.size() == pts.size())
		g_warn("To calculate CCF, points must have both veg (1) and ground classes (2).")

 * void fcLidarBLa(const std::list<std::unique_ptr<LiDARPoint> > &values, double *result) {
	double gnd = 0.0;
	double all = 0.0;
	for (const std::unique_ptr<LiDARPoint>& pt : values) {
		if (pt->isGround())
			gnd += pt->intensity();
		if (pt->cls() < 2) // TODO: This should perhaps be filtered by class to remove bogus points.
			all += pt->intensity();
	}
	result[0] = all != 0.0 ? 1.0 - std::sqrt(gnd / all) : -9999.0;
}

void fcLidarBLb(const std::list<std::unique_ptr<LiDARPoint> > &values, double *result) {
	double gndSingle = 0.0, gndLast = 0.0, first = 0.0, single = 0.0,
			intermediate = 0.0, last = 0.0, total = 0.0;
	for (const std::unique_ptr<LiDARPoint>& pt : values) {
		if (pt->isGround()) {
			if (pt->isSingle())
				gndSingle += pt->intensity();
			if (pt->isLast())
				gndLast += pt->intensity();
		}
		if (pt->isFirst())
			first += pt->intensity();
		if (pt->isSingle())
			single += pt->intensity();
		if (pt->isIntermediate())
			intermediate += pt->intensity();
		if (pt->isLast())
			last += pt->intensity();
		total += pt->intensity(); // TODO: This should perhaps be filtered by class to remove bogus points.
	}
	if (total == 0.0) {
		result[0] = -9999.0;
	} else {
		double denom = (first + single) / total
				+ std::sqrt((intermediate + last) / total);
		if (denom == 0.0) {
			result[0] = -9999.;
		} else {
			result[0] = (gndSingle / total + std::sqrt(gndLast / total))
					/ denom;
		}
	}
}

void fcLidarIR(const std::list<std::unique_ptr<LiDARPoint> > &values, double *result) {
	double canopy = 0.0, total = 0.0;
	for (const std::unique_ptr<LiDARPoint>& pt : values) {
		if (!pt->isGround())
			canopy += pt->intensity();
		total += pt->intensity();
	}
	result[0] = total != 0.0 ? canopy / total : -9999.0;
}

void fcLidarRR(const std::list<std::unique_ptr<LiDARPoint> > &values, double *result) {
	unsigned int canopy = 0, total = 0;
	for (const std::unique_ptr<LiDARPoint>& pt : values) {
		if (!pt->isGround())
			++canopy;
		++total;
	}
	result[0] = total != 0.0 ? (double) canopy / total : -9999.0;
}

void fcLidarFR(const std::list<std::unique_ptr<LiDARPoint> > &values, double *result) {
	unsigned int canopy = 0, total = 0;
	for (const std::unique_ptr<LiDARPoint>& pt : values) {
		if (pt->isFirst()) {
			if (!pt->isGround())
				++canopy;
			++total;
		}
	}
	result[0] = total != 0.0 ? (double) canopy / total : -9999.0;
}

void ccf(const std::list<std::unique_ptr<LiDARPoint> >& values, double *result, double threshold) {
	if (values.size() < 75) {
		result[0] = -9999.0;
	} else {
		double maxZ = -9999.0;
		for (const std::unique_ptr<LiDARPoint>& pt : values)
			maxZ = g_max(maxZ, pt->z());
		double htIncrement = (maxZ - threshold) / 20.0;
		double curHeight = threshold;
		for (int band = 0; band <= 20; ++band) {
			double count = 0;
			for (const std::unique_ptr<LiDARPoint>& pt : values) {
				if (pt->z() > curHeight)
					++count;
			}
			result[band] = (double) count / values.size();
			curHeight += htIncrement;
		}
	}
}

void gap(const std::list<std::unique_ptr<LiDARPoint> > &values, double *result, double threshold) {
	if (!values.size()) {
		result[0] = -9999.0;
	} else {
		int cnt = 0;
		for (const std::unique_ptr<LiDARPoint>& p : values) {
			if (p->z() > threshold)
				++cnt;
		}
		result[0] = 1.0 - ((double) cnt / values.size());
	}
}

CellGapFraction::CellGapFraction(unsigned char type, double threshold) :
		CellStats(), m_type(type), m_threshold(threshold) {
}

void CellGapFraction::threshold(double t) {
	m_threshold = t;
}

int CellGapFraction::bands() const {
	switch (m_type) {
	case GAP_CCF:
		return 21;
	default:
		return 1;
	}
}

void CellGapFraction::compute(double x, double y,
		const std::list<std::unique_ptr<LiDARPoint> > &values, double *result) {
	if (!values.size()) {
		result[0] = -9999.0;
	} else {
		switch (m_type) {
		case GAP_BLA:
			fcLidarBLa(values, result);
			break;
		case GAP_BLB:
			fcLidarBLb(values, result);
			break;
		case GAP_IR:
			fcLidarIR(values, result);
			break;
		case GAP_RR:
			fcLidarRR(values, result);
			break;
		case GAP_FR:
			fcLidarFR(values, result);
			break;
		case GAP_CCF:
			ccf(values, result, m_threshold);
			break;
		case GAP_GAP:
			gap(values, result, m_threshold);
			break;
		default:
			g_argerr("Unknown Gap Fraction method: " << m_type);
		}
	}
}
 *
 */



/*
 *         namespace pointstats_config {

            std::map<std::string, uint8_t> types = {
                {"Minimum", TYPE_MIN},
                {"Maximum", TYPE_MAX},
                {"Mean", TYPE_MEAN},
                {"Density", TYPE_DENSITY},
                {"Sample Variance", TYPE_VARIANCE},
                {"Sample Std. Dev.", TYPE_STDDEV},
                {"Population Variance", TYPE_PVARIANCE},
                {"Population Std. Dev.", TYPE_PSTDDEV},
                {"Count", TYPE_COUNT},
                {"Quantile", TYPE_QUANTILE},
                {"Median", TYPE_MEDIAN},
                {"Rugosity", TYPE_RUGOSITY},
                {"Kurtosis", TYPE_KURTOSIS},
                {"Skewness", TYPE_SKEW},
                {"Gap Fraction", TYPE_GAP_FRACTION},
                {"CoV", TYPE_COV}
            };
            std::map<std::string, uint8_t> attributes = {
                {"Height", ATT_HEIGHT},
                {"Intensity", ATT_INTENSITY}
            };
            std::map<std::string, uint8_t> gapFractionTypes = {
                {"IR", GAP_IR},
                {"BLa", GAP_BLA},
                {"BLb", GAP_BLB},
                {"RR", GAP_RR},
                {"FR", GAP_FR},
                {"CCF", GAP_CCF},
                {"GAP", GAP_GAP}
            };

            std::map<std::string, uint8_t> snapModes = {
                {"None" , SNAP_NONE},
                {"Grid" , SNAP_GRID},
                {"Origin" , SNAP_ORIGIN}
            };

            std::map<std::string, uint8_t> areaModes = {
                {"Full Cell", AREA_CELL},
                {"Radius", AREA_RADIUS}
            };

        } // config
 *
 */
