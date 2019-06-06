/*
 * hlrg_biometrics.cpp
 *
 *  Created on: Apr 3, 2018
 *      Author: rob
 */

#include <vector>

#include "pointcloud.hpp"
#include "pc_computer.hpp"

using namespace geo::pc::compute;

namespace {

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

	int hlrgCCF(const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filteredPts, int bands, double threshold, std::vector<double>& out) {

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

		size_t n = pts.size();

		double l1 = 0;
		double l2 = 0;
		double l3 = 0;
		double l4 = 0;

		for(size_t i = 0; i < n; ++i) {
			double v = i + 1;
			double cl1 = v - 1;
			double cl2 = cl1 * (v - 2) / 2.0;
			double cl3 = cl2 * (v - 3) / 3.0;
			double cr1 = n - v;
			double cr2 = cr1 * (n - v - 1) / 2.0;
			double cr3 = cr2 * (n - v - 2) / 3.0;
			double z = pts[i].value();
			l1 +=  z;
			l2 += (cl1 - cr1) * z;
			l3 += (cl2 - 2 * cl1 * cr1 + cr2) * z;
			l4 += (cl3 - 3 * cl2 * cr1 + 3 * cl1 * cr2 - cr3) * z;
		}

		double c2 = n * (n - 1) / 2.0;
		double c3 = c2 * (n - 2) / 3.0;
		double c4 = c3 * (n - 3) / 4.0;
		l1 /= n;
		l2 /= c2 / 2.0;
		l3 /= c3 / 3.0;
		l4 /= c4 / 4.0;

		double lmean = l1;
		double lcov = l2 / l1;
		double lskew = l3 / l2;
		double lkurt = l4 / l2; // TODO: Should be l3?

		out.push_back(lmean);
		out.push_back(lcov);
		out.push_back(lskew);
		out.push_back(lkurt);

		return 4;
	}

}

HLRGBiometricsComputer::HLRGBiometricsComputer(int bands, int minCount) :
		m_bands(bands),
		m_minCount(minCount) {
	m_stdDev.setBias(-1);
	m_perc.setPercentile(.85);
}

int HLRGBiometricsComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int HLRGBiometricsComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {

	std::vector<geo::pc::Point> _pts(filtered);
	int count = _pts.size();

	std::sort(_pts.begin(), _pts.end(), pointSort);

	double threshold = 0;

	if(rasterizer()) {
		if(rasterizer()->filter())
			threshold = rasterizer()->filter()->minZRange();
	}

	// "Rugosity"
	m_stdDev.compute(x, y, pts, filtered, radius, out);

	// "Gap fraction."
	if(count) {
		out.push_back(1.0 - ((double) filtered.size() / pts.size()));
	} else {
		out.push_back(1.0);
	}

	// "85th percentile"
	m_perc.compute(x, y, pts, filtered, radius, out);

	// "L-Moments"
	if(count) {
		hlrgLMoments(_pts, out);
	} else {
		for(int i = 0; i < 4; ++i)
			out.push_back(std::nan(""));
	}

	// "LHQ"
	if(count) {
		hlrgLHQ(_pts, m_bands, out);
	} else {
		for(int i = 0; i < m_bands; ++i)
			out.push_back(std::nan(""));
	}

	// "Canopy closure fraction."
	if(count) {
		hlrgCCF(pts, _pts, m_bands, threshold, out);
	} else {
		for(int i = 0; i < m_bands; ++i)
			out.push_back(0);
	}

	return bandCount();
}

int HLRGBiometricsComputer::bandCount() const {
	return 49;
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
