/*
 * pc_computer.hpp
 *
 *  Created on: Jan 27, 2018
 *      Author: rob
 */

#ifndef INCLUDE_PC_COMPUTER_HPP_
#define INCLUDE_PC_COMPUTER_HPP_

#include "pointcloud.hpp"

class IDWComputer : public geo::pc::Computer {
public:
	double exponent;
	IDWComputer(double exponent = 2);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

class MeanComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

class VarianceComputer : public MeanComputer {
private:
	double m_bias;
public:
	VarianceComputer(double bias = -1);
	void setBias(double bias);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

class StdDevComputer : public VarianceComputer {
public:
	StdDevComputer(double bias = -1);
	void setBias(double bias);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

class PercentileComputer : public geo::pc::Computer {
private:
	double m_percentile;
public:
	PercentileComputer(double percentile = .5);
	void setPercentile(double percentile);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

class MaxComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

class MinComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

class DensityComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

class CountComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

class SkewComputer : public geo::pc::Computer {
public:
	MeanComputer meanComp;
	StdDevComputer stdDevComp;
	SkewComputer(double bias = -1);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

class KurtosisComputer : public geo::pc::Computer {
public:
	MeanComputer meanComp;
	StdDevComputer stdDevComp;
	KurtosisComputer(double bias = -1);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

class CoVComputer : public geo::pc::Computer {
public:
	MeanComputer meanComp;
	StdDevComputer sampStdDevComp;
	double bias;
	CoVComputer(double bias = -1);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

class CanopyGapFractionComputer : public geo::pc::Computer {
public:
	MeanComputer meanComp;
	StdDevComputer sampStdDevComp;
	CanopyGapFractionComputer(double bias = -1);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

class HLRGBiometricsComputer : public geo::pc::Computer {
private:
	StdDevComputer m_stdDev;
	PercentileComputer m_perc;
	int m_bands;
	int m_minCount;
	double m_threshold;
public:
	HLRGBiometricsComputer(int bands, int minCount, double threshold);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

/*
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

class RugosityComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out);
	int bandCount() const;
};

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



#endif /* INCLUDE_PC_COMPUTER_HPP_ */
