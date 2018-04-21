/*
 * pc_computer.hpp
 *
 *  Created on: Jan 27, 2018
 *      Author: rob
 */

#ifndef INCLUDE_PC_COMPUTER_HPP_
#define INCLUDE_PC_COMPUTER_HPP_

#include "pointcloud.hpp"
#include "rbf.hpp"

namespace geo {
namespace pc {
namespace compute {

bool pointSort(const geo::pc::Point& a, const geo::pc::Point& b);

template <class T, class U>
int pointFilter(U begin, U end, T insert, geo::pc::PCPointFilter* filter) {
	int count = 0;
	for(auto& it = begin; it < end; ++it) {
		if(filter && !filter->keep(*it)) continue;
		insert = *it;
		++count;
	}
	return count;
}

class IDWComputer : public geo::pc::Computer {
public:
	double exponent;
	IDWComputer(double exponent = 2);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class MeanComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class VarianceComputer : public MeanComputer {
private:
	double m_bias;
public:
	VarianceComputer(double bias = -1);
	void setBias(double bias);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class StdDevComputer : public VarianceComputer {
public:
	StdDevComputer(double bias = -1);
	void setBias(double bias);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class PercentileComputer : public geo::pc::Computer {
private:
	double m_percentile;
public:
	PercentileComputer(double percentile = .5);
	void setPercentile(double percentile);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class MaxComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class MinComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class DensityComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class CountComputer : public geo::pc::Computer {
private:
	std::list<geo::pc::Point> m_pts;
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class SkewComputer : public geo::pc::Computer {
public:
	MeanComputer meanComp;
	StdDevComputer stdDevComp;
	SkewComputer(double bias = -1);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class KurtosisComputer : public geo::pc::Computer {
public:
	MeanComputer meanComp;
	StdDevComputer stdDevComp;
	KurtosisComputer(double bias = -1);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class CoVComputer : public geo::pc::Computer {
public:
	MeanComputer meanComp;
	StdDevComputer sampStdDevComp;
	double bias;
	CoVComputer(double bias = -1);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class CanopyGapFractionComputer : public geo::pc::Computer {
public:
	MeanComputer meanComp;
	StdDevComputer sampStdDevComp;
	CanopyGapFractionComputer(double bias = -1);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
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
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class RugosityComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

class RBFComputer : public geo::pc::Computer {
public:
	RBFComputer(geo::interp::RBF<geo::pc::Point>::Type type, double smoothing);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
};

} // compute
} // pc
} // geo

#endif /* INCLUDE_PC_COMPUTER_HPP_ */
