/*
 * pc_computer.hpp
 *
 *  Created on: Jan 27, 2018
 *      Author: rob
 */

#ifndef INCLUDE_PC_COMPUTER_HPP_
#define INCLUDE_PC_COMPUTER_HPP_

#include "pointcloud.hpp"
//#include "rbf.hpp"

namespace geo {
namespace pc {
namespace compute {


class IDWComputer : public geo::pc::Computer {
public:
	double exponent;
	IDWComputer(double exponent = 2);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

class MeanComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

class VarianceComputer : public MeanComputer {
private:
	double m_bias;
public:
	VarianceComputer(double bias = -1);
	void setBias(double bias);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

class StdDevComputer : public VarianceComputer {
public:
	StdDevComputer(double bias = -1);
	void setBias(double bias);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

class PercentileComputer : public geo::pc::Computer {
private:
	double m_percentile;
public:
	PercentileComputer(double percentile = .5);
	void setPercentile(double percentile);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

class MaxComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

class MinComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

class DensityComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

class CountComputer : public geo::pc::Computer {
private:
	std::list<geo::pc::Point> m_pts;
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

class SkewComputer : public geo::pc::Computer {
public:
	MeanComputer meanComp;
	StdDevComputer stdDevComp;
	SkewComputer(double bias = -1);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

class KurtosisComputer : public geo::pc::Computer {
public:
	MeanComputer meanComp;
	StdDevComputer stdDevComp;
	KurtosisComputer(double bias = -1);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

class CoVComputer : public geo::pc::Computer {
public:
	MeanComputer meanComp;
	StdDevComputer sampStdDevComp;
	double bias;
	CoVComputer(double bias = -1);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

class CanopyGapFractionComputer : public geo::pc::Computer {
public:
	MeanComputer meanComp;
	StdDevComputer sampStdDevComp;
	CanopyGapFractionComputer(double bias = -1);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

/**
 *
 */
class HLRGBiometricsComputer : public geo::pc::Computer {
private:
	StdDevComputer m_stdDev;
	PercentileComputer m_perc;
	int m_bands;
	int m_minCount;
public:
	/**
	 * \param bands The number of levels; with 20 bands, each represents a 5% slice.
	 * \param minCount Minimum number of points required for acomputation.
	 */
	HLRGBiometricsComputer(int bands, int minCount);
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};

class RugosityComputer : public geo::pc::Computer {
public:
	int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out);
	int bandCount() const;
	std::vector<std::pair<std::string,std::string>> bandMeta() const;
};


} // compute
} // pc
} // geo

#endif /* INCLUDE_PC_COMPUTER_HPP_ */
