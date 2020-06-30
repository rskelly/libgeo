/*
 * stats.cpp
 *
 *  Created on: Apr 3, 2018
 *      Author: rob
 */


#include <vector>

#include "pointcloud.hpp"
#include "pc_computer.hpp"

using namespace geo::pc::compute;

int MeanComputer::compute(double, double, const std::vector<geo::pc::Point>&, const std::vector<geo::pc::Point>& filtered, double, std::vector<double>& out) {
	if(!filtered.empty()) {
		double sum = 0;
		for(const geo::pc::Point& pt : filtered)
			sum += pt.value();
		out.push_back(sum / filtered.size());
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int MeanComputer::bandCount() const {
	return 1;
}

std::vector<std::pair<std::string, std::string>> MeanComputer::bandMeta() const {
	return {{"mean", "Mean"}};
}

VarianceComputer::VarianceComputer(double bias) :
	m_bias(bias) {}

void VarianceComputer::setBias(double bias) {
	m_bias = bias;
}

int VarianceComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	if(filtered.size() > 1) {
		std::vector<double> _out; // TODO: Instance variable.
		MeanComputer::compute(x, y, pts, filtered, radius, _out);
		double mean = _out[0];
		if(!std::isnan(mean)){
			double sum = 0;
			for(const geo::pc::Point& pt : filtered)
				sum += std::pow(pt.value() - mean, 2.0);
			out.push_back(sum / (filtered.size() + m_bias));
			return 1;
		}
	}
	out.push_back(std::nan(""));
	return 1;
}

int VarianceComputer::bandCount() const {
	return 1;
}

std::vector<std::pair<std::string, std::string>> VarianceComputer::bandMeta() const {
	return {{"variance", "Variance"}};
}

StdDevComputer::StdDevComputer(double bias) :
		VarianceComputer(bias) {
}

void StdDevComputer::setBias(double bias) {
	VarianceComputer::setBias(bias);
}

int StdDevComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	if(filtered.size() > 1) {
		std::vector<double> _out; // TODO: Instance variable.
		VarianceComputer::compute(x, y, pts, filtered, radius, _out);
		double variance = _out[0];
		if(!std::isnan(variance)) {
			out.push_back(std::sqrt(variance));
			return 1;
		}
	}
	out.push_back(std::nan(""));
	return 1;
}

int StdDevComputer::bandCount() const {
	return 1;
}

std::vector<std::pair<std::string, std::string>> StdDevComputer::bandMeta() const {
	return {{"stddev", "Standard Deviation"}};
}

int MaxComputer::compute(double, double, const std::vector<geo::pc::Point>&, const std::vector<geo::pc::Point>& filtered, double, std::vector<double>& out) {
	if(!filtered.empty()) {
		double max = DBL_MIN;
		for(const geo::pc::Point& pt : filtered) {
			if(pt.value() > max)
				max = pt.value();
		}
		out.push_back(max);
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int MaxComputer::bandCount() const {
	return 1;
}

std::vector<std::pair<std::string, std::string>> MaxComputer::bandMeta() const {
	return {{"max", "Maximum"}};
}

int MinComputer::compute(double, double, const std::vector<geo::pc::Point>&, const std::vector<geo::pc::Point>& filtered, double, std::vector<double>& out) {
	if(!filtered.empty()) {
		double min = DBL_MAX;
		for(const geo::pc::Point& pt : filtered) {
			if(pt.value() < min)
				min = pt.value();
		}
		out.push_back(min);
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int MinComputer::bandCount() const {
	return 1;
}

std::vector<std::pair<std::string, std::string>> MinComputer::bandMeta() const {
	return {{"min", "Minimum"}};
}

int CountComputer::compute(double, double, const std::vector<geo::pc::Point>&, const std::vector<geo::pc::Point>& filtered, double, std::vector<double>& out) {
	out.push_back(filtered.size());
	return 1;
}

int CountComputer::bandCount() const {
	return 1;
}

std::vector<std::pair<std::string, std::string>> CountComputer::bandMeta() const {
	return {{"count", "Count"}};
}

SkewComputer::SkewComputer(double bias) {
	stdDevComp.setBias(bias);
}

int SkewComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	if(!filtered.empty()) {
		// Fisher-Pearson
		std::vector<double> _out; // TODO: Instance variable.
		meanComp.compute(x, y, pts, filtered, radius, _out);
		double mean = _out[0];
		_out.clear();

		if(std::isnan(mean))
			return mean;

		double sum = 0.0;
		for (const geo::pc::Point& pt : filtered)
			sum += std::pow(pt.value() - mean, 3.0) / filtered.size();
		stdDevComp.compute(x, y, pts, filtered, radius, _out);
		double sd = _out[0];
		double skew = sum / std::pow(sd, 3.0);
		out.push_back(skew);
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int SkewComputer::bandCount() const {
	return 1;
}

std::vector<std::pair<std::string, std::string>> SkewComputer::bandMeta() const {
	return {{"skew", "Skewness"}};
}

KurtosisComputer::KurtosisComputer(double bias) {
	stdDevComp.setBias(bias);
}

int KurtosisComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	if(!filtered.empty()) {
		// Fisher-Pearson
		std::vector<double> _out;
		meanComp.compute(x, y, pts, filtered, radius, _out);
		double mean = _out[0];
		_out.clear();

		if(std::isnan(mean))
			return mean;

		double sum = 0.0;
		for (const geo::pc::Point& pt : filtered)
			sum += std::pow(pt.value() - mean, 4.0) / filtered.size();
		stdDevComp.compute(x, y, pts, filtered, radius, _out);
		double sd = _out[0];
		double kurt = sum / std::pow(sd, 4.0) - 3.0;
		out.push_back(kurt);
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int KurtosisComputer::bandCount() const {
	return 1;
}

std::vector<std::pair<std::string, std::string>> KurtosisComputer::bandMeta() const {
	return {{"kurtosis", "Kurtosis"}};
}

int CoVComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	if(!filtered.empty()) {
		// Fisher-Pearson
		std::vector<double> _out;
		meanComp.compute(x, y, pts, filtered, radius, _out);
		double mean = _out[0];
		_out.clear();

		if(std::isnan(mean))
			return mean;

		sampStdDevComp.compute(x, y, pts, filtered, radius, _out);
		double sd = _out[0];
		double cov = mean != 0 ? sd / mean : std::nan("");
		out.push_back(cov);
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int CoVComputer::bandCount() const {
	return 1;
}

std::vector<std::pair<std::string, std::string>> CoVComputer::bandMeta() const {
	return {{"cov", "Coefficient of Variance"}};
}
