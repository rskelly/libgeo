#include <vector>

#include <CGAL/Plane_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_2_algorithms.h>

#include "pointcloud.hpp"
#include "pc_computer.hpp"

using namespace geo::pc;

bool pointSort(const geo::pc::Point& a, const geo::pc::Point& b) {
	return a.value() < b.value();
}

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

int Computer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

IDWComputer::IDWComputer(double exponent) : exponent(exponent) {}

int IDWComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int IDWComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	if(!filtered.empty()) {
		double result = 0;
		double div = 0;
		for(const geo::pc::Point& pt : filtered) {
			double d = std::pow(x - pt.x(), 2.0) + std::pow(y - pt.y(), 2.0);
			if(d == 0) {
				result = pt.z();
				div = 1;
				break;
			} else {
				double w = std::pow(1 / d, exponent);
				result += w * pt.z();
				div += w;
			}
		}
		out.push_back(result / div);
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int IDWComputer::bandCount() const {
	return 1;
}

int MeanComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int MeanComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
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

VarianceComputer::VarianceComputer(double bias) : 
	m_bias(bias) {}

void VarianceComputer::setBias(double bias) {
	m_bias = bias;
}

int VarianceComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int VarianceComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	std::vector<double> _out; // TODO: Istance variable.
	MeanComputer::compute(x, y, pts, filtered, radius, _out);
	double mean = _out[0];
	if(!std::isnan(mean) && !filtered.empty()){
		double sum = 0;
		for(const geo::pc::Point& pt : filtered)
			sum += std::pow(pt.value() - mean, 2.0);
		out.push_back(sum / (filtered.size() + m_bias));
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int VarianceComputer::bandCount() const {
	return 1;
}

StdDevComputer::StdDevComputer(double bias) :
		VarianceComputer(bias) {
}

void StdDevComputer::setBias(double bias) {
	VarianceComputer::setBias(bias);
}

int StdDevComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int StdDevComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	std::vector<double> _out; // TODO: Instance variable.
	VarianceComputer::compute(x, y, pts, filtered, radius, _out);
	double variance = _out[0];
	if(std::isnan(variance) || filtered.empty()) {
		out.push_back(std::nan(""));
	} else {
		out.push_back(std::sqrt(variance));
	}
	return 1;
}

int StdDevComputer::bandCount() const {
	return 1;
}

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

int PercentileComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
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

int MaxComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int MaxComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
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

int MinComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int MinComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
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

int DensityComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int DensityComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	if(!filtered.empty()) {
		double area = radius * radius * M_PI;
		out.push_back(pts.size() / area);
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int DensityComputer::bandCount() const {
	return 1;
}

int CountComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int CountComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	out.push_back(filtered.size());
	return 1;
}

int CountComputer::bandCount() const {
	return 1;
}

SkewComputer::SkewComputer(double bias) {
	stdDevComp.setBias(bias);
}

int SkewComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
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

KurtosisComputer::KurtosisComputer(double bias) {
	stdDevComp.setBias(bias);
}

int KurtosisComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
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

int CoVComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
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

HLRGBiometricsComputer::HLRGBiometricsComputer(int bands, int minCount, double threshold) :
		m_bands(bands),
		m_minCount(minCount),
		m_threshold(threshold) {
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
		hlrgCCF(pts, _pts, m_bands, m_threshold, out);
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

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K> Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3 Point_3;
typedef K::Plane_3 Plane_3;
typedef Delaunay::Finite_faces_iterator Finite_faces_iterator;
typedef Delaunay::Face Face;

/**
 * Compute planar area.
 */
double computePArea(double x1, double y1, double z1, double x2, double y2,
		double z2, double x3, double y3, double z3) {
	double side0 = std::sqrt(std::pow(x1 - x2, 2.0) + std::pow(y1 - y2, 2.0) + std::pow(z1 - z2, 2.0));
	double side1 = std::sqrt(std::pow(x2 - x3, 2.0) + std::pow(y2 - y3, 2.0) + std::pow(z2 - z3, 2.0));
	double side2 = std::sqrt(std::pow(x3 - x1, 2.0) + std::pow(y3 - y1, 2.0) + std::pow(z3 - z1, 2.0));
	double s = (side0 + side1 + side2) / 2.0;
	return std::sqrt(s * (s - side0) * (s - side1) * (s - side2));
}

/**
 * Compute the area of a face.
 */
double computeFArea(const Face &face) {
	Point_3 p1 = face.vertex(0)->point();
	Point_3 p2 = face.vertex(1)->point();
	Point_3 p3 = face.vertex(2)->point();
	return computePArea(p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z(), p3.x(),
			p3.y(), p3.z());
}

double toPlane(const Point_3 &p, const Plane_3 &plane, const Point_3 &centroid) {
	return (p.x() * plane.a() + p.y() * plane.b() + plane.d()) / -plane.c();
}

double polyArea(const std::list<Point_3> &hull, const Plane_3 &plane,
		const Point_3 &centroid) {
	double area = 0.0;
	auto it0 = hull.begin();
	auto it1 = hull.begin();
	it1++;
	do {
		double z0 = toPlane(*it0, plane, centroid);
		double z1 = toPlane(*it1, plane, centroid);
		area += computePArea(it0->x(), it0->y(), z0, it1->x(), it1->y(), z1,
				centroid.x(), centroid.y(), centroid.z());
		it0++;
		it1++;
		if (it1 == hull.end())
			it1 = hull.begin();
	} while (it0 != hull.end());
	return area;
}

int RugosityComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter) {
	if(filter) {
		std::vector<geo::pc::Point> filtered;
		pointFilter(pts.begin(), pts.end(), std::back_inserter(filtered), filter);
		return compute(x, y, pts, filtered, radius, out);
	} else {
		return compute(x, y, pts, pts, radius, out);
	}
}

int RugosityComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	
	if(filtered.size() >= 3) {

		double area = radius * radius * M_PI;
		double density = filtered.size() / area;

		std::list<Point_3> verts;
		for (const geo::pc::Point& pt : filtered)
			verts.emplace_back(pt.x(), pt.y(), pt.value());

		// Convex hull and POBF
		std::list<Point_3> hull;
		Plane_3 plane;
		Point_3 centroid;
		CGAL::convex_hull_2(verts.begin(), verts.end(), std::back_inserter(hull), Gt());
		CGAL::linear_least_squares_fitting_3(hull.begin(), hull.end(), plane, centroid, CGAL::Dimension_tag<0>());

		// POBF surface area.
		double parea = polyArea(hull, plane, centroid);

		// If the poly area is zero, quit.
		if(parea <= 0) {
			out.push_back(std::nan(""));
		} else {

			// Delaunay 3D surface area.
			double tarea = 0.0;
			Delaunay dt(verts.begin(), verts.end());
			for (Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it)
				tarea += computeFArea(*it);

			// TODO: This is an attempt at modelling the relationship between the ACR and
			// density. The fractal dimension is involved. Should be redone and documented.
			double densityFactor = 1.0 / (2.49127261 + 9.01659384 * std::sqrt(density * 32.65748276));

			double acr = (tarea / parea) * densityFactor;

			out.push_back(acr);
		}
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int RugosityComputer::bandCount() const {
	return 1;
}

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
