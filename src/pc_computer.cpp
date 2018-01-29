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

IDWComputer::IDWComputer(double exponent) : exponent(exponent) {}

double IDWComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	if(pts.empty())
		return std::nan("");
	double result = 0;
	double div = 0;
	for(const geo::pc::Point& pt : pts) {
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
	return result / div;
}

double MeanComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	double sum = 0;
	int count = 0;
	for(const geo::pc::Point& pt : pts) {
		sum += pt.value();
		++count;
	}
	return count > 0 ? sum / count : std::nan("");
}

VarianceComputer::VarianceComputer(double bias) : bias(bias) {}

double VarianceComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	double mean = MeanComputer::compute(x, y, pts, radius);
	if(std::isnan(mean))
		return mean;
	double sum = 0;
	for(const geo::pc::Point& pt : pts)
		sum += std::pow(pt.value() - mean, 2.0);
	return sum / pts.size();
}

StdDevComputer::StdDevComputer(double bias) :
		VarianceComputer(bias) {
}

double StdDevComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	double variance = VarianceComputer::compute(x, y, pts, radius);
	if(std::isnan(variance))
		return variance;
	return std::sqrt(variance);
}

bool pointSort(const geo::pc::Point& a, const geo::pc::Point& b) {
	return a.value() < b.value();
}

PercentileComputer::PercentileComputer(double percentile) :
	percentile(percentile) {
	if(percentile <= 0 || percentile >= 100)
		g_runerr("Percentile must be between 0 and 100.");
}

double PercentileComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	if(pts.empty())
		return std::nan("");
	std::vector<geo::pc::Point> _pts(pts.begin(), pts.end());
	std::sort(_pts.begin(), _pts.end(), pointSort);
	size_t size = _pts.size();
	if(size % 2 == 0) {
		size_t idx = (size_t) (size * percentile) - 1;
		return (_pts[idx + 1].value() + _pts[idx].value()) / 2.0;
	} else {
		return _pts[(size_t) (size * percentile)].value();
	}
}

double MaxComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	double max = std::numeric_limits<double>::lowest();
	for(const geo::pc::Point& pt : pts) {
		if(pt.value() > max)
			max = pt.value();
	}
	return max;
}

double MinComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	double min = std::numeric_limits<double>::max();
	for(const geo::pc::Point& pt : pts) {
		if(pt.value() < min)
			min = pt.value();
	}
	return min;
}

double DensityComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	double area = radius * radius * M_PI;
	if(pts.empty())
		return std::nan("");
	return pts.size() / area;
}

double CountComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	return (double) pts.size();
}

SkewComputer::SkewComputer(double bias) {
	stdDevComp.bias = bias;
}

double SkewComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	if(pts.empty())
		return std::nan("");
	// Fisher-Pearson
	double mean = meanComp.compute(x, y, pts, radius);
	if(std::isnan(mean))
		return mean;
	size_t count = pts.size();
	double sum = 0.0;
	for (const geo::pc::Point& pt : pts)
		sum += std::pow(pt.value() - mean, 3.0) / count;
	double sd = stdDevComp.compute(x, y, pts, radius);
	double skew = sum / std::pow(sd, 3.0);
	return skew;
}

KurtosisComputer::KurtosisComputer(double bias) {
	stdDevComp.bias = bias;
}

double KurtosisComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	if(pts.empty())
		return std::nan("");
	// Fisher-Pearson
	double mean = meanComp.compute(x, y, pts, radius);
	if(std::isnan(mean))
		return mean;
	size_t count = pts.size();
	double sum = 0.0;
	for (const geo::pc::Point& pt : pts)
		sum += std::pow(pt.value() - mean, 4.0) / count;
	double sd = stdDevComp.compute(x, y, pts, radius);
	double kurt = sum / std::pow(sd, 4.0) - 3.0;
	return kurt;
}

double CoVComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	if(pts.empty())
		return std::nan("");
	// Fisher-Pearson
	double mean = meanComp.compute(x, y, pts, radius);
	if(std::isnan(mean))
		return mean;
	double sd = sampStdDevComp.compute(x, y, pts, radius);
	double cov = mean != 0 ? sd / mean : std::nan("");
	return cov;
}

double CanopyGapFractionComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	g_runerr("Not implemented.");
}

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


double RugosityComputer::compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius) {
	if(pts.size() < 3)
		return std::nan("");

	double area = radius * radius * M_PI;
	double density = pts.size() / area;

	std::list<Point_3> verts;
	for (const geo::pc::Point& pt : pts)
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
	if(parea <= 0)
		return std::nan("");

	// Delaunay 3D surface area.
	double tarea = 0.0;
	Delaunay dt(verts.begin(), verts.end());
	for (Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it)
		tarea += computeFArea(*it);

	// TODO: This is an attempt at modelling the relationship between the ACR and
	// density. The fractal dimension is involved. Should be redone and documented.
	double densityFactor = 1.0 / (2.49127261 + 9.01659384 * std::sqrt(density * 32.65748276));

	double acr = (tarea / parea) * densityFactor;

	return acr;
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
