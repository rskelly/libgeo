/*
 * rbf.hpp
 *
 *  Created on: Apr 19, 2018
 *      Author: rob
 */

#ifndef INCLUDE_RBF_HPP_
#define INCLUDE_RBF_HPP_

#include <vector>
#include <Eigen/Core>
#include <Eigen/LU>

#define EPS std::numeric_limits<double>::epsilon()

/**
 * Calculate the thin plate spline per
 * Arcgis (http://www.spatialanalysisonline.com/HTML/index.html?radial_basis_and_spline_functi.htm)
 * @param r The norm.
 * @param c The smoothing parameter.
 */
double tps0(double r, double c) {
	return c * c * r * r * std::log(c * r);
}

/**
 * Calculate the thin plate spline per
 * Surfer (http://www.spatialanalysisonline.com/HTML/index.html?radial_basis_and_spline_functi.htm)
 * @param r The norm.
 * @param c The smoothing parameter.
 */
double tps1(double r, double c) {
	return c * c * r * r * std::log(c * c + r * r);
}

/**
 * Calculate the multiquadratic per
 * http://www.spatialanalysisonline.com/HTML/index.html?radial_basis_and_spline_functi.htm
 * @param r The norm.
 * @param c The smoothing parameter.
 */
double multiquad(double r, double c) {
	return std::sqrt(r * r + c * c);
}

/**
 * Calculate the inverse multiquadratic per
 * http://www.spatialanalysisonline.com/HTML/index.html?radial_basis_and_spline_functi.htm
 * @param r The norm.
 * @param c The smoothing parameter.
 */
double invmultiquad(double r, double c) {
	return 1.0 / std::sqrt(r * r + c * c);
}

/**
 * Calculate the distance between points represented by T.
 * Implementations of T must returns coordinates from the []
 * operator: 0 for x, 1 for y and 2 for z.
 * @param a A point.
 * @param b A point.
 */
template <class T>
double dist(const T& a, const T& b) {
	return std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2) + EPS;
}

namespace geo {
namespace interp {

/**
 * The radial basis function calculator.
 */
template <class T>
class RBF {
private:
	bool m_built;
	double m_smoothing;
	double (*m_rbf)(double, double);
	std::vector<T> m_pts;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_Ai;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_c;

public:

	enum Type {
		MultiQuadratic,
		InvMultiQuadratic,
		ThinPlateSpline,
		ThinPlateSplineAlt,
		//MultiLog,
		//NaturalCubicSpline
	};

	/**
	 * Build a radial basis function caclulator of the given type
	 * with the given smoothing. Smoothing is specific to the chosen method:
	 * for some, 0 implies no smoothing, for others it causes a discontinuity.
	 * @param type The basis function.
	 * @param smoothing The smoothing parameter.
	 */
	RBF(Type type, double smoothing) :
		m_built(false),
		m_smoothing(smoothing),
		m_rbf(nullptr) {

		switch(type) {
		case MultiQuadratic:
			m_rbf = &multiquad;
			break;
		case InvMultiQuadratic:
			m_rbf = &invmultiquad;
			break;
		case ThinPlateSpline:
			m_rbf = &tps0;
			break;
		case ThinPlateSplineAlt:
			m_rbf = &tps1;
			break;
		default:
			g_argerr("Unknown type: " << type);
		}
	}

	/**
	 * Add points to the calculator.
	 * @param begin A beginning iterator into a container of points.
	 * @param end An ending iterator.
	 */
	template <class I>
	void add(I begin, I end) {
		m_built = false;
		m_pts.insert(m_pts.begin(), begin, end);
	}

	/**
	 * Add a point to the calculator.
	 * @param pt A point.
	 */
	void add(const T& pt) {
		m_built = false;
		m_pts.push_back(pt);
	}

	/**
	 * Build the matrices required to calculate the RBF per
	 * http://www.spatialanalysisonline.com/HTML/index.html?radial_basis_and_spline_functi.htm.
	 */
	void build() {

		size_t size = m_pts.size();

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A(size + 1, size + 1);

		m_c.resize(size + 1, 1);

		for(size_t i = 0; i < size; ++i) {
			const T& a = m_pts[i];
			for(size_t j = 0; j < size; ++j) {
				const T& b = m_pts[j];
				double z = (*m_rbf)(dist(a, b), m_smoothing);
				A(i, j) = z;
			}
		}

		for(size_t i = 0; i < size; ++i) {
			A(size, i) = 1;
			A(i, size) = 1;
		}

		A(size, size) = 0;
		m_Ai = A.inverse();
		m_built = true;
	}

	/**
	 * Clear the points list and reset. A rebuild is required.
	 */
	void clear() {
		m_pts.clear();
		m_built = false;
	}

	/**
	 * Calculate the radial basis function at the given point.
	 * @param pt A point.
	 * @return The value of the function at this point.
	 */
	double compute(const T& pt) {
		if(!m_built)
			build();

		size_t size = m_pts.size();

		for(size_t i = 0; i < size; ++i)
			m_c(i, 0) = (*m_rbf)(dist(pt, m_pts[i]), m_smoothing);
		m_c(size, 0) = 1;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> d = m_Ai * m_c;

		double z = 0;
		for(size_t i = 0; i < size; ++i)
			z += d(i, 0) * m_pts[i][2];

		return z;
	}



};

} // interp
} // geo

#endif /* INCLUDE_RBF_HPP_ */
