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

// Arcgis (http://www.spatialanalysisonline.com/HTML/index.html?radial_basis_and_spline_functi.htm)
double tps0(double r, double c) {
	return c * c * r * r * std::log(c * r);
}

double tps1(double r, double c) {
	return c * c * r * r * std::log(c * c + r * r);
}

double multiquad(double r, double c) {
	return std::sqrt(r * r + c * c);
}

double invmultiquad(double r, double c) {
	return 1.0 / std::sqrt(r * r + c * c);
}

template <class T>
double dist(const T& a, const T& b) {
	return std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2) + 0.00000001;
}

template <class T>
class RBF {
private:
	bool m_built;
	double m_smoothing;
	double (*m_rbf)(double, double);
	std::vector<T> m_pts;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_A;
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

	template <class I>
	void add(I begin, I end) {
		m_built = false;
		m_pts.insert(m_pts.begin(), begin, end);
	}

	void add(const T& pt) {
		m_built = false;
		m_pts.push_back(pt);
	}

	// http://www.spatialanalysisonline.com/HTML/index.html?radial_basis_and_spline_functi.htm
	void build() {

		size_t size = m_pts.size();

		m_A.resize(size + 1, size + 1);
		m_c.resize(size + 1, 1);

		for(size_t i = 0; i < size; ++i) {
			const T& a = m_pts[i];
			for(size_t j = 0; j < size; ++j) {
				const T& b = m_pts[j];
				double z = (*m_rbf)(dist(a, b), m_smoothing);
				m_A(i, j) = z;
			}
		}

		for(size_t i = 0; i < size; ++i) {
			m_A(size, i) = 1;
			m_A(i, size) = 1;
		}

		m_A(size, size) = 0;
		m_Ai = m_A.inverse();
		m_built = true;
	}

	void clear() {
		m_pts.clear();
		m_built = false;
	}

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


#endif /* INCLUDE_RBF_HPP_ */
