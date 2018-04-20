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

template <class T>
class RBF {
private:
	bool m_built;
	std::vector<T> m_pts;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_A;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_c;

	double tps(double val, double c = 1) {
		return c * c * val * val * std::log(c * val);
	}

	double dist(const T& a, const T& b) {
		return std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2) + 0.00000001;
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
				double z = tps(dist(a, b));
				m_A(i, j) = z;
			}
		}

		for(size_t i = 0; i < size; ++i) {
			m_A(size, i) = 1;
			m_A(i, size) = 1;
		}

		m_A(size, size) = 0;

		m_built = true;
	}

public:

	RBF() :
		m_built(false) {}

	template <class I>
	void add(I begin, I end) {
		m_built = false;
		m_pts.push_back(begin, end);
	}

	void add(const T& pt) {
		m_built = false;
		m_pts.push_back(pt);
	}

	double compute(const T& pt) {
		if(!m_built)
			build();

		size_t size = m_pts.size();

		for(size_t i = 0; i < size; ++i)
			m_c(i, 0) = tps(dist(pt, m_pts[i]));
		m_c(size, 0) = 1;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> d = m_A.inverse() * m_c;

		double z = 0;
		double m = d(size, 0);
		for(size_t i = 0; i < size; ++i)
			z += d(i, 0) * m_pts[i][2];

		return z;
	}



};


#endif /* INCLUDE_RBF_HPP_ */
