/*
 * loess.hpp
 *
 *  Created on: May 23, 2018
 *      Author: rob
 */

#ifndef INCLUDE_LOESS_HPP_
#define INCLUDE_LOESS_HPP_

#include <mutex>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/Core>

#include "ds/kdtree.hpp"

#define MAX_INT std::numeric_limits<int>::max()

namespace geo {
namespace interp {

template <class T>
inline double __dist(const T& a, const T& b) {
	return std::sqrt(std::pow(a.x - b.x, 2.0) + std::pow(a.y - b.y, 2.0));
}

inline double __tricube(double dist, double bandwidth) {
	if(dist > bandwidth * 0.5) {
		return 0;
	} else {
		return std::pow(1 - std::pow(dist / bandwidth * 0.5, 3.0), 3.0);
	}
}

template <class T>
class Loess {
private:
	geo::ds::KDTree<T> m_tree;
	double m_bandwidth;
	std::mutex m_treeMtx;

public:

	/**
	 * @param data A list of T, where T is indexable; each element is a coordinate.
	 * @param bandwidth The radius of the neighbourhood of an estimate.
	 */
	Loess(const std::vector<T>& data, double bandwidth) :
		m_bandwidth(bandwidth) {
		std::lock_guard<std::mutex> lk(m_treeMtx);
		m_tree.add(data.begin(), data.end());
		m_tree.build();
	}

	/**
	 * Estimate the value at the given location.
	 *
	 * @param pt A query point.
	 */
	double estimate(const T& pt) {
		std::vector<T> found;
		std::vector<double> dist;
		int count;
		{
			std::lock_guard<std::mutex> lk(m_treeMtx);
			count = m_tree.radSearch(pt, m_bandwidth, 1024, std::back_inserter(found), std::back_inserter(dist));
		}
		if(count >= 3) {
			// TODO: Check for colinearity.

			Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mtx;
			Eigen::VectorXd vec;

			mtx.resize(count, 3);
			vec.resize(count);

			for(size_t i = 0; i < count; ++i) {
				const T& f = found[i];
				double d = __dist(f, pt);
				mtx(i, 0) = f[0];
				mtx(i, 1) = f[1];
				mtx(i, 2) = 1;
				vec[i] = __tricube(d, m_bandwidth) * f[2];
			}

			std::cerr << mtx.bdcSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(vec) << "\n";
		}
		return 0;
	}

};

} // interp
} // geo



#endif /* INCLUDE_LOESS_HPP_ */
