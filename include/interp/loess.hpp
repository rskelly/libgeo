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

std::mutex _mtx;

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
			Eigen::Matrix<double, 3, Eigen::Dynamic> mtx(3, count);

			for(size_t i = 0; i < count; ++i) {
				const T& f = found[i];
				double d = __dist(f, pt);
				mtx(0, i) = f[0];
				mtx(1, i) = f[1];
				mtx(2, i) = f[2];
			}

			// Recenter on zero.
			Eigen::Vector3d cent = mtx.rowwise().mean();
			mtx.colwise() -= cent;

			Eigen::JacobiSVD<Eigen::MatrixXd> svd(mtx, Eigen::ComputeFullU);
			//Eigen::Vector3d norm = mtx.bdcSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(vec).matrix().col(2);
			Eigen::Vector3d norm = svd.matrixU().col(2);
			double z = -norm.dot(cent);
			std::cerr << z << "\n";
			return z;
		}
		return 0;
	}

};

} // interp
} // geo



#endif /* INCLUDE_LOESS_HPP_ */
