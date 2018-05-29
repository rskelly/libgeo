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

inline double __dist(double xa, double ya, double xb, double yb) {
	return std::sqrt(std::pow(xa - xb, 2.0) + std::pow(ya - yb, 2.0));
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
		std::vector<double> found;
		std::vector<double> dist;
		int count;
		{
			//std::lock_guard<std::mutex> lk(m_treeMtx);
			count = m_tree.radSearch(pt, m_bandwidth, 99999, std::back_inserter(found), std::back_inserter(dist));
		}
		if(count >= 3) {
			Eigen::Matrix<double, 3, Eigen::Dynamic> mtx(3, count);

			int equal = 0;
			int dims = m_tree.dims();
			double z = std::numeric_limits<double>::quiet_NaN();
			for(size_t i = 0; i < count; ++i) {
				double x0 = found[i * dims];
				double y0 = found[i * dims + 1];
				double z0 = found[i * dims + 2];
				double d = __dist(x0, y0, pt.x, pt.y);
				double t = __tricube(d, m_bandwidth);
				mtx(0, i) = x0;
				mtx(1, i) = y0;
				mtx(2, i) = t * z0;
				if(z == z0)
					++equal;
				z = z0;
			}

			if(equal == count - 1)
				return z;

			// Recenter on zero.
			Eigen::Vector3d cent = mtx.rowwise().mean();
			mtx.colwise() -= cent;
			Eigen::JacobiSVD<Eigen::MatrixXd> svd(mtx, Eigen::ComputeThinU | Eigen::ComputeThinV);
			Eigen::Vector3d norm = svd.matrixU().col(2);
			z = norm.dot(Eigen::Vector3d(0, 0, cent[2]));
			z = cent[2] > 0 ? std::abs(z) : -std::abs(z);
			/*
			{
				std::lock_guard<std::mutex> lk(_mtx);
				//std::cerr << "solve " << slv << "\n";
				std::cerr << "cent " << cent << "\n";
				std::cerr << "mtx " << mtx << "\n";
				std::cerr << "norm " << norm << "\n";
				std::cerr << "z " << z << "\n";
			}
			*/
			return z;
		}
		return 0;
	}

};

} // interp
} // geo



#endif /* INCLUDE_LOESS_HPP_ */
