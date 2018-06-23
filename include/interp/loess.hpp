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

#define MAX_INT std::numeric_limits<int>::max()

namespace geo {
namespace interp {

inline double __sqdist(double x0, double y0, double x1, double y1) {
	return (x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1);
}

template <class T>
class QNode {
private:
	int m_binSize;
	bool m_split;
	const std::vector<T>* m_pts;
	double m_minx;
	double m_maxx;
	double m_miny;
	double m_maxy;
	QNode* m_a;
	QNode* m_b;
	QNode* m_c;
	QNode* m_d;
	std::vector<size_t> m_indices;

	QNode* getNode(int idx) {
		switch(idx) {
		case 0: return m_a ? m_a : (m_a = new QNode(m_binSize, m_pts, m_minx, m_miny, m_minx + (m_maxx - m_minx) / 2.0, m_miny + (m_maxy - m_miny) / 2.0));
		case 1: return m_b ? m_b : (m_b = new QNode(m_binSize, m_pts, m_minx + (m_maxx - m_minx) / 2.0, m_miny, m_maxx, m_miny + (m_maxy - m_miny) / 2.0));
		case 2: return m_d ? m_d : (m_d = new QNode(m_binSize, m_pts, m_minx, m_miny + (m_maxy - m_miny) / 2.0, m_minx + (m_maxx - m_minx) / 2.0, m_maxy));
		case 3: return m_c ? m_c : (m_c = new QNode(m_binSize, m_pts, m_minx + (m_maxx - m_minx) / 2.0, m_miny + (m_maxy - m_miny) / 2.0, m_maxx, m_maxy));
		default:
			throw std::runtime_error("Unknown node index.");
		}
	}

	void split() {
		double midx = m_minx + (m_maxx - m_minx) / 2.0;
		double midy = m_miny + (m_maxy - m_miny) / 2.0;
		for(size_t idx : m_indices) {
			const T& pt = m_pts->at(idx);
			int nidx = (((int) (pt.y >= midy)) << 1) | ((int) (pt.x >= midx));
			getNode(nidx)->m_indices.push_back(idx);
		}
		m_indices.clear();
		m_split = true;
	}


public:
	QNode(int binSize, const std::vector<T>* pts, double minx, double miny, double maxx, double maxy) :
		m_binSize(binSize), m_split(false), m_pts(pts),
		m_minx(minx), m_maxx(maxx), m_miny(miny), m_maxy(maxy),
		m_a(nullptr), m_b(nullptr), m_c(nullptr), m_d(nullptr) {
	}

	void build() {
		for(size_t i = 0; i < m_pts->size(); ++i) {
			const T& pt = m_pts->at(i);
			if(pt.x >= m_minx && pt.x < m_maxx && pt.y >= m_miny && pt.y < m_maxy)
				m_indices.push_back(i);
		}
		if(m_binSize < m_indices.size())
			split();
	}

	template <class Iter>
	int find(const T& pt, double radius, Iter out) {
		if((pt.x + radius) < m_minx || (pt.x - radius) > m_maxx || (pt.y + radius) < m_miny || (pt.y - radius) > m_miny)
			return 0;
		int count = 0;
		if(!m_split) {
			radius *= radius;
			for(size_t idx : m_indices) {
				const T& pt0 = m_pts->at(idx);
				double dist = __sqdist(pt0.x, pt0.y, pt.x, pt.y);
				if(dist <= radius) {
					*out = idx;
					++count;
				}
			}
		} else {
			if(m_a)
				count += m_a->find(pt, radius, out);
			if(m_b)
				count += m_b->find(pt, radius, out);
			if(m_c)
				count += m_c->find(pt, radius, out);
			if(m_d)
				count += m_d->find(pt, radius, out);
		}
		return count;
	}

	~QNode() {
		if(m_a) delete m_a;
		if(m_b) delete m_b;
		if(m_c) delete m_c;
		if(m_d) delete m_d;
	}

};

template <class T>
class QTree {
private:
	QNode<T>* m_node;

public:
	QTree() : m_node(nullptr) {}

	QTree(const std::vector<T>& pts, int binSize = 100) : QTree() {
		build(pts, binSize);
	}

	void build(const std::vector<T>& pts, int binSize = 100) {
		double xmin = 99999999, ymin = 99999999;
		double xmax = -99999999, ymax = -99999999;
		for(const T& pt : pts) {
			if(pt.x < xmin) xmin = pt.x;
			if(pt.x > xmax) xmax = pt.x;
			if(pt.y < ymin) ymin = pt.y;
			if(pt.y > ymax) ymax = pt.y;
		}
		m_node = new QNode<T>(binSize, &pts, xmin - 0.0001, ymin - 0.0001, xmax + 0.0001, ymax + 0.0001);
		m_node->build();
	}

	template <class Iter>
	int find(const T& pt, double radius, Iter out) {
		return m_node->find(pt, radius, out);
	}

	~QTree() {
		if(m_node) delete m_node;
	}

};

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
	QTree<T> m_tree;
	std::vector<T> m_data;
	double m_bandwidth;

public:

	/**
	 * @param data A list of T, where T is indexable; each element is a coordinate.
	 * @param bandwidth The radius of the neighbourhood of an estimate.
	 */
	Loess(const std::vector<T>& data, double bandwidth) :
		m_bandwidth(bandwidth) {
		m_data = data;
		m_tree.build(m_data, 100);
	}

	/**
	 * Estimate the value at the given location.
	 *
	 * @param pt A query point.
	 */
	double estimate(const T& pt) {
		std::vector<size_t> found;
		int count = m_tree.find(pt, m_bandwidth, std::back_inserter(found));
		if(count >= 3) {
			Eigen::Matrix<double, 3, Eigen::Dynamic> mtx(3, count);

			bool equal = true;
			double z = 0;
			for(size_t i = 0; i < count; ++i) {
				const T& f = m_data[found[i]];
				double d = __dist(f, pt);
				double t = __tricube(d, m_bandwidth);
				z = f[2];
				mtx(0, i) = f.x;
				mtx(1, i) = f.y;
				mtx(2, i) = t * z;
				if(i > 0 && z != m_data[found[i - 1]].z)
					equal = false;
			}

			// All the elevations are equal, just return the elevation.
			if(equal) {
				//std::cerr << "skip\n";
				return z;
			}

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
