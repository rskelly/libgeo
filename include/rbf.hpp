/*
 * rbf.hpp
 *
 *  Created on: Apr 19, 2018
 *      Author: rob
 */

#ifndef INCLUDE_RBF_HPP_
#define INCLUDE_RBF_HPP_

#include <vector>
#include <mutex>

#include <Eigen/Core>
#include <Eigen/LU>

#include "ds/kdtree.hpp"
#include "kmeans.hpp"

#define PI 3.141592653589793
#define INF std::numeric_limits<double>::infinity()
#define EPS std::numeric_limits<double>::epsilon()


using namespace geo::ds;


std::mutex __mtx;


// Arcgis (http://www.spatialanalysisonline.com/HTML/index.html?radial_basis_and_spline_functi.htm)
double tps0(double r, double sigma, double smoothing) {
	return smoothing * smoothing * r * r * std::log(smoothing * r);
}

double tps1(double r, double sigma, double smoothing) {
	return (smoothing * smoothing + r * r) * std::log(smoothing * smoothing + r * r);
}

double multiquad(double r, double sigma, double smoothing) {
	return std::sqrt(r * r + smoothing * smoothing);
}

double invmultiquad(double r, double sigma, double smoothing) {
	return 1.0 / std::sqrt(r * r + smoothing * smoothing);
}

double gaussian(double r, double sigma, double smoothing) {
	return (1.0 / (sigma * std::sqrt(2 * PI))) * std::exp(-0.5 * std::pow(r / sigma, 2));
	//return std::exp(-r*r);
}

template <class T>
double dist(const T& a, const T& b) {
	return std::sqrt(std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2));
}


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> RBFMatrix;

template <class T>
class RBF {
private:
	bool m_built;
	double m_range;
	double m_sigma;
	double m_smoothing;
	int m_clusters;
	int m_samples;

	double (*m_rbf)(double, double, double);
	
	std::vector<T> m_pts;
	// Cluster centres.
	std::vector<T> m_clusterMeans;

	RBFMatrix m_Ai;
	RBFMatrix m_c;
	KDTree<T> m_tree;


	double computeR(const T& a, const T& b) {
		double r = dist(a, b) / m_range;
		return r <=0 ? EPS : r;
	}

public:

	enum Type {
		MultiQuadratic,
		InvMultiQuadratic,
		ThinPlateSpline,
		ThinPlateSplineAlt,
		Gaussian
		//MultiLog,
		//NaturalCubicSpline
	};

	RBF(Type type) :
		m_built(false),
		m_range(0),
		m_sigma(0),
		m_smoothing(0),
		m_clusters(0),
		m_samples(0),
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
		case Gaussian:
			m_rbf = &gaussian;
			break;
		default:
			g_argerr("Unknown type: " << type);
		}
	}

	void setRange(double range) {
		m_range = range;
	}

	double range() const {
		return m_range;
	}

	void setSigma(double sigma) {
		m_sigma = sigma;
	}

	double sigma() const {
		return m_sigma;
	}

	void setSmoothing(double smoothing) {
		m_smoothing = smoothing;
	}

	double smoothing() const {
		return m_smoothing;
	}

	void setSamples(int samples) {
		m_samples = samples;
	}

	int samples() const {
		return m_samples;
	}

	void setClusters(int clusters) {
		m_clusters = clusters;
	}

	int clusters() const {
		return m_clusters;
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

		m_clusterMeans.clear();

		if(m_samples && m_clusters)
			g_runerr("Only use one of samples or clusters.")
		if(m_clusters < 0)
			g_runerr("Clusters must be greater than zero.")
		if(m_samples < 0)
			g_runerr("Samples must be greater than zero.")

		if(m_clusters) {
			std::unordered_map<size_t, std::list<T> > clusters;
			kmeans(m_pts, m_clusters, m_clusterMeans, clusters);
		} else if(m_samples) {
			m_clusterMeans.assign(m_pts.begin(), m_pts.begin() + std::min(m_pts.size(), (size_t) m_samples));
		} else {
			m_clusterMeans.assign(m_pts.begin(), m_pts.end());
		}

		size_t size = m_clusterMeans.size();

		m_tree.destroy();
		m_tree.add(m_clusterMeans.begin(), m_clusterMeans.end());
		m_tree.build();

		g_debug("building matrices")
		{
			Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A(size + 1, size + 1);
			m_c.resize(size + 1, 1);

			for(size_t i = 0; i < size; ++i) {
				const T& a = m_clusterMeans[i];
				for(size_t j = 0; j < size; ++j) {
					const T& b = m_clusterMeans[j];
					double r = computeR(a, b);
					double z = (*m_rbf)(r, m_sigma, m_smoothing);
					A(i, j) = z;
				}
			}

			for(size_t i = 0; i < size; ++i) {
				A(size, i) = 1;
				A(i, size) = 1;
			}

			A(size, size) = 0;
			m_Ai = A.inverse();
		}
		g_debug("done building")

		m_built = true;
	}

	void clear() {
		m_pts.clear();
		m_clusterMeans.clear();
		m_built = false;
	}

	double compute(const T& pt) {
		if(!m_built)
			build();

		size_t size = m_clusterMeans.size();

		int count = 0;
		double dst = INF;
		{
			std::vector<T> pts;
			std::vector<double> dsts;
			std::lock_guard<std::mutex> lk(__mtx);
			if((count = m_tree.knn(pt, 1, std::back_inserter(pts), std::back_inserter(dsts))))
				dst = dsts[0];
		}

		if(!count || dst > m_range * 3)
			return 0;

		RBFMatrix c(size + 1, 1);
		c(size, 0) = 1;

		double r, z;
		size_t i = 0;
		for(const T& p : m_clusterMeans) {
			r = computeR(p, pt);
			z = (*m_rbf)(r, m_sigma, m_smoothing);
			c(i++, 0) = z;
		}

		RBFMatrix b = m_Ai * c;

		z = 0;
		i = 0;
		for(const T& p : m_clusterMeans)
			z += b(i++, 0) * p.z;

		return z;
	}



};


#endif /* INCLUDE_RBF_HPP_ */
