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

#define PI 3.141592653589793
#define INF std::numeric_limits<double>::infinity()

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
}

template <class T>
double dist(const T& a, const T& b) {
	return std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2) + 0.00000001;
}

template <class T>
void kmeans(std::vector<T>& pts, int n, 
	std::vector<T>& means, std::unordered_map<size_t, std::list<T> >& clusters) {

	if(n > pts.size())
		g_runerr("Too few points for meaningful clustering.")

	// Random shuffle and take the first n as cluster means.
	std::random_shuffle(pts.begin(), pts.end());
	means.assign(pts.begin(), pts.begin() + n);

	// Assign each point to the correct initial cluster.
	for(const T& p : pts) {
		double d0, d = INF;
		size_t idx = 0;
		for(size_t i = 0; i < means.size(); ++i) {
			if((d0 = dist(means[i], p)) < d) {
				d = d0;
				idx = i;
			}
		}
		clusters[idx].emplace_back(p);
	}

	// Iterate over the entire point set, assigning each point to 
	// the nearest mean.
	int changed;
	std::unordered_map<size_t, std::list<T> > clusters0;
	do {
		changed = 0;
		// Re-compute the means of the cluters and update the means list.
		g_debug("means")
		for(auto& pair : clusters) {
			double x = 0;
			double y = 0;
			for(const T& p : pair.second) {
				x += p[0];
				y += p[1];
			}
			// Update the 2d position in the means list.
			means[pair.first][0] = x / pair.second.size();
			means[pair.first][1] = y / pair.second.size();
		}
		// Iterate over the cluster member points and reassign to
		// the new cluster. Count each reassignment.
		g_debug("reassign")
		for(auto& pair : clusters) {
			for(const T& p : pair.second) {
				double d0, d = INF;
				size_t i = 0, idx = 0;
				for(const T& p0 : means) {
					if((d0 = dist(p0, p)) < d) {
						idx = i;
						d = d0;
					}
					++i;
				}
				// Transfer the point to the correct cluster.
				clusters0[idx].emplace_back(p);
				if(idx != pair.first)
					++changed;
			}
		}
		clusters.swap(clusters0);
		for(auto& pair : clusters0)
			pair.second.clear();
		g_debug("kmeans changed " << changed)
	} while(changed);
}

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> RBFMatrix;

template <class T>
class RBF {
private:
	bool m_built;
	double m_sigma;
	double m_smoothing;

	double (*m_rbf)(double, double, double);
	
	std::vector<T> m_pts;
	RBFMatrix m_Ai;
	RBFMatrix m_c;

	// Map for weights matrices for each cluster.
	std::unordered_map<size_t, RBFMatrix> m_clusterWeights;
	// Cluster centres.
	std::vector<T> m_clusterMeans;

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

	RBF(Type type, double sigma, double smoothing) :
		m_built(false),
		m_sigma(sigma),
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
		case Gaussian:
			m_rbf = &gaussian;
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

		g_debug("building matrices")
		{
			Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A(size + 1, size + 1);
			m_c.resize(size + 1, 1);

			for(size_t i = 0; i < size; ++i) {
				const T& a = m_pts[i];
				for(size_t j = 0; j < size; ++j) {
					const T& b = m_pts[j];
					double d = dist(a, b);
					double z = (*m_rbf)(d, m_sigma, m_smoothing);
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

		/*
		g_debug("clustering")
		{
			// Perform kmeans to get cluster centres.
			m_clusterMeans.clear();
			std::unordered_map<size_t, std::list<T> > clusters;
			kmeans(m_pts, 20, m_clusterMeans, clusters);

			// Calculate the weights for each cluster.
			for(size_t i = 0; i < m_clusterMeans.size(); ++i) {
				const T& p0 = m_clusterMeans[i];
				RBFMatrix c(size + 1, 1);
				size_t j = 0;
				for(const T& p : clusters[i]) {
					double d = dist(p0, p);
					double z = (*m_rbf)(d, m_sigma, m_smoothing);
					c(j++, 0) = z;
				}
				c(size, 0) = 1;
				// Save the cluster weights.
				m_clusterWeights[i] = m_Ai * c;
			}
		}
		g_debug("done clustering")
		*/

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
		/*
		size_t idx;
		double d0, d = INF;
		for(size_t i = 0; i < m_clusterMeans.size(); ++i) {
			if((d0 = dist(pt, m_clusterMeans[i])) < d) {
				d = d0;
				idx = i;
			}
		}

		RBFMatrix D = m_clusterWeights[idx];

		double z = 0;
		for(size_t i = 0; i < size; ++i)
			z += D(i, 0) * pt.z;
		*/

		RBFMatrix c(m_pts.size() + 1, 1);
		c(m_pts.size(), 0) = 1;

		double d, z;
		size_t i = 0;
		for(const T& p : m_pts) {
			d = dist(p, pt);
			z = (*m_rbf)(d, m_sigma, m_smoothing);
			c(i++, 0) = z;
		}

		RBFMatrix b = m_Ai * c;

		z = 0;
		i = 0;
		for(const T& p : m_pts)
			z += b(i++, 0) * p.z;

		return z;
	}



};


#endif /* INCLUDE_RBF_HPP_ */
