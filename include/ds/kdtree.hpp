/*
 * kdtree.hpp
 *
 *  Created on: Nov 17, 2017
 *      Author: rob
 */

#ifndef INCLUDE_DS_KDTREE_HPP_
#define INCLUDE_DS_KDTREE_HPP_

#include <vector>
#include <set>
#include <algorithm>
#include <iterator>

#include "ANN/ANN.h"

#include "geo.hpp"
#include "util.hpp"

using namespace geo::util;

namespace geo {
namespace ds {


template <class T>
class KDTree {
private:
	std::vector<T> m_items;
	std::vector<ANNpoint> m_pts;
	Buffer m_buf;
	ANNkd_tree* m_tree;
	size_t m_dims;

public:

	KDTree(size_t dims) :
		m_tree(nullptr),
		m_dims(dims) {
	}

	void add(T& item) {
		m_items.push_back(item);
	}

	template <class Iter>
	void add(Iter begin, Iter end) {
		m_items.insert(m_items.begin(), begin, end);
	}

	void build() {

		if(m_items.size() < 1)
			g_runerr("Not enough items.");

		// Clean up existing tree, etc.
		if(m_tree)
			destroy();

		// Set up a buffer for the point data, and an array of pointers into the buffer.
		m_buf.resize(m_items.size() * sizeof(double) * m_dims);
		m_pts.resize(m_items.size());

		// Write points into the buffer and save pointers.
		double* data = static_cast<double*>(m_buf.buf);
		for(size_t i = 0; i < m_items.size(); ++i) {
			for(size_t j = 0; j < m_dims; ++j)
				data[i + j] = m_items[i][j];
			m_pts[i] = data + i;
		}
		// Set up the tree.
		m_tree = new ANNkd_tree(static_cast<ANNpointArray>(m_pts.data()), m_items.size(), m_dims);
	}

	void destroy() {
		delete m_tree;
		m_tree = nullptr;
		m_buf.resize(0);
		m_items.clear();
		m_pts.clear();
	}

	template <class TIter, class DIter>
	void knn(const T& item, size_t count, TIter titer, DIter diter, double eps = 0.0) {

		if(!m_tree)
			g_runerr("Tree not build. Forget to call build?");
		if(count > m_items.size())
			count = m_items.size();
		if(count < 1)
			g_runerr("Count too small: " << count);

		// Turn the search item into an array of doubles castable to ANNpoint.
		std::vector<double> pt(m_dims);
		for(size_t i = 0; i < m_dims; ++i)
			pt[i] = item[i];

		// Create arrays for indices and distances.
		std::vector<ANNidx> idx(count);
		std::vector<ANNdist> dist(count);

		// Perform search.
		m_tree->annkSearch(static_cast<ANNpoint>(pt.data()), count, static_cast<ANNidxArray>(idx.data()),
				static_cast<ANNdistArray>(dist.data()), eps);

		// Populate output iterators.
		for(size_t i = 0; i < count; ++i) {
			*titer = m_items[idx[i]];
			++titer;
			*diter = std::sqrt(dist[i]);
			++diter;
		}
	}

	~KDTree() {
		destroy();
	}
};


}
}



#endif /* INCLUDE_DS_KDTREE_HPP_ */
