/*
 * This is wrapper around the ANN KDTree.
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

/**
 * A 3D point class that can be used in the KDTree.
 */
class KDPoint {
private:
	double m_x, m_y, m_z;

public:

	KDPoint(double x, double y, double z = 0) :
		m_x(x), m_y(y), m_z(z) {
	}

	double x() const {
		return m_x;
	}

	double y() const {
		return m_y;
	}

	double z() const {
		return m_z;
	}

	double operator[](int idx) const {
		switch(idx % 3) {
		case 0: return x();
		case 1: return y();
		default: return z(); // Who cares.
		}
	}
};

/**
 * A KDTree.
 *
 * It is assumed that the template is a class with a [] operator
 * which will return a value corresponding to the dimension,
 * which is the modulus of the given index.
 */
template <class T>
class KDTree {
private:
	std::vector<T> m_items;
	ANNpointArray m_pts;
	ANNkd_tree* m_tree;
	size_t m_dims;

public:

	/**
	 * Construct the KDTree with the given number of dimensions.
	 */
	KDTree(size_t dims) :
		m_pts(nullptr),
		m_tree(nullptr),
		m_dims(dims) {
	}

	/**
	 * Destroy the tree.
	 */
	void destroy() {
		if(m_tree) {
			delete m_tree;
			m_tree = nullptr;
		}
		m_items.clear();
		if(m_pts) {
			annDeallocPts(m_pts);
			m_pts = nullptr;
		}
	}

	/**
	 * Add an item to the tree.
	 * @param item An item.
	 */
	void add(T& item) {
		m_items.push_back(item);
	}

	/**
	 * Add the items to the tree.
	 * @param begin The start iterator.
	 * @param end The end iterator.
	 */
	template <class Iter>
	void add(Iter begin, Iter end) {
		m_items.insert(m_items.begin(), begin, end);
	}

	/**
	 * Build the tree. This destroys the existing tree and
	 * attempts to rebuild it using the existing list of items.
	 * If there are no items, just destroys the existing tree
	 * and does nothing.
	 */
	void build() {

		// Clean up existing tree, etc.
		if(m_tree)
			destroy();

		if(m_items.empty())
			g_runerr("Not enough items.");

		// Set up the points buffer.
		m_pts = annAllocPts(m_items.size(), m_dims);

		// Write points into the buffer and save pointers.
		for(size_t i = 0; i < m_items.size(); ++i) {
			for(size_t j = 0; j < m_dims; ++j)
				m_pts[i][j] = m_items[i][j];
		}

		// Set up the tree.
		m_tree = new ANNkd_tree(m_pts, m_items.size(), m_dims);
	}

	/**
	 * Returns the number of items in the tree.
	 */
	int size() const {
		return m_items.size();
	}

	/**
	 * Returns a reference to the vector containing all
	 * items added to the tree.
	 */
	const std::vector<T>& items() const {
		return m_items;
	}

	/**
	 * Perform a k-nearest neighbour search on the items. Returns
	 * the number of items that were found and appends the found items
	 * and distances to the given iterators.
	 * @param item The search item; of the same type as the tree items.
	 * @param count The number of items to return.
	 * @param titer A back_inserter for the found items.
	 * @param diter A back_inserter for the item distances.
	 * @param eps The allowable error bound.
	 */
	template <class TIter, class DIter>
	int knn(const T& item, size_t count, TIter titer, DIter diter, double eps = 0.0) {

		if(!m_tree) {
			g_warn("Tree not built. Forget to call build?");
			return 0;
		}

		if(m_items.empty())
			g_runerr("Count too small: " << count);

		// If the count is higher than the number of items, shrink it.
		if(count > m_items.size())
			count = m_items.size();

		// Turn the search item into an array of doubles castable to ANNpoint.
		std::vector<ANNcoord> pt(m_dims);
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
			if(idx[i] == -1)
				return i;
			*titer = m_items[idx[i]];
			++titer;
			*diter = std::sqrt(dist[i]);
			++diter;
		}
		return count;
	}

	/**
	 * Search the tree by radius with an upper bound on the number of elements returned. If
	 * The number of returned elements is equal to the bound, then there are probably
	 * more points than were returned.
	 * @param item The search point.
	 * @param radius The search radius.
	 * @param maxCount The maximum number of items to return.
	 * @param titer The output point iterator.
	 * @param diter The distance iterator.
	 * @param eps The error bound.
	 */
	template <class TIter, class DIter>
	int radSearch(const T& item, double radius, int maxCount, TIter titer, DIter diter, double eps = 0.0) {

		if(!m_tree) {
			g_warn("Tree not built. Forget to call build?");
			return 0;
		}

		if(radius < 0)
			g_runerr("Radius too small: " << radius);
		if(maxCount < 0)
			g_runerr("Max count too small: " << maxCount);

		// Turn the search item into an array of doubles castable to ANNpoint.
		std::vector<ANNcoord> pt(m_dims);
		for(size_t i = 0; i < m_dims; ++i)
			pt[i] = item[i];

		// Create arrays for indices and distances.
		std::vector<ANNidx> idx(maxCount);
		std::vector<ANNdist> dist(maxCount);

		// Perform search.
		m_tree->annkFRSearch(static_cast<ANNpoint>(pt.data()), radius * radius, maxCount,
				static_cast<ANNidxArray>(idx.data()), static_cast<ANNdistArray>(dist.data()), eps);

		// Populate output iterators.
		for(size_t i = 0; i < maxCount; ++i) {
			if(idx[i] == -1)
				return i;
			*titer = m_items[idx[i]];
			++titer;
			*diter = std::sqrt(dist[i]);
			++diter;
		}

		return maxCount;
	}

	~KDTree() {
		destroy();
	}
};


}
}



#endif /* INCLUDE_DS_KDTREE_HPP_ */
