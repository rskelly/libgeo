/*
 * This is wrapper around the ANN KDTree.
 *
 *  Created on: Nov 17, 2017
 *      Author: rob
 */

#ifndef INCLUDE_DS_KDTREE_HPP_
#define INCLUDE_DS_KDTREE_HPP_

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>

#include "geo.hpp"
#include "util.hpp"
#include "ds/mvector.hpp"

using namespace geo::util;

namespace {

	using namespace geo::ds;

	template <class T>
	class itemsort {
	public:
		int scale;
		itemsort(int scale) : scale(scale) {
		}

		bool operator()(const T& a, const T& b) {
			return morton(a.x(), a.y(), scale) < morton(b.x(), b.y(), scale);
		}
	};

	double NULL_BOUNDS[] = {
			std::numeric_limits<double>::max(),
			std::numeric_limits<double>::max(),
			std::numeric_limits<double>::lowest(),
			std::numeric_limits<double>::lowest()
	};
} // anon


namespace geo {
namespace ds {


/**
 * A file-backed qtree.
 *
 * It is assumed that the template is a class with a [] operator
 * which will return a value corresponding to the dimension,
 * which is the modulus of the given index.
 */
template <class T>
class mqtree {
private:

	mvector<T> m_items;
	std::map<size_t, size_t> m_index;
	size_t m_minMort;
	size_t m_maxMort;
	int m_scale;

public:

	mqtree<T>(int scale = 1) :
		m_minMort(-1), m_maxMort(0), m_scale(scale) {
	}

	mvector<T>& items() {
		return m_items;
	}

	void build() {
		itemsort<T> sorter(m_scale);
		m_items.sort(sorter);
		m_index.clear();
		m_minMort = -1;
		m_maxMort = 0;
		T item;
		size_t mort, lastMort = -1;
		for(size_t i = 0; i < m_items.size(); ++i) {
			m_items.get(i, item);
			mort = morton(item.x(), item.y(), m_scale);
			if(mort != lastMort) {
				m_index[mort] = i;
				lastMort = mort;
				if(mort < m_minMort) m_minMort = mort;
				if(mort > m_maxMort) m_maxMort = mort;
			}
		}
	}

	void add(const T& item) {
		m_items.push(item);
	}

	void add(std::vector<T>& items) {
		m_items.push(items);
	}

	void clear() {
		m_items.clear();
		m_index.clear();
	}

	bool contains(double x, double y) {
		return m_index.find(morton(x, y, m_scale)) != m_index.end();
 	}

	size_t size() const {
		return m_items.size();
	}

	template <class TIter>
	size_t knn(const T& pt, size_t n, TIter iter) {

		return 0;
	}

	inline double dist(const T& a, const T& b) {
		return std::pow(a[0] - b[0], 2.0) + std::pow(a[1] - b[1], 2.0);
	}

	template <class TIter>
	size_t search(const T& pt, double radius, TIter iter) {

		size_t count = 0;

		std::map<size_t, size_t>::iterator sit = m_index.upper_bound(morton(pt[0] - radius, pt[1] - radius, m_scale) - 1);
		size_t start = sit->second;
		std::map<size_t, size_t>::iterator eit = m_index.lower_bound(morton(pt[0] + radius, pt[1] + radius, m_scale) + 1);
		if(eit == m_index.end())
			--eit;
		size_t end = eit->second;
		if(end <= start)
			end = start + 1;
		static std::vector<T> items;
		items.resize(end - start);
		m_items.get(start, items, end - start);
		for(const T& item : items) {
			if(dist(item, pt) < radius * radius) {
				*iter = item;
				++iter;
				++count;
			}
		}
		return count;
	}

	~mqtree() {
	}


};

}
}



#endif /* INCLUDE_DS_KDTREE_HPP_ */
