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
	private:
		int m_scale;
	public:
		itemsort(int scale) : m_scale(scale) {
		}

		bool operator()(const T& a, const T& b) {
			// Scale each coordinate so it's between 0 and maxint.
			static uint32_t max = -1;
			uint32_t ax = (uint32_t) (a.x() / m_scale);
			uint32_t ay = (uint32_t) (a.y() / m_scale);
			uint32_t bx = (uint32_t) (b.x() / m_scale);
			uint32_t by = (uint32_t) (b.y() / m_scale);
			return morton(ax, ay) < morton(bx, by);
		}
	};

	template <class T>
	inline double dist(const T& a, const T& b) {
		return std::pow(a[0] - b[0], 2.0) + std::pow(a[1] - b[1], 2.0);
	}

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
	mvector<size_t> m_index;
	size_t m_minMort;
	size_t m_maxMort;
	int m_scale;

public:

	mqtree<T>(int scale) :
		m_minMort(-1), m_maxMort(0), m_scale(scale) {
	}

	mvector<T>& items() {
		return m_items;
	}

	void build() {
		g_trace("Building tree...")
		g_trace("Sorting points...")
		itemsort<T> sorter(m_scale);
		m_items.sort(sorter);
		g_trace('Sorted.')

		g_trace("Building index...")
		m_minMort = -1;
		m_maxMort = 0;
		T item;
		size_t mort, lastMort = -1;
		uint32_t sx, sy;
		std::vector<size_t> items;
		items.reserve(MEM_LIMIT / sizeof(size_t));
		for(size_t i = 0; i < m_items.size(); ++i) {
			m_items.get(i, item);
			sx = (uint32_t) (item.x() / m_scale);
			sy = (uint32_t) (item.y() / m_scale);
			mort = morton(sx, sy);
			if(mort != lastMort) {
				items.push_back(mort);
				items.push_back(i);
				lastMort = mort;
				if(mort < m_minMort) m_minMort = mort;
				if(mort > m_maxMort) m_maxMort = mort;
				if(items.size() >= MEM_LIMIT / sizeof(size_t)) {
					m_index.push(items, items.size());
					items.clear();
				}
			}
		}
		if(!items.empty())
			m_index.push(items, items.size());
		g_trace("Index built.")
		g_trace("Finished building tree.")
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
		uint32_t a = (uint32_t) (x / m_scale);
		uint32_t b = (uint32_t) (y / m_scale);
		uint64_t mort = morton(a, b);
		return mort >= m_minMort && mort <= m_maxMort;
 	}

	size_t size() const {
		return m_items.size();
	}

	template <class TIter>
	size_t knn(const T& pt, size_t n, TIter iter) {

		return 0;
	}

	size_t findMortIdxHigh(size_t mort, size_t start = -1, size_t end = -1) {
		if(start == (size_t) -1) {
			start = 0;
			end = m_items.size();
		}
		start -= start % 2;
		end -= end % 2;
		size_t tmp, mtmp;
		if(end - start <= 2) {
			m_index.get(end, mtmp);
			m_index.get(end + 1, tmp);
			return tmp;
		} else {
			size_t mid = (start + end) / 2;
			m_index.get(mid, tmp);
			if(tmp == mort) {
				m_index.get(mid + 1, tmp);
				return tmp;
			} else if(tmp > mort) {
				return findMortIdxHigh(mort, start, mid);
			} else {
				return findMortIdxHigh(mort, mid + 2, end);
			}
		}
	}

	size_t findMortIdxLow(size_t mort, size_t start = -1, size_t end = -1) {
		if(start == (size_t) -1) {
			start = 0;
			end = m_items.size();
		}
		start -= start % 2;
		end -= end % 2;
		size_t tmp, mtmp;
		if(end - start <= 2) {
			m_index.get(start, mtmp);
			m_index.get(start + 1, tmp);
			return tmp;
		} else {
			size_t mid = (start + end) / 2;
			mid -= mid % 2;
			m_index.get(mid, tmp);
			if(tmp == mort) {
				m_index.get(mid + 1, tmp);
				return tmp;
			} else if(tmp > mort) {
				return findMortIdxLow(mort, start, mid);
			} else {
				return findMortIdxLow(mort, mid + 2, end);
			}
		}
	}


	template <class TIter, class DIter>
	size_t search(const T& pt, double radius, TIter piter, DIter diter) {

		size_t count = 0;
		uint32_t x1 = (uint32_t) ((pt[0] - radius) / m_scale);
		uint32_t y1 = (uint32_t) ((pt[1] - radius) / m_scale);
		uint32_t x2 = (uint32_t) ((pt[0] + radius) / m_scale);
		uint32_t y2 = (uint32_t) ((pt[1] + radius) / m_scale);
		size_t m1 = morton(x1, y1);
		size_t m2 = morton(x2, y2);

		// The orders might be reversed if it's like UTM and the coords are upside down.
		if(m1 > m2) {
			size_t tmp = m1;
			m1 = m2;
			m2 = tmp;
		}

		size_t start = findMortIdxLow(m1);
		size_t end = findMortIdxHigh(m2);
		/*
		// Get the indices one before and one after the Morton code for thecoordinate.
		auto sit = m_index.lower_bound(m1);
		auto eit = m_index.upper_bound(m2);

		// If the end iterator is at the end, backit up.
		if(eit == m_index.end())
			--eit;

		// Get the start and end indices. If they overlap, create a gap of 1.
		size_t start = sit->second;
		size_t end = eit->second;
		*/

		// Retrieve the items from the mvector.
		static std::vector<T> items;
		items.resize(end - start + 1);
		m_items.get(start, items, items.size());

		// Check the distances from the query point and add to iterators.
		double d;
		for(const T& item : items) {
			if((d = dist(item, pt)) < radius * radius) {
				*piter = item;
				*diter = d;
				++piter;
				++diter;
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
