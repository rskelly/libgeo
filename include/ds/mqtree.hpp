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

	/**
	 * Instances of this class are used to sort points accorting to the Morton
	 * (z-curve) order. The bits of each cordinate are spread out, shifted and
	 * ored together to give a 1-dimensional order. Floating-point
	 * coordinates are scaled first to preserve some fixed amount of precision.
	 */
	template <class T>
	class itemsort {
	private:
		int m_scale;				///<! The scale for converting floating-point coordinates to ints.
		double m_minx, m_miny;		///<! The minimum x and y coordinate values, for shifting.

	public:
		/**
		 * Construct a sorter with the given corner coordinate and scale.
		 *
		 * \param minx The minimum x coordinate.
		 * \param miny The minimum y coordinate.
		 * \param scale A scale value for coorinates.
		 */
		itemsort(double minx, double miny, int scale) :
			m_minx(minx), m_miny(miny), m_scale(scale) {
			if(scale <= 0)
				g_runerr("Scale must be larger than zero: " << scale)
		}

		/**
		 * Return true if point a is less than point b in the Morton order.
		 *
		 * \param a A point.
		 * \param b A point.
		 * \return True if point a is less than point b in the Morton order.
		 */
		bool operator()(const T& a, const T& b) {
			return morton((a.x() - m_minx) * m_scale, (a.y() - m_miny) * m_scale) < morton((b.x() - m_minx) * m_scale, (b.y() - m_miny) * m_scale);
		}
	};

	/**
	 * Returns the squared distance between two points.
	 *
	 * \param a A point.
	 * \param b A point.
	 * \return The squared distance between two points.
	 */
	template <class T>
	inline double dist(const T& a, const T& b) {
		return std::pow(a[0] - b[0], 2.0) + std::pow(a[1] - b[1], 2.0);
	}

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

	mvector<T> m_items;				///<! A list of the stored items.
	mvector<size_t> m_index;		///<! An index that relates the Morton index to the first position in the items vector where points with that index are found.
	size_t m_minMort, m_maxMort;	///<! The minimum and maximum Morton indices found in the structure.
	double m_minx, m_miny;			///<! The minimum x and y coordinates. Used for shifting coordinates.
	int m_scale;					///<! The scale for converting floating-point coordinates to ints while preserving a fixed amount of precision.
	bool m_needsBuild;				///<! If true, the tree must be rebuilt before querying.

	/**
	 * Find the index associated with the Morton index or the next greater one.
	 *
	 * \param mort The Morton index.
	 * \param high Set true to get the next index the same or above the query, false to get the one the same or below.
	 * \param start A starting position for search. Defaults to 0.
	 * \param end An ending position for search. Defaults to the last index.
	 * \return The index into the items list.
	 */
	size_t findIndex(size_t mort, bool high, size_t start = -1, size_t end = -1) {
		if(start == (size_t) -1) {
			start = 0;
			end = m_index.size() - 1;
		}
		start -= start % 2;
		end -= end % 2;
		size_t tmp, mtmp;
		if(end - start <= 2) {
			if(high) {
				m_index.get(end, mtmp);
				m_index.get(end + 1, tmp);
			} else {
				m_index.get(start, mtmp);
				m_index.get(start + 1, tmp);
			}
			return tmp;
		} else {
			size_t mid = start + (end - start) / 2;
			mid -= mid % 2;
			m_index.get(mid, tmp);
			if(tmp == mort) {
				m_index.get(mid + 1, tmp);
				return tmp;
			} else if(tmp > mort) {
				return findIndex(mort, high, start, mid);
			} else {
				return findIndex(mort, high, mid, end);
			}
		}
	}

	uint32_t toX(double x) {
		return std::max(0.0, x - m_minx) * m_scale;
	}

	uint32_t toY(double y) {
		return std::max(0.0, y - m_miny) * m_scale;
	}

public:

	/**
	 * Construct an mqtree.
	 *
	 * \param scale The scale factor for decreasing the precision of the data.
	 * \param limit The memory threshold (bytes) that triggers the use of file-backed storage.
	 * \param minx The minimum x coordinate.
	 * \param miny The minimum y coordinate.
	 */
	mqtree<T>(int scale, size_t limit = 0, double minx = 0, double miny = 0) :
		m_minMort(-1), m_maxMort(0),
		m_minx(minx), m_miny(miny),
		m_scale(scale),
		m_needsBuild(true) {

		// The maximum number of index items is 2*the number of stored objects.
		// Get the proportion that will allow enough space in the index in the
		// worst case.
		float prop = (float) sizeof(T) / (2 * sizeof(size_t));
		if(prop > 1)
			prop = 1 / prop;
		m_items.setMemLimit((size_t) limit * prop);
		m_index.setMemLimit((size_t) limit * (1.0 - prop));
	}

	/**
	 * Return a reference to the item storage.
	 *
	 * \return A reference to the item storage.
	 */
	mvector<T>& items() {
		return m_items;
	}

	/**
	 * Build or rebuild the tree from the unsorted list of items.
	 */
	void build() {
		if(!m_needsBuild)
			return;
		g_trace("Building tree...")
		g_trace("Sorting points...")
		itemsort<T> sorter(m_minx, m_miny, m_scale);
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
			sx = toX(item.x());
			sy = toY(item.y());
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
		m_needsBuild = false;
	}

	/**
	 * Add a single item to the tree. Must (re)build before using.
	 *
	 * \param item An item.
	 */
	void add(const T& item) {
		m_items.push(item);
		m_needsBuild = true;
	}

	/**
	 * Add a vector of items to the tree. Must (re)build before using.
	 *
	 * \param items A vector of items.
	 */
	void add(std::vector<T>& items) {
		m_items.push(items);
		m_needsBuild = true;
	}

	/**
	 * Clear the tree and remove all items.
	 */
	void clear() {
		m_items.clear();
		m_index.clear();
		m_needsBuild = true;
	}

	/**
	 * Returns true if the tree contains the point within its bounds.
	 * If the radius is given, points within that distance are considered
	 * contained.
	 *
	 * \param x The x coordinate.
	 * \param y The y coordinate.
	 * \return True if the tree contains the point within the radius.
	 */
	bool contains(double x, double y, double radius = 0) {
		if(radius > 0)
			g_warn("Radius contains not implemented.")
		uint32_t a = toX(x);
		uint32_t b = toY(y);
		uint64_t mort = morton(a, b);
		return mort >= m_minMort && mort <= m_maxMort;
 	}

	/**
	 * Return the number of items in the tree.
	 *
	 * \return The number of items in the tree.
	 */
	size_t size() const {
		return m_items.size();
	}

	/**
	 * Return the n neares neighbours to the given point.
	 *
	 * \param pt A search point.
	 * \param n The number results.
	 * \param iter An insert iterator.
	 * \return The number of items found.
	 */
	template <class TIter>
	size_t knn(const T& pt, size_t n, TIter iter) {
		g_runerr("Not implemented.")
		return 0;
	}

	/**
	 * Search for points within the given radius.
	 *
	 * \param pt The search point.
	 * \param radius The search radius.
	 * \param piter The insert iterator for result points.
	 * \return The number of items found.
	 */
	template <class TIter>
	size_t search(const T& pt, double radius, TIter piter) {

		if(m_needsBuild)
			build();

		uint32_t x1 = toX(pt.x() - radius);
		uint32_t y1 = toY(pt.y() - radius);
		uint32_t x2 = toX(pt.x() + radius);
		uint32_t y2 = toY(pt.y() + radius);

		size_t count = 0;

		size_t m1 = morton(x1, y1);
		size_t m2 = morton(x2, y2);

		if(m1 > m2)
			g_runerr("Morton indices out of order. This is impossible.")

		//std::cout << m1 << ", " << m2 << ": " << (m2 - m1) << "\n";

		size_t start = findIndex(m1, false);
		size_t end = findIndex(m2, true);

		// Retrieve the items from the mvector.
		size_t maxSize = (1024 * 1024 * 10) / sizeof(T);
		if(end - start < maxSize) {
			std::vector<T> items(end - start + 1);
			if(!m_items.get(start, items, items.size()))
				return 0;

			// Check the distances from the query point and add to iterators.
			double d;
			for(const T& item : items) {
				if((d = dist(item, pt)) <= radius * radius) {
					*piter = item;
					++piter;
					++count;
				}
			}
		} else {
			std::vector<T> items(maxSize);
			for(size_t i = start; i <= end; i += maxSize) {
				size_t len = std::min(maxSize, end - i + 1);
				len = m_items.get(i, items, len);
				// Check the distances from the query point and add to iterators.
				double d;
				for(size_t j = 0; j < len; ++j) {
					const T& item = items[j];
					if((d = dist(item, pt)) < radius * radius) {
						*piter = item;
						++piter;
						++count;
					}
				}
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
