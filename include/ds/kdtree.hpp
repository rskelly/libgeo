/*
 * This is wrapper around the ANN KDTree.
 *
 *  Created on: Nov 17, 2017
 *      Author: rob
 */

#ifndef INCLUDE_DS_KDTREE_HPP_
#define INCLUDE_DS_KDTREE_HPP_

#include <sys/mman.h>

#include <vector>
#include <set>
#include <algorithm>
#include <iterator>

#include "ANN/ANN.h"

#include "geo.hpp"
#include "util.hpp"

#define EPS std::numeric_limits<double>::epsilon()

using namespace geo::util;

namespace geo {
namespace ds {

template <class T>
class mvector {
private:
	T* m_data;
	std::unique_ptr<TmpFile> m_file;
	size_t m_idx;
public:
	mvector(size_t size = 1) :
		m_data(nullptr), m_idx(0) {
		resize(size);
	}

	void resize(size_t size) {
		if(size < 1) size = 1;
		size *= sizeof(T);
		size_t oldSize = 0;
		if(!m_file.get()) {
			m_file.reset(new TmpFile(size));
		} else if(size > m_file->size) {
			oldSize = m_file->size;
			m_file->resize(size);
		} else { return; }
		if(!m_data) {
			m_data = (T*) mmap(0, size, PROT_READ|PROT_WRITE, MAP_SHARED, m_file->fd, 0);
		} else {
			m_data = (T*) mremap(m_data, oldSize, size, MREMAP_MAYMOVE, 0);
		}
	}

	bool push(const T& item) {
		if(m_idx * sizeof(T) >= m_file->size)
			resize((m_idx + 1) *  2);
		std::memcpy(m_data + m_idx, &item, sizeof(T)) ;
		++m_idx;
		return true;
	}

	bool insert(size_t idx, const T& item) {
		if(idx >= m_idx)
			return false;
		std::memcpy(m_data + idx, &item, sizeof(T)) ;
		return true;
	}

	T operator[](size_t idx) const {
		if(idx < m_idx) {
			T item;
			std::memcpy(&item, m_data + idx, sizeof(T));
			return item;
		} else {
			throw std::runtime_error("Index out of range.");
		}
	}

	void clear() {
		m_idx = 0;
	}

	T* data() {
		return m_data;
	}

	bool empty() const {
		return m_idx == 0;
	}

	size_t size() const {
		return m_idx;
	}

	~mvector() {
		munmap(m_data, m_file->size);
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
	// Props for in-core use.
	std::vector<T*> m_items;
	std::vector<ANNcoord> m_vcoords;

	// Props for out-of-core use.
	mvector<T> m_mitems;
	mvector<ANNcoord> m_mcoords;

	ANNkd_tree* m_tree;
	ANNpoint* m_pts;
	ANNcoord* m_coords;
	size_t m_dims;
	bool m_mapped;

public:

	/**
	 * Construct the KDTree with the given number of dimensions.
	 */
	KDTree(size_t dims = 3, bool mapped = false) :
		m_tree(nullptr),
		m_pts(nullptr),
		m_coords(nullptr),
		m_dims(dims),
		m_mapped(mapped) {
	}

	/**
	 * Destroy the tree.
	 *
	 * \param destroyItems If true, also destroys the stored items.
	 */
	void destroy(bool destroyItems = true) {
		if(m_tree) {
			delete m_tree;
			m_tree = nullptr;
		}
		if(m_mapped) {
			if(destroyItems)
				m_mitems.clear();
			m_mcoords.clear();
		} else {
			if(destroyItems) {
				for(T* item : m_items)
					delete item;
				m_items.clear();
			}
			m_vcoords.clear();
		}
	}

	/**
	 * Add an item to the tree.
	 * @param item An item.
	 */
	void add(T* item) {
		if(m_mapped) {
			add(*item);
		} else {
			m_items.push_back(item);
		}
	}

	void add(const T& item) {
		m_mitems.push(item);
	}

	/**
	 * Add the items to the tree.
	 * @param begin The start iterator.
	 * @param end The end iterator.
	 */
	template <class Iter>
	void add(Iter begin, Iter end) {
		if(m_mapped) {
			while(begin != end) {
				m_mitems.push(*(*begin));
				++begin;
			}
		} else {
			while(begin != end) {
				m_items.push_back(*begin);
				++begin;
			}
		}
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
			destroy(false);

		if(empty())
			g_runerr("Not enough items.");

		// The code below does the same as annAllocPts: the m_pts
		// array stores pointers to the coordinates in m_items
		// Set up the points buffer.
		if(m_mapped) {
			m_mcoords.resize(m_mitems.size() * m_dims);
			m_coords = (ANNcoord*) m_mcoords.data();
			m_pts = (ANNpoint*) m_coords;
		} else {
			m_vcoords.resize(m_items.size() * m_dims);
			m_coords = (ANNcoord*) m_vcoords.data();
			m_pts = (ANNpoint*) m_coords;
		}

		// Write points into the buffer and save pointers.
		if(m_mapped) {
			for(size_t i = 0, j = 0; i < m_mitems.size(); ++i, j += m_dims) {
				T item = m_mitems[i];
				for(size_t k = 0; k < m_dims; ++k)
					m_coords[j + k] = item[k];
			}
		} else {
			for(size_t i = 0, j = 0; i < m_items.size(); ++i, j += m_dims) {
				for(size_t k = 0; k < m_dims; ++k)
					m_coords[j + k] = (*m_items[i])[k];
			}
		}

		// Set up the tree.
		m_tree = new ANNkd_tree(static_cast<ANNpointArray>(m_pts), size(), m_dims, 128, ANN_KD_SUGGEST);
	}

	/**
	 * Returns the number of items in the tree.
	 */
	int size() const {
		return m_mapped ? m_mitems.size() : m_items.size();
	}

	bool empty() const {
		return size() == 0;
	}

	/**
	 * Returns a reference to the vector containing all
	 * coordinates added to the tree.
	 */
	const std::vector<T*>& items() const {
		return m_items;
	}

	/**
	 * Returns a reference to the vector containing all
	 * coordinates added to the tree.
	 */
	const size_t items(std::vector<T>& items) const {
		if(m_mapped) {
			for(size_t i = 0; i < m_mitems.size(); ++i)
				items.push_back(std::move(m_mitems[i]));
			return m_mitems.size();
		} else {
			for(size_t i = 0; i < m_items.size(); ++i)
				items.push_back(*m_items[i]);
			return m_items.size();
		}
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
	int knn(const T& item, size_t count, TIter titer, DIter diter, double eps = EPS) const {
		//std::lock_guard<std::mutex> lk(m_mtx);

		if(!m_tree) {
			g_warn("Tree not built. Forget to call build?");
			return 0;
		}

		if(m_items.empty() || m_mitems.empty())
			g_runerr("Count too small: " << count);

		// If the count is higher than the number of items, shrink it.
		if(m_mapped) {
			if(count > m_mitems.size())
				count = m_mitems.size();
		} else {
			if(count > m_items.size())
				count = m_items.size();
		}

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
		size_t ct = 0;
		for(size_t i = 0; i < count; ++i) {
			if(idx[i] > 0) {
				if(m_mapped) {
					*titer = m_items[idx[i]];
					*diter = dist[i];
					++ct;
				} else {
					*titer = *m_items[idx[i]];
					*diter = dist[i];
					++ct;
				}
			}
		}
		return ct;
	}

	/**
	 * Search the tree by radius with an upper bound on the number of elements returned. If
	 * The number of returned elements is equal to the bound, then there are probably
	 * more points than were returned.
	 * @param item The search point.
	 * @param radius The search radius.
	 * @param maxCount The number of items to return.
	 * @param titer The output point iterator.
	 * @param diter The distance iterator.
	 * @param eps The error bound.
	 */
	template <class TIter, class DIter>
	int radSearch(const T& item, double radius, int count, TIter titer = nullptr, DIter diter = nullptr, double eps = EPS) const {
		//std::lock_guard<std::mutex> lk(m_mtx);

		if(!m_tree) {
			g_warn("Tree not built. Forget to call build?");
			return 0;
		}

		if(radius < 0)
			g_runerr("Radius too small: " << radius);

		// Turn the search item into an array of doubles castable to ANNpoint.
		std::vector<ANNcoord> pt(m_dims);
		for(size_t i = 0; i < m_dims; ++i)
			pt[i] = item[i];

		// Create arrays for indices and distances.
		std::vector<ANNidx> idx(count);
		std::vector<ANNdist> dist(count);

		if(count == 0) {
			count = m_tree->annkFRSearch(static_cast<ANNpoint>(pt.data()), radius, 0);

			idx.resize(count);
			dist.resize(count);
		}

		// Perform search.
		int res = m_tree->annkFRSearch(static_cast<ANNpoint>(pt.data()), radius, count,
				static_cast<ANNidxArray>(idx.data()), static_cast<ANNdistArray>(dist.data()), eps);

		if(res < count)
			count = res;

		// Populate output iterators.
		for(size_t i = 0; i < count; ++i) {
			if(idx[i] != ANN_NULL_IDX) {
				if(m_mapped) {
					*titer = m_mitems[idx[i]];
					*diter = dist[i];
				} else {
					*titer = *m_items[idx[i]];
					*diter = dist[i];
				}
			} else {
				*diter = -1;
			}
		}

		return count;
	}

	int dims() const {
		return m_dims;
	}

	~KDTree() {
		destroy();
	}
};


}
}



#endif /* INCLUDE_DS_KDTREE_HPP_ */
