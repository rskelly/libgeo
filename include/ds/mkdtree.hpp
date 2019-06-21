/*
 * This is wrapper around the ANN KDTree.
 *
 *  Created on: Nov 17, 2017
 *      Author: rob
 */

#ifndef INCLUDE_DS_KDTREE_HPP_
#define INCLUDE_DS_KDTREE_HPP_

#include <sys/mman.h>

#include <cstdlib>
#include <vector>
#include <set>
#include <algorithm>
#include <iterator>
#include <list>

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

	bool get(size_t idx, T& item) const {
		if(idx < m_idx) {
			std::memcpy(&item, m_data + idx, sizeof(T));
			return true;
		}
		return false;
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

template <class T>
class ItemSort {
public:
	T pt;
	ItemSort(T& pt) : pt(pt) {
	}

	bool operator<(const T& a, const T& b) {
		double d1 = std::pow(a[0] - pt[0], 2.0) + std::pow(a[1] - pt[1], 2.0);
		double d2 = std::pow(b[0] - pt[0], 2.0) + std::pow(b[1] - pt[1], 2.0);
		return d1 < d2;
	}
};

template <class T>
class node {
public:
	constexpr static int maxDepth = 100;

	double bounds[4] { std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest() };
	double midx, midy;
	bool leaf;
	int depth;
	int dims;
	size_t start;
	size_t end;
	mvector<T>* items;

	node<T>* nodes[0] = { nullptr };
	node<T>* parent;

	node(mvector<T>* items, int dims) :
		node(0, 0, items->size(), items, dims) {}

	node(int depth, size_t start, size_t end, mvector<T>* items, int dims) :
		midx(0), midy(0), leaf(false), depth(depth), dims(dims), start(start), end(end), items(items), parent(nullptr) {
	}

	bool intersects(double box[4]) {
		return !(box[0] > bounds[2] && box[2] < bounds[0] && box[1] > bounds[3] && box[3] < bounds[1]);
	}

	bool contains(double x, double y) {
		return x >= bounds[0] && x <= bounds[2] && y >= bounds[1] && y <= bounds[3];
	}

	void build() {

		if(depth == maxDepth) {
			leaf = true;
			return;
		}

		T item;
		for(size_t i = start; i < end; ++i) {
			items->get(i, item);
			if(item[0] < bounds[0]) bounds[0] = item[0];
			if(item[0] > bounds[2]) bounds[2] = item[0];
			if(item[1] < bounds[1]) bounds[1] = item[1];
			if(item[1] > bounds[3]) bounds[3] = item[1];

		}

		double midx = (bounds[2] + bounds[0]) / 2.0;
		double midy = (bounds[3] + bounds[1]) / 2.0;

		mvector<size_t> indices(end - start);
		size_t idx;
		size_t idxCounts[4] = {0};
		for(size_t i = start; i < end; ++i) {
			items->get(i, item);
			idx = ((item[0] < midx) << 1) | (item[1] < midy);
			indices->set(i, idx);
			idxCounts[idx]++;
		}

		T tmpItem;
		size_t curIdx = 0, tmpIdx, count = 0;
		for(size_t i = start; i < end; ++i) {
			indices->get(i, idx);
			if(idx != curIdx) {
				for(size_t j = i + 1; j < end; ++j) {
					indices->get(j, tmpIdx);
					if(idx != tmpIdx) {
						items->get(j, tmpItem);
						items->get(i, item);
						items->insert(i, tmpItem);
						items->insert(j, item);
					}
				}
			}
			if(++count > idxCounts[curIdx]) {
				++curIdx;
				count = 0;
			}
		}

		size_t s = 0, e;
		for(int i = 0; i < 4; ++i) {
			if(idxCounts[i] > 0) {
				e = s + idxCounts[i];
				nodes[i] = new node<T>(depth + 1, s, e, items, dims);
				nodes[i]->parent = this;
				s = e;
			}
		}
	}

	size_t size() const {
		return end - start;
	}

	template <class TIter>
	size_t knn(const T& pt, size_t n, TIter iter) {

		if(!contains(pt[0], pt[1])) return 0;

		node<T>* tmp = this;
		bool found;

		do {
			found = false;
			for(int i = 0; i < 4; ++i) {
				if(tmp->nodes[i] && tmp->nodes[i]->contains(pt[0], pt[1])) {
					tmp = tmp->nodes[i];
					found = true;
					break;
				}
			}
		} while(found);

		if(tmp->parent)
			tmp = tmp->parent;

		std::list<T> lst;
		T item;
		for(size_t i = tmp->start; i < tmp->end; ++i) {
			items->get(i, item);
			lst.push_back(std::move(item));
		}

		size_t count = 0;
		std::sort(lst.begin(), lst.end(), geo::ds::ItemSort(pt));
		for(T& item : tmp) {
			*iter = item;
			++iter;
			++count;
			if(--n == 0)
				break;
		}

		return count;
	}

	template <class TIter>
	size_t search(const T& pt, double radius, TIter iter) {
		return search(pt, radius, iter, {pt[0] - radius, pt[1] - radius, pt[0] + radius, pt[1] + radius});
	}

	template <class TIter>
	size_t search(const T& pt, double radius, TIter iter, double bounds[4]) {

		if(!intersects(bounds))
			return 0;

		double c0 = pt[depth % dims];
		double d0 = std::pow(c0 - median, 2);

		std::list<node<T>*> q;
		q.push_back(this);

		while(!q.empty()) {

			node<T>* n = q.front();
			q.pop_front();

			if(d0 < best.front().first)
				best.emplace_front(d0, n);

			if(c0 < median) {
				if(n->left)
					q.push_back(n->left);
				if(d0 < radius * radius && n->right)
					q.push_back(n->right);
			} else {
				if(n->right)
					q.push_back(n->right);
				if(d0 < radius * radius && n->left)
					q.push_back(n->left);
			}
		}

		size_t count = 0;
		for(const auto& it : best) {
			node<T>* n = it.second;
			for(size_t i = n->start; i < n->end; ++i) {
				*iter = items[i];
				++iter;
				++count;
			}
		}

		return count;
	}

	~node() {
		if(left) delete left;
		if(right) delete right;
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
class mKDTree {
private:
	mvector<T> m_items;
	node<T>* m_root;
	size_t m_dims;

public:

	/**
	 * Construct the KDTree with the given number of dimensions.
	 */
	mKDTree(size_t dims = 3) :
		m_root(nullptr),
		m_dims(dims) {
	}

	/**
	 * Destroy the tree.
	 *
	 * \param destroyItems If true, also destroys the stored items.
	 */
	void destroy(bool destroyItems = true) {
		delete m_root;
		m_root = nullptr;
		if(destroyItems)
			m_items.clear();
	}

	/**
	 * Add an item to the tree.
	 * @param item An item.
	 */
	void add(const T& item) {
		m_items.push(item);
	}

	/**
	 * Add the items to the tree.
	 * @param begin The start iterator.
	 * @param end The end iterator.
	 */
	template <class Iter>
	void add(Iter begin, Iter end) {
		while(begin != end) {
			m_items.push(*(*begin));
			++begin;
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
		destroy(false);

		if(empty())
			g_runerr("Not enough items.");

		m_root = new node<T>(&m_items, m_dims);
	}

	/**
	 * Returns the number of items in the tree.
	 */
	int size() const {
		return m_items.size();
	}

	bool empty() const {
		return size() == 0;
	}

	/**
	 * Returns a reference to the vector containing all
	 * coordinates added to the tree.
	 */
	const mvector<T*>& items() const {
		return m_items;
	}

	/**
	 * Perform a k-nearest neighbour search on the items. Returns
	 * the number of items that were found and appends the found items
	 * and distances to the given iterators.
	 * @param item The search item; of the same type as the tree items.
	 * @param count The number of items to return.
	 * @param titer A back_inserter for the found items.
	 */
	template <class TIter>
	size_t knn(const T& item, size_t count, TIter titer) const {
		return m_root->knn(item, count, titer);
	}

	/**
	 * Search the tree by radius with an upper bound on the number of elements returned. If
	 * The number of returned elements is equal to the bound, then there are probably
	 * more points than were returned.
	 * @param item The search point.
	 * @param radius The search radius.
	 * @param titer The output point iterator.
	 */
	template <class TIter>
	size_t search(const T& item, double radius, TIter titer = nullptr) const {
		return m_root->search(item, radius, titer);
	}

	int dims() const {
		return m_dims;
	}

	~mKDTree() {
		destroy();
	}
};


}
}



#endif /* INCLUDE_DS_KDTREE_HPP_ */
