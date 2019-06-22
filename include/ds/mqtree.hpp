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
#include <list>

#include <fstream>

#include "geo.hpp"
#include "util.hpp"
#include "ds/mvector.hpp"

using namespace geo::util;

namespace {

	using namespace geo::ds;

	template <class T>
	class itemsort {
	public:
		T pt;
		itemsort(T& pt) : pt(pt) {
		}

		bool operator()(const T& a, const T& b) {
			double d1 = std::pow(a[0] - pt[0], 2.0) + std::pow(a[1] - pt[1], 2.0);
			double d2 = std::pow(b[0] - pt[0], 2.0) + std::pow(b[1] - pt[1], 2.0);
			return d1 < d2;
		}
	};

	class index {
	public:
		size_t from;
		size_t idx;
		index(size_t from, size_t idx) :
			from(from), idx(idx) {}
	};

	class indexsort {
	public:
		bool operator()(const index& a, const index& b) {
			if(a.idx == b.idx) {
				return a.from < b.from;	// Attempt to preserve locality.
			} else {
				return a.idx < b.idx;
			}
		}
	};

	int __oid = 0;


	template <class T>
	class node {
	private:
		constexpr static int maxDepth = 10;
		constexpr static int minSize = 256;

		double bounds[4];
		double midx, midy;
		bool leaf;
		int depth;
		int dims;
		size_t start;
		size_t end;
		mvector<T>* items;

		node<T>* nodes[4] = { nullptr };
		node<T>* parent;

		inline int nindex(double x, double y) {
			return ((x < midx) << 1) | (y < midy);
		}

		inline void nbounds(int idx, double b[4]) {
			b[0] = idx & 2 ? bounds[0] : midx;
			b[1] = idx & 1 ? bounds[1] : midy;
			b[2] = idx & 2 ? midx : bounds[2];
			b[3] = idx & 1 ? midy : bounds[3];
		}

		void build(double bbounds[4]) {

			for(int i = 0; i < 4; ++i)
				bounds[i] = bbounds[i];

			if(depth == maxDepth || end - start <= minSize) {
				std::ofstream out("test.csv", std::ios::app);
				out << std::setprecision(9);
				out << ++__oid << "," << bounds[0];
				for(int i = 1; i < 4; ++i)
					out << "," << bounds[i];
				out << "\n";
				leaf = true;
				return;
			}

			midx = (bounds[2] + bounds[0]) / 2.0;
			midy = (bounds[3] + bounds[1]) / 2.0;

			size_t idxCounts[4] = {0};
			{
				T item;
				size_t idx;
				indexsort idxSort;
				std::vector<index> indices;

				// Collect the node indices and the current position in the list.
				for(size_t i = start; i < end; ++i) {
					items->get(i, item);
					idx = nindex(item[0], item[1]);
					indices.emplace_back(i, idx);
					idxCounts[idx]++;
				}

				// Sort on node index.
				std::sort(indices.begin(), indices.end(), idxSort);

				if(indices.size() > 1000000) {
					mvector<T> tmp(end - start);
					for(size_t i = 0; i < indices.size(); ++i) {
						items->get(indices[i].from, item);
						tmp.insert(i, item);
					}
					for(size_t i = 0; i < indices.size(); ++i) {
						tmp.get(i, item);
						items->insert(i + start, item);
					}
				} else {
					std::vector<T> tmp(end - start);
					for(size_t i = 0; i < indices.size(); ++i)
						items->get(indices[i].from, tmp[i]);
					for(size_t i = 0; i < indices.size(); ++i)
						items->insert(i + start, tmp[i]);
				}

			}

			size_t s = start, e;
			double cbounds[4];
			for(int i = 0; i < 4; ++i) {
				if(idxCounts[i] > 0) {
					e = s + idxCounts[i];
					nbounds(i, cbounds);
					nodes[i] = new node<T>(depth + 1, s, e, items, dims);
					nodes[i]->parent = this;
					nodes[i]->build(cbounds);
					s = e;
				}
			}
		}

	public:

		node(mvector<T>* items, int dims) :
			node(0, 0, items->size(), items, dims) {}

		node(int depth, size_t start, size_t end, mvector<T>* items, int dims) :
			midx(0), midy(0), leaf(false), depth(depth), dims(dims),
			start(start), end(end), items(items), parent(nullptr) {
		}

		bool intersects(double box[4]) {
			return !(box[0] > bounds[2] || box[2] < bounds[0] || box[1] > bounds[3] || box[3] < bounds[1]);
		}

		bool contains(double x, double y) {
			return x >= bounds[0] && x <= bounds[2] && y >= bounds[1] && y <= bounds[3];
		}

		void build() {
			T item;
			double bbounds[] =  { std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest() };
			for(size_t i = start; i < end; ++i) {
				items->get(i, item);
				if(item[0] < bbounds[0]) bbounds[0] = item[0];
				if(item[0] > bbounds[2]) bbounds[2] = item[0];
				if(item[1] < bbounds[1]) bbounds[1] = item[1];
				if(item[1] > bbounds[3]) bbounds[3] = item[1];

			}
			double w = bbounds[2] - bbounds[0];
			double h = bbounds[3] - bbounds[1];
			if(w > h) {
				double dif = (w - h) / 2;
				bbounds[1] -= dif;
				bbounds[3] += dif;
			} else {
				double dif = (h - w) / 2;
				bbounds[0] -= dif;
				bbounds[2] += dif;
			}
			build(bbounds);
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
			std::sort(lst.begin(), lst.end(), itemsort<T>(pt));
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
			double sbounds[] = {pt[0] - radius, pt[1] - radius, pt[0] + radius, pt[1] + radius};
			return search(pt, radius, iter, sbounds);
		}

		template <class TIter>
		size_t search(const T& pt, double radius, TIter iter, double sbounds[4]) {

			if(!intersects(sbounds))
				return 0;

			size_t count = 0;

			if(leaf) {
				T item;
				for(size_t i = start; i < end; ++i) {
					items->get(i, item);
					double d = std::pow(pt[0] - item[0], 2.0) + std::pow(pt[1] - item[1], 2.0);
					if(d < radius * radius) {
						*iter = item;
						++iter;
						++count;
					}
				}
			} else {
				for(int i = 0; i < 4; ++i) {
					if(nodes[i])
						count += nodes[i]->search(pt, radius, iter, sbounds);
				}
			}

			return count;
		}

		~node() {
			for(int i = 0; i < 4; ++i) {
				if(nodes[i]) delete nodes[i];
			}
		}


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
	node<T>* m_root;
	size_t m_dims;

public:

	/**
	 * Construct the KDTree with the given number of dimensions.
	 */
	mqtree(size_t dims = 3) :
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
		m_root->build();
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

	~mqtree() {
		destroy();
	}
};


}
}



#endif /* INCLUDE_DS_KDTREE_HPP_ */
