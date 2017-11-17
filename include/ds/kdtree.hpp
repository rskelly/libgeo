/*
 * kdtree.hpp
 *
 *  Created on: Nov 17, 2017
 *      Author: rob
 */

#ifndef INCLUDE_DS_KDTREE_HPP_
#define INCLUDE_DS_KDTREE_HPP_

namespace rs {
namespace ds {

template <class T>
class KDTreeSorter {
private:
	int m_idx;
public:
	KDTreeSorter(int idx) :
		m_idx(idx) {
	}

	bool operator()(T* a, T* b) {
		return (*a)[m_idx] < (*b)[m_idx];
	}
};

template <class T>
class KDTreeNDSorter {
private:
	int m_dims;
	T* m_item;

public:
	KDTreeNDSorter(T* item, int dims) :
		m_item(item),
		m_dims(dims) {
	}

	bool operator()(T* a, T* b) {
		double dista = 0;
		double distb = 0;
		for(int i = 0; i < m_dims; ++i) {
			dista += g_sq((*a)[i] - (*m_item)[i]);
			distb += g_sq((*b)[i] - (*m_item)[i]);
		}
		return dista < distb;
	}
};

template <class T>
class KDTree {
private:
	KDTree* m_left;
	KDTree* m_right;
	KDTree* m_parent;
	int m_dims;
	int m_depth;
	T* m_item;

	KDTree(int dims, int depth, KDTree* parent) :
		m_left(nullptr), m_right(nullptr),
		m_parent(parent),
		m_dims(dims),
		m_depth(depth),
		m_item(nullptr) {
	}

	double dist(T* a,  T* b) {
		double dist = 0;
		for(int i = 0; i < m_dims; ++i)
			dist += g_sq((*a)[i] - (*b)[i]);
		return dist;
	}

	double dist(T* a,  T* b, int dim) {
		return g_sq((*a)[dim] - (*b)[dim]);
	}

public:

	KDTree(int dims) :
		KDTree(dims, 0, nullptr) {
	}

	template <class U>
	void add(U begin, U end) {
		std::vector<T*> items(begin, end);
		add(items);
	}

	void add(std::vector<T*>& items) {
		KDTreeSorter sorter(m_depth % m_dims);
		std::sort(items.begin(), items.end(), sorter);
		int size = items.size();
		m_item = items[size / 2];
		if(size > 1) {
			m_left = new KDTree(m_dims, m_depth + 1);
			m_right = new KDTree(m_dims, m_depth + 1);
			m_left->add(items.begin(), items.begin() + size / 2);
			m_right->add(items.begin + size / 2 + 1, items.end());
		}
	}

	template <class U>
	void knn(const T& item, int count, U iter) {

		std::vector<T*> items;

		// Crawl down the tree looking for the closest leaf.
		KDTree* n = this;
		while(n) {
			int idx = n->m_depth % m_dims;
			if(item[idx] < m_item[idx]) {
				if(!n->m_left)
					break;
				n = n->m_left;
			} else {
				if(!n->m_right)
					break;
				n = n->m_right;
			}
		}

		T* best = n->m_item;
		items.push_front(best);

		// Get the distance to the leaf item.
		double d = dist(&item, items.front());

		// Crawl back up the tree looking for closer items.
		n = n->m_parent;
		while(n) {
			// If the current point is closer to the hyperplane than the
			// current distance, search the other side of the tree.
			if(dist(&item, n->m_item, n->m_depth % m_dims) < d) {
				n->knn(item, count, std::back_inserter(items));
			}
			// Trim the array if it's too large.
			if(items.size() > count) {
				KDTreeNDSorter sorter(&item, m_dims);
				std::sort(items.begin(), items.end(), sorter);
				items.resize(count);
			}

			n = n->m_parent;
		}

		// Trim the array if it's too large.
		if(items.size() > count) {
			KDTreeNDSorter sorter(&item, m_dims);
			std::sort(items.begin(), items.end(), sorter);
			items.resize(count);
		}

		// Copy to output iterator.
		for(T* i : items) {
			*iter = i;
			++iter;
		}
	}

	~KDTree() {
		if(m_left)
			delete m_left;
		if(m_right)
			delete m_right;
	}
};


}
}



#endif /* INCLUDE_DS_KDTREE_HPP_ */
