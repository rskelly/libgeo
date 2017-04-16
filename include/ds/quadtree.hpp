/*
 * quadtree.hpp
 *
 *  Created on: Mar 29, 2017
 *      Author: rob
 */

#ifndef _QUADTREE_HPP_
#define _QUADTREE_HPP_

#include "util.hpp"

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <list>

using namespace geo::util;

namespace geo {
namespace ds {

template <class T>
class QuadTree {
private:
	std::list<std::tuple<double, double, T> > m_items;
	Bounds m_bounds;
	std::unordered_map<uint8_t, std::unique_ptr<QuadTree> > m_nodes;
	int m_maxDepth;
	int m_depth;
	QuadTree* m_parent;
	uint8_t m_idx;
	std::list<std::tuple<double, double, T> > m_allItems;

	uint8_t idx(double x, double y) {
		uint8_t i = 0;
		if(x >= m_bounds.midx())
			i |= 1;
		if(y >= m_bounds.midy())
			i |= 2;
		return i;
	}

	QuadTree* getNode(uint8_t idx) {
		if(m_nodes.find(idx) == m_nodes.end())
			return nullptr;
		return m_nodes[idx];
	}

	QuadTree* getNode(double x, double y) {
		uint8_t i = idx(x, y);
		if(m_nodes.find(i) == m_nodes.end()) {
			Bounds b;
			if(i & 1) {
				b.minx(m_bounds.midx());
				b.maxx(m_bounds.maxx());
			} else {
				b.minx(m_bounds.minx());
				b.maxx(m_bounds.midx());
			}
			if(i & 2) {
				b.miny(m_bounds.midy());
				b.maxy(m_bounds.maxy());
			} else {
				b.miny(m_bounds.miny());
				b.maxy(m_bounds.midy());
			}
			std::unique_ptr<QuadTree> t(new QuadTree(b, m_maxDepth, m_depth + 1, this, i));
			m_nodes[i] = std::move(t);
		}
		return m_nodes[i].get();
	}

	QuadTree(const Bounds& bounds, int maxDepth, int depth, QuadTree* parent, uint8_t idx) :
		m_bounds(bounds),
		m_maxDepth(maxDepth),
		m_depth(depth),
		m_parent(parent),
		m_idx(idx) {}

public:
	QuadTree(const Bounds& bounds, int maxDepth) :
		m_bounds(bounds),
		m_maxDepth(maxDepth),
		m_depth(0),
		m_parent(nullptr),
		m_idx(0) {

	}

	void add(double x, double y, T data) {
		if(m_depth == 0)
			m_allItems.push_back(std::make_tuple(x, y, data));
		if(m_depth == m_maxDepth) {
			m_items.push_back(std::make_tuple(x, y, data));
		} else {
			getNode(x, y)->add(x, y, data);
		}
	}

	template <class U>
	int search(const Bounds& bounds, U inserter) {
		int count = 0;
		if(m_bounds.intersects(bounds)) {
			if(m_depth == m_maxDepth) {
				for(const auto& it : m_items) {
					double x = std::get<0>(it);
					double y = std::get<1>(it);
					if(bounds.contains(x, y)) {
						*inserter = it;
						++inserter;
						++count;
					}
				}
			} else {
				for(auto& it : m_nodes)
					count += it.second->search(bounds, inserter);
			}
		}
		return count;
	}

	class Sorter {
	public:
		double x, y;
		Sorter(double x, double y) : x(x), y(y) {}
		bool operator()(std::tuple<double, double, T>& a, std::tuple<double, double, T>& b) {
			double da = g_sq(x - std::get<0>(a)) + g_sq(y - std::get<1>(a));
			double db = g_sq(x - std::get<0>(b)) + g_sq(y - std::get<1>(a));
			return da < db;
		}
	};

	template <class U>
	int nearest(double x, double y, size_t count, U inserter, double radius = 10.0) {
		Bounds bounds(x - radius, y - radius, x + radius, y + radius);
		size_t rcount = 0;
		do {
			std::list<std::tuple<double, double, T> > result;
			search(bounds, std::back_inserter(result));
			if(result.size() < count) {
				double minx = bounds.midx() - bounds.width();
				double maxx = bounds.midx() + bounds.width();
				double miny = bounds.midy() - bounds.height();
				double maxy = bounds.midy() + bounds.height();
				bounds.extend(minx, miny);
				bounds.extend(maxx, maxy);
				continue;
			}
			if(result.size() > count) {
				result.sort(Sorter(x, y));
				auto iter = result.begin();
				std::advance(iter, count);
				result.erase(iter, result.end());
			}
			rcount = 0;
			for(const auto& it : result) {
				*inserter = it;
				++inserter;
				++rcount;
			}
		} while(rcount < count && !bounds.contains(m_bounds));
		return rcount;
	}

	int size() const {
		return m_allItems.size();
	}

	const std::list<std::tuple<double, double, T> >& items() const {
		return m_allItems;
	}

};

static uint64_t s_id = 0;

static uint64_t s_nextId() {
	uint64_t id = s_id;
	++s_id;
	return id;
}

template <class T>
class FileQuadTree;

template <class T>
class FileQuadTreeNode {
	friend class FileQuadTree<T>;
protected:

	typedef std::tuple<uint64_t, double, double, T> Item;

	std::unordered_map<char, std::unique_ptr<FileQuadTreeNode<T> > > m_nodes;
	std::vector<size_t> m_offsets;
	geo::util::Bounds m_bounds;
	MappedFile& m_data;
	int m_depth;
	int m_maxDepth;
	size_t m_count;
	size_t m_cacheSize;

	uint8_t idx(double x, double y) {
		uint8_t i = 0;
		if(x >= m_bounds.midx())
			i |= 1;
		if(y >= m_bounds.midy())
			i |= 2;
		return i;
	}

	FileQuadTreeNode* getNode(double x, double y) {
		uint8_t i = idx(x, y);
		if(m_nodes.find(i) == m_nodes.end()) {
			Bounds b;
			if(i & 1) {
				b.minx(m_bounds.midx());
				b.maxx(m_bounds.maxx());
			} else {
				b.minx(m_bounds.minx());
				b.maxx(m_bounds.midx());
			}
			if(i & 2) {
				b.miny(m_bounds.midy());
				b.maxy(m_bounds.maxy());
			} else {
				b.miny(m_bounds.miny());
				b.maxy(m_bounds.midy());
			}
			std::unique_ptr<FileQuadTreeNode> t(new FileQuadTreeNode(m_data, b, m_depth + 1, m_maxDepth, m_cacheSize));
			m_nodes[i] = std::move(t);
		}
		return m_nodes[i].get();
	}

	FileQuadTreeNode(MappedFile& data, const geo::util::Bounds& bounds,
			int depth, int maxDepth, int cacheSize) :
		m_bounds(bounds),
		m_data(data),
		m_depth(depth), m_maxDepth(maxDepth),
		m_count(0),
		m_cacheSize(cacheSize) {
	}

public:

	// Adds an item to the tree at the given position; returns
	// the item's ID.
	uint64_t addItem(double x, double y, const T& item) {
		if(m_depth < m_maxDepth) {
			return getNode(x, y)->addItem(x, y, item);
		} else {
			// The size of a "row" in mapped memory.
			size_t size = m_cacheSize * sizeof(Item);
			// The offset into the mapped memory for the start of the current cache.
			size_t offset;

			// If the count has surpassed the size of the current block,
			// get a new offset and start a new block.
			if(m_count % m_cacheSize == 0) {
				offset = s_nextId();
				m_offsets.push_back(offset);
				// Resize the memory if necessary.
				m_data.reset(offset * size + size);
			} else {
				 offset = m_offsets[m_offsets.size() - 1];
			}

			// The position of the new item.
			size_t position = offset * size + (m_count % m_cacheSize) * sizeof(Item);

			// Write the item.
			Item i = std::make_tuple(position, x, y, item);
			char* fdata = (char*) m_data.data() + position;
			std::memcpy(fdata, &i, sizeof(Item));

			++m_count;

			return position;
		}
	}

	// Get an item from the tree.
	bool getItem(Item& item, uint64_t id) {
		if(id + sizeof(Item) > m_data.size())
			return false;
		char* data = (char*) m_data.data() + id;
		std::memcpy(&item, data, sizeof(Item));
		return true;
	}

	// Update the item in the tree.
	bool updateItem(const Item& item, uint64_t id) {
		if(id + sizeof(Item) > m_data.size())
			return false;
		char* data = (char*) m_data.data() + id;
		std::memcpy(data, &item, sizeof(Item));
		return true;
	}

	int count() const {
		if(m_depth < m_maxDepth) {
			int c = 0;
			for(const auto& it : m_nodes)
				c += it.second->count();
			return c;
		} else {
			return m_count;
		}
	}

	template <class U>
	int search(const Bounds& bounds, U inserter) {
		int count = 0;
		if(m_bounds.intersects(bounds)) {
			if(m_depth == m_maxDepth) {
				// Move on to check items in mapped memory.
				std::vector<Item> items(m_cacheSize);
				size_t size = m_cacheSize * sizeof(Item);
				for(size_t i = 0; i < m_offsets.size(); ++i) {
					size_t offset = m_offsets[i] * size;
					char* fdata = (char*) m_data.data() + offset;
					char* vdata = (char*) items.data();
					std::memcpy(vdata, fdata, size);
					for(size_t j = 0; j < g_min(m_cacheSize, m_count - i * m_cacheSize); ++j) {
						const Item& it = items[j];
						double x = std::get<1>(it);
						double y = std::get<2>(it);
						if(bounds.contains(x, y)) {
							*inserter = it;
							++inserter;
							++count;
						}
					}
				}
			} else {
				for(auto& it : m_nodes)
					count += it.second->search(bounds, inserter);
			}
		}
		return count;
	}

	class Sorter {
	public:
		double x, y;
		Sorter(double x, double y) : x(x), y(y) {}
		bool operator()(std::tuple<uint64_t, double, double, T>& a, std::tuple<uint64_t, double, double, T>& b) {
			double da = g_sq(x - std::get<1>(a)) + g_sq(y - std::get<2>(a));
			double db = g_sq(x - std::get<1>(b)) + g_sq(y - std::get<2>(a));
			return da < db;
		}
	};

	template <class U>
	int nearest(double x, double y, size_t count, U inserter, double radius = 10.0) {
		Bounds bounds(x - radius, y - radius, x + radius, y + radius);
		size_t rcount = 0;
		do {
			std::list<std::tuple<double, double, T> > result;
			search(bounds, std::back_inserter(result));
			if(result.size() < count) {
				double minx = bounds.midx() - bounds.width();
				double maxx = bounds.midx() + bounds.width();
				double miny = bounds.midy() - bounds.height();
				double maxy = bounds.midy() + bounds.height();
				bounds.extend(minx, miny);
				bounds.extend(maxx, maxy);
				continue;
			}
			if(result.size() > count) {
				result.sort(Sorter(x, y));
				auto iter = result.begin();
				std::advance(iter, count);
				result.erase(iter, result.end());
			}
			rcount = 0;
			for(const auto& it : result) {
				*inserter = it;
				++inserter;
				++rcount;
			}
		} while(rcount < count && !bounds.contains(m_bounds));
		return rcount;
	}

};

template <class T>
class FileQuadTree : public FileQuadTreeNode<T> {
private:
	MappedFile m_data;
	std::vector<FileQuadTreeNode<T>*> m_inodes;
	size_t m_idx;
	size_t m_nidx;
	FileQuadTreeNode<T>* m_iqt;

public:
	FileQuadTree(const geo::util::Bounds& bounds, int maxDepth, int cacheSize = 4096) :
		FileQuadTreeNode<T>(m_data, bounds, 0, maxDepth, cacheSize),
		m_idx(0), m_nidx(0), m_iqt(nullptr) {
	}


	void reset() {
		m_idx = 0;
		m_nidx = 0;
		m_inodes.clear();
		std::queue<FileQuadTreeNode<T>*> q;
		m_iqt = this;
		q.push(m_iqt);
		while(!q.empty()) {
			FileQuadTreeNode<T>* t = q.front();
			q.pop();
			if(t->m_depth == t->m_maxDepth) {
				m_inodes.push_back(t);
			} else {
				for(auto& p : t->m_nodes)
					q.push(p.second.get());
			}
		}
	}

	bool next(std::tuple<uint64_t, double, double, T>& item) {
		while(m_nidx < m_inodes.size()) {
			if(m_idx < m_inodes[m_nidx]->m_items.size()) {
				std::memcpy(&item, &(m_inodes[m_nidx]->m_items[m_idx]), sizeof(std::tuple<uint64_t, double, double, T>));
				++m_idx;
				return true;
			} else {
				++m_nidx;
			}
		}
		return false;
	}
};

}
}

#endif /* _QUADTREE_HPP_ */
