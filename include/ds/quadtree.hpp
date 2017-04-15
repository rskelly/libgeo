/*
 * quadtree.hpp
 *
 *  Created on: Mar 29, 2017
 *      Author: rob
 */

#ifndef _QUADTREE_HPP_
#define _QUADTREE_HPP_

#include "util.hpp"

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

static size_t s_fqtnOffset = 0;

template <class T>
class FileQuadTreeNode {
private:

	typedef std::tuple<double, double, T> Item;

	MappedFile& m_data;
	std::vector<size_t> m_offsets;
	std::vector<Item> m_items;
	geo::util::Bounds m_bounds;
	std::unordered_map<char, std::unique_ptr<FileQuadTreeNode<T> > > m_nodes;
	size_t m_count;
	size_t m_cacheSize;
	int m_depth;
	int m_maxDepth;

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

protected:
	FileQuadTreeNode(MappedFile& data, const geo::util::Bounds& bounds,
			int depth, int maxDepth, int cacheSize) :
		m_data(data),
		m_bounds(bounds),
		m_count(0),
		m_cacheSize(cacheSize),
		m_depth(depth), m_maxDepth(maxDepth) {
		m_items.resize(m_cacheSize);
	}

public:

	void add(double x, double y, const T& item) {
		if(m_depth < m_maxDepth) {
			getNode(x, y)->add(x, y, item);
		} else {
			m_items[m_count % m_cacheSize] = std::make_tuple(x, y, item);
			++m_count;
			if(m_count > 0 && m_count % m_cacheSize == 0) {
				size_t fOffset = s_fqtnOffset++;
				size_t size = m_cacheSize * sizeof(Item);
				size_t offset = fOffset * size;
				m_data.reset((uint64_t) size * (fOffset + 1), true);
				char* fdata = (char*) m_data.data();
				char* vdata = (char*) m_items.data();
				std::memcpy(fdata + offset, vdata, size);
				m_offsets.push_back(fOffset);
			}
		}
	}

	int count() const {
		if(m_depth < m_maxDepth) {
			int c = 0;
			for(const auto& it : m_nodes)
				c += it.second.count();
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
				// If there are items in memory, process them first.
				if(m_count % m_cacheSize != 0) {
					for(size_t j = 0; j < m_count % m_cacheSize; ++j) {
						const Item& it = m_items[j];
						double x = std::get<0>(it);
						double y = std::get<1>(it);
						if(bounds.contains(x, y)) {
							*inserter = it;
							++inserter;
							++count;
						}
					}
				}
				// Move on to check items in mapped memory.
				std::vector<Item> items(m_cacheSize);
				size_t size = m_cacheSize * sizeof(Item);
				for(size_t i = 0; i < m_offsets.size(); ++i) {
					size_t offset = m_offsets[i];
					char* fdata = (char*) m_data.data();
					char* vdata = (char*) items.data();
					std::memcpy(vdata, fdata + offset * size, size);
					for(size_t j = 0; j < m_count - i * m_cacheSize; ++j) {
						const Item& it = m_items[j];
						double x = std::get<0>(it);
						double y = std::get<1>(it);
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

};

template <class T>
class FileQuadTree : public FileQuadTreeNode<T> {
private:
	MappedFile m_data;

public:
	FileQuadTree(const geo::util::Bounds& bounds, int maxDepth, int cacheSize = 4096) :
		FileQuadTreeNode<T>(m_data, bounds, 0, maxDepth, cacheSize) {
	}
};

}
}

#endif /* _QUADTREE_HPP_ */
