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

	QuadTree* getNode(uint8_t pidx, uint8_t idx) {
		QuadTree* n = nullptr;
		if(m_parent) {
			n = m_parent->getNode(pidx);
			if(n)
				n = n->getNode(idx);
		}
		return n;
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
};


#endif /* _QUADTREE_HPP_ */
