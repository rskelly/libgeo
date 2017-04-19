#ifndef __QTREE_HPP__
#define __QTREE_HPP__

#include "util.hpp"
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace geo::util;

static size_t s_idx = 0;
static size_t s_bufSize = 64;

template <class T>
class CatItem {
private:
	size_t m_idx;
	size_t m_iidx;
	size_t m_iiidx;

public:
	size_t count;
	std::vector<std::pair<size_t, size_t> > items;

	CatItem() : count(0) {}

	size_t newPosition() {
		if(items.empty() || (count % s_bufSize) == 0)
			items.push_back(std::make_pair(0, s_idx++));
		std::pair<size_t, size_t>& item = items[items.size() - 1];
		size_t pos = item.second * s_bufSize * sizeof(T) + item.first * sizeof(T);
		item.first++;
		++count;
		return pos;
	}

	void reset() {
		m_idx = 0;
		m_iidx = 0;
		m_iiidx = 0;
	}

	bool nextPosition(size_t* idx) {
		while(m_idx < items.size()) {
			if(m_iidx < items[m_idx].first) {
				*idx = items[m_idx].second * s_bufSize * sizeof(T) + m_iidx * sizeof(T);
				++m_iidx;
				return true;
			} else {
				m_iidx = 0;
				++m_idx;
			}
		}
		return false;
	}
};

template <class T>
class QTree {
private:
	// The geographic boundary corresponding to this tree. Will be reshaped to a square 
	// using the longer side.
	Bounds m_bounds;
	// The number of bins along each side.
	size_t m_bins;
	// The length of one side in map units.
	double m_size; 
	// A catalog to keep pointers into the mapped memory.
	std::unordered_map<size_t, CatItem<T> > m_catalog;
	// The mapped memory.
	MappedFile m_mapped;	

public:
	QTree(const Bounds& bounds, size_t bins) :
		m_bounds(bounds), 
		m_bins(bins) {

		m_size = g_max(m_bounds.width(), m_bounds.height());
		m_bounds.set(m_bounds.midx() - m_size / 2, m_bounds.midy() - m_size / 2, 
			m_bounds.midx() + m_size / 2, m_bounds.midy() + m_size / 2);
	}

	const Bounds& bounds() const {
		return m_bounds;
	}

	size_t toBinX(double x) {
		return (size_t) (x - m_bounds.minx()) / m_size * m_bins;
	}

	size_t toBinY(double y) {
		return (size_t) (y - m_bounds.miny()) / m_size * m_bins;
	}

	size_t index(double x, double y) {
		return toBinY(y) * m_bins + toBinX(x);
	}

	std::list<size_t> indices(double x, double y, double radius, bool outer = false) {
		size_t b0 = g_max(0, toBinX(x - radius));
		size_t b2 = g_min(m_bins, toBinX(x + radius));
		size_t b1 = g_max(0, toBinY(y - radius));
		size_t b3 = g_min(m_bins, toBinY(y + radius));
		std::list<size_t> lst;
		if(outer) {
			for(size_t a = b0; a <= b2; ++a)
				lst.push_back(b1 * m_bins + a);
			for(size_t b = b1 + 1; b <= b3 - 1; ++b) {
				lst.push_back(b * m_bins + b0);
				lst.push_back(b * m_bins + b2);
			}
			for(size_t a = b0; a <= b2; ++a)
				lst.push_back(b3 * m_bins + a);
		} else {
			for(size_t b = b1; b <= b3; ++b)
				for(size_t a = b0; a <= b2; ++a) {
					lst.push_back(b * m_bins + a);
			}
		}
		return lst;
	}

	void addItem(double x, double y, const T& item) {
		size_t idx = index(x, y);
		size_t pos = m_catalog[idx].newPosition();
		if(m_mapped.size() < pos + sizeof(T))
			m_mapped.reset(pos + sizeof(T));
		char* data = (char*) m_mapped.data() + pos;
		std::memcpy(data, &item, sizeof(T));
	}

	template <class U>
	void search(double x, double y, double radius, U output) {
		if(!m_bounds.contains(x, y))
			return;
		std::list<size_t> idx = indices(x, y, radius);
		char* data = (char*) m_mapped.data();
		radius = g_sq(radius);
		for(const size_t& i : idx) {
			if(m_catalog.find(i) != m_catalog.end()) {
				CatItem<T>& ci = m_catalog[i];
				ci.reset();
				size_t pos = 0;
				T item;
				while(ci.nextPosition(&pos)) {
					if(pos + sizeof(T) > m_mapped.size())
						g_runerr("Requested position is out of bounds: " << pos);
					std::memcpy(&item, data + pos, sizeof(T));
					if(g_sq(item.x - x) + g_sq(item.y - y) <= radius) {
						*output = item;
						++output;
					}
				}
			}
		}
	}

	template <class U>
	void search(const Bounds& bounds, U output) {
		if(!m_bounds.intersects(bounds))
			return;
		std::list<size_t> idx = indices(bounds.midx(), bounds.midy(), g_max(bounds.width(), bounds.height()));
		char* data = (char*) m_mapped.data();
		for(const size_t& i : idx) {
			if(m_catalog.find(i) != m_catalog.end()) {
				CatItem<T>& ci = m_catalog[i];
				ci.reset();
				size_t pos = 0;
				T item;
				while(ci.nextPosition(&pos)) {
					if(pos + sizeof(T) > m_mapped.size())
						g_runerr("Requested position is out of bounds: " << pos);
					std::memcpy(&item, data + pos, sizeof(T));
					if(bounds.contains(item.x, item.y)) {
						*output = item;
						++output;
					}
				}
			}
		}
	}

	struct _sorter {
		double x, y;
		_sorter(double x, double y) : x(x), y(y) {}
		bool operator()(const T& a, const T& b) {
			return  (g_sq(x - a.x) + g_sq(y - a.y)) < (g_sq(x - b.x) + g_sq(y - b.y));
		}
	};

	template <class U>
	void nearest(double x, double y, size_t n, U output) {
		char* data = (char*) m_mapped.data();
		std::vector<T> tmp;
		size_t radius = 1;
		while(radius <= m_bins) {
			
			std::list<size_t> idx = indices(x, y, radius, true); // Only the outer ring; inner ring already searched.

			for(const size_t& i : idx) {
				
				if(m_catalog.find(i) != m_catalog.end()) {
					CatItem<T>& ci = m_catalog[i];
					ci.reset();
					size_t pos = 0;
					T item;
					while(ci.nextPosition(&pos)) {
						if(pos + sizeof(T) > m_mapped.size())
							g_runerr("Requested position is out of bounds: " << pos);
						std::memcpy(&item, data + pos, sizeof(T));
						tmp.push_back(std::move(item));
					}
				}
			}

			if(tmp.size() >= n) {
				if(tmp.size() > n) {
					struct _sorter sorter(x, y);
					std::sort(tmp.begin(), tmp.end(), sorter);
				}
				for(size_t i = 0; i < n; ++i) {
					*output = std::move(tmp[0]);
					++output;
				}
				break;
			} 

			++radius;
		}
	}

	~QTree() {
	}

};

#endif