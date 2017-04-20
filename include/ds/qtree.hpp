#ifndef __QTREE_HPP__
#define __QTREE_HPP__

#include "util.hpp"
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <thread>
#include <condition_variable>

using namespace geo::util;

static size_t s_position = 0;
static size_t s_bufSize = 1024;

class Writer {
public:
	virtual void doWrite() = 0;
	virtual ~Writer() {}
};

void writer(Writer* qt) {
	qt->doWrite();
}

template <class T>
class CatItem {
private:
	size_t m_idx;
	size_t m_iidx;

public:
	std::vector<std::pair<size_t, size_t> > items;

	CatItem() : m_idx(0), m_iidx(0) {}

	size_t update(size_t count) {
		size_t pos = s_position;
		s_position += count * sizeof(T);
		items.push_back(std::make_pair(count, pos));
		return pos;
	}

	void reset() {
		m_idx = 0;
		m_iidx = 0;
	}

	bool next(size_t* idx) {
		while(m_idx < items.size()) {
			if(m_iidx < items[m_idx].first) {
				*idx = items[m_idx].second;
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
class QTree : public Writer {
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
	// Cache points for writing.
	std::unordered_map<size_t, std::vector<T> > m_cache;
	// The mapped memory.
	MappedFile m_mapped;	
	std::queue<T> m_wq;
	std::thread m_writer;
	bool m_running;
	std::mutex m_mtx;
	std::condition_variable m_cdn;

public:
	QTree(const Bounds& bounds, size_t bins) :
		m_bounds(bounds), 
		m_bins(bins),
		m_running(false) {

		m_size = g_max(m_bounds.width(), m_bounds.height());
		m_bounds.set(m_bounds.midx() - m_size / 2, m_bounds.midy() - m_size / 2, 
			m_bounds.midx() + m_size / 2, m_bounds.midy() + m_size / 2);

	}

	void startWrite() {
		if(!m_running) {
			m_running = true;
			m_writer = std::thread(writer, this);
		}
	}

	void finishWrite() {
		if(m_running) {
			m_cdn.notify_one();
			m_running = false;
			m_writer.join();
		}
	}

	void doWrite() {
		while(m_running) {
			std::cerr << m_wq.size() << "\n";
			std::unique_lock<std::mutex> lock(m_mtx);
		    while(m_wq.empty()) {
			    m_cdn.wait(lock);
		    }
			T item = m_wq.front();
			m_wq.pop();
			m_cdn.notify_one();
			size_t idx = index(item.x, item.y);
			m_cache[idx].push_back(std::move(item));
			if(m_cache[idx].size() == s_bufSize)
				flush(idx);
		}
		flushAll();
	}

	const Bounds& bounds() const {
		return m_bounds;
	}

	int toBinX(double x) {
		return (int) (x - m_bounds.minx()) / m_size * m_bins;
	}

	int toBinY(double y) {
		return (int) (y - m_bounds.miny()) / m_size * m_bins;
	}

	size_t index(double x, double y) {
		return (size_t) toBinY(y) * m_bins + toBinX(x);
	}

	std::list<size_t> indices(double x, double y, double radius, bool outer = false) {
		int b0 = g_max(0, toBinX(x - radius));
		int b2 = g_max(1, g_min((int) m_bins, toBinX(x + radius)));
		int b1 = g_max(0, toBinY(y - radius));
		int b3 = g_max(1, g_min((int) m_bins, toBinY(y + radius)));
		std::list<size_t> lst;
		if(outer) {
			for(int a = b0; a <= b2; ++a)
				lst.push_back(b1 * m_bins + a);
			for(int b = b1 + 1; b <= b3 - 1; ++b) {
				lst.push_back(b * m_bins + b0);
				lst.push_back(b * m_bins + b2);
			}
			for(int a = b0; a <= b2; ++a)
				lst.push_back(b3 * m_bins + a);
		} else {
			for(int b = b1; b <= b3; ++b)
				for(int a = b0; a <= b2; ++a) {
					lst.push_back(b * m_bins + a);
			}
		}
		return lst;
	}

	void flushAll() {
		for(const auto& it : m_cache) {
			size_t idx = it.first;
			const std::vector<T>& items = it.second;
			if(!items.empty()) {
				size_t pos = m_catalog[idx].update(items.size());
				m_mapped.write((void*) items.data(), pos, items.size() * sizeof(T));
			}
		}
		m_cache.clear();
	}

	void flush(size_t idx) {
		std::vector<T>& items = m_cache[idx];
		if(!items.empty()) {
			size_t pos = m_catalog[idx].update(items.size());
			m_mapped.write((void*) items.data(), pos, items.size() * sizeof(T));
		}
		m_cache.erase(idx);
	}

	void addItem(const T& item) {
		startWrite();
		{
			std::lock_guard<std::mutex> lock(m_mtx);
			m_wq.push(item);
		    m_cdn.notify_one();
		}
	}

	template <class U>
	void search(double x, double y, double radius, U output) {
		if(!m_bounds.contains(x, y))
			return;
		finishWrite();
		std::list<size_t> idx = indices(x, y, radius);
		radius = g_sq(radius);
		for(const size_t& i : idx) {
			if(m_catalog.find(i) != m_catalog.end()) {
				CatItem<T>& ci = m_catalog[i];
				ci.reset();
				size_t pos = 0;
				T item;
				while(ci.next(&pos)) {
					if(!m_mapped.read(item, pos))
						g_runerr("Failed to read item.");
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
		finishWrite();
		std::list<size_t> idx = indices(bounds.midx(), bounds.midy(), g_max(bounds.width(), bounds.height()));
		for(const size_t& i : idx) {
			if(m_catalog.find(i) != m_catalog.end()) {
				CatItem<T>& ci = m_catalog[i];
				ci.reset();
				size_t pos = 0;
				T item;
				while(ci.next(&pos)) {
					if(!m_mapped.read(item, pos))
						g_runerr("Failed to read item.");
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
		if(!m_bounds.contains(x, y))
			return;
		finishWrite();
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
					while(ci.next(&pos)) {
						if(!m_mapped.read(item, pos))
							g_runerr("Failed to read item.");
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
