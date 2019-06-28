/*
 * mvector.hpp
 *
 *  Created on: Jun 21, 2019
 *      Author: rob
 */

#ifndef INCLUDE_DS_MVECTOR_HPP_
#define INCLUDE_DS_MVECTOR_HPP_


#include <sys/mman.h>

#include <cstdlib>
#include <memory>
#include <list>

#include "geo.hpp"
#include "util.hpp"

using namespace geo::util;

namespace {

template <class T>
class Sorter {
public:
	bool operator()(const T& a, const T& b) const {
		return a < b;
	}
};

}

#define MEM_LIMIT 1024*1024*10

namespace geo {
namespace ds {

template <class T>
class mvector {
private:
	T* m_data;
	std::unique_ptr<TmpFile> m_file;
	size_t m_idx;
	size_t m_count;

	template <class Sort>
	void sort(size_t s, size_t e, Sort sorter) {
		if(s < e) {
			size_t p = partition(s, e, sorter);
			sort(s, p, sorter);
			sort(p + 1, e, sorter);
		}
	}

	template <class Sort>
	size_t partition(size_t s, size_t e, Sort sorter) {
		// If the chunk size is less than the configured limit,  sort in memory.
		if(e - s <= MEM_LIMIT / sizeof(T)) {
			T pivot, a, b;
			size_t o = s; // offset.
			bool swapped = false;
			std::vector<T> buf(e - s + 1);
			if(!get(s, buf, buf.size()))
				g_runerr("Invalid index: " << s)
			pivot = buf[(e - s) / 2];
			--s;
			++e;
			while(true) {
				do {
					++s;
					a = buf[s - o];
				} while(sorter(a, pivot));
				do {
					--e;
					b = buf[e - o];
				} while(sorter(pivot, b));
				if(s >= e) {
					if(swapped)
						insert(o, buf, buf.size());
					return e;
				}
				buf[e - o] = a;
				buf[s - o] = b;
				swapped = true;
			}
		} else {
			T pivot, a, b;
			if(!get((s + e) / 2, pivot))
				g_runerr("Invalid index: " << (s + e) / 2)
			--s;
			++e;
			while(true) {
				do {
					++s;
					if(!get(s, a))
						g_runerr("Invalid index: " << s)
				} while(sorter(a, pivot));
				do {
					--e;
					if(!get(e, b))
						g_runerr("Invalid index: " << e)
				} while(sorter(pivot, b));
				if(s >= e)
					return e;
				insert(e, a);
				insert(s, b);
			}
		}
	}

public:
	mvector(size_t size = 1) :
		m_data(nullptr), m_idx(0), m_count(0) {
		resize(size);
	}

	void resize(size_t count) {
		size_t size = count * sizeof(T);
		size_t oldSize = 0;
		if(!m_file.get()) {
			m_file.reset(new TmpFile(size));
		} else if(size > m_file->size) {
			oldSize = m_file->size;
			m_file->resize(size);
		} else {
			return;
		}
		if(!m_data) {
			m_data = (T*) mmap(0, size, PROT_READ|PROT_WRITE, MAP_SHARED, m_file->fd, 0);
		} else {
			m_data = (T*) mremap(m_data, oldSize, size, MREMAP_MAYMOVE, 0);
		}
		m_count = count;
	}

	void swap(mvector& other) {
		m_file.swap(other.m_file);

		T* data = m_data;
		m_data = other.m_data;
		other.m_data = data;

		size_t tmp = m_idx;
		m_idx = other.m_idx;
		other.m_idx = tmp;

		tmp = m_count;
		m_count = other.m_count;
		other.m_count = tmp;
	}

	bool push(const T& item) {
		if(m_idx >= m_count)
			resize(m_count *  2);
		std::memcpy(m_data + m_idx, &item, sizeof(T)) ;
		++m_idx;
		return true;
	}

	bool push(const std::vector<T>& items, size_t len) {
		if(m_idx + len >= m_count) {
			while(m_idx + len >= m_count)
				m_count *=  2;
			resize(m_count);
		}
		std::memcpy(m_data + m_idx, items.data(), len * sizeof(T)) ;
		m_idx += items.size();
		return true;
	}

	bool insert(size_t idx, const T& item) {
		if(idx >= m_count)
			resize(m_count *  2);
		std::memcpy(m_data + idx, &item, sizeof(T)) ;
		if(idx > m_idx)
			m_idx = idx;
		return true;
	}

	bool insert(size_t idx, const std::vector<T>& item, size_t len) {
		if(idx + len >= m_count) {
			while(idx + len >= m_count)
				m_count *=  2;
			resize(m_count);
		}
		std::memcpy(m_data + idx, item.data(), len * sizeof(T)) ;
		if(idx > m_idx)
			m_idx = idx;
		return true;
	}

	bool get(size_t idx, T& item) const {
		if(idx < m_idx) {
			std::memcpy(&item, m_data + idx, sizeof(T));
			return true;
		}
		return false;
	}

	size_t get(size_t idx, std::vector<T>& items, size_t len) const {
		if(idx > m_idx)
			return 0;
		if(idx + len > m_idx)
			len = m_idx - idx;
		items.resize(len);
		std::memcpy(items.data(), m_data + idx, len * sizeof(T));
		return len;
	}

	T operator[](size_t idx) const {
		if(idx < m_count) {
			T item;
			std::memcpy(&item, m_data + idx, sizeof(T));
			return item;
		} else {
			throw std::runtime_error("Index out of range.");
		}
	}

	void clear() {
		m_idx = 0;
		m_count = 0;
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

	void sort() {
		Sorter<T> sorter;
		sort(sorter);
	}

	template <class Sort>
	void sort(Sort sorter) {
		g_trace("Sorting mvector...")
		sort(0, size() - 1, sorter);
		g_trace("Sorted.")
	}

	~mvector() {
		if(m_file.get())
			munmap(m_data, m_file->size);
	}
};

} // ds
} // geo


#endif /* INCLUDE_DS_MVECTOR_HPP_ */
