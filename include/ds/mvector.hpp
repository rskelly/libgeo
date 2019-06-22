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

#include "geo.hpp"
#include "util.hpp"

using namespace geo::util;

namespace geo {
namespace ds {

template <class T>
class mvector {
private:
	T* m_data;
	std::unique_ptr<TmpFile> m_file;
	size_t m_idx;
	size_t m_count;

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

	bool insert(size_t idx, const T& item) {
		if(idx >= m_count)
			return false;
		std::memcpy(m_data + idx, &item, sizeof(T)) ;
		return true;
	}

	bool get(size_t idx, T& item) const {
		if(idx < m_count) {
			std::memcpy(&item, m_data + idx, sizeof(T));
			return true;
		}
		return false;
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

	~mvector() {
		munmap(m_data, m_file->size);
	}
};

} // ds
} // geo


#endif /* INCLUDE_DS_MVECTOR_HPP_ */
