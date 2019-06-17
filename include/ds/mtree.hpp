
#include <sys/mman.h>

#include <cstring>
#include <vector>

#include "geo.hpp"

namespace geo {
namespace ds {

template <class T>
class mtree {
private:
	T* m_data;
	size_t m_size;
	size_t m_initSize;
	size_t m_idx;
	size_t m_iter;
	size_t m_maxCount;
	mtree<T>* m_nodes[4];
	double m_bounds[4];
	bool m_split;

	void init() {
		m_size = m_initSize;
		m_data = (T*) mmap(0, m_size * sizeof(T), PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_PRIVATE, 0, 0);
		if(!m_data)
			g_runerr("Failed to map memory in mtree.")
	}

	void expand() {
		if(m_size == 0) {
			init();
		} else {
			size_t size = m_size * 2;
			m_data = (T*) mremap(m_data,  m_size * sizeof(T), size * sizeof(T), MREMAP_MAYMOVE, 0);
			m_size = size;
			if(!m_data)
				g_runerr("Failed to remap memory in mtree.")
		}
	}

	void free() {
		munmap(m_data, m_size);
		m_data = nullptr;
		m_size = 0;
		m_idx = 0;
		m_iter = 0;
	}

	mtree<T>* getNode(const T& item) {
		double midx = (m_bounds[0] + m_bounds[2]) / 2.0;
		double midy = (m_bounds[1] + m_bounds[3]) / 2.0;
		int idx = 0;
		if(item.x() >= midx)
			idx |= 1;
		if(item.y() >= midy)
			idx |= 2;
		double bounds[] = {
			idx & 1 ? midx : m_bounds[0],
			idx & 2 ? midy : m_bounds[1],
			idx & 1 ? m_bounds[2] : midx,
			idx & 2 ? m_bounds[3] : midy,
		};
		if(!m_nodes[idx])
			m_nodes[idx] = new mtree<T>(bounds, m_maxCount, m_initSize);
		return m_nodes[idx];
	}

	void split() {
		reset();
		T item;
		while(next(item)) {
			mtree<T>* node = getNode(item);
			node->add(item);
		}
		free();
		m_split = true;
	}

	bool contains(const T& item, double bounds[4]) {
		return (item.x() >= bounds[0] && item.x() <= bounds[2] && item.y() >= bounds[1] && item.y() <= bounds[3]);
	}

public:

	mtree(double bounds[4], size_t maxCount = 256, size_t count = 8) :
		m_data(nullptr),
		m_size(0),
		m_initSize(count),
		m_idx(0),
		m_iter(0),
		m_maxCount(maxCount),
		m_nodes({nullptr}),
		m_split(false) {

		for(int i = 0; i < 4; ++i)
			m_bounds[i] = bounds[i];
	}

	void clear() {
		for(int i = 0; i < 4; ++i) {
			if(m_nodes[i])
				delete m_nodes[i];
			m_nodes[i] = nullptr;
		}
		free();
	}

	size_t size() {
		size_t s = 0;
		if(m_split) {
			for(mtree<T>* n : m_nodes)
				if(n) s += n->size();
		} else {
			s += m_idx;
		}
		return s;
	}

	void reset() {
		m_iter = 0;
		if(m_split) {
			for(mtree<T>* n : m_nodes) {
				if(n)
					n->reset();
			}
		}
	}

	bool next(T& item) {
		if(m_split) {
			for(mtree<T>* n : m_nodes) {
				if(n->next(item))
					return true;
			}
		} else {
			if(m_iter < m_idx) {
				std::memcpy(&item, m_data + m_iter, sizeof(T));
				++m_iter;
				return true;
			}
		}
		return false;
	}

	void add(const T& item) {
		if(!m_split && m_idx == m_maxCount)
			split();
		if(m_split) {
			getNode(item)->add(item);
		} else {
			if(m_idx >= m_size)
				expand();
			std::memcpy(m_data + m_idx, &item, sizeof(T));
			++m_idx;
		}
	}

	bool intersects(double bounds[4]) {
		return !(bounds[2] < m_bounds[0] && bounds[0] > m_bounds[2] && bounds[3] < m_bounds[1] && bounds[1] > m_bounds[3]);
	}

	int get(double bounds[4], std::vector<T>& items) {
		int count = 0;
		if(intersects(bounds)) {
			if(m_split) {
				for(mtree<T>* n : m_nodes) {
					if(n)
						count += n->get(bounds, items);
				}
			} else {
				reset();
				T item;
				while(next(item)) {
					if(contains(item, bounds)) {
						items.push_back(item);
						++count;
					}
				}
			}
		}
		return count;
	}

	~mtree() {
		clear();
	}

};

}
}
