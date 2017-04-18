#ifndef __QTREE_HPP__
#define __QTREE_HPP__

#include "util.hpp"
#include <unordered_map>

using namespace geo::util;

static uint64_t s_idx = 0;
static uint64_t s_bufSize = 64;

template <class T>
class CatItem {
private:
	uint64_t m_idx;
	uint64_t m_iidx;
	uint64_t m_iiidx;

public:
	uint64_t count;
	std::vector<std::pair<uint64_t, uint64_t> > items;

	CatItem() : count(0) {}

	uint64_t newPosition() {
		if(items.empty() || ++count % s_bufSize == 0)
			items.push_back(std::make_pair(0, s_idx++));
		std::pair<uint64_t, uint64_t>& item = items[items.size() - 1];
		item.first++;
		return item.second * s_bufSize * sizeof(T) + item.first * sizeof(T);
	}

	void reset() {
		m_idx = 0;
		m_iidx = 0;
		m_iiidx = 0;
	}

	bool nextPosition(uint64_t* idx) {
		while(m_idx < items.size()) {
			if(m_iidx < items[m_idx].first) {
				uint64_t pos = items[m_idx].second * s_bufSize * sizeof(T) + m_iidx * sizeof(T);
				++m_iidx;
				return pos;
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
	Bounds m_bounds;
	int m_bins;
	std::unordered_map<uint64_t, CatItem<T> > m_catalog;
	std::FILE* m_file;
	std::string m_filename;

public:
	QTree(const Bounds& bounds, int bins) :
		m_bounds(bounds), 
		m_bins(bins) {
		m_filename = Util::tmpFile();
		m_file = std::fopen(m_filename.c_str(), "wb+");
	}

	uint64_t toBinX(double x) {
		return (uint64_t) (x - m_bounds.minx()) / m_bounds.width() * m_bins;
	}

	uint64_t toBinY(double y) {
		return (uint64_t) (y - m_bounds.miny()) / m_bounds.height() * m_bins;
	}

	uint64_t morton(double x, double y) {
		uint64_t ix = toBinX(x);
		uint64_t iy = toBinY(y);
		uint64_t out = 0;
		int shift = 1;
		for(int i = 0; i < 31; ++i)
			out |= ((iy >> i) & 0xf) << shift++;
		return ix | iy;
	}

	void addItem(double x, double y, const T& item) {
		uint64_t bin = morton(x, y);
		uint64_t pos = m_catalog[bin].newPosition();
		std::fseek(m_file, pos, SEEK_SET);
		std::fwrite(&item, sizeof(T), 1, m_file);
	}

	template <class U>
	void search(double x, double y, double radius, U output) {
		uint64_t min = morton(x - radius, y - radius);
		uint64_t max = morton(x + radius, y + radius);
		radius *= radius;
		for(uint64_t i = min; i <= max; ++i) {
			if(m_catalog.find(i) != m_catalog.end()) {
				CatItem<T>& ci = m_catalog[i];
				ci.reset();
				uint64_t pos = 0;
				T item;
				while(ci.nextPosition(&pos)) {
					std::fseek(m_file, pos, SEEK_SET);
					std::fread(&item, sizeof(T), 1, m_file);
					if(std::pow(item.x - x, 2) + std::pow(item.y - y, 2) <= radius) {
						*output = item;
						++output;
					}
				}
			}
		}
	}

	~QTree() {
		std::fclose(m_file);
		Util::rm(m_filename);
	}

};

#endif