/*
 * This is wrapper around the ANN KDTree.
 *
 *  Created on: Nov 17, 2017
 *      Author: rob
 */

#ifndef INCLUDE_DS_MQTREE_HPP_
#define INCLUDE_DS_MQTREE_HPP_

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iterator>
#include <unordered_map>

#include "geo.hpp"
#include "util.hpp"

using namespace geo::util;

namespace geo {
namespace ds {

/**
 * A file-backed qtree.
 *
 * It is assumed that the template is a class with a [] operator
 * which will return a value corresponding to the dimension,
 * which is the modulus of the given index.
 */
template <class T>
class mqtree {
public:
	enum Mode {
		Memory,
		File
	};

private:

	double m_minx, m_miny;				///<! The minimum x and y coordinates. Used for shifting coordinates.
	double m_maxx, m_maxy;				///<! The minimum x and y coordinates. Used for shifting coordinates.
	int m_tileCount;					///<! The number of tiles, or files that can be created to store blocks of points.
	double m_tileSize;					///<! The tile size required to achieve the tile count given the data extent.
	int m_cols;
	int m_rows;
	std::unordered_map<size_t, int> m_files;			///<! File descriptors of temporary files. For file mode.
	std::unordered_map<size_t, std::vector<T>> m_mem;	///<! Vectors for items. For memory mode.
	std::unordered_map<size_t, int> m_counts;			///<! Number of points in each file.
	int m_size;
	Mode m_mode;										///<! Mode.

public:

	/**
	 * Construct an mqtree.
	 *
	 * \param scale The scale factor for decreasing the precision of the data.
	 * \param limit The memory threshold (bytes) that triggers the use of file-backed storage.
	 * \param minx The minimum x coordinate.
	 * \param miny The minimum y coordinate.
	 * \param minx The maximum x coordinate.
	 * \param miny The maximum y coordinate.
	 * \param maxTiles the maximum number of tiles, or files, created to store blocks of points.
	 */
	mqtree<T>(double minx, double miny, double maxx, double maxy, int tileCount, Mode mode = Memory) :
		m_minx(minx), m_miny(miny),
		m_maxx(maxx), m_maxy(maxy),
		m_tileCount(tileCount),
		m_size(0), m_mode(mode) {

		double w = m_maxx - m_minx;
		double h = m_maxy - m_miny;
		m_tileSize = std::max(w / std::sqrt(m_tileCount), h / std::sqrt(m_tileCount));
		m_cols = (int) std::ceil(w / m_tileSize);
		m_rows = (int) std::ceil(h / m_tileSize);
		m_tileCount = m_cols * m_rows;
	}

	/**
	 * Add a single item to the tree. Must (re)build before using.
	 *
	 * \param item An item.
	 */
	void add(const T& item) {

		if(!contains(item.x(), item.y()))
			return;

		int col = (int) (item.x() - m_minx) / (m_maxx - m_minx) * m_cols;
		int row = (int) (item.y() - m_miny) / (m_maxy - m_miny) * m_rows;
		size_t idx = row * m_cols + col;
		if(m_mode == File && m_files.find(idx) == m_files.end()) {
			if(!(m_files[idx] = open("/tmp", O_TMPFILE|O_RDWR|O_EXCL, S_IRWXU)))
				g_runerr("Failed to open temp file for mqtree.");
			m_counts[idx] = 0;
		}
		//lseek(m_files[idx], m_counts[idx], SEEK_SET);
		if(m_mode == File) {
			if(!write(m_files[idx], &item, sizeof(T)))
				g_runerr("Failed to write item in qtree.");
		} else {
			m_mem[idx].push_back(item);
		}
		++m_counts[idx];
		++m_size;
	}

	/**
	 * Add a vector of items to the tree. Must (re)build before using.
	 *
	 * \param items A vector of items.
	 */
	void add(std::vector<T>& items) {
		for(const T& item : items)
			add(item);
	}

	/**
	 * Clear the tree and remove all items.
	 */
	void clear() {
		if(m_mode == File) {
			for(auto& it : m_files) {
				close(it.second);
				it.second = 0;
			}
		} else {
			m_mem.clear();
		}
		m_size = 0;
	}

	/**
	 * Returns true if the tree contains the point within its bounds.
	 * If the radius is given, points within that distance are considered
	 * contained.
	 *
	 * \param x The x coordinate.
	 * \param y The y coordinate.
	 * \return True if the tree contains the point within the radius.
	 */
	bool contains(double x, double y, double radius = 0) {
		return (x >= m_minx - radius && x <= m_maxx + radius && y >= m_miny - radius && y <= m_maxy + radius);
 	}

	/**
	 * Return the number of items in the tree.
	 *
	 * \return The number of items in the tree.
	 */
	size_t size() const {
		return m_size;
	}

	/**
	 * Return the n neares neighbours to the given point.
	 *
	 * \param pt A search point.
	 * \param n The number results.
	 * \param iter An insert iterator.
	 * \return The number of items found.
	 */
	template <class TIter>
	size_t knn(const T& pt, size_t n, TIter iter) {
		g_runerr("Not implemented.")
		return 0;
	}

	/**
	 * Search for points within the given radius.
	 *
	 * \param pt The search point.
	 * \param radius The search radius.
	 * \param piter The insert iterator for result points.
	 * \return The number of items found.
	 */
	template <class TIter>
	size_t search(const T& pt, double radius, TIter piter) {

		int rr = (int) std::ceil(radius / m_tileSize);
		int col = (int) (pt.x() - m_minx) / (m_maxx - m_minx) * m_cols;
		int row = (int) (pt.y() - m_miny) / (m_maxy - m_miny) * m_rows;
		int col0 = std::max(0, col - rr);
		int col1 = std::min(m_cols - 1, col + rr);
		int row0 = std::max(0, row - rr);
		int row1 = std::min(m_rows - 1, row + rr);

		size_t count = 0;
		T pt0;
		for(int r = row0; r <= row1; ++r) {
			for(int c = col0; c <= col1; ++c) {
				size_t idx = r * m_cols + c;
				if(m_mode == File) {
					if(m_files.find(idx) != m_files.end()) {
						int ct = m_counts[idx];
						int fd = m_files[idx];
						lseek(fd, 0, SEEK_SET);
						for(int i = 0; i < ct; ++i) {
							read(fd, &pt0, sizeof(T));
							if(std::pow(pt.x() - pt0.x(), 2.0) + std::pow(pt.y() - pt0.y(), 2.0) <= radius * radius) {
								*piter = pt0;
								++piter;
								++count;
							}
						}
					}
				} else {
					for(const T& p : m_mem[idx]) {
						if(std::pow(pt.x() - p.x(), 2.0) + std::pow(pt.y() - p.y(), 2.0) <= radius * radius) {
							*piter = p;
							++piter;
							++count;
						}
					}
				}
			}
		}

		return count;
	}

	~mqtree() {
		clear();
	}

};

}
}



#endif /* INCLUDE_DS_MQTREE_HPP_ */
