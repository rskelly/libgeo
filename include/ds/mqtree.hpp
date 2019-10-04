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
#include <unordered_set>
#include <list>

#include "geo.hpp"
#include "util.hpp"

using namespace geo::util;

namespace geo {
namespace ds {

template <class T>
class lrunode {
public:

	std::vector<T> cache;
	lrunode* left;
	lrunode* right;
	std::string path;

	void load(const std::string& path) {
		if(path == this->path)
			return;

		this->path = path;

		int handle, pos;
		size_t size;

		if((handle = open(path.c_str(), O_RDONLY|O_EXCL, S_IRWXU)) > 0) {
			if((pos = lseek(handle, 0, SEEK_END)) > 0) {
				lseek(handle, 0, SEEK_SET);
				if(read(handle, &size, sizeof(size_t)) > 0) {
					cache.resize(size);
					read(handle, cache.data(), size * sizeof(T));
				}
			}
			close(handle);

		} else {

			cache.resize(0);

		}
	}

	void flush() {

		if(cache.empty())
			return;

		int handle, pos;
		size_t size = cache.size();

		if((handle = open(path.c_str(), O_WRONLY|O_EXCL|O_CREAT|O_TRUNC, S_IRWXU)) > -1) {

			if(!write(handle, &size, sizeof(size_t)))
				std::cerr << "Failed to write cache size.\n";
			if(!write(handle, cache.data(), size * sizeof(T)))
				std::cerr << "Failed to write cache.\n";
			close(handle);

		}

		cache.resize(0);

	}

};

template <class T>
class mqnode;

template <class T>
class lru {
private:
	lrunode<T>* m_first;
	lrunode<T>* m_last;
	std::unordered_map<std::string, lrunode<T>*> m_nodes;
	std::unordered_set<std::string> m_files;
	size_t m_size;

	void movelast(lrunode<T>* node) {
		if(!m_last) {
			m_first = m_last = node;
		} else if(node != m_last) {
			if(node->left)
				node->left->right = node->right;
			if(node->right)
				node->right->left = node->left;
			m_last->right = node;
			node->left = m_last;
			node->right = nullptr;
			m_last = node;
		}
	}

public:

	lru(size_t size = 10) :
		m_first(nullptr), m_last(nullptr),
		m_size(size) {
	}

	std::vector<T>& cache(mqnode<T>* node) {
		const std::string& path = node->path();
		lrunode<T>* n = nullptr;
		if(m_nodes.find(path) != m_nodes.end()) {
			n = m_nodes[path];
		} else if(m_nodes.size() < m_size) {
			n = new lrunode<T>();
			m_nodes[path] = n;
			n->load(path);
		} else {
			n = m_first;
			n->flush();
			m_nodes.erase(n->path);
			m_nodes[n->path] = n;
			n->load(path);
		}
		if(n != m_last)
			movelast(n);
		return n->cache;
	}

	void remove(mqnode<T>* node) {
		util::rem(node->path());
	}

	void clear() {
		for(const std::string& path: m_files)
			util::rem(path);
		m_files.clear();
		for(auto& it : m_nodes)
			delete it.second;
		m_nodes.clear();
	}

	~lru() {
		clear();
	}
};

template <class T>
class mqtree;

template <class T>
class mqnode {
public:
	lru<T>* m_cache;
	mqnode* m_parent;
	double m_bounds[4];
	double m_midx;
	double m_midy;
	mqnode* m_nodes[4];
	int m_depth;
	size_t m_maxSize;
	int m_maxDepth;
	char m_idx;
	bool m_split;
	std::string m_path;
	size_t m_size;

	mqnode(char idx, double minx, double miny, double maxx, double maxy,
			int depth, size_t maxSize, int maxDepth,
			mqnode<T>* parent, lru<T>* cache) :
		m_cache(cache), m_parent(parent),
		m_depth(depth), m_maxSize(maxSize), m_maxDepth(maxDepth),
		m_idx(idx), m_split(false),
		m_size(0) {

		m_midx = (minx + maxx) / 2.0;
		m_midy = (miny + maxy) / 2.0;
		m_bounds[0] = minx;
		m_bounds[1] = miny;
		m_bounds[2] = maxx;
		m_bounds[3] = maxy;
		for(int i = 0; i < 4; ++i)
			m_nodes[i] = nullptr;
	}

	const std::string& path() {
		if(m_path.empty()) {
			std::list<int> ids;
			mqnode* n = this;
			while(n) {
				ids.push_back(n->m_idx);
				n = n->m_parent;
			}
			std::reverse(ids.begin(), ids.end());
			std::stringstream ss;
			for(int id : ids)
				ss << "_" << id;
			m_path = "/tmp/tree" + ss.str();
		}
		return m_path;
	}

	void add(const T& item) {
		if(m_split) {
			node(idx(item))->add(item);
			++m_size;
		} else if(m_size >= m_maxSize && m_depth < m_maxDepth) {
			split();
			node(idx(item))->add(item);
			++m_size;
		} else {
			m_cache->cache(this).push_back(item);
			++m_size;
		}
	}

	size_t size() const {
		return m_size;
	}

	char idx(const T& item) {
		return ((char) (item.x() >= m_midx) << 1) | (char) (item.y() >= m_midy);
	}

	mqnode<T>* node(char i) {
		char ix = i >> 1;
		char iy = i & 1;
		mqnode* n;
		if(!(n = m_nodes[i])) {
			n = m_nodes[i] = new mqnode(i,
					ix ? m_midx : m_bounds[0],
					iy ? m_midy : m_bounds[1],
					ix ? m_bounds[2] : m_midx,
					iy ? m_bounds[3] : m_midy,
					m_depth + 1, m_maxSize, m_maxDepth,
					this, m_cache
			);
		}
		return n;
	}

	void split() {
		std::vector<T>& c = m_cache->cache(this);
		for(T& item : c)
			node(idx(item))->add(item);
		// m_cache->remove(this);
		m_split = true;
	}

	bool contains(const T& item, double radius = 0) {
		return (item.x() + radius) >= m_bounds[0] && (item.x() - radius) < m_bounds[2] &&
				(item.y() + radius) >= m_bounds[1] && (item.y() - radius) < m_bounds[3];
	}

	void clear() {
		for(char i = 0; i < 4; ++i) {
			if(m_nodes[i])
				m_nodes[i]->clear();
		}
		m_cache->remove(this);
	}

	template <class TIter>
	size_t search(const T& pt, double radius, TIter& piter) {
		size_t count = 0;
		if(contains(pt, radius)) {
			if(m_split) {
				for(char i = 0; i < 4; ++i) {
					if(m_nodes[i])
						count += m_nodes[i]->search(pt, radius, piter);
				}
			} else {
				for(const T& item : m_cache->cache(this)) {
					if(std::pow(item.x() - pt.x(), 2.0) + std::pow(item.y() - pt.y(), 2.0) <= radius * radius) {
						*piter = item;
						++piter;
						++count;
					}
				}
			}
		}
		return count;
	}

	~mqnode() {
		for(char i = 0; i < 4; ++i) {
			if(m_nodes[i])
				delete m_nodes[i];
		}
	}

};

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
	int m_maxDepth;						///<! The maximum tree depth.
	double m_tileSize;					///<! The tile size required to achieve the tile count given the data extent.
	std::unordered_map<size_t, std::vector<T>> m_mem;	///<! Vectors for items. For memory mode.
	std::unordered_map<size_t, int> m_counts;			///<! Number of points in each file.
	int m_size;
	Mode m_mode;										///<! Mode.
	lru<T> m_lru;
	mqnode<T>* m_root;

	/**
	 * Get the cache object associated with the index.
	 *
	 * Return file-backed or memory version depending on mode.
	 */
	std::vector<T>& cache(size_t idx) {
		if(m_mode == File) {
			return m_lru.cache(idx);
		} else {
			return m_mem[idx];
		}
	}

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
	 * \param maxDepth The maximum depth of the tree. Zero is no limit.
	 * \param mode Determines whether file-backed storage is used, or memory.
	 */
	mqtree<T>(double minx, double miny, double maxx, double maxy, int maxDepth = 0, Mode mode = Memory) :
		m_minx(minx), m_miny(miny),
		m_maxx(maxx), m_maxy(maxy),
		m_maxDepth(maxDepth),
		m_size(0), m_mode(mode) {

		double w = m_maxx - m_minx;
		double h = m_maxy - m_miny;
		if(w < h) {
			w = h;
		} else {
			h = w;
		}
		m_tileSize = std::max(m_maxx - m_minx, m_maxy - m_miny);
		m_root = new mqnode<T>(0, m_minx, m_miny, m_minx + m_tileSize, m_miny + m_tileSize, 0, 10000, 10, nullptr, &m_lru);
	}

	/**
	 * Add a single item to the tree. Must (re)build before using.
	 *
	 * \param item An item.
	 */
	void add(const T& item) {
		m_root->add(item);
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
		m_root->clear();
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
	bool contains(const T& item, double radius = 0) {
		return m_root->contains(item, radius);
 	}

	/**
	 * Return the number of items in the tree.
	 *
	 * \return The number of items in the tree.
	 */
	size_t size() const {
		return m_root->size();
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
	size_t search(const T& pt, double radius, TIter& piter) {
		return m_root->search(pt, radius, piter);
	}

	~mqtree() {
		m_lru.clear();
		delete m_root;
	}

};

}
}



#endif /* INCLUDE_DS_MQTREE_HPP_ */
