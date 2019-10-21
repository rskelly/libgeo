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
#include <list>

#include "geo.hpp"
#include "util.hpp"

#define BUF_SIZE (512 * 1024)

using namespace geo::util;

namespace {

int mkpath(const char* file_path, mode_t mode) {
    //assert(file_path && *file_path);
    for (char* p = strchr((char*) file_path + 1, '/'); p; p = strchr(p + 1, '/')) {
        *p = '\0';
        if (mkdir(file_path, mode) == -1) {
            if (errno != EEXIST) {
                *p = '/';
                return -1;
            }
        }
        *p = '/';
    }
    return 0;
}

}


namespace geo {
namespace ds {

template <class T>
class lrucache;

template <class T>
class lrunode {
private:
	lrucache<T>* m_lru;
	std::vector<T> m_cache;
	std::string m_path;

public:

	lrunode* left;
	lrunode* right;

	lrunode(lrucache<T>* lru) :
		m_lru(lru), left(nullptr), right(nullptr) {
	}

	std::vector<T>& cache() {
		return m_cache;
	}

	const std::string& path() const {
		return m_path;
	}

	void path(const std::string& path) {
		m_path = path;
	}

	void load(const std::string& path) {
		if(path != m_path) {

			flush();

			m_path = path;

			static std::vector<char> buf(BUF_SIZE);

			int handle;
			size_t size;

			if((handle = open(path.c_str(), O_RDWR, 0777)) > 0) {
				if(read(handle, buf.data(), BUF_SIZE) > sizeof(T)) {
					std::memcpy(&size, buf.data(), sizeof(size_t));
					m_cache.resize(size);
					std::memcpy(m_cache.data(), buf.data() + sizeof(size_t), size * sizeof(T));
				}
				close(handle);
			} else if(errno != ENOENT){
				std::cerr << "Failed to open file for load: " << path << " (" << strerror(errno) << ")\n";
			}
		}
	}

	void flush() {
		if(!m_cache.empty()) {

			static std::vector<char> buf(BUF_SIZE);

			int handle;
			size_t size = m_cache.size();

			if(!mkpath(m_path.c_str(), 0777)) {
				if((handle = open(m_path.c_str(), O_RDWR|O_CREAT|O_TRUNC, 0777)) > 0) {
					std::memcpy(buf.data(), &size, sizeof(size_t));
					std::memcpy(buf.data() + sizeof(size_t), m_cache.data(), size * sizeof(T));
					if(!write(handle, buf.data(), BUF_SIZE))
						std::cerr << "Failed to write cache.\n";
					close(handle);
				} else {
					std::cerr << "Failed to open file for flush: " << m_path << " (" << strerror(errno) << ")\n";
				}
			} else {
				std::cerr << "Failed to create dir for flush: " << m_path << " (" << strerror(errno) << ")\n";
			}

			m_cache.resize(0);
		}
	}

};

template <class T>
class mqnode;

template <class T>
class lrucache {
private:
	std::unordered_map<std::string, lrunode<T>*> m_nodes;	///<! The mapping from path to node instance.
	lrunode<T>* m_first;									///<! The first (oldest) node.
	lrunode<T>* m_last;										///<! The last (newest) node.
	size_t m_size;											///<! The number of nodes allowed before swapping.
	size_t m_blkSize;
	size_t m_stats[3];										///<! Statistics for cache hits, new nodes and swapped nodes.

	/**
	 * \brief Move the given node to the end of the list.
	 */
	void movelast(lrunode<T>* node) {
		if(!m_last) {
			m_first = m_last = node;
		} else if(node == m_first) {
			m_first = node->right;
			m_first->left = nullptr;
			m_last->right = node;
			node->left = m_last;
			node->right = nullptr;
			m_last = node;
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

	/**
	 * \brief Create the LRU cache with the given maximum number of elements.
	 */
	lrucache(size_t size = 100, size_t blkSize = 4096) :
		m_first(nullptr), m_last(nullptr),
		m_size(size), m_blkSize(blkSize) {

		for(int i = 0; i < 3; ++i)
			m_stats[i] = 0;
	}

	void blkSize(size_t blkSize) {
		m_blkSize = blkSize;
	}

	size_t blkSize() const {
		return m_blkSize;
	}

	/**
	 * \brief Get the cache array for the given node.
	 *
	 * Will be returned if it's in memory, created if it doesn't exist, loaded from file if it does.
	 * Node is moved to the end of the LRU list.
	 */
	std::vector<T>& cache(mqnode<T>* node) {
		const std::string& path = node->path();
		lrunode<T>* n = nullptr;
		if(m_nodes.find(path) != m_nodes.end()) {
			// The node is in the list. Select it.
			std::cout << path << "\n";
			n = m_nodes[path];
			m_stats[0]++;
		} else if(m_nodes.size() < m_size) {
			// The node is not in the list and there's room. Create one.
			n = new lrunode<T>(this);
			m_nodes[path] = n;
			n->path(path);
			m_stats[1]++;
		} else {
			// The node is not in the list and there's no room. Swap the first node.
			n = m_first;
			m_nodes.erase(n->path());
			m_nodes[path] = n;
			n->load(path);
			m_stats[2]++;
		}
		// If the node isn't last, move it up.
		if(n != m_last)
			movelast(n);
		return n->cache();
	}

	/**
	 * \brief Write all nodes' contents to file.
	 */
	void flush() {
		for(const auto& it : m_nodes)
			it.second->flush();
	}

	/**
	 * \brief Delete the node's file.
	 */
	void remove(mqnode<T>* node) {
		util::rem(node->path());
	}

	/**
	 * \brief Delete all nodes and remove their files.
	 */
	void clear() {
		for(auto& it : m_nodes)
			delete it.second;
		m_nodes.clear();
	}

	/**
	 * \brief Print statistics for cache hits, misses and switches.
	 *
	 * A hit is when the node is in memory. A miss or switch is when the cache hasn't been used
	 * and needs to be created. A switch happens when the LRU list is full and an element
	 * is reused for the new cache.
	 */
	void printStats() const {
		size_t sum = m_stats[0] + m_stats[1] + m_stats[2];
		if(sum) {
			std::cout << std::setprecision(2) << std::fixed;
			std::cout << "Stats: hits: " << ((float) m_stats[0] / sum * 100) << "%; loads: " << ((float) m_stats[1] / sum * 100) << "%; switches: " << ((float) m_stats[2] / sum * 100) << "%\n";
		}
	}

	~lrucache() {
		clear();
	}
};

template <class T>
class mqtree;

template <class T>
class mqnode {
public:
	mqtree<T>* m_tree;
	mqnode* m_parent;							///<! Pointer to this node's parent.
	mqnode* m_nodes[4];							///<! Pointers to this node's children (if instantiated).
	double m_bounds[4];							///<! Geographic bounds of node.
	double m_midx;								///<! Lateral midpoint of bounds.
	double m_midy;								///<! Vertical midpoint of bounds.
	int m_depth;								///<! Depth of this node. Root is zero.
	char m_idx;									///<! The index of this node relative to the parent (0-3).
	bool m_split;								///<! True if the node has been split.
	std::string m_path;							///<! The path to this node's temp file.
	size_t m_size;								///<! The current number of items in the node or its children.

	/**
	 * \brief Create an mqtree node.
	 *
	 * \param idx The index of the ndoe in the children list of the parent.
	 * \param minx The minimum x coordinate of the bounding box.
	 * \param miny The minimum y coordinate of the bounding box.
	 * \param maxx  The maximum x coordinate of the bounding box.
	 * \param maxy The maximum y coordinate of the bounding box.
	 * \param dept The depth of this node relative to the root.
	 * \param maxSize The maximum number of items allowed in a node before splitting.
	 * \param maxDepth The maximum depth of the tree.
	 * \param parent The parent node; nullptr if root.
	 * \param cache The pointer to the LRU cache.
	 */
	mqnode(char idx, double minx, double miny, double maxx, double maxy, int depth,
			mqnode<T>* parent, mqtree<T>* tree) :
		m_parent(parent), m_depth(depth), m_idx(idx), m_split(false),
		m_size(0), m_tree(tree) {

		m_midx = (minx + maxx) / 2.0;
		m_midy = (miny + maxy) / 2.0;
		m_bounds[0] = minx;
		m_bounds[1] = miny;
		m_bounds[2] = maxx;
		m_bounds[3] = maxy;
		for(int i = 0; i < 4; ++i)
			m_nodes[i] = nullptr;
	}

	/**
	 * \brief Compute the path from parent to this node as a file in /tmp.
	 *
	 * Only occurs once. Is stored thereafter.
	 */
	const std::string& path() {
		if(m_path.empty()) {

			// Generate the path list and reverse.
			std::list<int> ids;
			mqnode* n = this;
			while(n) {
				ids.push_back(n->m_idx);
				n = n->m_parent;
			}
			std::reverse(ids.begin(), ids.end());

			// Compile the root path.
			std::stringstream ss;
			ss << m_tree->rootPath();

			// Add the subdirs and filename.
			for(int id : ids)
				ss << "/" << id;
			ss << "/dat";

			m_path = ss.str();
		}
		return m_path;
	}

	/**
	 * \brief Add the item to the tree.
	 */
	void add(const T& item) {
		if(m_split) {
			node(idx(item))->add(item);
			++m_size;
		} else if(m_size >= m_tree->maxCount()) {
			split();
			node(idx(item))->add(item);
			++m_size;
		} else {
			m_tree->lru().cache(this).push_back(item);
			++m_size;
		}
	}

	/**
	 * \brief Return the size of the tree.
	 */
	size_t size() const {
		return m_size;
	}

	/**
	 * \brief Return the index (i.e. quadrant) of the given item.
	 */
	char idx(const T& item) {
		return ((char) (item.x() >= m_midx) << 1) | (char) (item.y() >= m_midy);
	}

	/**
	 * \brief Return the node from the given index (quadrant). Create if required.
	 */
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
					m_depth + 1, this, m_tree
			);
		}
		return n;
	}

	/**
	 * \brief When the node's capacity is exceeded, split it into quadrants.
	 */
	void split() {
		if(m_depth >= m_tree->maxDepth())
			std::cerr << "Max depth exceeded: " << m_depth << "\n";
		std::vector<T>& c = m_tree->lru().cache(this);
		m_tree->lru().remove(this); // The file is removed before being replaced with a directory.
		for(T& item : c)
			node(idx(item))->add(item);
		m_split = true;
	}

	/**
	 * \brief Returns true of the node contains the given point (considers the radius if given).
	 */
	bool contains(const T& item, double radius = 0) {
		return item.x() >= (m_bounds[0] - radius) && item.x() <= (m_bounds[2] + radius) &&
				item.y() >= (m_bounds[1] - radius) && item.y() <= (m_bounds[3] + radius);
	}

	/**
	 * \brief Clear the tree and remove its nodes.
	 */
	void clear() {
		for(char i = 0; i < 4; ++i) {
			if(m_nodes[i])
				m_nodes[i]->clear();
		}
		m_tree->lru().remove(this);
	}

	/**
	 * \brief Perform a radius search of the tree, add all results to the iterator. Return the count.
	 */
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
				double r2 = radius * radius;
				double d;
				const std::vector<T>& c = m_tree->lru().cache(this);
				for(const T& item : c) {
					d = std::pow(item.x() - pt.x(), 2.0) + std::pow(item.y() - pt.y(), 2.0);
					if(d <= r2) {
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
private:

	int m_maxDepth;						///<! The maximum tree depth.
	size_t m_maxCount;					///<! Max number of items in each cell.
	lrucache<T> m_lru;
	mqnode<T>* m_root;
	std::string m_rootPath;
	int m_blkSize;

public:

	mqtree<T>() :
		m_maxDepth(0),
		m_maxCount(0),
		m_root(nullptr),
		m_blkSize(0) {
	}

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
	mqtree<T>(double minx, double miny, double maxx, double maxy, int maxDepth = 100) {
		init(minx, miny, maxx, maxy, maxDepth);
	}

	void init(double minx, double miny, double maxx, double maxy, int maxDepth = 100) {
		m_maxDepth = maxDepth;

		// Create a root path.
		m_rootPath = util::tmpdir("/tmp/mqtreeXXXXXX");

		// Get the side length of the table region.
		double side = std::max(maxx - minx, maxy - miny);

		// Get the block size for IO.
		struct stat st;
		stat(m_rootPath.c_str(), &st);
		m_blkSize = st.st_blksize;

		// The max count takes the block size into consideration. It should
		// be a multipl of the ideal block size, as near as possible.
		// The size of a size_t is subtracted because that's the space
		// required for the count at the start of the file.
		m_maxCount = (BUF_SIZE - sizeof(size_t)) / sizeof(T);

		m_root = new mqnode<T>(0, minx, miny, minx + side, miny + side, 0, nullptr, this);
	}

	size_t maxCount() const {
		return m_maxCount;
	}

	int maxDepth() const {
		return m_maxDepth;
	}

	size_t blkSize() const {
		return m_blkSize;
	}

	lrucache<T>& lru() {
		return m_lru;
	}

	const std::string& rootPath() const {
		return m_rootPath;
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

	void printStats() {
		m_lru.printStats();
	}

	~mqtree() {
		m_lru.clear();
		delete m_root;
		util::rem(m_rootPath);
	}

};

}
}



#endif /* INCLUDE_DS_MQTREE_HPP_ */
