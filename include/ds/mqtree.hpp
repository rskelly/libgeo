/*
 * This is wrapper around the ANN KDTree.
 *
 *  Created on: Nov 17, 2017
 *      Author: rob
 */

#ifndef INCLUDE_DS_MQTREE_HPP_
#define INCLUDE_DS_MQTREE_HPP_

#define NOMINMAX	// To prevent windows redefining these.

#ifdef _WIN32
#include <Windows.h>
#include <io.h>
typedef int mode_t;

#else
#include <unistd.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <list>
#include <atomic>
#include <algorithm>
#include <queue>

#include "geo.hpp"
#include "util.hpp"

#define BUF_SIZE 1024

using namespace geo::util;

namespace geo {
namespace ds {

/**
 * \brief Loads/stores cached elements for the lrucache.
 */
template <class T>
class lrunode {
private:
	bool m_locked;				///<! If the node is in use is should not be expirable. This is a problem in tree traversals. Lock to prevent.
	size_t m_key;				///<! This node's key in the lookup table.
	std::string m_path;			///<! The path to the storage file.
	std::vector<T> m_cache;		///<! The list of elements stored in the leaf.

public:

	lrunode* left;				///<! The left (older) node.
	lrunode* right;				///<! The right (newer) node.

	lrunode() :
		m_locked(false), m_key(0),
		left(nullptr), right(nullptr) {
	}

	static size_t bufSize() {
		// The buffer size is proportional to the size of T + plus the length variable.
		return sizeof(T) * BUF_SIZE + sizeof(size_t);
	}

	void key(size_t key) {
		m_key = key;
	}

	size_t key() const {
		return m_key;
	}

	void lock() {
		m_locked = true;
	}

	void unlock() {
		m_locked = false;
	}

	bool locked() const {
		return m_locked;
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
			int handle;
			size_t size;
			m_path = path;
			std::vector<char> buf(bufSize());
			if((handle = open(path.c_str(), O_RDWR, 0777)) > 0) {
				if(read(handle, buf.data(), buf.size()) > (int) sizeof(T)) {
					std::memcpy(&size, buf.data(), sizeof(size_t));
					m_cache.resize(size);
					std::memcpy((void*) m_cache.data(), buf.data() + sizeof(size_t), size * sizeof(T));
					// Offset by sizeof(size_t) because the size is the first element in the file;
					// cast to char* because the offset is in bytes.
				}
				close(handle);
			} else if(errno != ENOENT){
				g_runerr("Failed to open file for load: " << path << " (" << strerror(errno) << ")");
			}
		}
	}

	void flush() {
		if(!m_cache.empty()) {
			int handle;
			size_t size = m_cache.size();
			std::vector<char> buf(bufSize());
			if(makedir(m_path)) {
				if((handle = open(m_path.c_str(), O_RDWR|O_CREAT|O_TRUNC, 0777)) > 0) {
					std::memcpy(buf.data(), &size, sizeof(size_t));
					std::memcpy(buf.data() + sizeof(size_t), m_cache.data(), size * sizeof(T));
					// Cast to char* because the offset is in bytes.
					if(!write(handle, buf.data(), buf.size()))
						g_runerr("Failed to write cache.");
					close(handle);
				} else {
					g_runerr("Failed to open file for flush: " << m_path << " (" << strerror(errno) << ")");
				}
			} else {
				g_runerr("Failed to create dir for flush: " << m_path << " (" << strerror(errno) << ")");
			}
			m_cache.resize(0);
		}
	}

	~lrunode() {
		util::rem(m_path);
	}

};

/**
 * \brief Maintains a cache of elements based on age of last use. Oldest elements are expired to make
 * way for new ones.
 */
template <class T>
class lrucache {
private:
	std::unordered_map<size_t, lrunode<T>*> m_nodes;	///<! The mapping from path to node instance.
	lrunode<T>* m_first;								///<! The first (oldest) node.
	lrunode<T>* m_last;									///<! The last (newest) node.
	size_t m_size;										///<! The number of nodes allowed before swapping.
	size_t m_blkSize;									///<! The size of data blocks used for storing to file.
	size_t m_stats[3];									///<! Statistics for cache hits, new nodes and swapped nodes.

	/**
	 * \brief Move the given node to the end of the list.
	 *
	 * This is done when a node is accessed to mark it as recently-used.
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
	lrucache(size_t size = 1000, size_t blkSize = 4096) :
		m_first(nullptr), m_last(nullptr),
		m_size(size), m_blkSize(blkSize) {

		for(int i = 0; i < 3; ++i)
			m_stats[i] = 0;
	}

	void size(size_t size) {
		m_size = size;
	}

	void blkSize(size_t blkSize) {
		m_blkSize = blkSize;
	}

	size_t blkSize() const {
		return m_blkSize;
	}

	/**
	 * \brief Return a copy of the cache and destroy the lrunode.
	 */
	std::vector<T> detach(size_t key, const std::string& path) {
		lrunode<T>* n = get(key, path);
		std::vector<T> cache;
		if(n) {
			cache = std::move(n->cache());
			if(n->left) {
				n->left->right = nullptr;
				m_last = n->left;
			} else {
				m_last = m_first = nullptr;
			}
			m_nodes.erase(key);
			remove(path);
			delete n;
		}
		return cache;
	}


	/**
	 * \brief Get the cache array for the given node.
	 *
	 * Will be returned if it's in memory, created if it doesn't exist, loaded from file if it does.
	 * Node is moved to the end of the LRU list.
	 */
	lrunode<T>* get(size_t key, const std::string& path) {
		lrunode<T>* n = nullptr;
		if(m_nodes.find(key) != m_nodes.end()) {
			// The node is in the list. Select it.
			n = m_nodes.at(key);
			n->lock();
			m_stats[0]++;
		} else if(m_nodes.size() < m_size) {
			// The node is not in the list and there's room. Create one.
			n = new lrunode<T>();
			n->path(path);
			n->key(key);
			n->lock();
			m_nodes.insert(std::make_pair(key, n));
			m_stats[1]++;
		} else {
			// The node is not in the list and there's no room. Swap the first node.
			n = m_first;
			while(n && n->locked())
				n = n->right;
			if(!n || n->locked()) {
				std::cerr << "No available LRU nodes. Expanding.\n";
				n = new lrunode<T>();
				n->path(path);
				n->key(key);
				n->lock();
				m_nodes.insert(std::make_pair(key, n));
				m_stats[1]++;
			} else {
				m_nodes.erase(n->key());
				n->load(path);
				n->key(key);
				n->lock();
				m_nodes.insert(std::make_pair(key, n));
				m_stats[2]++;
			}
		}
		// If the node isn't last, move it up.
		if(n != m_last)
			movelast(n);
		return n;
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
	void remove(const std::string& path) {
		util::rem(path);
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
	mqnode<T>* m_parent;						///<! Pointer to this node's parent.
	mqnode<T>* m_nodes[4];						///<! Pointers to this node's children (if instantiated).
	double m_bounds[4];							///<! Geographic bounds of node.
	double m_midx;								///<! Lateral midpoint of bounds.
	double m_midy;								///<! Vertical midpoint of bounds.
	int m_depth;								///<! Depth of this node. Root is zero.
	int m_idx;									///<! The index of this node relative to the parent (0-3).
	bool m_split;								///<! True if the node has been split.
	std::string m_path;							///<! The path to this node's temp file.
	size_t m_key;								///<! A numerical key for locating the node in the cache.
	size_t m_size;								///<! The current number of items in the node or its children.
	size_t m_iter;								///<! An index to keep track of iteration.

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
	 * \param tree A pointer to the tree.
	 */
	mqnode(int idx, double minx, double miny, double maxx, double maxy, int depth,
			mqnode<T>* parent, mqtree<T>* tree) :
		m_tree(tree), m_parent(parent), m_depth(depth), m_idx(idx), m_split(false),
		m_key(tree->nextKey()), m_size(0), m_iter(0) {

		m_midx = (minx + maxx) / 2.0;
		m_midy = (miny + maxy) / 2.0;
		m_bounds[0] = minx;
		m_bounds[1] = miny;
		m_bounds[2] = maxx;
		m_bounds[3] = maxy;
		for(int i = 0; i < 4; ++i)
			m_nodes[i] = nullptr;
	}

	uint64_t key() {
		return m_key;
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
			lrunode<T>* n = m_tree->lru().get(m_key, path());
			n->cache().push_back(item);
			n->unlock();
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
	int idx(const T& item) {
		return ((int) (item.x() >= m_midx) << 1) | (int) (item.y() >= m_midy);
	}

	/**
	 * \brief Return the node from the given index (quadrant). Create if required.
	 */
	mqnode<T>* node(int i) {
		int ix = i >> 1;
		int iy = i & 1;
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

	double minx() const {
		return m_bounds[0];
	}

	double miny() const {
		return m_bounds[1];
	}

	double maxx() const {
		return m_bounds[2];
	}

	double maxy() const {
		return m_bounds[3];
	}

	/**
	 * \brief When the node's capacity is exceeded, split it into quadrants.
	 */
	void split() {
		if(m_depth >= m_tree->maxDepth())
			std::cerr << "Max depth exceeded: " << m_depth << "\n";
		std::vector<T> c = m_tree->lru().detach(m_key, path());
		for(T& item : c)
			node(idx(item))->add(item);
		m_split = true;
		m_iter = 0;
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
		for(int i = 0; i < 4; ++i) {
			if(m_nodes[i])
				m_nodes[i]->clear();
		}
		m_tree->lru().remove(path());
	}

	/**
	 * \brief Perform a radius search of the tree, add all results to the iterator. Return the count.
	 */
	template <class TIter>
	size_t search(const T& pt, double radius, TIter& piter) {
		size_t count = 0;
		if(contains(pt, radius)) {
			if(m_split) {
				for(int i = 0; i < 4; ++i) {
					if(m_nodes[i])
						count += m_nodes[i]->search(pt, radius, piter);
				}
			} else {
				double r2 = radius * radius;
				double d;
				lrunode<T>* n = m_tree->lru().get(m_key, path());
				const std::vector<T>& c = n->cache();
				for(const T& item : c) {
					d = geo::sq(item.x() - pt.x()) + geo::sq(item.y() - pt.y());
					if(d <= r2) {
						piter = item;
						++count;
					}
				}
				n->unlock();
			}
		}
		return count;
	}

	/**
	 * \brief Container class for priority_queue in knn.
	 */
	class kqitem {
	public:
		T item;				// A leaf item.
		mqnode<T>* node;	// A node.

		/**
		 * Construct the instance with a node.
		 */
		kqitem(mqnode<T>* node) :
			node(node) {}

		/**
		 * Construct the instance with a leaf.
		 */
		kqitem(T& item) :
			item(item), node(nullptr) {}

		/**
		 * Return true if this is a leaf; false if a node.
		 */
		bool isItem() const {
			return node == nullptr;
		}

		/**
		 * Return the distance from the point to the contained object.
		 */
		double dist(const T& pt) const {
			if(node) {
				// TODO: Inefficient.
				double minx = node->minx();
				double miny = node->miny();
				double maxx = node->maxx();
				double maxy = node->maxy();
				double d0 = linedist(pt.x(), pt.y(), minx, miny, minx, maxy);
				double d1 = linedist(pt.x(), pt.y(), minx, miny, maxx, miny);
				double d2 = linedist(pt.x(), pt.y(), maxx, miny, maxx, miny);
				double d3 = linedist(pt.x(), pt.y(), minx, maxy, maxx, maxy);
				return std::min(std::min(d0, d1), std::min(d2, d3));
			} else {
				return std::sqrt(geo::sq(item.x() - pt.x()) + geo::sq(item.y() - pt.y()));
			}
		}

	};

	class kqcomp {
	public:
		T pt;
		kqcomp(const T& pt) :
			pt(pt) {}
		bool operator()(const kqitem& a, const kqitem& b) {
			return a.dist(pt) < b.dist(pt);
		}
	};

	template <class TIter>
	size_t knn(const T& pt, int count, TIter& piter) {
		// 1) Add nodes to queue.
		// 2) Unpack node at top of queue. Repeat until top of queue contains leaves.
		// 3) Remove leaves. If count satisfied, return, else goto 3.
		kqcomp comp(pt);
		std::priority_queue<kqitem, std::vector<kqitem>, kqcomp> q(comp);
		q.emplace(this);
		int found = 0;
		while(!q.empty() && found < count) {
			if(q.top().isItem()) {
				// If the front of the queue is an item, add to the output and continue.
				piter = q.top().item;
				q.pop();
				++found;
				continue;
			} else {
				// If the front of the queue is a node, expand it and continue.
				mqnode<T>* node = q.top().node;
				if(!node->m_split) {
					node->reset();
					T item;
					while(node->next(item))
						q.emplace(item);
				} else {
					for(int i = 0; i < 4; ++i) {
						if(node->m_nodes[i])
							q.emplace(node->m_nodes[i]);
					}
				}
			}
		}
		return found;
	}

	void reset() {
		m_iter = 0;
		for(int i = 0; i < 4; ++i) {
			if(m_nodes[i])
				m_nodes[i]->reset();
		}
	}

	/**
	 * \brief Return the next item in a depth-first traversal.
	 *
	 * \param The item to populate.
	 * \return True if an item was found.
	 */
	bool next(T& item) {
		if(!m_split) {
			if(m_iter < m_size) {
				// Get the list and return the ndex item in it.
				lrunode<T>* n = m_tree->lru().get(m_key, path());
				const std::vector<T>& c = n->cache();
				if(!c.empty()) {
					item = c[m_iter];
					++m_iter;
					return true;
				}
			}
		} else {
			// Check each node in turn, if it returns false,
			// try the next. If all return false, return false.
			while(m_iter < 4) {
				if(m_nodes[m_iter] && m_nodes[m_iter]->next(item))
					return true;
				++m_iter;
			}
		}
		return false;
	}

	~mqnode() {
		for(int i = 0; i < 4; ++i) {
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
	friend class mqnode<T>;
private:

	int m_maxDepth;						///<! The maximum tree depth.
	size_t m_maxCount;					///<! Max number of items in each cell.
	lrucache<T> m_lru;
	mqnode<T>* m_root;
	std::string m_rootPath;
	int m_blkSize;
	uint64_t m_nextKey;					///<! Retrieved by child nodes for identification.

protected:

	uint64_t nextKey() {
		return ++m_nextKey;
	}

public:

	mqtree<T>() :
		m_maxDepth(0),
		m_maxCount(0),
		m_root(nullptr),
		m_blkSize(0),
		m_nextKey(0) {
	}

	/**
	 * Construct an mqtree.
	 *
	 * \param minx The minimum x coordinate.
	 * \param miny The minimum y coordinate.
	 * \param minx The maximum x coordinate.
	 * \param miny The maximum y coordinate.
	 * \param cacheSize The LRU nodes stored in memory. Smaller for less memory/more disk use.
	 * \param maxDepth The maximum depth of the tree. Zero is no limit.
	 */
	mqtree<T>(double minx, double miny, double maxx, double maxy, int cacheSize, int maxDepth = 100) :
		mqtree<T>() {
		init(minx, miny, maxx, maxy, cacheSize, maxDepth);
	}

	void init(double minx, double miny, double maxx, double maxy, int cacheSize, int maxDepth = 100) {
		m_maxDepth = maxDepth;
		m_lru.size(cacheSize);

		// Create a root path.
		m_rootPath = util::tmpdir("mqtree");

		// Get the side length of the table region.
		double side = geo::max(maxx - minx, maxy - miny);

		// Get the block size for IO.
#ifdef _WIN32
		m_blkSize = 4096; // TODO: Hack. Find out how windows does this.
#else
		struct stat st;
		stat(m_rootPath.c_str(), &st);
		m_blkSize = st.st_blksize;
#endif

		// The max count takes the block size into consideration.
		m_maxCount = (lrunode<T>::bufSize() - sizeof(size_t)) / sizeof(T);
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
	 * \brief Reset the iterator.
	 */
	void reset() {
		m_root->reset();
	}

	/**
	 * \brief Return the next item in a depth-first traversal.
	 *
	 * \param The item to populate.
	 * \return True if an item was found.
	 */
	bool next(T& item) {
		return m_root->next(item);
	}

	/**
	 * Clear the tree and remove all items.
	 */
	void clear() {
		if(m_root)
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
		return m_root->knn(pt, n, iter);
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
