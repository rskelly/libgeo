#ifndef __QTREE_HPP__
#define __QTREE_HPP__

#include "util.hpp"

#include "geos/geom/Geometry.h"
#include "geos/geom/Envelope.h"
#include "geos/geom/Point.h"
#include "geos/geom/GeometryFactory.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <thread>
#include <condition_variable>

using namespace geo::util;
using namespace geos::geom;

namespace geo {
	namespace ds {
		namespace util {

			// Iterface for MQTree that allows writing
			// to disk.
			class Writer {
			public:
				virtual void doWrite() = 0;
				virtual ~Writer() {}
			};

			static void writer(Writer* qt) {
				qt->doWrite();
			}

			// A catalog item for MQTRee. Stores the number of
			// items in a disk location, and its offset from the
			// beginning of the file.
			template <class T>
			class CatItem {
			private:
				uint64_t m_idx;
				uint64_t m_iidx;

			public:
				std::vector<std::pair<uint64_t, uint64_t> > items;

				CatItem() : m_idx(0), m_iidx(0) {}

				uint64_t update(uint64_t count, uint64_t& position) {
					uint64_t pos = position;
					position += count * sizeof(T);
					items.push_back(std::make_pair(count, pos));
					return pos;
				}

				// Reset the position iterator.
				void reset() {
					m_idx = 0;
					m_iidx = 0;
				}

				// Get the position in memory of the next item in the catalog.
				bool next(uint64_t* pos) {
					while(m_idx < items.size()) {
						if(m_iidx < items[m_idx].first) {
							// The start position is stored; add to it the position
							// of the item at the given index.
							*pos = items[m_idx].second + m_iidx * (uint64_t) sizeof(T);
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
		}

		using namespace geo::ds::util;

		// An in-memory quadtree.
		template <class T>
		class QTree {
		protected:
			// The geographic boundary corresponding to this tree. Will be reshaped to a square
			// using the longer side.
			Bounds m_bounds;
			// The four sub-quads
			std::unique_ptr<QTree> m_nodes[4];
			// The list of items stored in the current node if not split.
			std::list<T> m_items;
			// The max tree depth. TODO: Should be small enough to prevent degenerate nodes.
			uint32_t m_maxDepth;
			// Max number of items before splitting a node.
			uint32_t m_maxCount;
			// The depth of the current node.
			uint32_t m_depth;
			// True if the current node has been split.
			bool m_split;

			// An iterator for traversing the current node's items.
			typename std::list<T>::iterator m_iter;
			// An index for traversing the current node's sub-nodes.
			uint8_t m_iterIdx;

			// Sorts the items.
			struct qtree_sorter {
				double x, y;
				qtree_sorter(double x, double y) : x(x), y(y) {}
				bool operator()(const T& a, const T& b) {
					return  (g_sq(x - a.x) + g_sq(y - a.y)) < (g_sq(x - b.x) + g_sq(y - b.y));
				}
			};

			// Returns the index of a sub-node given a point.
			int index(double x, double y) {
				return ((y >= m_bounds.midy()) << 1) | (x >= m_bounds.midx());
			}

			// Returns the node for the given index; creates it if necessary.
			QTree* node(uint8_t idx) {
				if(!m_nodes[idx].get()) {
					Bounds bounds;
					if(idx & 1) {
						bounds.minx(m_bounds.midx());
						bounds.maxx(m_bounds.maxx());
					} else {
						bounds.minx(m_bounds.minx());
						bounds.maxx(m_bounds.midx());
					}
					if(idx & 2) {
						bounds.miny(m_bounds.midy());
						bounds.maxy(m_bounds.maxy());
					} else {
						bounds.miny(m_bounds.miny());
						bounds.maxy(m_bounds.midy());
					}
					m_nodes[idx].reset(new QTree(bounds, m_maxDepth, m_maxCount, m_depth + 1));
				}
				return m_nodes[idx].get();
			}

			// Returns the node for the given point; creates it if necessary.
			QTree* node(double x, double y) {
				return node(index(x, y));
			}

			// Split the node; distribute items to subnodes.
			void split() {
				for(T& item : m_items)
					node(item.x, item.y)->addItem(std::move(item));
				m_items.clear();
				m_split = true;
			}

			// Search for all points inside nodes intersecting the bounding box.
			void findIntersecting(const Bounds& bounds, std::list<T>& output) {
				if(m_bounds.intersects(bounds)) {
					if(!m_split) {
						output.insert(output.end(), m_items.begin(), m_items.end());
					} else {
						for(int i = 0; i < 4; ++i) {
							if(m_nodes[i].get())
								m_nodes[i]->findIntersecting(bounds, output);
						}
					}
				}
			}

			// Construct a sub-node.
			QTree(const Bounds& bounds, int maxDepth, int maxCount, int depth) :
				m_bounds(bounds),
				m_maxDepth(maxDepth),
				m_maxCount(maxCount),
				m_depth(depth),
				m_split(false),
				m_iterIdx(0) {

			}

		public:

			// Construct a QTree with the given bounds, depth and count.
			QTree(const Bounds& bounds, int maxDepth, int maxCount) :
				m_bounds(bounds),
				m_maxDepth(maxDepth),
				m_maxCount(maxCount),
				m_depth(0),
				m_split(false),
				m_iterIdx(0) {

				m_bounds.cube();
			}

			// The total number of items in the node.
			uint64_t count() const {
				if(m_split) {
					return m_items.size();
				} else {
					uint64_t c = 0;
					for(int i = 0; i < 4; ++i) {
						if(m_nodes[i].get())
							c += m_nodes[i]->count();
					}
					return c;
				}
			}

			// The bounds of the node.
			const Bounds& bounds() const {
				return m_bounds;
			}

			// Add an item to the node.
			void addItem(const T& item) {
				if(!m_bounds.contains(item.x, item.y))
					g_runerr("Item is out of bounds for QTree.");
				if(m_depth < m_maxDepth && m_items.size() == m_maxCount)
					split();
				if(m_split) {
					node(item.x, item.y)->addItem(item);
				} else {
					m_items.push_back(item);
				}
			}

			void removeItem(const T& uitem) {
				if(m_bounds.contains(uitem.x, uitem.y)) {
					if(!m_split) {
						for(auto it = m_items.begin(); it != m_items.end(); ) {
							if(it->x == uitem.x && it->y == uitem.y) {
								it = m_items.erase(it);
							} else {
								++it;
							}
						}
					} else {
						node(uitem.x, uitem.y)->removeItem(uitem);
					}
				}
			}

			// Update the item in the tree. Updating the position is not allowed.
			// For that, remove the item and re-insert it.
			void updateItem(const T& uitem) {
				if(m_bounds.contains(uitem.x, uitem.y)) {
					if(!m_split) {
						for(T& item : m_items) {
							if(item.x == uitem.x && item.y == uitem.y)
								item = uitem;
						}
					} else {
						node(uitem.x, uitem.y)->updateItem(uitem);
					}
				}
			}

			// Search for points within [radius] of the coordinate.
			template <class U>
			int search(double x, double y, double radius, U output) {
				// Search with an inside radius of zero.
				return search(x, y, radius, 0, output);
			}

			// Search for points within [outside] of the coordinate, but not within [inside].
			template <class U>
			int search(double x, double y, double outside, double inside, U output) {
				int count = 0;
				// Search within the bounding box first.
				Bounds bounds(x - outside, y - outside, x + outside, y + outside);
				if(m_bounds.intersects(bounds)) {
					std::list<T> result;
					// Find all nodes that intersect the bounds and return their contents.
					findIntersecting(bounds, result);
					// Find all items inside the outside radius, and outside the inside radius.
					outside = g_sq(outside);
					inside = g_sq(inside);
					for(T& item : result) {
						double dist = g_sq(item.x - x) + g_sq(item.y - y);
						if(dist <= outside && dist > inside) {
							*output = item;
							++output;
							++count;
						}
					}
				}
				return count;
			}

			// Search for points inside the bounding box.
			template <class U>
			int search(const Bounds& bounds, U output) {
				int count = 0;
				if(m_bounds.intersects(bounds)) {
					std::list<T> result;
					// Find all nodes that intersect the bounds and return their contents.
					findIntersecting(bounds, result);
					// Find all items contained in the bounds.
					for(T& item : result) {
						if(bounds.contains(item.x, item.y)) {
							*output = item;
							++output;
							++count;
						}
					}
				}
				return count;
			}

			// Search for points inside the GEOS geometry.
			template <class U>
			int search(const Geometry& geom, U output) {
				int count = 0;
				// Create an envelope from the geometry.
				Envelope* env = (Envelope*) geom.getEnvelope();
				// Make sure the geom intersects with the tree.
				Bounds bounds(env->getMinX(), env->getMinY(), env->getMaxX(), env->getMaxY());
				if(m_bounds.intersects(bounds)) {
					std::list<T> result;
					// Find all nodes with relevant items.
					findIntersecting(bounds, result);
					GeometryFactory gf;
					// Find all items that are inside the geometry.
					for(T& item : result) {
						geos::geom::Point* pt = gf.createPoint(Coordinate(item.x, item.y));
						if(geom.contains(pt)) {
							*output = item;
							++output;
							++count;
						}
					}
				}
				return count;
			}

			// Find the [n] nearest items to the coordinate.
			// If [outside] and [inside] are given, they're the outer and inner radii to search within.
			// and will be doubled on each iteration.
			template <class U>
			int nearest(double x, double y, uint64_t n, U output, double outside = 1, double inside = 0) {
				int count = 0;
				std::list<T> result;
				// Search while the outside radius is smaller than the bounds.
				while(outside <= m_bounds.width()) {
					search(x, y, outside, inside, std::back_inserter(result));
					if(result.size() >= n) {
						if(result.size() > n) {
							struct qtree_sorter sorter(x, y);
							result.sort(sorter);
							//std::sort(result.begin(), result.end(), sorter);
						}
						uint64_t i = 0;
						for(T& item : result) {
							*output = item;
							++output;
							++count;
							if(++i == n)
								break;
						}
						break;
					}
					inside = outside;
					outside *= 2; // If we don't find enough, try again with a larger radius.
				}
				return count;
			}

			// Reset the iterator on this node and children.
			void reset() {
				if(!m_split) {
					m_iter = m_items.begin();
				} else {
					for(int i = 3; i >= 0; --i) {
						if(m_nodes[i].get()) {
							m_nodes[i]->reset();
							m_iterIdx = i;
						}
					}
				}
			}

			// Get the next item in this subtree.
			bool next(T& item) {
				if(m_split) {
					while(m_iterIdx < 4) {
						if(m_nodes[m_iterIdx].get()) {
							if(m_nodes[m_iterIdx]->next(item))
								return true;
						}
						++m_iterIdx;
					}
				} else if(m_iter != m_items.end()) {
					item = *m_iter;
					++m_iter;
					return true;
				}
				return false;
			}

			~QTree() {
			}

		};

		// An file-backed quadtree.
		template <class T>
		class FQTree {
		protected:
						// The geographic boundary corresponding to this tree. Will be reshaped to a square
			// using the longer side.
			Bounds m_bounds;
			// The four sub-quads
			std::unique_ptr<FQTree> m_nodes[4];
			// The list of items stored in the current node if not split.
			std::list<T> m_items;
			// Items used by the iterator.
			std::list<T> m_iterItems;
			// The max tree depth. TODO: Should be small enough to prevent degenerate nodes.
			uint32_t m_maxDepth;
			// Max number of items before splitting a node.
			uint32_t m_maxCount;
			// The depth of the current node.
			uint32_t m_depth;
			// The number of items in the node.
			uint32_t m_count;
			// True if the current node has been split.
			bool m_split;

			// An iterator for traversing the current node's items.
			typename std::list<T>::iterator m_iter;
			// An index for traversing the current node's sub-nodes.
			uint8_t m_iterIdx;

			std::string m_fdir;
			std::string m_fpath;
			size_t m_fpos;

			std::mutex m_fmtx;

			// Sorts the items.
			struct qtree_sorter {
				double x, y;
				qtree_sorter(double x, double y) : x(x), y(y) {}
				bool operator()(const T& a, const T& b) {
					return  (g_sq(x - a.x) + g_sq(y - a.y)) < (g_sq(x - b.x) + g_sq(y - b.y));
				}
			};

			// Returns the index of a sub-node given a point.
			int index(double x, double y) {
				return ((y >= m_bounds.midy()) << 1) | (x >= m_bounds.midx());
			}

			// Returns the node for the given index; creates it if necessary.
			QTree* node(uint8_t idx) {
				if(!m_nodes[idx].get()) {
					Bounds bounds;
					if(idx & 1) {
						bounds.minx(m_bounds.midx());
						bounds.maxx(m_bounds.maxx());
					} else {
						bounds.minx(m_bounds.minx());
						bounds.maxx(m_bounds.midx());
					}
					if(idx & 2) {
						bounds.miny(m_bounds.midy());
						bounds.maxy(m_bounds.maxy());
					} else {
						bounds.miny(m_bounds.miny());
						bounds.maxy(m_bounds.midy());
					}
					std::string dir = Util::pathJoin(m_fdir, std::to_string(idx));
					m_nodes[idx].reset(new FQTree(dir, bounds, m_maxDepth, m_maxCount, m_depth + 1));
				}
				return m_nodes[idx].get();
			}

			// Returns the node for the given point; creates it if necessary.
			QTree* node(double x, double y) {
				return node(index(x, y));
			}

			void load(std::vector<T>& items) {
				std::FILE* f = std::fopen(m_fpath.c_str(), "rb");
				if(!f)
					g_runerr("Failed to open file for split.");
				std::fread(items.data(), sizeof(T), items.size(), f);
			}

			// Split the node; distribute items to subnodes.
			void split() {
				flush();
				std::lock_guard<std::mutex> lk(m_fmtx);
				std::vector<T> items(m_fpos / sizeof(T));
				for(T& item : items)
					node(item.x, item.y)->addItem(std::move(item));
				m_split = true;
				m_count = 0;
				Util::rm(m_fpath);
			}

			// Search for all points inside nodes intersecting the bounding box.
			void findIntersecting(const Bounds& bounds, std::list<T>& output) {
				if(m_bounds.intersects(bounds)) {
					flush();
					if(!m_split) {
						std::vector<T> items;
						load(items);
						output.insert(output.end(), items.begin(), items.end());
					} else {
						for(int i = 0; i < 4; ++i) {
							if(m_nodes[i].get())
								m_nodes[i]->findIntersecting(bounds, output);
						}
					}
				}
			}

			// Construct a sub-node.
			FQTree(const std::string& dir, const Bounds& bounds, int maxDepth, int maxCount, int depth) :
				m_bounds(bounds),
				m_maxDepth(maxDepth),
				m_maxCount(maxCount),
				m_depth(depth),
				m_count(0),
				m_split(false),
				m_iterIdx(0),
				m_fdir(dir),
				m_fpos(0) {

				Util::mkdir(dir);
				m_fpath = Util::pathJoin(dir, "items");
			}

		public:

			// Construct a QTree with the given bounds, depth and count.
			FQTree(const std::string& dir, const Bounds& bounds, int maxDepth, int maxCount) :
				m_bounds(bounds),
				m_maxDepth(maxDepth),
				m_maxCount(maxCount),
				m_depth(0),
				m_count(0),
				m_split(false),
				m_iterIdx(0),
				m_fdir(dir),
				m_fpos(0) {

				m_bounds.cube();
				Util::mkdir(dir);
				m_fpath = Util::pathJoin(dir, "items");
			}

			// The total number of items in the node.
			uint64_t count() const {
				if(m_split) {
					return m_count;
				} else {
					uint64_t c = 0;
					for(int i = 0; i < 4; ++i) {
						if(m_nodes[i].get())
							c += m_nodes[i]->count();
					}
					return c;
				}
			}

			// The bounds of the node.
			const Bounds& bounds() const {
				return m_bounds;
			}

			// Add an item to the node.
			void addItem(const T& item) {
				if(!m_bounds.contains(item.x, item.y))
					g_runerr("Item is out of bounds for FQTree.");
				if(m_depth < m_maxDepth && m_count == m_maxCount)
					split();
				if(m_split) {
					node(item.x, item.y)->addItem(item);
				} else {
					m_items.push_back(item);
				}
			}

			void removeItem(const T& uitem) {
				if(m_bounds.contains(uitem.x, uitem.y)) {
					flush();
					if(!m_split) {
						load(m_items);
						for(auto it = m_items.begin(); it != m_items.end(); ) {
							if(it->x == uitem.x && it->y == uitem.y) {
								it = m_items.erase(it);
							} else {
								++it;
							}
						}
						flush();
					} else {
						node(uitem.x, uitem.y)->removeItem(uitem);
					}
				}
			}

			// Update the item in the tree. Updating the position is not allowed.
			// For that, remove the item and re-insert it.
			void updateItem(const T& uitem) {
				flush();
				if(m_bounds.contains(uitem.x, uitem.y)) {
					if(!m_split) {
						load(m_items);
						for(T& item : m_items) {
							if(item.x == uitem.x && item.y == uitem.y)
								item = uitem;
						}
						flush();
					} else {
						node(uitem.x, uitem.y)->updateItem(uitem);
					}
				}
			}

			// Search for points within [radius] of the coordinate.
			template <class U>
			int search(double x, double y, double radius, U output) {
				// Search with an inside radius of zero.
				return search(x, y, radius, 0, output);
			}

			// Search for points within [outside] of the coordinate, but not within [inside].
			template <class U>
			int search(double x, double y, double outside, double inside, U output) {
				int count = 0;
				// Search within the bounding box first.
				Bounds bounds(x - outside, y - outside, x + outside, y + outside);
				if(m_bounds.intersects(bounds)) {
					std::list<T> result;
					// Find all nodes that intersect the bounds and return their contents.
					findIntersecting(bounds, result);
					// Find all items inside the outside radius, and outside the inside radius.
					outside = g_sq(outside);
					inside = g_sq(inside);
					for(T& item : result) {
						double dist = g_sq(item.x - x) + g_sq(item.y - y);
						if(dist <= outside && dist > inside) {
							*output = item;
							++output;
							++count;
						}
					}
				}
				return count;
			}

			// Search for points inside the bounding box.
			template <class U>
			int search(const Bounds& bounds, U output) {
				int count = 0;
				if(m_bounds.intersects(bounds)) {
					std::list<T> result;
					// Find all nodes that intersect the bounds and return their contents.
					findIntersecting(bounds, result);
					// Find all items contained in the bounds.
					for(T& item : result) {
						if(bounds.contains(item.x, item.y)) {
							*output = item;
							++output;
							++count;
						}
					}
				}
				return count;
			}

			// Search for points inside the GEOS geometry.
			template <class U>
			int search(const Geometry& geom, U output) {
				int count = 0;
				// Create an envelope from the geometry.
				Envelope* env = (Envelope*) geom.getEnvelope();
				// Make sure the geom intersects with the tree.
				Bounds bounds(env->getMinX(), env->getMinY(), env->getMaxX(), env->getMaxY());
				if(m_bounds.intersects(bounds)) {
					std::list<T> result;
					// Find all nodes with relevant items.
					findIntersecting(bounds, result);
					GeometryFactory gf;
					// Find all items that are inside the geometry.
					for(T& item : result) {
						geos::geom::Point* pt = gf.createPoint(Coordinate(item.x, item.y));
						if(geom.contains(pt)) {
							*output = item;
							++output;
							++count;
						}
					}
				}
				return count;
			}

			// Find the [n] nearest items to the coordinate.
			// If [outside] and [inside] are given, they're the outer and inner radii to search within.
			// and will be doubled on each iteration.
			template <class U>
			int nearest(double x, double y, uint64_t n, U output, double outside = 1, double inside = 0) {
				int count = 0;
				std::list<T> result;
				// Search while the outside radius is smaller than the bounds.
				while(outside <= m_bounds.width()) {
					search(x, y, outside, inside, std::back_inserter(result));
					if(result.size() >= n) {
						if(result.size() > n) {
							struct qtree_sorter sorter(x, y);
							result.sort(sorter);
						}
						uint64_t i = 0;
						for(T& item : result) {
							*output = item;
							++output;
							++count;
							if(++i == n)
								break;
						}
						break;
					}
					inside = outside;
					outside *= 2; // If we don't find enough, try again with a larger radius.
				}
				return count;
			}

			// Reset the iterator on this node and children.
			void reset() {
				flush();
				if(!m_split) {
					m_iterItems.clear();
					load(m_iterItems);
					m_iter = m_iterItems.begin();
				} else {
					for(int i = 3; i >= 0; --i) {
						if(m_nodes[i].get()) {
							m_nodes[i]->reset();
							m_iterIdx = i;
						}
					}
				}
			}

			// Get the next item in this subtree.
			bool next(T& item) {
				if(m_split) {
					while(m_iterIdx < 4) {
						if(m_nodes[m_iterIdx].get()) {
							if(m_nodes[m_iterIdx]->next(item))
								return true;
						}
						++m_iterIdx;
					}
				} else if(m_iter != m_iterItems.end()) {
					item = *m_iter;
					++m_iter;
					return true;
				}
				m_iterItems.clear();
				return false;
			}

			void flush() {
				if(!m_split) {
					if(!m_items.empty()) {
						std::lock_guard<std::mutex> lk(m_fmtx);
						std::FILE* f = std::fopen(m_fpath.c_str(), "rb+");
						if(!f) {
							f = std::fopen(m_fpath.c_str(), "wb");
							if(!f)
								g_runerr("Can't open file: " << m_fpath);
						}
						std::fseek(f, m_fpos, SEEK_SET);
						int written = 0;
						if(!(written = std::fwrite(m_items, sizeof(T), m_items.size(), f)))
							g_runerr("Failed to write items.");
						std::fclose(f);
						m_fpos += sizeof(T) * written;
						m_items.clear();
					}
				} else {
					for(auto& t : m_nodes)
						t->flush();
				}
			}

			~FQTree() {
			}

		};

	}
}
#endif
