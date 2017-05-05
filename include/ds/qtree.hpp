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
			void search(double x, double y, double radius, U output) {
				// Search with an inside radius of zero.
				search(x, y, radius, 0, output);
			}

			// Search for points within [outside] of the coordinate, but not within [inside].
			template <class U>
			void search(double x, double y, double outside, double inside, U output) {
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
						}
					}
				}
			}

			// Search for points inside the bounding box.
			template <class U>
			void search(const Bounds& bounds, U output) {
				if(m_bounds.intersects(bounds)) {
					std::list<T> result;
					// Find all nodes that intersect the bounds and return their contents.
					findIntersecting(bounds, result);
					// Find all items contained in the bounds.
					for(T& item : result) {
						if(bounds.contains(item.x, item.y)) {
							*output = item;
							++output;
						}
					}
				}

			}

			// Search for points inside the GEOS geometry.
			template <class U>
			void search(const Geometry& geom, U output) {
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
						}
					}
				}
			}

			// Find the [n] nearest items to the coordinate.
			// If [outside] and [inside] are given, they're the outer and inner radii to search within.
			// and will be doubled on each iteration.
			template <class U>
			void nearest(double x, double y, uint64_t n, U output, double outside = 1, double inside = 0) {

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
							if(++i == n)
								break;
						}
						break;
					}
					inside = outside;
					outside *= 2; // If we don't find enough, try again with a larger radius.
				}

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

		template <class T>
		class MQTree : public Writer {
		private:
			// The geographic boundary corresponding to this tree. Will be reshaped to a square
			// using the longer side.
			Bounds m_bounds;
			// The number of bins along each side.
			uint64_t m_bins;
			// The numer of cached points before writing.
			uint64_t m_bufSize;
			// The current write position in the memory.
			uint64_t m_position;
			// The length of one side in map units.
			double m_size;
			// The number of elements contained in the tree.
			uint64_t m_count;
			// The number of writer threads to start.
			size_t m_cores;
			// A catalog to keep pointers into the mapped memory.
			std::unordered_map<uint64_t, CatItem<T> > m_catalog;
			// Cache points for writing.
			std::unordered_map<uint64_t, std::vector<T> > m_cache;
			// The mapped memory.
			std::unique_ptr<MappedFile> m_mapped;
			std::queue<T> m_wq;
			std::vector<std::thread> m_writers;
			bool m_running;
			std::mutex m_mtx;
			std::condition_variable m_cdn;

			int toBinX(double x) {
				return (int) (x - m_bounds.minx()) / m_size * m_bins;
			}

			int toBinY(double y) {
				return (int) (y - m_bounds.miny()) / m_size * m_bins;
			}

			int index(double x, double y) {
				return toBinY(y) * m_bins + toBinX(x);
			}

			std::list<int> indices(double x, double y, double radius, bool outer = false) {
				int b0 = g_max(0, toBinX(x - radius));
				int b2 = g_max(1, g_min((int) m_bins - 1, toBinX(x + radius)));
				int b1 = g_max(0, toBinY(y - radius));
				int b3 = g_max(1, g_min((int) m_bins - 1, toBinY(y + radius)));
				std::list<int> lst;
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
					for(int b = b1; b <= b3; ++b) {
						for(int a = b0; a <= b2; ++a)
							lst.push_back(b * m_bins + a);
					}
				}
				return lst;
			}

			void flushAll() {
				//std::cerr << "flush all" << "\n";
				for(const auto& it : m_cache) {
					uint64_t idx = it.first;
					const std::vector<T>& items = it.second;
					if(!items.empty()) {
						uint64_t pos = m_catalog[idx].update(items.size(), m_position);
						m_mapped->write((void*) items.data(), pos, items.size() * sizeof(T));
					}
				}
				m_cache.clear();
			}

			void flush(uint64_t idx) {
				//std::cerr << "flush " << idx << "\n";
				std::vector<T>& items = m_cache[idx];
				if(!items.empty()) {
					uint64_t pos = m_catalog[idx].update(items.size(), m_position);
					m_mapped->write((void*) items.data(), pos, items.size() * sizeof(T));
				}
				m_cache.erase(idx);
			}

			void startWrite() {
				if(!m_running) {
					m_running = true;
					for(int i = 0; i < m_cores; ++i)
						m_writers.push_back(std::thread(writer, this));
				}
			}

			void finishWrite() {
				if(m_running) {
					m_running = false;
					m_cdn.notify_all();
					for(std::thread& w : m_writers)
						w.join();
				}
			}

			void doWrite() {
				T item;
				while(m_running || !m_wq.empty()) {
					{
						std::unique_lock<std::mutex> lock(m_mtx);
						while(m_running && m_wq.empty())
							m_cdn.wait(lock);
						if(m_wq.empty())
							continue;
						item = m_wq.front();
						m_wq.pop();
					}
					m_cdn.notify_one();
					int idx = index(item.x, item.y);
					if(idx < 0 || idx >= (int) g_sq(m_bins))
						g_runerr("Illegal index in doWrite: " << idx << " max: " << g_sq(m_bins));
					m_cache[idx].push_back(std::move(item));
					if(m_cache[idx].size() == m_bufSize)
						flush(idx);
				}
				flushAll();
			}

			struct _sorter {
				double x, y;
				_sorter(double x, double y) : x(x), y(y) {}
				bool operator()(const T& a, const T& b) {
					return  (g_sq(x - a.x) + g_sq(y - a.y)) < (g_sq(x - b.x) + g_sq(y - b.y));
				}
			};


		public:

			MQTree(const Bounds& bounds, uint64_t bins = 512, uint64_t bufSize = 64, bool mapped = false) :
				m_bounds(bounds),
				m_bins(bins),
				m_bufSize(bufSize),
				m_position(0),
				m_count(0),
				m_running(false),
				m_cores(4) {

				m_size = g_max(m_bounds.width(), m_bounds.height());
				m_bounds.set(m_bounds.midx() - m_size / 2, m_bounds.midy() - m_size / 2,
					m_bounds.midx() + m_size / 2, m_bounds.midy() + m_size / 2);

				m_mapped.reset(new MappedFile(false));

				startWrite();
			}

			void finalize() {
				finishWrite();
			}

			uint64_t count() const {
				return m_count;
			}

			const Bounds& bounds() const {
				return m_bounds;
			}

			// Add an item to the store.
			void addItem(const T& item) {
				if(!m_bounds.contains(item.x, item.y))
					g_runerr("Item is out of bounds for MQTree.");
				{
					std::lock_guard<std::mutex> lock(m_mtx);
					m_wq.push(item);
				}
				++m_count;
				m_cdn.notify_one();
			}

			// Search for points within [radius] of the coordinate.
			template <class U>
			void search(double x, double y, double radius, U output) {
				if(!m_bounds.contains(x, y))
					return;
				finishWrite();
				std::list<int> idx = indices(x, y, radius);
				radius = g_sq(radius);
				for(const int& i : idx) {
					if(i < 0 || i >= (int) g_sq(m_bins))
						g_runerr("Illegal index in search (1): " << i << " max: " << g_sq(m_bins));
					if(m_catalog.find(i) != m_catalog.end()) {
						CatItem<T>& ci = m_catalog[i];
						ci.reset();
						uint64_t pos = 0;
						T item;
						while(ci.next(&pos)) {
							if(!m_mapped->read(item, pos))
								g_runerr("Failed to read item.");
							if(g_sq(item.x - x) + g_sq(item.y - y) <= radius) {
								item.pos = pos;
								*output = item;
								++output;
							}
						}
					}
				}
			}

			// Search for points inside the bounding box.
			template <class U>
			void search(const Bounds& bounds, U output) {
				if(!m_bounds.intersects(bounds))
					return;
				finishWrite();
				std::list<int> idx = indices(bounds.midx(), bounds.midy(), g_max(bounds.width(), bounds.height()));
				for(const int& i : idx) {
					if(i < 0 || (uint64_t) i >= g_sq(m_bins))
						g_runerr("Illegal index in search (2): " << i << " max: " << g_sq(m_bins));
					if(m_catalog.find(i) != m_catalog.end()) {
						CatItem<T>& ci = m_catalog[i];
						ci.reset();
						uint64_t pos = 0;
						T item;
						while(ci.next(&pos)) {
							if(!m_mapped->read(item, pos))
								g_runerr("Failed to read item.");
							if(bounds.contains(item.x, item.y)) {
								item.pos = pos;
								*output = item;
								++output;
							}
						}
					}
				}
			}

			// Search for points inside the GEOS geometry.
			template <class U>
			void search(const Geometry& geom, U output) {
				Envelope* env = (Envelope*) geom.getEnvelope();
				Bounds bounds(env->getMinX(), env->getMinY(), env->getMaxX(), env->getMaxY());
				if(!m_bounds.intersects(bounds))
					return;
				finishWrite();
				GeometryFactory gf;
				std::list<int> idx = indices(bounds.midx(), bounds.midy(), g_max(bounds.width(), bounds.height()));
				for(const int& i : idx) {
					if(i < 0 || i >= g_sq(m_bins))
						g_runerr("Illegal index in search (2): " << i << " max: " << g_sq(m_bins));
					if(m_catalog.find(i) != m_catalog.end()) {
						CatItem<T>& ci = m_catalog[i];
						ci.reset();
						uint64_t pos = 0;
						T item;
						while(ci.next(&pos)) {
							if(!m_mapped->read(item, pos))
								g_runerr("Failed to read item.");
							geos::geom::Point* pt = gf.createPoint(Coordinate(item.x, item.y));
							if(geom.contains(pt)) {
								item.pos = pos;
								*output = item;
								++output;
							}
							delete pt;
						}
					}
				}
			}

			template <class U>
			void nearest(double x, double y, uint64_t n, U output) {
				if(!m_bounds.contains(x, y))
					return;
				finishWrite();
				std::vector<T> tmp;
				uint64_t radius = 1;
				while(radius <= m_bins) {
					std::list<int> idx = indices(x, y, radius, true); // Only the outer ring; inner ring already searched.
					for(const int& i : idx) {
						if(i < 0 || i >=(int)  g_sq(m_bins))
							g_runerr("Illegal index in nearest: " << i << " max: " << g_sq(m_bins));
						if(m_catalog.find(i) != m_catalog.end()) {
							CatItem<T>& ci = m_catalog[i];
							ci.reset();
							uint64_t pos = 0;
							T item;
							while(ci.next(&pos)) {
								if(!m_mapped->read(item, pos))
									g_runerr("Failed to read item.");
								item.pos = pos;
								tmp.push_back(std::move(item));
							}
						}
					}

					if(tmp.size() >= n) {
						if(tmp.size() > n) {
							struct _sorter sorter(x, y);
							std::sort(tmp.begin(), tmp.end(), sorter);
						}
						for(uint64_t i = 0; i < n; ++i) {
							*output = std::move(tmp[0]);
							++output;
						}
						break;
					}

					++radius;
				}
			}

			std::vector<uint64_t> m_iter;
			uint64_t m_iterIdx;

			void reset() {
				for(const auto& it : m_catalog) {
					m_catalog[it.first].reset();
					m_iter.push_back(it.first);
				}
				m_iterIdx = 0;
			}

			bool next(T& item) {
				while(m_iterIdx < m_iter.size()) {
					uint64_t pos = 0;
					if(!m_catalog[m_iter[m_iterIdx]].next(&pos)) {
						++m_iterIdx;
						continue;
					}
					m_mapped->read(item, pos);
					item.pos = pos;
					return true;
				}
				return false;
			}

			void updateItem(const T& item) {
				m_mapped->write(item, item.pos);
			}

			~MQTree() {
			}

		};
	}
}

#endif
