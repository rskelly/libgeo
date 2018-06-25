#ifndef __QTREE_HPP__
#define __QTREE_HPP__
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <thread>
#include <condition_variable>
#include <fstream>

#include "geos/geom/Geometry.h"
#include "geos/geom/Envelope.h"
#include "geos/geom/Point.h"
#include "geos/geom/GeometryFactory.h"

#include "util.hpp"

using namespace geo::util;
using namespace geos::geom;

namespace geo {
	namespace ds {

		// Sorts the items.
		template <class T>
		struct qtree_sorter {
			double x, y;
			qtree_sorter(double x, double y) : 
				x(x), 
				y(y) {
			}

			bool operator()(const T& a, const T& b) {
				return  (g_sq(x - a.x()) + g_sq(y - a.y())) < (g_sq(x - b.x()) + g_sq(y - b.y()));
			}
		};

		// An in-memory quadtree.
		template <class T>
		class QTree {
		protected:
			Bounds m_bounds; 		///<! The geographic boundary corresponding to this tree. Will be reshaped to a square using the longer side.
			std::list<T> m_items;	///<! The list of items stored in the current node if not split.
			QTree* m_nodes[4];		///<! The four sub-quads
			int m_maxDepth;			///<! The max tree depth. TODO: Should be small enough to prevent degenerate nodes.
			int m_maxCount;			///<! Max number of items before splitting a node.
			int m_depth;			///<! The depth of the current node.
			int m_iterIdx; 			///<! An index for traversing the current node's sub-nodes.
			bool m_split;			///<! True if the current node has been split.

			typename std::list<T>::iterator m_iter; ///<! An iterator for traversing the current node's items.


			/**
			 * Returns the index of a sub-node given a point.
			 *
			 * @param x X-coordinate.
			 * @param y Y-coordinate.
			 */
			int index(double x, double y) {
				return ((y >= m_bounds.midy()) << 1) | (x >= m_bounds.midx());
			}

			/**
			 * Returns the node for the given index; creates it if necessary.
			 *
			 * @param idx The node index.
			 */
			QTree* node(int idx) {
				if(!m_nodes[idx]) {
					Bounds bounds;
					if(idx & 1) {
						bounds.minx(m_bounds.midx());
					} else {
						bounds.maxx(m_bounds.midx());
					}
					if(idx & 2) {
						bounds.miny(m_bounds.midy());
					} else {
						bounds.maxy(m_bounds.midy());
					}
					m_nodes[idx] = new QTree(bounds, m_maxDepth, m_maxCount, m_depth + 1);
				}
				return m_nodes[idx];
			}

			/**
			 * Returns the node for the given point; creates it if necessary.
			 *
			 * @param x X-coordinate.
			 * @param y Y-coordinate.
			 */
			QTree* node(double x, double y) {
				return node(index(x, y));
			}

			/**
			 * Split the node; distribute items to subnodes.
			 */
			void split() {
				for(T& item : m_items)
					node(item.x(), item.y())->addItem(std::move(item));
				m_items.clear();
				m_split = true;
			}

			/**
			 * Search for all points inside nodes intersecting the bounding box.
			 *
			 * @param bounds The bounds to search in.
			 * @param output The output list.
			 */
			void findIntersecting(const Bounds& bounds, std::list<T>& output) {
				if(m_bounds.intersects(bounds)) {
					if(!m_split) {
						output.insert(output.end(), m_items.begin(), m_items.end());
					} else {
						for(int i = 0; i < 4; ++i) {
							if(m_nodes[i])
								m_nodes[i]->findIntersecting(bounds, output);
						}
					}
				}
			}

			/**
			 * Construct a sub-node.
			 *
			 * @param bounds The 3D bounding box of the tree.
			 * @param maxDepth The maximum depth of a leaf.
			 * @param maxCount The maximum number of items in a leaf.
			 * @param depth The depth of the new node.
			 */
			QTree(const Bounds& bounds, int maxDepth, int maxCount, int depth) :
				m_bounds(bounds),
				m_maxDepth(maxDepth),
				m_maxCount(maxCount),
				m_depth(depth),
				m_iterIdx(0),
				m_split(false) {
				for(int i = 0; i < 4; ++i)
					m_nodes[i] = nullptr;
			}

		public:

			/**
			 * Construct a QTree with the given bounds, depth and count.
			 *
			 * @param bounds The 3D bounding box of the tree.
			 * @param maxDepth The maximum depth of a leaf.
			 * @param maxCount The maximum number of items in a leaf.
			 */
			QTree(const Bounds& bounds, int maxDepth, int maxCount) :
				QTree(bounds, maxDepth, maxCount, 0) {
			}

			/**
			 * Initialize a QTree with the given bounds, depth and count.
			 *
			 * @param bounds The 3D bounding box of the tree.
			 * @param maxDepth The maximum depth of a leaf.
			 * @param maxCount The maximum number of items in a leaf.
			 */
			void init(const Bounds& bounds, int maxDepth, int maxCount) {
				for(int i = 0; i < 4; ++i) {
					if(m_nodes[i]) delete m_nodes[i];
					m_nodes[i] = nullptr;
				}
				m_bounds = bounds;
				m_maxDepth = maxDepth;
				m_maxCount = maxCount;
				m_depth = 0;
				m_split = false;
				m_iterIdx = 0;
				m_bounds.cube();
			}

			/**
			 * The total number of items in the node.
			 *
			 * @return The total number of items in the node.
			 */
			size_t count() const {
				if(!m_split) {
					return m_items.size();
				} else {
					size_t c = 0;
					for(int i = 0; i < 4; ++i) {
						if(m_nodes[i])
							c += m_nodes[i]->count();
					}
					return c;
				}
			}

			/**
			 * The bounds of the node.
			 *
			 * @return The bounds of the node.
			 */
			const Bounds& bounds() const {
				return m_bounds;
			}

			/**
			 * Add an item to the node.
			 *
			 * @param item An item.
			 */
			void addItem(const T& item) {
				if(!m_bounds.contains(item.x(), item.y())) {
					g_warn("Item is out of bounds for QTree.");
					return;
				}
				if(m_depth < m_maxDepth && m_items.size() == m_maxCount)
					split();
				if(m_split) {
					node(item.x(), item.y())->addItem(item);
				} else {
					m_items.push_back(item);
				}
			}

			/**
			 * Remove the item from the tree.
			 *
			 * @param uitem The item to remove.
			 */
			void removeItem(const T& uitem) {
				if(m_bounds.contains(uitem.x(), uitem.y())) {
					if(!m_split) {
						for(auto it = m_items.begin(); it != m_items.end(); ) {
							if(it->x == uitem.x() && it->y == uitem.y()) {
								it = m_items.erase(it);
							} else {
								++it;
							}
						}
					} else {
						node(uitem.x(), uitem.y())->removeItem(uitem);
					}
				}
			}

			/**
			 * Update the item in the tree. Updating the position is not allowed.
			 * For that, remove the item and re-insert it.
			 *
			 * @param uitem An item to update.
			 */
			void updateItem(const T& uitem) {
				if(m_bounds.contains(uitem.x(), uitem.y())) {
					if(!m_split) {
						for(T& item : m_items) {
							if(item.x() == uitem.x() && item.y() == uitem.y())
								item = uitem;
						}
					} else {
						node(uitem.x(), uitem.y())->updateItem(uitem);
					}
				}
			}

			/**
			 * Search for points within [radius] of the coordinate.
			 *
			 * @param x The x-coordinate.
			 * @param y The y-coordinate.
			 * @param radius The search radius.
			 * @param output An iterator for output items.
			 * @return The number of found items.
			 */
			template <class U>
			int search(double x, double y, double radius, U output) {
				// Search with an inside radius of zero.
				return search(x, y, radius, 0, output);
			}

			/**
			 * Search for points within [outside] of the coordinate, but not within [inside].
			 *
			 * @param x The x-coordinate.
			 * @param y The y-coordinate.
			 * @param inside The inside search radius.
			 * @param outside The outside search radius.
			 * @param output An iterator for output items.
			 * @return The number of found items.
			 */
			template <class U>
			size_t search(double x, double y, double outside, double inside, U output) {
				size_t count = 0;
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
						double dist = g_sq(item.x() - x) + g_sq(item.y() - y);
						if(dist <= outside && dist > inside) {
							*output = item;
							++output;
							++count;
						}
					}
				}
				return count;
			}

			/**
			 * Search for points inside the bounding box.
			 *
			 * @param bounds The search bounds.
			 * @param output An iterator for output items.
			 * @return The number of found items.
			 */
			template <class U>
			size_t search(const Bounds& bounds, U output) {
				size_t count = 0;
				if(m_bounds.intersects(bounds)) {
					std::list<T> result;
					// Find all nodes that intersect the bounds and return their contents.
					findIntersecting(bounds, result);
					// Find all items contained in the bounds.
					for(T& item : result) {
						if(bounds.contains(item.x(), item.y())) {
							*output = item;
							++output;
							++count;
						}
					}
				}
				return count;
			}

			/**
			 * Search for points inside the GEOS geometry.
			 *
			 * @param geom A Geometry.
			 * @param output An iterator for output items.
			 * @return The number of found items.
			 */
			template <class U>
			size_t search(const Geometry& geom, U output) {
				size_t count = 0;
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
						geos::geom::Point* pt = gf.createPoint(Coordinate(item.x(), item.y()));
						if(geom.contains(pt)) {
							*output = item;
							++output;
							++count;
						}
						delete pt;
					}
				}
				return count;
			}

			/**
			 * Find the [n] nearest items to the coordinate.
			 * If [outside] and [inside] are given, they're the outer and inner radii to search within.
			 * and will be doubled on each iteration.
			 *
			 * @param x The x-coordinate.
			 * @param y The y-coordinate.
			 * @param n The number of items to find.
			 * @param output An iterator for output items.
			 * @param outside The outside search radius.
			 * @param inside The inside search radius.
			 * @return The number of found items.
			 */
			template <class U>
			size_t nearest(double x, double y, size_t n, U output, double outside = 1, double inside = 0) {
				size_t count = 0;
				std::list<T> result;
				// Search while the outside radius is smaller than the bounds.
				while(outside <= m_bounds.width()) {
					search(x, y, outside, inside, std::back_inserter(result));
					if(result.size() >= n) {
						if(result.size() > n) {
							qtree_sorter<T> sorter(x, y);
							result.sort(sorter);
						}
						size_t i = 0;
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

			/**
			 * Reset the iterator on this node and children.
			 */
			void reset() {
				if(!m_split) {
					m_iter = m_items.begin();
				} else {
					for(int i = 3; i >= 0; --i) {
						if(m_nodes[i]) {
							m_nodes[i]->reset();
							m_iterIdx = i;
						}
					}
				}
			}

			/**
			 * Get the next item in this subtree.
			 *
			 * @param item An item to set (out).
			 * @return True if an item was found.
			 */
			bool next(T& item) {
				if(m_split) {
					while(m_iterIdx < 4) {
						if(m_nodes[m_iterIdx]) {
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
				for(int i = 0; i < 4; ++i)
					if(m_nodes[i]) delete m_nodes[i];
			}

		};
	}
}
#endif

