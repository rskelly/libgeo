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

			class Writer {
			public:
				virtual void doWrite() = 0;
				virtual ~Writer() {}
			};

			static void writer(Writer* qt) {
				qt->doWrite();
			}

			template <class T>
			class CatItem {
			private:
				size_t m_idx;
				size_t m_iidx;

			public:
				std::vector<std::pair<size_t, size_t> > items;

				CatItem() : m_idx(0), m_iidx(0) {}

				size_t update(size_t count, size_t& position) {
					size_t pos = position;
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
				bool next(size_t* pos) {
					while(m_idx < items.size()) {
						if(m_iidx < items[m_idx].first) {
							// The start position is stored; add to it the position
							// of the item at the given index.
							*pos = items[m_idx].second + m_iidx * sizeof(T);
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

		template <class T>
		class QTree : public Writer {
		private:
			// The geographic boundary corresponding to this tree. Will be reshaped to a square
			// using the longer side.
			Bounds m_bounds;
			// The number of bins along each side.
			size_t m_bins;
			// The numer of cached points before writing.
			size_t m_bufSize;
			// The current write position in the memory.
			size_t m_position;
			// The length of one side in map units.
			double m_size;
			// The number of elements contained in the tree.
			size_t m_count;
			// The number of writer threads to start.
			size_t m_cores;
			// A catalog to keep pointers into the mapped memory.
			std::unordered_map<size_t, CatItem<T> > m_catalog;
			// Cache points for writing.
			std::unordered_map<size_t, std::vector<T> > m_cache;
			// The mapped memory.
			MappedFile m_mapped;
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
					size_t idx = it.first;
					const std::vector<T>& items = it.second;
					if(!items.empty()) {
						size_t pos = m_catalog[idx].update(items.size(), m_position);
						m_mapped.write((void*) items.data(), pos, items.size() * sizeof(T));
					}
				}
				m_cache.clear();
			}

			void flush(size_t idx) {
				//std::cerr << "flush " << idx << "\n";
				std::vector<T>& items = m_cache[idx];
				if(!items.empty()) {
					size_t pos = m_catalog[idx].update(items.size(), m_position);
					m_mapped.write((void*) items.data(), pos, items.size() * sizeof(T));
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
					for(auto& writer : m_writers)
						writer.join();
				}
			}

			void doWrite() {
				T item;
				while(m_running) {
					{
						std::unique_lock<std::mutex> lock(m_mtx);
						while(m_running && m_wq.empty())
							m_cdn.wait(lock);
						item = std::move(m_wq.front());
						m_wq.pop();
						int idx = index(item.x, item.y);
						if(idx < 0 || idx >= (int) g_sq(m_bins))
							g_runerr("Illegal index in doWrite: " << idx << " max: " << g_sq(m_bins));
						m_cache[idx].push_back(std::move(item));
						if(m_cache[idx].size() == m_bufSize)
							flush(idx);
					}
					m_cdn.notify_one();
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

			QTree(const Bounds& bounds, size_t bins = 512, size_t bufSize = 64) :
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

				startWrite();
			}

			void finalize() {
				finishWrite();
			}

			size_t count() const {
				return m_count;
			}

			const Bounds& bounds() const {
				return m_bounds;
			}

			// Add an item to the store.
			void addItem(const T& item) {
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
						size_t pos = 0;
						T item;
						while(ci.next(&pos)) {
							if(!m_mapped.read(item, pos))
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
					if(i < 0 || (size_t) i >= g_sq(m_bins))
						g_runerr("Illegal index in search (2): " << i << " max: " << g_sq(m_bins));
					if(m_catalog.find(i) != m_catalog.end()) {
						CatItem<T>& ci = m_catalog[i];
						ci.reset();
						size_t pos = 0;
						T item;
						while(ci.next(&pos)) {
							if(!m_mapped.read(item, pos))
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
						size_t pos = 0;
						T item;
						while(ci.next(&pos)) {
							if(!m_mapped.read(item, pos))
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
			void nearest(double x, double y, size_t n, U output) {
				if(!m_bounds.contains(x, y))
					return;
				finishWrite();
				std::vector<T> tmp;
				size_t radius = 1;
				while(radius <= m_bins) {
					std::list<int> idx = indices(x, y, radius, true); // Only the outer ring; inner ring already searched.
					for(const int& i : idx) {
						if(i < 0 || i >=(int)  g_sq(m_bins))
							g_runerr("Illegal index in nearest: " << i << " max: " << g_sq(m_bins));
						if(m_catalog.find(i) != m_catalog.end()) {
							CatItem<T>& ci = m_catalog[i];
							ci.reset();
							size_t pos = 0;
							T item;
							while(ci.next(&pos)) {
								if(!m_mapped.read(item, pos))
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
						for(size_t i = 0; i < n; ++i) {
							*output = std::move(tmp[0]);
							++output;
						}
						break;
					}

					++radius;
				}
			}

			std::vector<size_t> m_iter;
			size_t m_iterIdx;

			void reset() {
				for(const auto& it : m_catalog) {
					m_catalog[it.first].reset();
					m_iter.push_back(it.first);
				}
				m_iterIdx = 0;
			}

			bool next(T& item) {
				while(m_iterIdx < m_iter.size()) {
					size_t pos = 0;
					if(!m_catalog[m_iter[m_iterIdx]].next(&pos)) {
						++m_iterIdx;
						continue;
					}
					m_mapped.read(item, pos);
					item.pos = pos;
					return true;
				}
				return false;
			}

			void updateItem(const T& item) {
				m_mapped.write(item, item.pos);
			}

			~QTree() {
			}

		};

	}
}

#endif
