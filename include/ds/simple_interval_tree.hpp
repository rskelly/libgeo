/*
 * range_tree.hpp
 *
 *  Created on: Feb 18, 2017
 *      Author: rob
 */

#ifndef _SIMPLE_INTERVAL_TREE_HPP_
#define _SIMPLE_INTERVAL_TREE_HPP_

#include <memory>
#include <list>
#include <cmath>

namespace geo {

	namespace ds {

		template <class T, class U>
		class SimpleIntervalTree {
		protected:
			T m_key;
			std::list<U> m_values;
			std::unique_ptr<SimpleIntervalTree> m_lt;
			std::unique_ptr<SimpleIntervalTree> m_gt;

			SimpleIntervalTree(T key, U value) :
					m_key(key) {
				m_values.push_back(value);
			}

		public:
			SimpleIntervalTree() :
					m_key(nan("")) {}

			// Add an item to the tree.
			void add(T key, U value) {
				if(key == m_key || std::isnan(m_key)) {
					m_key = key;
					m_values.push_back(value);
				} else if(key < m_key) {
					if(!m_lt.get()) {
						m_lt.reset(new SimpleIntervalTree(key, value));
					} else {
						m_lt->add(key, value);
					}
				} else {
					if(!m_gt.get()) {
						m_gt.reset(new SimpleIntervalTree(key, value));
					} else {
						m_gt->add(key, value);
					}
				}
			}

			// Find the items in the interval that contains the query.
			template <class I>
			bool find(T query, I iter) const {
				if(query == m_key) {
					for(const U& v : m_values) {
						*iter = v;
						++iter;
					}
					return !m_values.empty();
				} else if(query < m_key) {
					if(m_lt.get())
						return m_lt->find(query, iter);
				} else {
					if(m_gt.get()) {
						return m_gt->find(query, iter);
					} else {
						for(const U& v : m_values) {
							*iter = v;
							++iter;
						}
						return !m_values.empty();
					}
				}
				return false;
			}

			bool find(T query, U* value) const {
				if(query == m_key) {
					if(!m_values.empty()) {
						*value = *(m_values.begin());
						return true;
					}
				} else if(query < m_key) {
					if(m_lt.get())
						return m_lt->find(query, value);
				} else {
					if(m_gt.get() && m_gt->m_key <= query) {
						return m_gt->find(query, value);
					} else {
						*value = *(m_values.begin());
						return true;
					}
				}
				return false;
			}

		};


	} // ds

} // geo

#endif /* _SIMPLE_INTERVAL_TREE_HPP_ */
