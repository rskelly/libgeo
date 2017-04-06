/*
 * interval_tree.hpp
 *
 *  Created on: Apr 5, 2017
 *      Author: rob
 */

#ifndef _INTERVAL_TREE_HPP_
#define _INTERVAL_TREE_HPP_

template <class T>
class IntervalTree {
private:
	double m_min, m_max, m_mid;
	IntervalTree* m_left, m_right;
	std::vector<std::tuple<double, double, T> > m_items;

public:
	IntervalTree(double min, double max) :
		m_min(min), m_max(max), m_mid(min + (max - min) / 2.0),
		m_left(nullptr), m_right(nullptr) {}

	void add(double min, double max, T value) {
		if(max < m_mid) {
			if(!m_left)
				m_left = new IntervalTree(m_min, m_mid);
			m_left->add(min, max, value);
		} else if(min > m_mid) {
			if(!m_right)
				m_right = new IntervalTree(m_mid, m_max);
		} else {
			m_items.push_back(std::make_tuple(min, max, value));
		}
	}

	template <class U>
	void find(double x, U iter) {

	}

	template <class U>
	void find(double x0, double x1, U iter) {

	}

	~IntervalTree() {
		if(m_left) delete m_left;
		if(m_right) delete m_right;
	}
};



#endif /* _INTERVAL_TREE_HPP_ */
