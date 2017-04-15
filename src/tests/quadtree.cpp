/*
 * quadtree.cpp
 *
 *  Created on: Apr 14, 2017
 *      Author: rob
 */

#include "ds/quadtree.hpp"
#include "util.hpp"

	int main(int argc, char** argv) {

	using namespace geo::util;
	using namespace geo::ds;

	Bounds b(0, 0, 100, 100);

	FileQuadTree<char> qt(b, 1, 100);

	for(int i = 0; i < 1000; ++i) {
		double x = Util::random(0.0, 100.0);
		double y = Util::random(0.0, 100.0);
		char a = (char) Util::random(0, 65);
		std::cerr << x << ", " << y << "\n";
		qt.add(x, y, a);
	}

	Bounds q(4, 4, 6, 6);
	std::list<std::tuple<double, double, char> > lst;
	qt.search(q, std::back_inserter(lst));
	std::cerr << "search\n";
	for(const auto& it : lst)
		std::cerr << std::get<0>(it) << ", " << std::get<1>(it) << ", " << std::get<2>(it) << "\n";
}


