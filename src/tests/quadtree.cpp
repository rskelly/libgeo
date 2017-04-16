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

	std::vector<uint64_t> ids;
	for(int i = 0; i < 10; ++i) {
		double x = Util::random(0.0, 10.0);
		double y = Util::random(0.0, 10.0);
		char a = (char) Util::random(65, 91);
		std::cerr << x << ", " << y << "\n";
		ids.push_back(qt.addItem(x, y, a));
	}

	std::cerr << "ids\n";
	std::tuple<uint64_t, double, double, char> item = std::make_tuple(0, 0, 0, 0);
	for(uint64_t id : ids) {
		if(!qt.getItem(item, id)) {
			std::cerr << "failed " << id << "\n";
		} else {
			std::cerr << "found: " << id << ", " << std::get<0>(item) << ", " << std::get<1>(item) << ", " << std::get<2>(item) << ", " << std::get<3>(item) << "\n";
			qt.updateItem(std::make_tuple(id, std::get<1>(item), std::get<2>(item), 'x'), id);
			qt.getItem(item, id);
			std::cerr << "updated: " << id << ", " << std::get<0>(item) << ", " << std::get<1>(item) << ", " << std::get<2>(item) << ", " << std::get<3>(item) << "\n";
		}
	}
	Bounds q(4, 4, 6, 6);
	std::list<std::tuple<uint64_t, double, double, char> > lst;
	qt.search(q, std::back_inserter(lst));
	std::cerr << "search\n";
	for(const auto& it : lst)
		std::cerr << std::get<0>(it) << ", " << std::get<1>(it) << ", " << std::get<2>(it) << "\n";
}


