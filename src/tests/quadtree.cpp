/*
 * quadtree.cpp
 *
 *  Created on: Apr 14, 2017
 *      Author: rob
 */

#include "ds/quadtree.hpp"
#include "util.hpp"

class Boop {
public:
	double x;
	double y;
	Boop() : x(0), y(0) {}
	Boop(double x, double y) : x(x), y(y) {}
	std::string str() const {
		std::stringstream ss;
		ss << x << ", " << y;
		return ss.str();
	}
};

int main(int argc, char** argv) {

	using namespace geo::util;
	using namespace geo::ds;

	Bounds b(0, 0, 100, 100);

	FileQuadTree<Boop> qt(b, 1, 100);

	std::vector<uint64_t> ids;
	for(int i = 0; i < 10; ++i) {
		double x = Util::random(0.0, 10.0);
		double y = Util::random(0.0, 10.0);
		Boop a(x, y);
		std::cerr << x << ", " << y << "\n";
		ids.push_back(qt.addItem(x, y, a));
	}

	std::cerr << "ids\n";
	FileQuadTreeItem<Boop> item;
	for(uint64_t id : ids) {
		if(!qt.getItem(item, id)) {
			std::cerr << "failed " << id << "\n";
		} else {
			std::cerr << "found: " << id << ", " << item.id << ", " << item.x << ", " << item.y << ", " << item.item.str() << "\n";
			item.x = -1;
			qt.updateItem(item, id);
			qt.getItem(item, id);
			std::cerr << "updated: " << id << ", " << item.id << ", " << item.x << ", " << item.y << ", " << item.item.str() << "\n";
		}
	}
	Bounds q(-2, -2, 6, 6);
	std::list<FileQuadTreeItem<Boop> > lst;
	qt.search(q, std::back_inserter(lst));
	std::cerr << "search\n";
	for(const auto& item : lst)
		std::cerr << "found: " << item.id << ", " << item.x << ", " << item.y << ", " << item.item.str() << "\n";
}


