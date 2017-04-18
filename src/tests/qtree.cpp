#include "ds/qtree.hpp"
#include <list>

class PPoint {
public:
	double x;
	double y;
	PPoint(double x, double y) : x(x), y(y) {}
	PPoint() : x(0), y(0) {}
};

int main(int argc, char** argv) {
	Bounds b(0, 0, 100, 100);
	QTree<PPoint> t(b, 10);
	t.addItem(1, 1, PPoint(1, 1));
	t.addItem(2, 2, PPoint(2, 2));
	t.addItem(3, 3, PPoint(3, 3));
	t.addItem(99, 99, PPoint(99, 99));
	std::list<PPoint> res;
	t.search(99, 99, 1, std::back_inserter(res));
	for(const PPoint& tt : res)
		std::cerr << tt.x << ", " << tt.y << "\n";
}