#include "ds/qtree.hpp"
#include "util.hpp"
#include <list>

using namespace geo::util;

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
	for(int i = 0; i < 1000000; ++i) {
		double x = Util::random(0, 100);
		double y = Util::random(0, 100);
		t.addItem(x, y,  PPoint(x, y));
		//std::cerr << "added: " << x << ", " << y << "\n";
	}

	std::list<PPoint> res;
	t.search(50, 50, 50, std::back_inserter(res));
	for(const PPoint& tt : res)
		std::cerr << "found: " << tt.x << ", " << tt.y << "\n";
}