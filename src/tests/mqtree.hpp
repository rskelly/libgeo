#include <list>

#include "util.hpp"
#include "ds/mqtree.hpp"

#include "test.hpp"

using namespace geo::util;
using namespace geo::ds;

namespace geo {
namespace test {

class PPoint {
public:
	double _x;
	double _y;
	size_t pos;
	PPoint(double x, double y) : _x(x), _y(y), pos(0) {}
	PPoint() : _x(0), _y(0), pos(0) {}
	double x() const {
		return _x;
	}
	double y() const {
		return _y;
	}
	double operator[](int i) const {
		switch(i % 2) {
		case 0: return _x;
		default: return _y;
		}
	}
};

class mqtree : public test {
public:

	bool run(int argc, char** argv) {

		Bounds b(0, 0, 100, 100);

		geo::ds::mqtree<PPoint> t(100, 0, 100, 100);

		for(double y = 100; y < 200; y += 1) {
			for(double x = 100; x < 200; x += 1) {
				t.add(PPoint(x, y));
			}
		}
		t.build();

		std::list<PPoint> res;

		for(double y = 100; y < 200; y += 1) {
			for(double x = 100; x < 200; x += 1) {
				t.search(PPoint(x, y), 10, std::back_inserter(res));
				std::cout << res.size() << "\n";
				res.clear();
			}
		}

		return true;
	}
};

}
}
