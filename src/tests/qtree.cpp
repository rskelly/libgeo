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

	int count = 100;
	if(argc == 2)
		count = atoi(argv[1]);

	Bounds b(0, 0, 100, 100);

	double sx = Util::random(0, 100);
	double sy = Util::random(0, 100);
	double sr = 10.0;
	int sc = 0;

	QTree<PPoint> t(b, 10);
	Stopwatch sw;
	sw.start();
	for(int i = 0; i < count; ++i) {
		double x = Util::random(0, 100);
		double y = Util::random(0, 100);
		t.addItem(x, y,  PPoint(x, y));
		//std::cerr << "added: " << x << ", " << y << "\n";
		if(g_sq(x - sx) + g_sq(y - sy) <= g_sq(sr)) {
			//std::cerr << "1," << x << ", " << y << "\n";
			++sc;
		}
	}

	sw.stop();
	std::cerr << "added " << count << " in " << sw.millis() << "\n";

	std::list<PPoint> res;
	t.search(sx, sy, sr, std::back_inserter(res));
	std::cerr << "found " << res.size() << " items out of " << sc << "\n";
	//for(const PPoint& p : res)
		//std::cerr << "2," << p.x << ", " << p.y << "\n";

	for(int i = 0; i < 10; ++i) {
		res.clear();
		sw.reset();
		sw.start();
		t.search(Util::random(0, 100), Util::random(0, 100), 5, std::back_inserter(res));
		sw.stop();
		std::cerr << "found by search " << res.size() << " in " << sw.millis() << "\n";	
		//for(const PPoint& tt : res)
		//	std::cerr << "found: " << tt.x << ", " << tt.y << "\n";
	}

	for(int i = 0; i < 10; ++i) {
		double x = Util::random(0, 100);
		double y = Util::random(0, 100);
		res.clear();
		sw.reset();
		sw.start();
		t.nearest(x, y, 1, std::back_inserter(res));
		sw.stop();
		PPoint& d = res.front();
		double dist = std::sqrt(g_sq(x - d.x) + g_sq(y - d.y));
		std::cerr << "found nearest " << res.size() << " in " << sw.millis() << " at " << dist << "\n";	
		//for(const PPoint& tt : res)
		//	std::cerr << "found: " << tt.x << ", " << tt.y << "\n";
	}
}