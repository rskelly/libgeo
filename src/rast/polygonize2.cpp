#include <vector>

#include "raster.hpp"

#define NODE 1
#define EDGE 2
#define A 4
#define B 8
#define C 16
#define D 32
#define E 64
#define F 128
#define G 256
#define H 512
#define I 1024
#define J 2048
#define K 4096
#define L 8192
#define M 16384
#define N 32768

int getState(const std::vector<int>& row1, const std::vector<int>& row2, int col) {
	int p0 = row1[col - 1];
	int p1 = row1[col];
	int p2 = row2[col - 1];
	int p3 = row2[col];
	if(p0 == p1 && p2 == p3 && p0 != p2) {
		return EDGE|A;
	} else if(p0 == p2 && p1 == p3 && p0 != p1) {
		return EDGE|B;
	} else if(p0 == p1 && p0 == p2 && p0 != p3) {
		return EDGE|C;
	} else if(p0 == p1 && p1 == p3 && p0 != p2) {
		return EDGE|D;
	} else if(p1 == p2 && p2 == p3 && p0 != p1) {
		return EDGE|E;
	} else if(p0 == p2 && p0 == p3 && p0 != p1) {
		return EDGE|F;
	} else if(p0 != p1 && p0 != p2 && p0 != p3 && p1 != p3 && p3 != p2) {
		return NODE|G;
	} else if(p0 == p1 && p0 != p2 && p0 != p3 && p2 != p3) {
		return NODE|H;
	} else if(p0 != p1 && p0 != p2 && p2 == p3 && p1 != p3) {
		return NODE|I;
	} else if(p0 == p2 && p0 != p1 && p0 != p2 && p1 != p3) {
		return NODE|J;
	} else if(p0 != p1 && p0 != p2 && p1 == p2 && p2 != p3) {
		return NODE|K;
	} else if(p0 == p3 && p0 != p1 && p0 != p2 && p1 != p2) {
		return NODE|L;
	} else if(p0 != p1 && p1 == p2 && p0 != p2 && p0 != p3) {
		return NODE|M;
	} else if(p0 == p3 && p1 == p2 && p0 != p2) {
		return NODE|N;
	} else {
		return 0;
	}
}

using namespace geo::raster;

class Arc {
public:
	int leftId;
	int rightId;
	int col, row;
	Arc() : Arc(col, row, leftId, rightId) {}
	Arc(int col, int row, int leftId, int rightId) :
		col(col), row(row),
		leftId(leftId), rightId(rightId) {}
};

class

void Raster::polygonize(const std::string &filename, const std::string &layerName,
        const std::string &driver, int srid, int band, bool removeHoles,
		bool removeDangles,	geo::util::Status *status = nullptr, bool *cancel,
		const std::string& mask, int threads) {

	int cols = props().cols();
	int rows = props().rows();

	std::vector<int> buf1(cols + 2);
	std::vector<int> buf2(cols + 2);

	for(int row = 0; row <= rows; ++row) {

		std::fill(buf1.begin(), buf1.end(), 0);
		std::fill(buf2.begin(), buf2.end(), 0);
		if(row < rows)
			getIntRow(row, band, buf2.data() + sizeof(int));
		if(row > 0)
			getIntRow(row - 1, band, buf1.data() + sizeof(int));

		std::vector<Arc> arcs;
		int lastArc = -1;

		for(int col = 0; col < cols; ++col) {

			int state = getState(buf1, buf2, col + 1);

			if(state & NODE) {
				arcs.emplace_back(col, row, 0, 0);
				++lastArc;
			} else (state & EDGE) {

			}
		}
	}
}
