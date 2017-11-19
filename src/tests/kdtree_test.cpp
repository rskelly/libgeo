/*
 * kdtree.cpp
 *
 *  Created on: Nov 18, 2017
 *      Author: rob
 */

#include <list>
#include <iostream>

#include "ds/kdtree.hpp"

class Pt {
public:
	double x, y, z;
	Pt(double x, double y, double z) :
		x(x), y(y), z(z) {}
	double operator[](size_t idx) const {
		switch(idx % 3) {
		case 0: return x;
		case 1: return y;
		case 2: return z;
		default: return x; // never happen
		}
	}
};

int main(int argc, char** argv) {

	using namespace geo::ds;

	std::vector<Pt*> pts;
	for(int r = 0; r < 10; ++r) {
		for(int c = 0; c < 10; ++c) {
			Pt* p = new Pt((double) c, (double) r, (double) (c * r % 11));
			pts.push_back(p);
		}
	}

	Pt pt(3, 3, 9);
	std::vector<Pt*> output;
	std::vector<double> distances;

	KDTree<Pt> kt(3);
	kt.add(pts.begin(), pts.end());
	kt.build();
	kt.knn(pt, 4, std::back_inserter(output), std::back_inserter(distances));

	std::cerr << "Found items:\n";
	for(size_t i = 0; i < output.size(); ++i)
		std::cerr << (*output[i])[0] << "," << (*output[i])[1] << "," << (*output[i])[2] << "; " << distances[i] << "\n";

}

