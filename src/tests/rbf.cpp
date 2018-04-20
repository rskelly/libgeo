/*
 * rbf.cpp
 *
 *  Created on: Apr 19, 2018
 *      Author: rob
 */

#include <iostream>

#include "rbf.hpp"

class Pt{
public:
	double x, y, z;
	Pt(double x, double y, double z) :
		x(x), y(y), z(z) {}
	double operator[](size_t idx) const {
		switch(idx % 3) {
		case 1: return y;
		case 2: return z;
		default: return x;
		}
	}
};

int main() {

	RBF<Pt> rb;

	rb.add(Pt(0, 0, 10));
	rb.add(Pt(10, 10, 0));

	for(int r = 0; r < 10; ++r) {
		for(int c = 0; c < 10; ++c) {
			std::cerr << "r: " << r << "; c: " << c << "; z: " << rb.compute(Pt(c, r, 0)) << "\n";
		}
	}
	return 0;
}

