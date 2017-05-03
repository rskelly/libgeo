/*
 * smooth.cpp
 *
 *  Created on: May 3, 2017
 *      Author: rob
 */

#include "raster.hpp"

int main(int argc, char** argv) {

	using namespace geo::raster;

	Raster r("/home/rob/Documents/git/treetops/_data/carmanah.tif");
	GridProps p(r.props());
	p.setWritable(true);
	Raster s("../data/tests/output/smoothed.tif", p);

	r.smooth(s, 0.8, 5);
}


