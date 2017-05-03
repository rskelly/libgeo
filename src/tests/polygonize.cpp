/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include "raster.hpp"

int main(int argc, char** argv) {

	using namespace geo::raster;

	Raster test("/home/rob/Documents/git/libgeo/data/tests/polygons.tif");
	test.polygonize("/home/rob/Documents/git/libgeo/data/tests/output/polygons.sqlite", "polygons", "sqlite", 26910, 1, false, false, nullptr, nullptr);

	return 0;
}


