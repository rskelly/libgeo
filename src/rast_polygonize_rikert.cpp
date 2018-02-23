/*
 * rast_polygonize_rikert.cpp
 *
 *  Created on: Feb 22, 2018
 *      Author: rob
 */

#include "raster.hpp"

using namespace geo::raster;

enum State {

};
void readBlock(Raster& rast, int cols, int row, int* block) {
	rast.getIntRow(0, row, block);
	rast.getIntRow(0, row + 1, block + cols);
}

const int BORDER = 1;
const int NOBORDER = 2;
const int NODE = 4;

int getState(int col, int cols, std::vector<int>& block) {
	int a = block[col];
	int b = block[col + 1];
	int c = block[col + cols];
	int d = block[col + cols + 1];
	if((a == c && b == d)
			|| (b == d && c == d && a != b)
			|| (a == b && b == d && a != c)) {
		return NOBORDER;
	} else if((a != b && a != c && b == d)) {
		return NOBORDER | NODE;
	} else if((a == c && c == d && a != b)
			|| ( a == b && a == c && a != d)) {
		return BORDER;
	} else if((a != b && a != c && c == d)
			|| (a != b && a != c && a != d)
			|| (a == b && a != c && a != d)) {
		return BORDER | NODE;
	} else {
		g_runerr("Invalid event.");
	}
}

// Vectorize the raster.
void Raster::polygonize(const std::string &filename, const std::string &layerName,
    const std::string &driver, uint16_t srid = 0, uint16_t band = 1, uint16_t threads = 1,
	int bufSize = 0, bool removeHoles = false, bool removeDangles = false,
	geo::util::Status *status = nullptr, bool *cancel = nullptr, const std::string& mask = "") {

	Raster rast(filename);
	const GridProps& props = rast.props();
	int cols = props.cols();
	int rows = props.rows();

	// Block containing two rows of pixels.
	std::vector<int> block(cols * 2);

	int state;

	for(int row = 0; row < rows - 1; ++row) {

		readBlock(rast, cols, row, block.data());

		for(int col = 0; col < cols - 1; ++col) {

			state = getState(col, cols, block);



		}
	}
}


