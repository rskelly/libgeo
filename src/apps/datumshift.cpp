/*
 * datumshift.cpp
 *
 *  Created on: Feb 14, 2018
 *      Author: rob
 */

#include <grid.hpp>
#include <proj_api.h>


using namespace geo::raster;

double interp(Raster& grid, double x, double y) {

	const GridProps& props = grid.props();
	int col = props.toCol(x);
	int row = props.toRow(y);
	double cx = props.toCentroidX(col);
	double cy = props.toCentroidY(row);
	double resX = props.resolutionX();
	double resY = props.resolutionY();

	// Determine which four cells capture the point.
	int col0, col1;
	int row0, row1;
	if((x < cx && resX > 0) || (x > cx && resX < 0)) {
		col0 = col - 1;
		col1 = col;
	} else {
		col0 = col;
		col1 = col + 1;
	}
	if((y < cy && resY > 0) || (y > cy && resY < 0)) {
		row0 = row - 1;
		row1 = row;
	} else {
		row0 = row;
		row1 = row + 1;
	}

	// If out of bounds, quit.
	if(row0 < 0 || row1 >= props.rows() || col0 < 0 || col1 >= props.cols())
		return std::nan("");

	// Interpolate...
	double cx0 = props.toCentroidX(col0);
	double cx1 = props.toCentroidX(col1);
	double cy0 = props.toCentroidY(row0);
	double cy1 = props.toCentroidY(row1);

	double z0 = grid.getFloat(col0, row0); // TL
	double z1 = grid.getFloat(col1, row0); // TR
	double z2 = grid.getFloat(col0, row1); // BL
	double z3 = grid.getFloat(col1, row1); // BR

	double iz0 = z0 + (z1 - z0) * ((x - cx0) / (cx1 - cx0)); // Top horizontal interp.
	double iz1 = z2 + (z3 - z2) * ((x - cx0) / (cx1 - cx0)); // Bottom horizontal interp.

	double iz3 = iz0 + (iz1 - iz0) * ((y - cy0) / (cy1 - cy0)); // Vertical interp.

	return iz3;
}

double toDeg(double c) {
	return c * 180.0 / M_PI;
}

void usage() {
	std::cerr << "Usage: datumshift <input file> <output file> <datum file> <source srid> <grid srid>\n";
}

int main(int argc, char** argv) {

	if(argc < 4) {
		usage();
		return 1;
	}

	std::string inputfile(argv[1]);
	std::string outputfile(argv[2]);
	std::string gridfile(argv[3]);
	int ssrid = atoi(argv[4]);
	int gsrid = atoi(argv[5]);

	Raster input(inputfile);
	const GridProps& iprops = input.props();
	double nodata = iprops.nodata();

	GridProps oprops(iprops);
	oprops.setWritable(true);
	Raster output(outputfile, oprops);

	Raster grid(gridfile);

	int cols = iprops.cols();
	int rows = iprops.rows();

	std::stringstream ss;
	ss << "+init=epsg:" << ssrid;
	projPJ sp = pj_init_plus(ss.str().c_str());

	std::stringstream gs;
	gs << "+init=epsg:" << gsrid;
	projPJ gp = pj_init_plus(gs.str().c_str());

	double iz;

	for(int row = 0; row < rows; ++row) {
		for(int col = 0; col < cols; ++col) {

			if((iz = input.getFloat(col, row)) == nodata)
				continue;

			double x = iprops.toCentroidX(col);
			double y = iprops.toCentroidY(row);
			if(pj_transform(sp, gp, 1, 1, &x, &y, NULL))
				g_runerr("Failed to transform point.");

			double z = interp(grid, toDeg(x), toDeg(y));

			output.setFloat(col, row, std::isnan(z) ? nodata : iz - z);

		}
	}
}
