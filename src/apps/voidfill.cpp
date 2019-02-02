/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include "raster.hpp"

using namespace geo::raster;

void usage() {
	std::cerr << "Usage: voidfill [options] <input raster> <output raster>\n"
			<< " -b  <band>       The band. Default 1.\n"
			<< " -m  <area>       Maximum area to fill. Square map units.\n"
			<< " -e               Fill voids on edges (otherwise don't).\n";
}

int fillVoid(Grid& mask, Grid& rast, int col, int row) {
	if(mask.getInt(col, row, 1) != 1)
		return 0;
	const GridProps& props = mask.props();
	double nd = rast.props().nodata();
	double v, s = 0;
	int ct = 0;
	for(int r = std::max(0, row - 1); r < std::min(props.rows(), row + 2); ++r) {
		for(int c = std::max(0, col - 1); c < std::min(props.cols(), col + 2); ++c) {
			if((v = rast.getFloat(c, r, 1)) != nd) {
				s += v;
				++ct;
			}
		}
	}
	if(ct) {
		rast.setFloat(col, row, s / ct - std::numeric_limits<double>::min(), 1);
		mask.setInt(col, row, 2, 1);
		return 1;
	}
	return 0;
}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	std::string infile;
	std::string outfile;
	int band = 1;
	double maxarea = 0;
	bool edges = false;
	int mode = 0;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-e") {
			edges = true;
		} else if(v == "-b") {
			band = atoi(argv[++i]);
		} else if(v == "-m") {
			maxarea = atof(argv[++i]);
		} else if(mode == 0) {
			infile = v;
			++mode;
		} else if(mode == 1) {
			outfile = v;
			++mode;
		}
	}

	if(band < 1) {
		std::cerr << "Illegal band number: " << band << "\n";
		usage();
		return 1;
	}

	if(infile.empty() || outfile.empty()) {
		std::cerr << "Input and output filenames required.\n";
		usage();
		return 1;
	}

	MemRaster mask;
	MemRaster rast;
	int pxarea;
	{
		Raster input(infile);
		GridProps props(input.props());
		pxarea = (int) std::ceil(maxarea / std::abs(props.resolutionX()) * std::abs(props.resolutionY()));
		props.setWritable(true);
		props.setBands(1);
		rast.init(props, true);
		props.setDataType(DataType::Byte);
		mask.init(props, true);
		input.writeTo(rast, props.cols(), props.rows(), 0, 0, 0, 0, band, 1);
		mask.fillInt(0, 1);
	}

	const GridProps& props = mask.props();
	int cols = props.cols();
	int rows = props.rows();
	double nd = rast.props().nodata();

	geo::raster::TargetFillOperator<double, int> op1(&rast, 1, &mask, 1, nd, 1); // Mark for filling (nd --> 1)
	geo::raster::TargetFillOperator<int, int> op2(&mask, 1, &mask, 1, 1, 2); // Mark for filling (1 --> 2)
	int cmin = 0, cmax = 0, rmin = 0, rmax = 0, area = 0;
	int v;

	for(int row = 0; row < rows; ++row) {
		if(row % 100 == 0)
			std::cerr << "Row " << row << " of " << rows << "\n";
		for(int col = 0; col < cols; ++col) {

			if(mask.getInt(col, row, 1) != 0 || (v = rast.getFloat(col, row, 1)) != nd)
				continue;

			Grid::floodFill(col, row, op1, false, &cmin, &rmin, &cmax, &rmax, &area);

			// Skip if the area is too large, or if it's on the edge and
			// edges are not desired.
			if(area == 0 || (pxarea > 0 && area > pxarea) || (!edges && (cmin == 0 || rmin == 0 || cmax >= cols - 1 || rmax >= rows - 1))) {
				Grid::floodFill(col, row, op2, false);
				continue;
			}

			int count = 0;
			do {
				count = 0;
				for(int r = rmin; r <= rmax; ++r) {
					for(int c = cmin; c <= cmax; ++c) {
						count += fillVoid(mask, rast, c, r);
					}
				}
			} while(count > 0);

			//Grid::floodFill(col, row, op2, false);

		}
	}

	{
		GridProps props(rast.props());
		props.setWritable(true);
		Raster output(outfile, props);
		rast.writeTo(output, props.cols(), props.rows(), 0, 0, 0, 0, 1, 1);
	}

	return 0;
}


