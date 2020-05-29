/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include <geos_c.h>

#include "grid.hpp"

using namespace geo::grid;

void usage() {
	std::cerr << "Usage: voidfill [options] <input raster> <output raster>\n"
			<< " -b  <band>       The band. Default 1.\n"
			<< " -m  <area>       Maximum area to fill. Square map units.\n"
			<< " -e               Fill voids on edges (otherwise don't).\n";
}

int fillVoid(Band<int>& mask, Band<float>& rast, int col, int row) {
	if(mask.get(col, row) != 1)
		return 0;
	const GridProps& props = mask.props();
	double nd = rast.props().nodata();
	double v, s = 0;
	int ct = 0;
	for(int r = std::max(0, row - 1); r < std::min(props.rows(), row + 2); ++r) {
		for(int c = std::max(0, col - 1); c < std::min(props.cols(), col + 2); ++c) {
			if((v = rast.get(c, r)) != nd) {
				s += v;
				++ct;
			}
		}
	}
	if(ct) {
		rast.set(col, row, s / ct - std::numeric_limits<double>::min());
		mask.set(col, row, 2);
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
	bool useGeomMask = true;

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

	Band<int> mask;
	Band<float> rast;
	int pxarea;
	{
		rast.init(infile, 0, false, true, nullptr);
		GridProps props(rast.props());
		pxarea = (int) std::ceil(maxarea / std::abs(props.resX()) * std::abs(props.resY()));
		props.setWritable(true);
		props.setBands(1);
		rast.init(props, true);
		props.setDataType(DataType::Byte);
		mask.init(props, true);
		mask.fill(0);
	}


	const GridProps& mprops = mask.props();
	int cols = mprops.cols();
	int rows = mprops.rows();
	const GridProps& rprops = rast.props();
	double nd = rprops.nodata();

	geo::grid::TargetFillOperator<float, int> op1(&rast, 1, &mask, 1, nd, 1); // Mark for filling (nd --> 1)
	geo::grid::TargetFillOperator<int, int> op2(&mask, 1, &mask, 1, 1, 2); // Mark for filling (1 --> 2)
	int cmin = 0, cmax = 0, rmin = 0, rmax = 0, area = 0;
	int v;

	// Temp nodata for mask.
	double tnd = nd - 1;

	if(useGeomMask){
		// Build concave hull to produce mask.
		std::cerr << "Building concave hull mask.\n";
		std::vector<std::pair<int, int>> chull;
		double v;
		for(int row = 0; row < rows; ++row) {
			if(row % 100 == 0)
				std::cerr << "Row " << row << " of " << rows << "\n";
			for(int col = 0; col < cols; ++col) {

				if((v = rast.get(col, row)) == nd)
					continue;

				bool isEdge = false;
				for(int rr = row - 1; !isEdge && rr < row + 2; ++rr) {
					for(int cc = col - 1; !isEdge && cc < col + 2; ++cc) {
						if(rprops.hasCell(cc, rr))
							isEdge = (rast.get(cc, rr) == nd);
					}
				}

				if(isEdge)
					chull.emplace_back(col, row);
			}
		}

		GEOSContextHandle_t gctx = initGEOS_r(0, 0);
		GEOSCoordSequence* seq = GEOSCoordSeq_create_r(gctx, chull.size(), 2);
		GEOSGeometry* mp = GEOSGeom_createPoint_r(gctx, seq);
		GEOSGeometry* hull = GEOSConvexHull_r(gctx, mp);
		GEOSGeometry* shull = GEOSSimplify_r(gctx, hull, 500.0);
		GEOSGeom_destroy_r(gctx, mp);
		GEOSGeom_destroy_r(gctx, hull);

		GEOSCoordSequence* cseq = GEOSCoordSeq_create_r(gctx, 1, 2);
		for(int row = 0; row < rows; ++row) {
			if(row % 100 == 0)
				std::cerr << "Row " << row << " of " << rows << "\n";
			for(int col = 0; col < cols; ++col) {
				GEOSCoordSeq_setX_r(gctx, cseq, 0, (double) col);
				GEOSCoordSeq_setY_r(gctx, cseq, 0, (double) row);
				GEOSGeometry* pt = GEOSGeom_createPoint_r(gctx, cseq);
				if(GEOSContains_r(gctx, shull, pt) == 0 && rast.get(col, row) == nd)
					rast.set(col, row, tnd);
				GEOSGeom_destroy_r(gctx, pt);
			}
		}
		GEOSCoordSeq_destroy_r(gctx, cseq);
	}

	for(int row = 0; row < rows; ++row) {
		if(row % 100 == 0)
			std::cerr << "Row " << row << " of " << rows << "\n";
		for(int col = 0; col < cols; ++col) {

			if(mask.get(col, row) != 0 || (v = rast.get(col, row)) != nd)
				continue;

			Band<float>::floodFill(col, row, op1, false, &cmin, &rmin, &cmax, &rmax, &area);

			// Skip if the area is too large, or if it's on the edge and
			// edges are not desired.
			if(area == 0 || (pxarea > 0 && area > pxarea) || (!edges && (cmin == 0 || rmin == 0 || cmax >= cols - 1 || rmax >= rows - 1))) {
				Band<int>::floodFill(col, row, op2, false);
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
		}
	}

	for(int row = 0; row < rows; ++row) {
		if(row % 100 == 0)
			std::cerr << "Row " << row << " of " << rows << "\n";
		for(int col = 0; col < cols; ++col) {
			if(rast.get(col, row) == tnd)
				rast.set(col, row, nd);
		}
	}

	{
		GridProps props(rast.props());
		props.setWritable(true);
		Band<float> output(outfile, props);
		rast.writeTo(output);
	}

	return 0;
}


