/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include "grid.hpp"

using namespace geo::grid;

void usage() {
	std::cerr << "Usage: voidfill [options] <input raster> <output raster>\n"
			<< " -b  <band>       The band. Default 1.\n"
			<< " -m  <area>       Maximum area to fill. Square map units.\n"
			<< " -e               Fill voids on edges (otherwise don't).\n"
			<< " -d  <mode>       Mode: 0=min, 1=mean, 2=median, 3=max. Default 0.\n"
			<< " -n  <n>          The radius of the alpha disc.\n";
}

int fillVoids(Band<int>& mask, Band<float>& inrast, Band<float>& outrast, int col0, int row0, int col1, int row1, int mode) {

	const GridProps& props = inrast.props();
	float nd = props.nodata();
	float v, s;
	int ct = 0;
	std::vector<float> v0;

	// Initialize the values.
	switch(mode) {
	case 0: s = geo::maxvalue<float>(); break;
	case 3: s = geo::minvalue<float>(); break;
	default: s = 0; break;
	}

	// Find the edge pixel values.
	for(int r = row0; r <= row1; ++r) {
		for(int c = col0; c <= col1; ++c) {
			if(mask.get(c, r) != 2)
				continue;
			// Check the pixels around the masked pixel for values.
			for(int rr = r - 1; rr < r + 2; ++rr) {
				for(int cc = c - 1; cc < c + 2; ++cc) {
					if(props.hasCell(cc, rr) && (v = inrast.get(cc, rr)) != nd) {
						switch(mode) {
						case 0:
							if(v < s) s = v;
							break;
						case 3:
							if(v > s) s = v;
							break;
						case 2:
							v0.push_back(v);
							break;
						default:
							s += v;
							break;
						}
						++ct;
					}
				}
			}
		}
	}

	float m = nd;

	// Calculate the output value.
	// TODO: Add a spline, IDW (etc.) interp method.
	if(ct) {
		switch(mode) {
		case 0:
		case 3:
			m = s;
			break;
		case 2:
			if(v0.size() % 2 == 0) {
				std::sort(v0.begin(), v0.end());
				m = (v0[v0.size() / 2 - 1] + v0[v0.size() / 2]) / 2.0;
			} else {
				m = v0[v0.size() / 2];
			}
			break;
		default:
			m = s / ct;
			break;
		}
	}

	// Write the new value to the null region, set mask to 3.
	for(int r = row0; r <= row1; ++r) {
		for(int c = col0; c <= col1; ++c) {
			if(mask.get(c, r) == 2) {
				outrast.set(c, r, m);
				mask.set(c, r, 3);
			}
		}
	}

	return 0;
}

/**
 * Set all edge-contacting null pixel regions to 2 in the mask.
 */
void edgeMask(Band<float>& dem, Band<int>& mask) {
	float nd = dem.props().nodata();
	int cols = dem.props().cols();
	int rows = dem.props().rows();
	TargetFillOperator<float, int> op1(&dem, 1, &mask, 1, nd, 2); // Mark for filling (nd --> 2)
	int minc, minr, maxc, maxr, area;
	mask.fill(0);
	for(int c = 0; c < cols; ++c) {
		if(mask.get(c, 0) != 2)
			dem.floodFill(c, 0, op1, false, &minc, &minr, &maxc, &maxr, &area);
		if(mask.get(c, rows - 1) != 2)
			dem.floodFill(c, rows - 1, op1, false, &minc, &minr, &maxc, &maxr, &area);
	}
	for(int r = 0; r < rows; ++r) {
		if(mask.get(0, r) != 2)
			dem.floodFill(0, r, op1, false, &minc, &minr, &maxc, &maxr, &area);
		if(mask.get(cols - 1, r) != 2)
			dem.floodFill(cols - 1, r, op1, false, &minc, &minr, &maxc, &maxr, &area);
	}
}

/**
 * Set all non-edge, null pixel regions to 1 in the mask; zero everywhere else.
 */
void ndMask(Band<float>& dem, Band<int>& mask) {

	edgeMask(dem, mask);

	float nd = dem.props().nodata();
	int cols = dem.props().cols();
	int rows = dem.props().rows();
	for(int r = 0; r < rows; ++r) {
		for(int c = 0; c < cols; ++c) {
			if(mask.get(c, r) != 2 && dem.get(c, r) == nd) {
				mask.set(c, r, 1);
			} else {
				mask.set(c, r, 0);
			}
		}
	}
}

/**
 * Set all non-edge, null pixel regions to 1 in the mask, according to the alpha shape parameter; zero everywhere else.
 */
void hullMask(Band<float>& dem, Band<int>& mask, double alpha) {

	edgeMask(dem, mask);

	float nd = dem.props().nodata();
	int cols = dem.props().cols();
	int rows = dem.props().rows();

	// Build the round kernel.
	int rad = (int) std::ceil(std::abs(alpha / dem.props().resX()));
	std::vector<std::pair<int, int>> k;
	for(int rr = -rad; rr < rad + 1; ++rr) {
		for(int cc = -rad; cc < rad + 1; ++cc) {
			if(geo::sq(rr) + geo::sq(cc) <= rad)
				k.emplace_back(cc, rr);
		}
	}

	// Fill the mask with 3 outside the alpha shape.
	const GridProps& props = mask.props();
	int step = rows / 100;
	int cc, rr;
	for(int r = 0; r < rows; ++r) {
		if(r % step == 0)
			g_debug("Row: " << r << " of " << rows);
		for(int c = 0; c < cols; ++c) {
			if(mask.get(c, r) == 2) {
				bool fill = true;
				for(const auto& it : k) {
					cc = c + it.first;
					rr = r + it.second;
					if(props.hasCell(cc, rr) && mask.get(cc, rr) == 2 && dem.get(cc, rr) != nd) {
						fill = false;
						break;
					}
				}
				if(fill) {
					for(const auto& it : k) {
						cc = c + it.first;
						rr = r + it.second;
						if(props.hasCell(cc, rr))
							mask.set(cc, rr, 3);
					}
				}
			}
		}
	}

	// Set
	for(int r = 0; r < rows; ++r) {
		for(int c = 0; c < cols; ++c) {
			if(mask.get(c, r) != 3 && dem.get(c, r) == nd) {
				mask.set(c, r, 1);
			} else {
				mask.set(c, r, 0);
			}
		}
	}

}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	std::string infile;
	std::string outfile;
	int band = 1;
	float maxarea = geo::maxvalue<float>();
	int mode = 0;
	bool useGeomMask = true;
	int state = 0;
	float n = 100;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-b") {
			band = atoi(argv[++i]);
		} else if(v == "-m") {
			maxarea = atof(argv[++i]);
		} else if(v == "-d") {
			mode = atoi(argv[++i]);
		} else if(v == "-n") {
			n = atof(argv[++i]);
		} else if(state == 0) {
			infile = v;
			++state;
		} else if(state == 1) {
			outfile = v;
			++state;
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

	Band<float> inrast(infile, 0, false, true);

	Band<float> outrast(outfile, inrast.props(), true);
	inrast.writeTo(outrast);

	GridProps mprops(inrast.props());
	mprops.setWritable(true);
	mprops.setBands(1);
	mprops.setDataType(DataType::Byte);
	Band<int> mask("/tmp/mask.tif", mprops, true);
	mask.fill(0);

	int cols = mprops.cols();
	int rows = mprops.rows();

	geo::grid::TargetFillOperator<int, int> op1(&mask, 1, &mask, 1, 1, 2); // Mark for filling (1 --> 2)
	geo::grid::TargetFillOperator<int, int> op2(&mask, 1, &mask, 1, 2, 3); // Mark for filling (2 --> 3)
	int cmin = 0, cmax = 0, rmin = 0, rmax = 0, area = 0;

	if(useGeomMask){
		// Build concave hull to produce mask.
		std::cerr << "Building concave hull mask.\n";
		hullMask(inrast, mask, n);
	} else {
		std::cerr << "Building edge mask.\n";
		ndMask(inrast, mask);
	}

	int maxpxarea;
	double q = maxarea / geo::sq(mprops.resX());
	if(q > (double) geo::maxvalue<int>()) {
		maxpxarea = geo::maxvalue<int>();
	} else {
		maxpxarea = (int) q;
	}
	if(maxpxarea < 1) maxpxarea = 1;

	for(int row = 0; row < rows; ++row) {
		if(row % 100 == 0)
			std::cerr << "Filling. Row " << row << " of " << rows << "\n";
		for(int col = 0; col < cols; ++col) {

			// Only hit pixels marked with 1.
			if(mask.get(col, row) != 1)
				continue;

			// Fill the target region with 2.
			Band<int>::floodFill(col, row, op1, false, &cmin, &rmin, &cmax, &rmax, &area);

			if(area >= maxpxarea) {
				// The filled area was too large. Ignore it by setting it to 3.
				Band<int>::floodFill(col, row, op2, false, &cmin, &rmin, &cmax, &rmax, &area);
			} else {
				// Fill the voids.
				fillVoids(mask, inrast, outrast, cmin, rmin, cmax, rmax, mode);
			}

		}
	}

	return 0;
}


