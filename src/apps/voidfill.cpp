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
			<< " -e               Fill voids on edges (otherwise don't).\n"
			<< " -d  <mode>       Mode: 0=min, 1=mean, 2=median, 3=max. Default 0.\n"
			<< " -n  <N>          N for concave hull. Usually 0-5. Default 2.5.\n";
}

int fillVoid(Band<int>& mask, Band<float>& rast, int col, int row, int mode) {
	if(mask.get(col, row) != 1)
		return 0;
	const GridProps& props = mask.props();
	double nd = rast.props().nodata();
	double v, s = 0;
	int ct = 0;
	std::vector<float> v0;
	for(int r = std::max(0, row - 1); r < std::min(props.rows(), row + 2); ++r) {
		for(int c = std::max(0, col - 1); c < std::min(props.cols(), col + 2); ++c) {
			if((v = rast.get(c, r)) != nd) {
				if(mode == 2) {
					v0.push_back(v);
				} else {
					s += v;
				}
				++ct;
			}
		}
	}
	if(ct) {
		float m;
		if(mode == 2) {
			if(v0.size() % 2 == 0) {
				std::sort(v0.begin(), v0.end());
				m = (v0[v0.size() / 2 - 1] + v0[v0.size() / 2]) / 2.0;
			} else {
				m = v0[v0.size() / 2];
			}
		} else {
			m = s / ct;
		}
		rast.set(col, row, m);
		mask.set(col, row, 2);
		return 1;
	}
	return 0;
}

GEOSGeometry* makeLine(double x1, double y1, double x2, double y2, GEOSContextHandle_t gctx) {
	GEOSCoordSequence* seq = GEOSCoordSeq_create_r(gctx, 2, 2);
	GEOSCoordSeq_setX_r(gctx, seq, 0, x1);
	GEOSCoordSeq_setY_r(gctx, seq, 0, y1);
	GEOSCoordSeq_setX_r(gctx, seq, 1, x2);
	GEOSCoordSeq_setY_r(gctx, seq, 1, y2);
	return GEOSGeom_createLineString_r(gctx, seq);
}

GEOSGeometry* makePoint(double x, double y, GEOSContextHandle_t gctx) {
	GEOSCoordSequence* seq = GEOSCoordSeq_create_r(gctx, 1, 2);
	GEOSCoordSeq_setX_r(gctx, seq, 0, x);
	GEOSCoordSeq_setY_r(gctx, seq, 0, y);
	return GEOSGeom_createPoint_r(gctx, seq);
}

bool getXY(const GEOSGeometry* geom, int i, double& x, double& y, GEOSContextHandle_t gctx) {
	const GEOSGeometry* pt = GEOSGeomGetPointN_r(gctx, geom, i);
	if(GEOSGeomGetX_r(gctx, pt, &x)) {
		if(GEOSGeomGetY_r(gctx, pt, &y)) {
			return true;
		}
	}
	return false;
}

// https://www.iis.sinica.edu.tw/page/jise/2012/201205_10.pdf
GEOSGeometry* concaveHull(const std::vector<std::pair<double, double>>& pts, double N, GEOSContextHandle_t gctx) {

	// Build the tree and list of points.
	std::vector<GEOSGeometry*> treePts;
	for(size_t i = 0; i < pts.size(); ++i)
		treePts.push_back(makePoint(pts[i].first, pts[i].second, gctx));

	// Build the hull and get the exterior ring.
	GEOSGeometry* chull;
	const GEOSGeometry* er;
	{
		GEOSGeometry* mp = GEOSGeom_createCollection_r(gctx, GEOS_MULTIPOINT, treePts.data(), treePts.size());
		chull = GEOSConvexHull_r(gctx, mp);
		er = GEOSGetExteriorRing_r(gctx, chull);
	}

	// Make the list of line segments.
	double x1, y1, x2, y2;
	int numPts = GEOSGeomGetNumPoints_r(gctx, er);
	std::vector<GEOSGeometry*> lines;
	GEOSSTRtree* tree = GEOSSTRtree_create_r(gctx, pts.size());

	for(int i = 0; i < numPts; ++i) {
		getXY(er, i, x1, y1, gctx);
		getXY(er, (i + 1) % numPts, x2, y2, gctx);
		GEOSGeometry* line = makeLine(x1, y1, x2, y2, gctx);
		lines.push_back(line);
		GEOSGeometry* pt = GEOSGeomGetPointN_r(gctx, er, i);
		if(1 == GEOSRelatePattern_r(gctx, chull, pt, "T**FF*FF*"))
			GEOSSTRtree_insert_r(gctx, tree, pt, pt);
	}

	// Do the "digging" operation.
	double eh, dd, dd1, dd2, x, y;
	//bool found = false;
	//do {
		std::vector<GEOSGeometry*> tmp;
		std::vector<size_t> destid;
		//found = false;
		for(int i = 0; i < (int) lines.size(); ++i) {
			GEOSGeometry* line = lines[i];
			const GEOSGeometry* l1 = lines[(i - 1) % lines.size()];
			const GEOSGeometry* l2 = lines[(i + 1) % lines.size()];
			const GEOSGeometry* pt = GEOSSTRtree_nearest_r(gctx, tree, line);
			GEOSDistance_r(gctx, pt, line, &dd);
			GEOSDistance_r(gctx, pt, l1, &dd1);
			GEOSDistance_r(gctx, pt, l2, &dd2);
			if(dd1 < dd || dd2 < dd)
				continue;
			GEOSGeomGetLength_r(gctx, line, &eh);
			if(dd > 0 && eh > 0 && eh / dd > N) {
				getXY(pt, 0, x, y, gctx);
				tmp.push_back(makeLine(x1, y1, x, y, gctx));
				tmp.push_back(makeLine(x, y, x2, y2, gctx));
				destid.push_back(i);
				//found = true;
				continue;
			} else {
				tmp.push_back(line);
			}
		}
		for(size_t i : destid)
			GEOSGeom_destroy_r(gctx, lines[i]);
		tmp.swap(lines);
		tmp.clear();
	//} while(found);

	//
	GEOSGeometry* hull = GEOSPolygonize_r(gctx, lines.data(), lines.size());

	GEOSGeom_destroy_r(gctx, chull);

	return hull;
}

#include <stdarg.h>

void msg(const char* fmt, ...) {
	va_list v;
	va_start(v, 1);
	char* a = va_arg(v, char*);
	char buf[256];
	sprintf(buf, fmt, a);
	g_warn(buf);
}

void hullMask(Band<int>& mask, const std::vector<std::pair<double, double>>& pts, double N) {
	GEOSContextHandle_t gctx = initGEOS_r(msg, msg);
	GEOSGeometry* hull = concaveHull(pts, N, gctx);
	for(int r = 0; r < mask.props().rows(); ++r) {
		for(int c = 0; c < mask.props().cols(); ++c) {
			GEOSGeometry* pt = makePoint((double) c, (double) r, gctx);
			if(GEOSContains_r(gctx, hull, pt))
				mask.set(c, r, 1);
			GEOSGeom_destroy_r(gctx, pt);
		}
	}
	GEOSGeom_destroy_r(gctx, hull);
	finishGEOS_r(gctx);
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
	int state = 0;
	double n = 2.5;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-e") {
			edges = true;
		} else if(v == "-b") {
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

	Band<int> mask;
	Band<float> rast;
	int pxarea;
	{
		rast.init(infile, 0, false, true, nullptr);
		GridProps props(rast.props());
		pxarea = (int) std::ceil(maxarea / std::abs(props.resX()) * std::abs(props.resY()));
		props.setWritable(true);
		props.setBands(1);
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

	// Temp nodata for mask.
	double tnd = rast.stats().min - 1;

	if(useGeomMask){
		// Build concave hull to produce mask.
		std::cerr << "Building concave hull mask.\n";
		std::vector<std::pair<double, double>> chull;
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

		hullMask(mask, chull, n);

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
							count += fillVoid(mask, rast, c, r, mode);
						}
					}
				} while(count > 0);
			}
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


