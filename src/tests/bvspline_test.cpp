/*
 * bvspline_test.cpp
 *
 *  Created on: Oct 26, 2019
 *      Author: rob
 */

#include <vector>
#include <fstream>

#include "util.hpp"

using namespace geo::util;
using namespace geo::util::csv;

int main(int argc, char** argv) {

	CSV csv("/home/rob/Desktop/ec/tmp/points.csv", true);

	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
	std::vector<double> w;

	std::vector<CSVValue> cx = csv.column("X");
	std::vector<CSVValue> cy = csv.column("Y");
	std::vector<CSVValue> cz = csv.column("top_z");

	double minx = 99999999.0;
	double miny = 99999999.0;
	double maxx = -99999999.0;
	double maxy = -99999999.0;
	double xx, yy;
	for(size_t i = 0; i < cx.size(); ++i) {
		x.push_back((xx = cx[i].asDouble()));
		y.push_back((yy = cy[i].asDouble()));
		z.push_back(cz[i].asDouble());
		if(xx < minx) minx = xx;
		if(yy < miny) miny = yy;
		if(xx > maxx) maxx = xx;
		if(yy > maxy) maxy = yy;
	}

	minx -= 100;
	miny -= 100;
	maxx += 100;
	maxy += 100;

	double res = 0.5;
	int cols = (int) (maxx - minx) / res + 1;
	int rows = (int) (maxy - miny) / res + 1;

	double smooth = x.size()  / 8;//- std::sqrt(x.size() * 2);
	geo::util::BivariateSpline bv;
	std::cout << bv.init(smooth, x, y, z, w, minx, miny, maxx, maxy) << "\n";

	x.clear();
	y.clear();

	for(int c = 0; c < cols; ++c)
		x.push_back(minx + c * res + res * 0.5);
	for(int r = 0; r < rows; ++r)
		y.push_back(miny + r * res + res * 0.5);

	std::cout << bv.evaluate(x, y, z) << "\n";

	std::ofstream out("/home/rob/Desktop/ec/tmp/out.csv");
	size_t ys = y.size();
	size_t xs = x.size();
	std::vector<double> grid(x.size() * y.size());
	for(size_t j = 0; j < y.size(); ++j) {
		for(size_t i = 0; i < x.size(); ++i) {
			out << x[i] << "," << y[j] << "," << z[i *  xs + j] << "\n";
			grid[(ys - j - 1) * xs + i] = z[i * ys + j];
		}
	}

	geo::util::saveGrid("/home/rob/Desktop/ec/tmp/out3.tif", grid, cols, rows, minx, maxy, 0.5, -0.5, "PROJCS[\"WGS 84 / UTM zone 10N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.01745329251994328,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",-123],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],AUTHORITY[\"EPSG\",\"32610\"],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH]]");
}
