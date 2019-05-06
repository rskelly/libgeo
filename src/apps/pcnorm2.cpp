/*
 * pcnorm2.cpp
 *
 *  Created on: May 5, 2019
 *      Author: rob
 */

#include <limits>
#include <fstream>

#include <liblas/liblas.hpp>

void getBounds(liblas::Reader& rdr, double* bounds) {

	bounds[0] = bounds[2] = bounds[4] = std::numeric_limits<double>::max();
	bounds[1] = bounds[3] = bounds[5] = std::numeric_limits<double>::lowest();

	rdr.Reset();
	while(rdr.ReadNextPoint()) {
		const liblas::Point& pt = rdr.GetPoint();
		if(pt.GetX() < bounds[0]) bounds[0] = pt.GetX();
		if(pt.GetX() > bounds[1]) bounds[1] = pt.GetX();
		if(pt.GetY() < bounds[2]) bounds[2] = pt.GetY();
		if(pt.GetY() > bounds[3]) bounds[3] = pt.GetY();
		if(pt.GetZ() < bounds[4]) bounds[4] = pt.GetZ();
		if(pt.GetZ() > bounds[5]) bounds[5] = pt.GetZ();
	}
}

void buildGrid(liblas::Reader& rdr, double* bounds, double resolution, int& cols, int& rows, std::vector<float>& grid) {

	bounds[0] -= resolution;
	bounds[1] += resolution;
	bounds[2] -= resolution;
	bounds[3] += resolution;

	cols = (int) std::ceil((bounds[1] - bounds[0]) / resolution);
	rows = (int) std::ceil((bounds[3] - bounds[2]) / resolution);

	std::vector<float> weights(cols * rows);
	grid.resize(cols * rows);
	std::fill(grid.begin(), grid.end(), 0);
	std::fill(weights.begin(), weights.end(), 0);

	rdr.Reset();
	while(rdr.ReadNextPoint()) {
		const liblas::Point& pt = rdr.GetPoint();
		if(pt.GetClassification().GetClass() != 2)
			continue;
		double px = pt.GetX();
		double py = pt.GetY();
		int col = (int) (px - bounds[0]) / resolution;
		int row = (int) (py - bounds[2]) / resolution;
		double rad = resolution * resolution;
		for(int r = row - 2; r < row + 3; ++r) {
			for(int c = col - 2; c < col + 3; ++c) {
				if(c >= 0 && r >= 0 && c < cols && r < rows) {
					double x = bounds[0] + (c * resolution) + resolution * 0.5;
					double y = bounds[2] + (r * resolution) + resolution * 0.5;
					double d = std::pow(x - px, 2.0) + std::pow(y - py, 2.0);
					if(d > rad)
						continue;
					double w = 1.0 - d / rad;
					grid[r * cols + c] += pt.GetZ() * w;
					weights[r * cols + c] += w;
				}
			}
		}
	}

	for(size_t i = 0; i < grid.size(); ++i) {
		if(weights[i] != 0)
			grid[i] /= weights[i];
	}
}

void normalize(liblas::Reader& rdr, liblas::Writer& wtr, double* bounds, double resolution, int cols, int rows, std::vector<float>& grid) {

	bounds[4] = std::numeric_limits<double>::max();
	bounds[5] = std::numeric_limits<double>::lowest();

	int col, row, col0, row0;
	double x, y, z, x0, y0, z0, w0, z1, w1, z2, w2, zn;
	rdr.Reset();
	while(rdr.ReadNextPoint()) {
		const liblas::Point& pt = rdr.GetPoint();
		// Point's home cell.
		col = (int) (pt.GetX() - bounds[0]) / resolution;
		row = (int) (pt.GetY() - bounds[2]) / resolution;
		// Center of home cell.
		x = bounds[0] + (col * resolution) + resolution * 0.5;
		y = bounds[2] + (row * resolution) + resolution * 0.5;
		// Cell offsets.
		col0 = pt.GetX() < x ? col - 1 : col + 1;
		row0 = pt.GetY() < y ? row - 1 : row + 1;
		// Centers of offset cells.
		x0 = bounds[0] + (col0 * resolution) + resolution * 0.5;
		y0 = bounds[2] + (row0 * resolution) + resolution * 0.5;
		// Z and weight for each cell.
		z0 = grid[row * cols + col];
		w0 = (std::pow(pt.GetX() - x, 2.0) + std::pow(pt.GetY() - y, 2.0)) / resolution;
		z1 = grid[row0 * cols + col];
		w1 = (std::pow(pt.GetX() - x, 2.0) + std::pow(pt.GetY() - y0, 2.0)) / resolution;
		z2 = grid[row * cols + col0];
		w2 = (std::pow(pt.GetX() - x0, 2.0) + std::pow(pt.GetY() - y, 2.0)) / resolution;
		// Weighted mean at point.
		zn = (z0 * w0 + z1 * w1 + z2 * w2) / (w0 + w1 + w2);
		// Make point and write it.
		liblas::Point pt0(pt);
		pt0.SetZ((z = pt.GetZ() - zn));
		wtr.WritePoint(pt0);
		// Adjust bounds.
		if(z < bounds[4]) bounds[4] = z;
		if(z > bounds[5]) bounds[5] = z;
	}
}

int main(int argc, char** argv) {

	std::string infile = argv[1];
	std::string outfile = argv[2];
	double resolution = atof(argv[3]);

	std::ifstream input(infile);
	liblas::ReaderFactory rf;
	liblas::Reader rdr = rf.CreateWithStream(input);
	const liblas::Header& rhdr = rdr.GetHeader();

	liblas::Header whdr(rhdr);

	double bounds[6];
	getBounds(rdr, bounds);

	double gbounds[6];
	for(int i = 0; i < 6; ++i)
		gbounds[i] = bounds[i];

	int cols, rows;
	std::vector<float> grid;
	buildGrid(rdr, gbounds, resolution, cols, rows, grid);

	std::ofstream output;
	liblas::WriterFactory wf;
	liblas::Create(output, outfile);
	liblas::Writer wtr(output, whdr);

	normalize(rdr, wtr, gbounds, resolution, cols, rows, grid);

	whdr.SetMin(bounds[0], bounds[2], gbounds[4]);
	whdr.SetMax(bounds[1], bounds[3], gbounds[5]);

	wtr.SetHeader(whdr);

}
