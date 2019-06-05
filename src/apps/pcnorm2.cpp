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

	double x, y, z;
	rdr.Reset();
	while(rdr.ReadNextPoint()) {
		const liblas::Point& pt = rdr.GetPoint();
		x = pt.GetX();
		y = pt.GetY();
		z = pt.GetZ();
		if(x < bounds[0]) bounds[0] = x;
		if(x > bounds[1]) bounds[1] = x;
		if(y < bounds[2]) bounds[2] = y;
		if(y > bounds[3]) bounds[3] = y;
		if(z < bounds[4]) bounds[4] = z;
		if(z > bounds[5]) bounds[5] = z;
	}
}

void buildGrid2(liblas::Reader& rdr, double* bounds, double resolution, int& cols, int& rows, std::vector<float>& grid) {

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

	double px, py, pz, w, cx, cy;
	int col, row;

	rdr.Reset();
	while(rdr.ReadNextPoint()) {
		const liblas::Point& pt = rdr.GetPoint();
		if(pt.GetClassification().GetClass() != 2)
			continue;
		px = pt.GetX();
		py = pt.GetY();
		pz = pt.GetZ();
		col = (int) ((px - bounds[0]) / resolution);
		row = (int) ((py - bounds[2]) / resolution);
		for(int r = row - 1; r < row + 2; ++r) {
			for(int c = col - 1; c < col + 2; ++c) {
				if(c < 0 || r < 0 || c >= cols || r >= rows) continue;
				cx = bounds[0] + c * resolution + resolution * 0.5;
				cy = bounds[2] + r * resolution + resolution * 0.5;
				w = (std::pow(cx - px, 2.0) + std::pow(cy - py, 2.0)) / (resolution * resolution);
				if(w > 1) continue;
				grid[r * cols + c] += pz * w;
				weights[r * cols + c] += w;

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

	double z, zn, px, py, cx, cy, w, d, gz;
	int col, row;
	rdr.Reset();
	while(rdr.ReadNextPoint()) {
		const liblas::Point& pt = rdr.GetPoint();
		px = pt.GetX();
		py = pt.GetY();
		col = (int) ((px - bounds[0]) / resolution);
		row = (int) ((py - bounds[2]) / resolution);
		cx = bounds[0] + col * resolution + resolution * 0.5;
		cy = bounds[2] + row * resolution + resolution * 0.5;
		if(px == cx) {
			col0 = col;
		} else if(px < cx) {
			col0 = std::max(0, col - 1);
		} else if(px > cx) {
			col0 = std::min(cols - 1, col + 1);
		}
		if(py == cy) {
			row0 = row;
		} else if(py < cy) {
			row0 = std::max(0, row - 1);
		} else {if(py > cy) {
			row0 = std::min(rows - 1, row + 1);
		}
		cx0 = bounds[0] + col0 * resolution + resolution * 0.5;
		cy0 = bounds[2] + row0 * resolution + resolution * 0.5;
		z0 = grid[row * cols + col];
		z1 = grid[row0 * cols + col];
		z2 = grid[row0 * cols + col0];
		1 -> cx, cy
		2 -> cx0, cy
		3 -> cx, cy0
		w0 = (cy - cy0) * (px - cx) + (cx - )
		d0 = std::sqrt(std::pow(cx - px, 2.0) + std::pow(cy - py, 2.0));
		d1 = std::sqrt(std::pow(cx - px, 2.0) + std::pow(cy0 - py, 2.0));
		d2 = std::sqrt(std::pow(cx0 - px, 2.0) + std::pow(cy0 - py, 2.0));
		zn =
		d = std::sqrt(std::pow(cx - px, 2.0) + std::pow(cy - py, 2.0));
				if(d == 0) {
					zn = gz;
					w = 0;
				} else {
					zn += gz / d;
					w += 1.0 / d;
				}
			}
		}
		if(w > 0)
			zn /= w;
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
