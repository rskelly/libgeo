/*
 * pc_rasterizer.cpp
 *
 *  Created on: Mar 17, 2018
 *      Author: rob
 */

#include <memory>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdint>
#include <list>
#include <unordered_map>
#include <condition_variable>
#include <thread>

#include "util.hpp"
#include "raster.hpp"
#include "pointcloud.hpp"
#include "pc_computer.hpp"
#include "externalmergesort.hpp"

using namespace geo::raster;
using namespace geo::pc;
using namespace geo::pc::compute;

const std::unordered_map<std::string, std::string> computerNames = {
		{"min", "The minimum value"},
		{"max", "The maximum value"},
		{"percentile-5", "The 5th percentile"},
		{"decile-1", "The 1st decile"},
		{"decile-2", "The 2nd decile"},
		{"quartile-1", "The 1st quartile"},
		{"decile-3", "The 3rd decile"},
		{"decile-4", "The 4th decile"},
		{"decile-5", "The 4th decile"},
		{"quartile-2", "The 2nd quartile"},
		{"median", "The median value"},
		{"decile-6", "The 6th decile"},
		{"decile-7", "The 7th decile"},
		{"quantile-3", "The 3rd quantile"},
		{"decile-8", "The 8th decile"},
		{"decile-9", "The 9th decile"},
		{"percentile-95", "The 95th percentile"},
		{"mean", "The mean value"},
		{"variance", "The variance with n-1"},
		{"std-dev", "The standard deviation with n-1"},
		{"rugosity-acr", "The arc-chord rugosity (DuPreez, 2004)"},
		{"idw-2", "Inverse distance weighting; coefficient 2"},
		{"hlrg-bio", "HLRG biometrics set"}
};

Computer* getComputer(const std::string& name) {
	if(name == "min") { 						return new MinComputer();
	} else if(name == "max") { 					return new MaxComputer();
	} else if(name == "percentile-5") { 		return new PercentileComputer(0.05);
	} else if(name == "decile-1") { 			return new PercentileComputer(0.1);
	} else if(name == "decile-2") { 			return new PercentileComputer(0.2);
	} else if(name == "quartile-1") { 			return new PercentileComputer(0.25);
	} else if(name == "decile-3") { 			return new PercentileComputer(0.3);
	} else if(name == "decile-4") { 			return new PercentileComputer(0.4);
	} else if(name == "decile-5") { 			return new PercentileComputer(0.5);
	} else if(name == "quartile-2") { 			return new PercentileComputer(0.5);
	} else if(name == "median") { 				return new PercentileComputer(0.5);
	} else if(name == "decile-6") { 			return new PercentileComputer(0.6);
	} else if(name == "decile-7") { 			return new PercentileComputer(0.7);
	} else if(name == "quartile-3") { 			return new PercentileComputer(0.75);
	} else if(name == "decile-8") { 			return new PercentileComputer(0.8);
	} else if(name == "decile-9") { 			return new PercentileComputer(0.9);
	} else if(name == "percentile-95") { 		return new PercentileComputer(0.95);
	} else if(name == "mean") { 				return new MeanComputer();
	} else if(name == "variance") { 			return new VarianceComputer();
	} else if(name == "std-dev") { 				return new StdDevComputer();
	} else if(name == "rugosity-acr") { 		return new RugosityComputer();
	} else if(name == "idw-2") {				return new IDWComputer();
	} else if(name == "hlrg-bio") {				return new HLRGBiometricsComputer(20, 75, 2.0);
	}
	g_runerr("Unknown computer name (" << name << ")");
}


Rasterizer::Rasterizer(const std::vector<std::string> filenames) :
	m_filter(nullptr) {
	for(const std::string& filename : filenames)
		m_files.emplace_back(filename);
}

const std::unordered_map<std::string, std::string>& Rasterizer::availableComputers() {
	return computerNames;
}

double Rasterizer::density(double resolution, double radius) {
	size_t count = 0;
	double w, h, cells, sum = 0;
	double fBounds[6];
	for(PCFile& f: m_files) {
		f.init();
		f.fileBounds(fBounds);
		w = fBounds[2] - fBounds[0];
		h = fBounds[3] - fBounds[1];
		cells = (w * h) / (resolution * resolution);
		if(cells > 0) {
			sum += f.pointCount() / cells;
			++count;
		}
	}

	double cell = radius > 0 ? (M_PI * radius * radius) / (resolution * resolution) : 1;
	return (sum / count) * 1.5 * cell;
}

void fixBounds(double* bounds, double resX, double resY, double& easting, double& northing) {

	double rx = std::abs(resX);
	double ry = std::abs(resY);

	double xmin = std::floor(std::min(bounds[0], bounds[2]) / rx) * rx;
	double ymin = std::floor(std::min(bounds[1], bounds[3]) / ry) * ry;
	double xmax = std::ceil(std::max(bounds[0], bounds[2]) / rx) * rx;
	double ymax = std::ceil(std::max(bounds[1], bounds[3]) / ry) * ry;

	if(std::isnan(easting)) {
		xmin -= rx;
		xmax += rx;
		easting = resX > 0 ? xmin : xmax;
	} else {
		if(resX > 0) {
			if(easting > xmin)
				g_argerr("Easting must be to the West of the boundary.")
			xmin = easting;
			xmax += rx;
		} else {
			if(easting < xmax)
				g_argerr("Easting must be to the West of the boundary.")
			xmax = easting;
			xmin -= rx;
		}
	}
	if(std::isnan(northing)) {
		ymin -= ry;
		ymax += ry;
		northing = resY > 0 ? ymin : ymax;
	} else {
		if(resY > 0) {
			if(northing > ymin)
				g_argerr("Easting must be to the North of the boundary.")
			ymin = northing;
			ymax += ry;
		} else {
			if(northing < ymax)
				g_argerr("Northing must be to the North of the boundary.")
			ymax = northing;
			ymin -= ry;
		}
	}

	bounds[resX > 0 ? 0 : 2] = xmin;
	bounds[resY > 0 ? 1 : 3] = ymin;
	bounds[resX > 0 ? 2 : 1] = xmax;
	bounds[resY > 0 ? 3 : 1] = ymax;

}

void Rasterizer::setFilter(PCPointFilter* filter) {
	//if(m_filter)
	//	delete m_filter;
	m_filter = filter;
}

Rasterizer::~Rasterizer() {
	//if(m_filter)
	//	delete m_filter;
}


void Rasterizer::rasterize(const std::string& filename, const std::vector<std::string>& _types,
		double resX, double resY, double easting, double northing, double radius, int srid, int memory, bool useHeader) {

	if(std::isnan(resX) || std::isnan(resY))
		g_runerr("Resolution not valid");

	if(radius < 0 || std::isnan(radius)) {
		radius = std::sqrt(std::pow(resX / 2, 2) * 2);
		g_warn("Invalid radius; using " << radius);
	}

	std::vector<std::string> types(_types);
	if(types.empty())
		g_argerr("No methods given; defaulting to mean");

	std::vector<std::unique_ptr<Computer> > computers;
	for(const std::string& name : types)
		computers.emplace_back(getComputer(name));

	g_trace("Checking file bounds");
	std::vector<std::string> filenames;
	double bounds[4] = {G_DBL_MAX_POS, G_DBL_MAX_POS, G_DBL_MIN_POS, G_DBL_MIN_POS};
	{
		double fBounds[6];
		for(PCFile& f: m_files) {
			for(const std::string& fn : f.filenames())
				filenames.push_back(fn);
			f.init(useHeader);
			f.fileBounds(fBounds);
			if(fBounds[0] < bounds[0]) bounds[0] = fBounds[0];
			if(fBounds[1] < bounds[1]) bounds[1] = fBounds[1];
			if(fBounds[2] > bounds[2]) bounds[2] = fBounds[2];
			if(fBounds[3] > bounds[3]) bounds[3] = fBounds[3];
		}
	}

	g_trace("Fixing bounds ")
	fixBounds(bounds, resX, resY, easting, northing);
	g_trace(" bounds: " << bounds[0] << ", " << bounds[1] << "; " << bounds[2] << ", " << bounds[3])

	int cols = (int) std::ceil((bounds[2] - bounds[0]) / resX);
	int rows = (int) std::ceil((bounds[3] - bounds[1]) / resY);
	g_trace(" cols: " << cols << ", rows: " << rows)

	int bandCount = 1;
	for(const std::unique_ptr<Computer>& comp : computers)
		bandCount += comp->bandCount();
	g_trace(" bands: " << bandCount)

	g_trace("Preparing raster")
	GridProps props;
	props.setTrans(easting, resX, northing, resY);
	props.setSize(cols, rows);
	props.setNoData(NODATA);
	props.setDataType(DataType::Float32);
	props.setSrid(srid);
	props.setWritable(true);
	props.setBands(bandCount);

	std::vector<std::unique_ptr<MemRaster> > rasters;
	for(int band = 1; band <= bandCount; ++band)
		rasters.emplace_back(new MemRaster(props, true));

	rasters[0]->fillFloat(0);
	for(int i= 1; i < bandCount; ++i)
		rasters[i]->fillFloat(NODATA);

	ExternalMergeSort es(filename, filenames, "/tmp", easting, northing, resX, resY);
	CountComputer countComp;

	es.sort(true);

	// The last finalized row.
	int lastR = 0;
	// The squared radius for distance comparisons without sqrt.
	double radiusSq = radius * radius;
	// Radius in pixel coords.
	int rad = (int) std::ceil(radius / std::max(std::abs(resX), std::abs(resY)));
	// Dictionary of cells/points.
	std::unordered_map<size_t, std::vector<geo::pc::Point> > cells;
	// List of populated cells.
	std::vector<bool> populated(cols * rows);
	// List of cells to remove from the dictionary.
	std::list<size_t> remove;

	geo::pc::Point pt;

	while(es.next(pt)) {

		double px = pt.x();
		double py = pt.y();
		int col = props.toCol(px);
		int row = props.toRow(py);
		size_t idx = row * cols + col;

		if(rad) {
			for(int r = std::max(0, row - rad); r < std::min(rows, row + rad + 1); ++r) {
				for(int c = std::max(0, col - rad); c < std::min(cols, col + rad + 1); ++c) {
					double cx = props.toCentroidX(c);
					double cy = props.toCentroidY(r);
					if(std::pow(px - cx, 2) + std::pow(py - cy, 2) <= radiusSq) {
						size_t i = r * cols + c;
						cells[i].push_back(std::move(pt));
						populated[i] = true;
					}
				}
			}
		} else {
			cells[idx].push_back(std::move(pt));
			populated[idx] = true;
		}

		if(row - rad > lastR) {
			for(int r = lastR; r < row - rad; ++r) {
				for(int c = 0; c < cols; ++c) {
					size_t i = r * cols + c;
					if(populated[i]) {
						const std::vector<geo::pc::Point>& pts = cells[i];
						finalize(c, r, props.toCentroidX(c), props.toCentroidY(r), radius, pts, rasters, computers);
						remove.push_back(i);
					}
				}
			}
			lastR = row - rad;
			for(const size_t& i : remove) {
				cells.erase(i);
				populated[i] = false;
			}
			remove.clear();
		}
	}

	g_debug("Finalizing the rest " << cells.size())
	for(auto& pair : cells) {
		int col = pair.first % cols;
		int row = pair.first / cols;
		finalize(col, row, props.toCentroidX(col), props.toCentroidY(row), radius, pair.second, rasters, computers);
	}

	Raster outrast(filename, props);
	for(int band = 1; band <= bandCount; ++band)
		rasters[band - 1]->writeTo(outrast, props.cols(), props.rows(), 0, 0, 0, 0, 1, band);

	g_debug("Done")
}

std::vector<geo::pc::Point> _filtered;
std::vector<double> _write;
std::vector<double> _out;

void Rasterizer::finalize(int col, int row, double x, double y, double radius,
		const std::vector<geo::pc::Point>& points, std::vector<std::unique_ptr<MemRaster> >& rasters,
		std::vector<std::unique_ptr<Computer> >& computers) {

	size_t count = points.size();
	if(count)
		count = m_filter->filter(points.begin(), points.end(), std::back_inserter(_filtered));
	_write.push_back(count);
	if(count) {
		for(size_t i = 0; i < computers.size(); ++i) {
			computers[i]->compute(x, y, points, _filtered, radius, _out);
			for(double val : _out)
				_write.push_back(std::isnan(val) ? NODATA : val);
		}
	} else {
		for(size_t i = 0; i < computers.size(); ++i) {
			for(int j = 0; j < computers[i]->bandCount(); ++j)
				_write.push_back(NODATA);
		}
	}

	for(size_t i = 0; i < _write.size(); ++i)
		rasters.at(i)->setFloat(col, row, _write[i]);

	_write.clear();
	_filtered.clear();
	_out.clear();
}
