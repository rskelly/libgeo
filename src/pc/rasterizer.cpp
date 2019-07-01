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
#include <algorithm>

#include "pointcloud.hpp"
#include "pc_computer.hpp"
#include "externalmergesort.hpp"
#include "grid.hpp"
#include "ds/kdtree.hpp"
#include "util.hpp"

using namespace geo::grid;
using namespace geo::pc;
using namespace geo::ds;
using namespace geo::pc::compute;
using namespace geo::pc::sort;

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
	} else if(name == "hlrg-bio") {				return new HLRGBiometricsComputer(20, 75);
	}
	g_runerr("Unknown computer name (" << name << ")");
}


Rasterizer::Rasterizer(const std::vector<std::string> filenames) :
	m_filter(nullptr),
	m_thin(0),
	m_nodata(-9999),
	m_limit(0) {
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
				g_argerr("Northing must be to the North of the boundary.")
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

void Rasterizer::setThin(int thin) {
	m_thin = thin;
}

void Rasterizer::setNoData(double nodata) {
	m_nodata = nodata;
}

void Rasterizer::setMemLimit(size_t limit) {
	m_limit = limit;
}

void Rasterizer::setFilter(PCPointFilter* filter) {
	m_filter = filter;
}

PCPointFilter* Rasterizer::filter() const {
	return m_filter;
}

Rasterizer::~Rasterizer() {
}


void Rasterizer::rasterize(const std::string& filename, const std::vector<std::string>& types,
		double resX, double resY, double easting, double northing, double radius, 
		int srid, bool useHeader, bool voids) {

	if(std::isnan(resX) || std::isnan(resY))
		g_runerr("Resolution not valid");

	double initrad = radius;
	radius = -1;

	if(radius < 0 || std::isnan(radius)) {
		radius = std::sqrt(std::pow(resX / 2, 2) * 2);
		g_warn("Invalid radius; using " << radius);
	}

	if(types.empty())
		g_argerr("No methods given; defaulting to mean");

	// Configure the computers.
	m_computers.clear();
	for(const std::string& name : types) {
		m_computers.emplace_back(getComputer(name));
		m_computers.back()->setRasterizer(this);
	}

	// Calculate the overall boundaries of the point cloud.
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

	mqtree<geo::pc::Point> tree(1);

	// "Fix" the bounds so they align with the resolution, and are oriented correctly.
	g_trace("Fixing bounds ")
	fixBounds(bounds, resX, resY, easting, northing);
	g_trace(" bounds: " << bounds[0] << ", " << bounds[1] << "; " << bounds[2] << ", " << bounds[3])

	// Compute the grid dimensions.
	int cols = (int) std::ceil((bounds[2] - bounds[0]) / resX);
	int rows = (int) std::ceil((bounds[3] - bounds[1]) / resY);
	g_trace(" cols: " << cols << ", rows: " << rows)

	// Work out the number of bands; each computer knows how many bands it will produce.
	int bandCount = 1;
	for(const std::unique_ptr<Computer>& comp : m_computers)
		bandCount += comp->bandCount();
	g_trace(" bands: " << bandCount)

	// Configure the raster properties.
	g_trace("Preparing raster")
	GridProps props;
	props.setTrans(easting, resX, northing, resY);
	props.setSize(cols, rows);
	props.setDataType(DataType::Float32);
	props.setSrid(srid);
	props.setWritable(true);
	props.setBands(bandCount);
	props.setNoData(m_nodata);

	// Add the points to the qtree. Note the scale must be chosen carefully.
	g_trace("Adding files to tree");
	{
		geo::pc::Point pt;
		for(PCFile& f: m_files) {
			for(const std::string& fn : f.filenames())
				filenames.push_back(fn);
			f.init(useHeader);
			while(f.next(pt)) {
				tree.add(pt);
			}
		}
	}
	g_trace(" " << tree.size() << " points added to tree.");
	tree.build();


	// Prepare the final output raster.
	Grid<double> outrast(filename, props);

	g_trace("Running...")
	{
		// The first layer is just the point count.
		CountComputer countComp;
		size_t count = 0;
		int band;

		// Temporary storage.
		std::vector<geo::pc::Point> pts;
		std::vector<geo::pc::Point> cpts;
		std::vector<geo::pc::Point> filtered;
		std::vector<double> out;
		std::vector<double> dist;
		geo::pc::Point spt;

		bool voidCells = false;

		// Prepare insertion iterators.
		auto citer = std::back_inserter(cpts);
		auto diter = std::back_inserter(dist);
		auto fiter = std::back_inserter(filtered);

		// Fill the raster with nodata first.
		outrast.fill(props.nodata());

		std::vector<geo::pc::Point> tmp;
		double aresX = std::abs(resX) / 2;
		double aresY = std::abs(resY) / 2;

		do {

			// If this is a void filling loop, increase the radius.
			if(voidCells)
				radius *= 2;

			for(int r = 1; r < rows - 1; ++r) {
				std::cout << "Row " << r << " of " << rows << "\n";
				for(int c = 1; c < cols - 1; ++c) {

					// Prepare a query point based on the grid location.
					spt.x(props.toX(c));
					spt.y(props.toY(r));

					// If this is a void filing loop, skip when the value is valid.
					if(voidCells && outrast.get(c, r, 1) != props.nodata())
						continue;

					// Search for points within the radius of the cell centre.
					if(tree.search(spt, radius, citer, diter)) {

						// If the radius was 0, filter the points to the cell bounds.
						if(initrad == 0) {
							tmp.swap(cpts);
							cpts.clear();
							for(geo::pc::Point& pt : tmp) {
								if(pt.x() < spt.x() - aresX || pt.x() > spt.x() + aresX || pt.y() < spt.y() - aresY || pt.y() > spt.y() + aresY)
									continue;
								cpts.push_back(std::move(pt));
							}
						}

						// Filter the points according to the configured filter.
						if(!cpts.empty()) {
							count = m_filter->filter(cpts.begin(), cpts.end(), fiter);

							if(m_thin > 0) {
								// Thin the points to the desired density, if required.
								if(filtered.size() < (size_t) m_thin) {
									filtered.clear();
									count = 0;
								} else {
									// Randomize, then take the first n points.
									std::random_shuffle(filtered.begin(), filtered.end());
									filtered.resize(m_thin);
									count = m_thin;
								}
							}
						}

						band = 0;

						// Write the point count.
						outrast.set(c, r, count, band++);

						if(count) {
							for(size_t i = 0; i < m_computers.size(); ++i) {
								out.clear();
								// Calculate the values in the computer and append to the raster.
								m_computers[i]->compute(spt.x(), spt.y(), cpts, filtered, radius, out);
								for(double val : out) {
									if(!std::isnan(val))
										outrast.set(c, r, val, band++);
								}
							}
						}

						filtered.clear();
						out.clear();
						cpts.clear();
						dist.clear();
						pts.clear();
					}
				}
			}

			// If void filling is enabled, check whether there are void cells and
			// enable the loop.
			if(voids) {
				voidCells = false;
				for(int r = 0; !voidCells && r < rows; ++r) {
					for(int c = 0; !voidCells && c < cols; ++c) {
						if(outrast.get(c, r, 1) == props.nodata())
							voidCells = true;
					}
				}
			}
		} while(voidCells);
	}

	g_debug("Done")
}

