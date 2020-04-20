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
#include "grid.hpp"
#include "util.hpp"

using namespace geo::grid;
using namespace geo::pc;
using namespace geo::ds;
using namespace geo::pc::compute;

const std::unordered_map<std::string, std::string> computerNames = {
		{"min", "The minimum value"},
		{"max", "The maximum value"},
		{"percentile-n", "The percentile"},
		{"decile-n", "The decile"},
		{"quartile-n", "The quartile"},
		{"median", "The median value"},
		{"mean", "The mean value"},
		{"variance", "The variance with n-1"},
		{"std-dev", "The standard deviation with n-1"},
		{"rugosity-acr", "The arc-chord rugosity (DuPreez, 2004)"},
		{"idw-2", "Inverse distance weighting; coefficient 2"},
		{"hlrg-bio", "HLRG biometrics set"}
};

Computer* getComputer(const std::string& name) {
	if(name.find("percentile") == 0) {
		size_t len = std::string("percentile-").size();
		int p  = std::stoi(name.substr(len));
		std::cout << "Percentile: " << p << "\n";
		return new PercentileComputer((double) p / 100.0);
	} else if(name.find("quartile") == 0) {
		size_t len = std::string("percentile-").size();
		int p = std::stoi(name.substr(len));
		std::cout << "Decile: " << p << "\n";
		return new PercentileComputer((double) (p * 25) / 100.0);
	} else if(name.find("decile") == 0) {
		size_t len = std::string("decile-").size();
		int p = std::stoi(name.substr(len));
		std::cout << "Decile: " << p << "\n";
		return new PercentileComputer((double) (p * 10) / 100.0);
	} else {
		if(name == "min") { 				return new MinComputer();
		} else if(name == "max") { 			return new MaxComputer();
		} else if(name == "median") { 			return new PercentileComputer(0.5);
		} else if(name == "mean") { 			return new MeanComputer();
		} else if(name == "variance") { 		return new VarianceComputer();
		} else if(name == "std-dev") { 			return new StdDevComputer();
		} else if(name == "rugosity-acr") { 		return new RugosityComputer();
		} else if(name == "idw-2") {			return new IDWComputer();
		} else if(name == "hlrg-bio") {			return new HLRGBiometricsComputer(20, 75);
		}
	}
	g_runerr("Unknown computer name (" << name << ")");
}


Rasterizer::Rasterizer(const std::vector<std::string> filenames) :
	m_filter(nullptr),
	m_thin(0),
	m_nodata(-9999),
	m_prefilter(false) {

	for(size_t i = 0; i < 4; ++i)
		m_bounds[i] = std::nan("");

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

void fixBounds(double* bounds, double resX, double resY, double& easting, double& northing, double* rasterBounds) {

	double rx = std::abs(resX);
	double ry = std::abs(resY);

	// If the raster bounds are given, use 'em.
	if(!std::isnan(rasterBounds[0])) {
		g_debug("Using given raster bounds");
		for(size_t i = 0; i < 4; ++i)
			bounds[i] = rasterBounds[i];
	}

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
			while(easting < xmin)
				easting += rx;
			while(easting > xmin)
				easting -= rx;
			xmin = easting;
			xmax += rx;
		} else {
			while(easting > xmax)
				easting -= rx;
			while(easting < xmax)
				easting += rx;
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
			while(northing < ymin)
				northing += ry;
			while(northing > ymin)
				northing -= ry;
			ymin = northing;
			ymax += ry;
		} else {
			while(northing > ymax)
				northing -= ry;
			while(northing < ymax)
				northing += ry;
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

void Rasterizer::setFilter(PCPointFilter* filter) {
	m_filter = filter;
}

PCPointFilter* Rasterizer::filter() const {
	return m_filter;
}

void Rasterizer::setPrefilter(bool prefilter) {
	m_prefilter = prefilter;
}

void Rasterizer::setBounds(double* bounds) {
	for(size_t i = 0; i < 4; ++i)
		m_bounds[i] = bounds[i];
}

Rasterizer::~Rasterizer() {
}


void Rasterizer::rasterize(const std::string& filename, const std::vector<std::string>& types,
		double resX, double resY, double easting, double northing, double radius, 
		const std::string& projection, bool useHeader, bool voids, double maxRadius) {

	if(std::isnan(resX) || std::isnan(resY))
		g_runerr("Resolution not valid");

	double initrad = radius;
	radius = -1;

	if(radius < 0 || std::isnan(radius)) {
		radius = std::sqrt(std::pow(resX / 2, 2) * 2);
		g_warn("Invalid radius; using " << radius);
	}

	if(types.empty())
		g_argerr("No methods given.");

	// Configure the computers.
	m_computers.clear();
	for(const std::string& name : types) {
		std::unique_ptr<Computer> comp(getComputer(name));
		comp->setFilters(m_filter->computerFilters(name));
		comp->setRasterizer(this);
		m_computers.push_back(std::move(comp));
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

	// "Fix" the bounds so they align with the resolution, and are oriented correctly.
	// Constrain to the given bounds if set.
	g_trace("Fixing bounds ")
	fixBounds(bounds, resX, resY, easting, northing, m_bounds);
	g_trace(" bounds: " << bounds[0] << ", " << bounds[1] << "; " << bounds[2] << ", " << bounds[3])

	mqtree<geo::pc::Point> tree(std::min(bounds[0], bounds[2]), std::min(bounds[1], bounds[3]),
			std::max(bounds[0], bounds[2]), std::max(bounds[1], bounds[3]), 100);

	// Compute the grid dimensions.
	int cols = (int) std::ceil((bounds[2] - bounds[0]) / resX);
	int rows = (int) std::ceil((bounds[3] - bounds[1]) / resY);
	g_trace(" cols: " << cols << ", rows: " << rows)

	// Work out the number of bands; each computer knows how many bands it will produce.
	int bandCount = 1;
	std::vector<std::string> bandMeta = {"count"};
	for(const std::unique_ptr<Computer>& comp : m_computers) {
		bandCount += comp->bandCount();
		std::vector<std::string> meta = comp->bandMeta();
		bandMeta.insert(bandMeta.end(), meta.begin(), meta.end());
	}
	g_trace(" bands: " << bandCount)

	// Configure the raster properties.
	g_trace("Preparing raster")
	GridProps props;
	props.setTrans(easting, resX, northing, resY);
	props.setSize(cols, rows);
	props.setDataType(DataType::Float32);
	props.setProjection(projection);
	props.setWritable(true);
	props.setBands(bandCount);
	props.setBandMetadata(bandMeta);
	props.setNoData(m_nodata);

	// Add the points to the qtree. Note the scale must be chosen carefully.
	g_trace("Adding files to tree");
	{
		size_t pts = 0;
		geo::pc::Point pt;
		for(PCFile& f: m_files) {
			for(const std::string& fn : f.filenames())
				filenames.push_back(fn);
			f.init(useHeader);
			while(f.next(pt)) {
				if(!m_prefilter || m_filter->keep(pt)) {
					tree.add(pt);
					if(++pts % 10000000 == 0)
						std::cout << pts << "pts\n";
				}
			}
		}
	}
	g_trace(" " << tree.size() << " points added to tree.");

	// Prepare the final output raster.
	Grid<float> outrast(filename, props);

	// The number of cells found to be void.
	size_t voidCount = 0;

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

		// Prepare insertion iterators.
		auto citer = std::back_inserter(cpts);
		auto fiter = std::back_inserter(filtered);

		// Fill the raster with nodata first.
		outrast.fill(props.nodata());

		std::vector<geo::pc::Point> tmp;
		double aresX = std::abs(resX) / 2;
		double aresY = std::abs(resY) / 2;

		std::cout << std::setprecision(2) << std::fixed;

		std::cout << "Rows: " << rows << "\n";
		for(int r = 0; r < rows; ++r) {
			if(r % 100 == 0)
				std::cout << "Row " << r << " of " << rows << " ("<< ((float) r / rows * 100) << "%)\n";
			for(int c = 0; c < cols; ++c) {

				// Prepare a query point based on the grid location.
				spt.x(props.toX(c));
				spt.y(props.toY(r));

				// Search for points within the radius of the cell centre.
				if(tree.search(spt, radius, citer)) {

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

					outrast.set(c, r, count, band++);

					if(count) {
						for(size_t i = 0; i < m_computers.size(); ++i) {
							out.clear();
							// Calculate the values in the computer and append to the raster.
							m_computers[i]->compute(spt.x(), spt.y(), cpts, filtered, radius, out);
							for(double val : out) {
								if(!std::isnan(val)) {
									outrast.set(c, r, val, band++);
								} else {
									outrast.set(c, r, m_nodata, band++);
								}
							}
						}
					} else {
						++voidCount;
					}

					filtered.clear();
					out.clear();
					cpts.clear();
					dist.clear();
					pts.clear();
				} else {
					++voidCount;
				}
			}
			//tree.printStats();
		}
	}

	if(voids && voidCount) {

		GridProps maskProps(props);
		maskProps.setBands(1);
		maskProps.setDataType(DataType::Byte);
		maskProps.setWritable(true);
		Grid<char> mask(maskProps);
		mask.fill(0);

		int fx0, fx1, fy0, fy1, fa;

		// Make a mask based on the count layer.
		TargetFillOperator<float, char> op1(&outrast, 0, &mask, 0, props.nodata(), 1);
		TargetFillOperator<char, char> op2(&mask, 0, &mask, 0, 1, 2);

		// Skip the first band which is always count.
		for(int b = 1; b < props.bands(); ++b) {
			std::cout << "Filling voids in band " << b << "\n";
			for(int r = 0; r < rows; ++r) {
				if(r % 10 == 0)
					std::cout << "Row " << r << " of " << rows << " (" << (float) r / rows * 100 << "%)\n";
				for(int c = 0; c < cols; ++c) {

					// The value is good. Skip.
					if(outrast.get(c, r, b) != props.nodata())
						continue;

					// The mask says this is an edge-connected pixel. Skip.
					if(mask.get(c, r, 0) == 2)
						continue;

					// Flood to find edge-connected pixels.
					Grid<float>::floodFill(c, r, op1, false, &fx0, &fy0, &fx1, &fy1, &fa);

					// If the fill touches an edge, set it to 2 and ignore it.
					if(fx0 == 0 || fx1 == cols - 1 || fy0 == 0 || fy1 == rows - 1) {
						Grid<char>::floodFill(c, r, op2, false, &fx0, &fy0, &fx1, &fy1, &fa);
						continue;
					}

					// Fill the pixel.
					int rad = 0;
					double v0, v = 0, w0, w = 0;
					bool found = false;
					do {
						++rad;
						for(int rr = r - rad; rr < r + rad + 1; ++rr) {
							for(int cc = c - rad; cc < c + rad + 1; ++cc) {
								if(cc < 0 || cc >= cols || rr < 0 || rr >= rows || (cc == c && rr == r))
									continue;
								if((v0 = outrast.get(cc, rr, b)) != props.nodata()) {
									w0 = 1.0 / (std::pow(props.toX(c) - props.toX(cc), 2.0) + std::pow(props.toY(r) - props.toY(rr), 2.0));	// not possible for d to be zero.
									w += w0;
									v += w0 * v0;
									found = true;
								}
							}
						}
					} while(!found && rad * std::abs(resY) < maxRadius);

					if(found)
						outrast.set(c, r, v / w, b);

				}
			}
		}
	}

	g_debug("Done")
}

