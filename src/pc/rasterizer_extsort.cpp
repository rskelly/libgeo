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
	} else if(name == "min") { 					return new MaxComputer();
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

void fixBounds(double* bounds, double resX, double resY, double* easting, double* northing) {

	double aresX = std::abs(resX);
	double aresY = std::abs(resY);

	{
		int a = 0, b = 2;
		if(resX < 0)
			a = 2, b = 0;
		bounds[a] = std::floor(bounds[a] / resX) * resX;
		bounds[b] = std::ceil(bounds[b] / resX) * resX;
		a = 1, b = 3;
		if(resY < 0)
			a = 3, b = 1;
		bounds[a] = std::floor(bounds[a] / resY) * resY;
		bounds[b] = std::ceil(bounds[b] / resY) * resY;
	}

	if(!std::isnan(*easting)) {
		if((resX > 0 && *easting < bounds[0]) || (resX < 0 && *easting > bounds[2]))
			g_argerr("The easting is within the data boundary.");
		double w = bounds[2] - bounds[0];
		if(resX > 0) {
			while(*easting + w < bounds[2])
				w += resX;
			bounds[0] = *easting;
			bounds[2] = *easting + w;
		} else {
			while(*easting - w > bounds[1])
				w += aresX;
			bounds[2] = *easting;
			bounds[0] = *easting - w;
		}
	} else {
		*easting = bounds[resX > 0 ? 0 : 2];
	}

	if(!std::isnan(*northing)) {
		if((resY > 0 && *northing < bounds[1]) || (resY < 0 && *northing > bounds[3]))
			g_argerr("The *northing is within the data boundary.");
		double h = bounds[3] - bounds[1];
		if(resY > 0) {
			while(*northing + h < bounds[3])
				h += resY;
			bounds[1] = *northing;
			bounds[3] = *northing + h;
		} else {
			while(*northing - h > bounds[1])
				h += aresY;
			bounds[3] = *northing;
			bounds[1] = *northing - h;
		}
	} else {
		*northing = bounds[resY > 0 ? 1 : 3];
	}
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

	g_trace("Fixing bounds")
	fixBounds(bounds, resX, resY, &easting, &northing);
	g_trace(" bounds: " << bounds[0] << ", " << bounds[1] << "; " << bounds[2] << ", " << bounds[3])

	int cols = (int) ((bounds[2] - bounds[0]) / std::abs(resX)) + 1;
	int rows = (int) ((bounds[3] - bounds[1]) / std::abs(resY)) + 1;
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

	ExternalMergeSort es(filename, filenames, "/tmp");
	CountComputer countComp;

	es.sort(true);

	std::unordered_map<size_t, std::vector<geo::pc::Point> > cells;
	geo::pc::Point pt;

	while(es.next(pt)) {

		size_t index = es.index(pt);
		cells[index].push_back(std::move(pt));

		if(cells.size() > 1000) {
			for(auto& pair : cells) {
				double x, y;
				es.position(pair.first, x, y);
				int col = props.toCol(x);
				int row = props.toRow(y);
				x = props.toCentroidX(col);
				y = props.toCentroidY(row);
				std::vector<geo::pc::Point> points;
				bool found = true;
				for(double y0 = y - radius; y0 <= y + radius; y0 += 1.0) {
					for(double x0 = x - radius; x0 <= x + radius; x0 += 1.0) {
						size_t idx0 = es.index(x0, y0);
						if(cells.find(idx0) == cells.end()) {
							found = false;
							break;
						} else {
							points.insert(points.end(), cells[pair.first].begin(), cells[pair.first].end());
						}
					}
					if(!found)
						break;
				}
				if(found) {
					std::vector<geo::pc::Point> filtered;
					std::vector<double> write;
					std::vector<double> out;
					size_t count = points.size();
					if(count)
						count = m_filter->filter(points.begin(), points.end(), std::back_inserter(filtered));
					write.push_back(count);
					if(count) {
						for(size_t i = 0; i < computers.size(); ++i) {
							computers[i]->compute(x, y, points, filtered, radius, out);
							for(double val : out)
								write.push_back(std::isnan(val) ? NODATA : val);
						}
					}
					int i = 0;
					for(double v : write)
						rasters.at(i++)->setFloat(col, row, v);

					//std::cerr << "computing " << pair.first << "\n";
				}
			}

			// TODO: This removal regime doesn't work.
			double x, y;
			es.position(index, x, y);
			x = props.toCentroidX(props.toCol(x));
			y = props.toCentroidY(props.toRow(y));
			std::list<size_t> remove;
			for(auto& pair : cells) {
				double x0, y0;
				es.position(pair.first, x0, y0);
				x0 = props.toCentroidX(props.toCol(x0));
				y0 = props.toCentroidY(props.toRow(y0));
				if(std::pow(x - x0, 2) + std::pow(y - y0, 2) > radius * radius)
					remove.push_back(pair.first);
			}
			for(size_t i : remove)
				cells.erase(i);
			g_debug("removed " << remove.size())
		}
	}

	for(auto& pair : cells) {
		double x, y;
		es.position(pair.first, x, y);
		int col = props.toCol(x);
		int row = props.toRow(y);
		x = props.toCentroidX(col);
		y = props.toCentroidY(row);
		std::vector<geo::pc::Point> points;
		bool found = true;
		for(double y0 = y - radius; y0 <= y + radius; y0 += 1.0) {
			for(double x0 = x - radius; x0 <= x + radius; x0 += 1.0) {
				size_t idx0 = es.index(x0, y0);
				if(cells.find(idx0) != cells.end())
					points.insert(points.end(), cells[pair.first].begin(), cells[pair.first].end());
			}
		}
		std::vector<geo::pc::Point> filtered;
		std::vector<double> write;
		std::vector<double> out;
		size_t count = points.size();
		if(count)
			count = m_filter->filter(points.begin(), points.end(), std::back_inserter(filtered));
		write.push_back(count);
		if(count) {
			for(size_t i = 0; i < computers.size(); ++i) {
				computers[i]->compute(x, y, points, filtered, radius, out);
				for(double val : out)
					write.push_back(std::isnan(val) ? NODATA : val);
			}
		}
		int i = 0;
		for(double v : write)
			rasters.at(i++)->setFloat(col, row, v);
	}

	g_debug("Finalizing the rest " << cells.size())

	Raster outrast(filename, props);
	for(int band = 1; band <= bandCount; ++band)
		rasters[band - 1]->writeTo(outrast, props.cols(), props.rows(), 0, 0, 0, 0, 1, band);

	g_debug("Done")
}

