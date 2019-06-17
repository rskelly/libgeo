/*
 * pc_rasterizer.cpp
 *
 *  Created on: Mar 17, 2018
 *      Author: rob
 */

#include <grid.hpp>
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
#include "pointcloud.hpp"
#include "pc_computer.hpp"
#include "pc_index.hpp"

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


// 1) Create byte grid (map) corresponding to output raster.
// 2) Iterate over bounds. For each one, increment the cells in the grid that are covered by the bounds plus the radius.
// 4) Iterate over bounds.
//	-- Create a cell if required and add values to it.
//  -- Decrement the count for that cell in the grid.
//  -- If the count is zero, finalize the cell.

class RCell {
public:
	std::vector<geo::pc::Point> values;
};

class RWrite {
public:
	int col;
	int row;
	std::vector<double> values;
	RWrite(int col, int row, const std::vector<double>& values) :
		col(col), row(row),
		values(values) {
	}
};

class RFinalize {
public:
	int col, row;
	double x, y;
	double radius;
	RCell cell;
	PCPointFilter* filter;
	std::vector<std::unique_ptr<Computer> >* computers;

	RFinalize(int col, int row, double x, double y, double radius, const RCell& cell,
			PCPointFilter* filter, std::vector<std::unique_ptr<Computer> >* computers) :
		col(col), row(row),
		x(x), y(y), radius(radius),
		cell(cell),
		filter(filter),
		computers(computers) {
	}

	void finalize(std::vector<double>& write) {
		//g_debug(idx << ", " << counts.count(idx) << ", " << counts[idx])
		std::vector<geo::pc::Point> filtered;
		std::vector<double> out;
		size_t count = cell.values.size();
		if(count)
			count = filter->filter(cell.values.begin(), cell.values.end(), std::back_inserter(filtered));
		write.push_back(count);
		if(count) {
			for(size_t i = 0; i < computers->size(); ++i) {
				(*computers)[i]->compute(x, y, cell.values, filtered, radius, out);
				for(double val : out)
					write.push_back(std::isnan(val) ? NODATA : val);
				out.clear();
			}
		}
	}
};

void rprocess(std::queue<RFinalize>* finalizeQ, std::queue<RWrite>* writeQ,
		std::condition_variable* fcond, std::condition_variable* wcond,
		std::mutex* fmtx, std::mutex* wmtx, bool* running) {
	std::vector<double> write;
	std::vector<RFinalize> finals;
	while(*running || !finalizeQ->empty()) {
		{
			std::unique_lock<std::mutex> lk(*fmtx);
			while((*running && finalizeQ->empty()) || writeQ->size() > 5000)
				fcond->wait(lk);
			while(!finalizeQ->empty()) {
				finals.emplace_back(finalizeQ->front());
				finalizeQ->pop();
			}
		}
		{
			std::unique_lock<std::mutex> lk(*wmtx);
			for(RFinalize& f : finals) {
				f.finalize(write);
				writeQ->emplace(f.col, f.row, write);
				write.clear();
			}
		}
		finals.clear();
		wcond->notify_one();
	}
}

void rwrite(std::queue<RWrite>* writeQ, std::vector<std::unique_ptr<MemRaster> >* rasters,
	std::condition_variable* wcond, std::condition_variable* fcond, std::mutex* wmtx, bool* running) {
	std::vector<RWrite> writes;
	while(*running || !writeQ->empty()) {
		{
			std::unique_lock<std::mutex> lk(*wmtx);
			while(*running && writeQ->empty()) 
				wcond->wait(lk);
			while(!writeQ->empty()) {
				writes.emplace_back(writeQ->front());
				writeQ->pop();
			}
			fcond->notify_all();
		}
		for(RWrite& w : writes) {
			int i = 0;
			for(double v : w.values)
				rasters->at(i++)->setFloat(w.col, w.row, v);
		}
		writes.clear();
	}
}

Rasterizer::Rasterizer(const std::vector<std::string> filenames) :
	m_filter(nullptr),
	m_files(filenames) {
	//for(const std::string& filename : filenames)
	//	m_files.emplace_back(filename);
}

const std::unordered_map<std::string, std::string>& Rasterizer::availableComputers() {
	return computerNames;
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
		double resX, double resY, double easting, double northing, double radius,
		int srid, int memory, bool useHeader) {

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

	PointCloud pc(m_files);
	pc.init(radius, resX, resY, &easting, &northing, useHeader);

	int bandCount = 1;
	for(const std::unique_ptr<Computer>& comp : computers)
		bandCount += comp->bandCount();
	g_trace(" bands: " << bandCount)

	g_trace("Preparing raster")
	GridProps props;
	props.setTrans(easting, resX, northing, resY);
	props.setSize(pc.cols(), pc.rows());
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

	geo::pc::Point pt;
	std::vector<size_t> indices;
	std::vector<size_t> finals;
	std::unordered_map<size_t, RCell> cells;
	std::vector<double> output;

	while(pc.next(pt, finals)) {
		double x, y;
		int col, row;
		indices.clear();
		pc.getIndices(pt, indices);
		for(size_t idx : indices)
			cells[idx].values.push_back(pt);
		for(size_t idx : finals) {
			const RCell& cell = cells[idx];
			pc.fromIndex(idx, col, row);
			pc.fromIndex(idx, x, y);
			RFinalize rf(col, row, x, y, radius, cell, m_filter, &computers);
			rf.finalize(output);
			cells.erase(idx);
			output.clear();
		}

	}

	Raster outrast(filename, props);
	for(int band = 1; band <= bandCount; ++band)
		rasters[band - 1]->writeTo(outrast, props.cols(), props.rows(), 0, 0, 0, 0, 1, band);

	g_debug("Done")
}

