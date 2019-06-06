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

using namespace geo::raster;
using namespace geo::pc;
using namespace geo::pc::compute;

namespace {

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
		} else if(name == "hlrg-bio") {				return new HLRGBiometricsComputer(20, 75);
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
					finals.push_back(std::move(finalizeQ->front()));
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
			wcond->notify_all();
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
					writes.push_back(std::move(writeQ->front()));
					writeQ->pop();
				}
				fcond->notify_all();
			}
			for(RWrite& w : writes) {
				int i = 0;
				for(double v : w.values)
					rasters->at(i++)->setFloat(w.col, w.row, v, 1);
			}
			writes.clear();
		}
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

}

Rasterizer::Rasterizer(const std::vector<std::string> filenames) :
	m_filter(nullptr),
	m_thin(0) {
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

void Rasterizer::setFilter(PCPointFilter* filter) {
	//if(m_filter)
	//	delete m_filter;
	m_filter = filter;
}

PCPointFilter* Rasterizer::filter() const {
	return m_filter;
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
	for(const std::string& name : types) {
		std::unique_ptr<geo::pc::Computer> comp(getComputer(name));
		comp->setRasterizer(this);
		computers.push_back(std::move(comp));
	}

	g_trace("Checking file bounds");
	double bounds[4] = {G_DBL_MAX_POS, G_DBL_MAX_POS, G_DBL_MIN_POS, G_DBL_MIN_POS};
	{
		double fBounds[6];
		for(PCFile& f: m_files) {
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

	rasters[0]->fillFloat(0, 1);
	for(int i= 1; i < bandCount; ++i)
		rasters[i]->fillFloat(NODATA, 1);

	std::unordered_map<size_t, int> counts;
	{
		double fBounds[6];
		int radCols = std::ceil(radius / std::abs(resX));
		for(PCFile& f: m_files) {
			f.fileBounds(fBounds);
			int c0 = std::min(props.toCol(fBounds[0]), props.toCol(fBounds[2])) - radCols;
			int r0 = std::min(props.toRow(fBounds[1]), props.toRow(fBounds[3])) - radCols;
			int c1 = std::max(props.toCol(fBounds[0]), props.toCol(fBounds[2])) + radCols;
			int r1 = std::max(props.toRow(fBounds[1]), props.toRow(fBounds[3])) + radCols;
			//g_debug("box " << c0 << ", " << r0 << ", " << c1 << ", " << r1 << ", " << radCols)
			for(int r = r0; r <= r1; ++r) {
				if(r < 0 || r >= rows) continue;
				for(int c = c0; c <= c1; ++c) {
					if(c < 0 || c >= cols) continue;
					counts[r * cols + c]++;
				}
			}
		}
	}

	liblas::ReaderFactory fact;
	CountComputer countComp;

	// The squared radius for comparison.
	double rad0 = radius * radius;

	// The radius of the "box" of pixels to check for a claim on the current point.
	int radpx = (int) std::ceil(radius / std::abs(props.resolutionX()));
	int i = 0;

	std::unordered_map<size_t, RCell> cells;
	std::queue<RFinalize> finalizeQ;
	std::queue<RWrite> writeQ;
	std::mutex fmtx;
	std::mutex wmtx;
	std::condition_variable fcond;
	std::condition_variable wcond;
	bool frunning = true;
	bool wrunning = true;

	std::thread processT(rprocess, &finalizeQ, &writeQ, &fcond, &wcond, &fmtx, &wmtx, &frunning);
	std::thread writeT(rwrite, &writeQ, &rasters, &wcond, &fcond, &wmtx, &wrunning);

	for(PCFile& file : m_files) {
		while(finalizeQ.size() > 5000 || writeQ.size() > 5000)
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
		g_debug("Reading file " << i++ << " of " << m_files.size());
		std::unordered_set<size_t> indexes;
		for(const std::string& filename : file.filenames()) {
			std::ifstream str(filename);
			liblas::Reader rdr = fact.CreateWithStream(str);
			while(rdr.ReadNextPoint()) {
				const geo::pc::Point pt(rdr.GetPoint());
				double x = pt.x();
				double y = pt.y();
				int col = props.toCol(x);
				int row = props.toRow(y);
				if(radius == 0) {
					size_t idx = row * cols + col;
					cells[idx].values.push_back(pt);
					indexes.insert(idx);
				} else {
					//g_trace(g_max(0, col - radpx) << ", " << g_min(cols, col + radpx + 1) << ", " << g_max(0, row - radpx) << ", " << g_min(rows, row + radpx + 1))
					for(int r = g_max(0, row - radpx); r < g_min(rows, row + radpx + 1); ++r) {
						for(int c = g_max(0, col - radpx); c < g_min(cols, col + radpx + 1); ++c) {
							double cx = props.toCentroidX(c);
							double cy = props.toCentroidY(r);
							double dist = std::pow(cx - x, 2.0) + std::pow(cy - y, 2.0);
							if(dist <= rad0) {
								size_t idx = r * cols + c;
								cells[idx].values.push_back(pt);
								indexes.insert(idx);
							}
						}
					}
				}
			}
		}

		for(size_t idx : indexes) {
			if(counts.find(idx) == counts.end())
				g_runerr("Col " << (idx % cols) << ", row " << (idx / cols) << " has already been processed.")
			if(--counts[idx] == 0) {
				int col = idx % cols;
				int row = idx / cols;
				double x = props.toCentroidX(col);
				double y = props.toCentroidY(row);
				finalizeQ.emplace(col, row, x, y, radius, cells[idx], m_filter, &computers);
				counts.erase(idx);
				cells.erase(idx);
			}
		}
		fcond.notify_all();
	}

	g_debug("Finalizing the rest " << cells.size())
	for(auto& item : cells) {
		size_t idx = item.first;
		int col = idx % cols;
		int row = idx / cols;
		double x = props.toCentroidX(col);
		double y = props.toCentroidY(row);
		finalizeQ.emplace(col, row, x, y, radius, cells[idx], m_filter, &computers);
	}
	fcond.notify_all();
	counts.clear();
	cells.clear();

	g_debug("Shutting down threads")
	frunning = false;
	fcond.notify_all();
	processT.join();
	g_debug("Finalizer done")

	wrunning = false;
	wcond.notify_all();
	writeT.join();
	g_debug("Writer done")

	Raster outrast(filename, props);
	for(int band = 1; band <= bandCount; ++band)
		rasters[band - 1]->writeTo(outrast, props.cols(), props.rows(), 0, 0, 0, 0, 1, band);

	g_debug("Done")
}

