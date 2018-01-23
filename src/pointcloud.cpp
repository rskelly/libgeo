#define _USE_MATH_DEFINES

#include <fstream>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <memory>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdint>
#include <list>

#include <math.h>

#include <CGAL/Plane_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <liblas/liblas.hpp>

#include "util.hpp"
#include "raster.hpp"
#include "pointcloud.hpp"
#include "ds/kdtree.hpp"


using namespace geo::raster;
using namespace geo::util;
using namespace geo::pc;

PCFile::PCFile(const std::string& filename, double x, double y) :
	m_x(x), m_y(y),
	m_bounds{99999999., 99999999., -999999999., -99999999., 99999999., -99999999.} {

	m_filenames.push_back(filename);
}

PCFile::PCFile(const std::vector<std::string>& filenames, double x, double y) :
	m_x(x), m_y(y),
	m_bounds{99999999., 99999999., -999999999., -99999999., 99999999., -99999999.},
	m_filenames(filenames) {
}

double PCFile::x() const {
	return m_x;
}

double PCFile::y() const {
	return m_y;
}

void PCFile::bounds(double* bounds) const {
	for(int i = 0; i < 6; ++i)
		bounds[i] = m_bounds[i];
}

const std::vector<std::string>& PCFile::filenames() const {
	return m_filenames;
}

void PCFile::init(bool useHeader) {
	liblas::ReaderFactory f;
	for(const std::string& filename : m_filenames) {
		std::ifstream str(filename, std::ios::in | std::ios::binary);
		liblas::Reader reader = f.CreateWithStream(str);
		if(useHeader) {
			const liblas::Header& hdr = reader.GetHeader();
			double minx = hdr.GetMinX();
			double miny = hdr.GetMinY();
			double minz = hdr.GetMinZ();
			double maxx = hdr.GetMaxX();
			double maxy = hdr.GetMaxY();
			double maxz = hdr.GetMaxZ();
			if(minx < m_bounds[0]) m_bounds[0] = minx;
			if(miny < m_bounds[1]) m_bounds[1] = miny;
			if(maxx > m_bounds[2]) m_bounds[2] = maxx;
			if(maxy > m_bounds[3]) m_bounds[3] = maxy;
			if(minz < m_bounds[4]) m_bounds[4] = minz;
			if(maxz > m_bounds[5]) m_bounds[5] = maxz;
		} else {
			while(reader.ReadNextPoint()) {
				const liblas::Point& pt = reader.GetPoint();
				double x = pt.GetX();
				double y = pt.GetY();
				double z = pt.GetZ();
				if(x < m_bounds[0]) m_bounds[0] = x;
				if(y < m_bounds[1]) m_bounds[1] = y;
				if(x > m_bounds[2]) m_bounds[2] = x;
				if(y > m_bounds[3]) m_bounds[3] = y;
				if(z < m_bounds[4]) m_bounds[4] = z;
				if(z > m_bounds[5]) m_bounds[5] = z;
			}
		}
	}
}

PCFile::~PCFile() {}


const long maxPoints = 20000000;

PCWriter::PCWriter(const std::string& filename, const liblas::Header& hdr, double x, double y) :
	m_fileIdx(0),
	m_returns(0), m_retNum{0,0,0,0,0},
	m_totalReturns(0),
	m_outBounds{99999999.0,99999999.0,-99999999.0,-99999999.0,99999999.0,-99999999.0},
	m_x(x), m_y(y),
	m_filename(filename),
	m_writer(nullptr), m_header(nullptr),
	m_dod(true) {

	m_header = new liblas::Header(hdr);

	open();
}

double PCWriter::x() const {
	return m_x;
}

double PCWriter::y() const {
	return m_y;
}

void PCWriter::bounds(double* bounds) const {
	for(int i = 0; i < 6; ++i)
		bounds[i] = m_outBounds[i];
}

const std::vector<std::string>& PCWriter::filenames() const {
	return m_filenames;
}

std::string PCWriter::nextFile() {
	std::stringstream ss;
	ss << m_filename << "_" << ++m_fileIdx << ".las";
	std::string filename = ss.str();
	m_filenames.push_back(filename);
	return filename;
}

void PCWriter::open() {
	close();
	std::string filename = nextFile();
	m_str.open(filename, std::ios::out | std::ios::binary);
	m_writer = new liblas::Writer(m_str, *m_header);
}

void PCWriter::deleteOnDestruct(bool dod) {
	m_dod = dod;
}

void PCWriter::close() {
	if(m_writer) {
		m_header->SetMin(m_outBounds[0], m_outBounds[2], m_outBounds[4]);
		m_header->SetMax(m_outBounds[1], m_outBounds[3], m_outBounds[5]);
		for(int i = 0; i < 5; ++i)
			m_header->SetPointRecordsByReturnCount(i, m_retNum[i]);
		m_header->SetPointRecordsCount(m_returns);
		m_writer->SetHeader(*m_header);
		m_writer->WriteHeader();
		m_str.close();
		delete m_writer;
		m_writer = nullptr;

		m_returns = 0;
		for(int i = 0; i < 5; ++i)
			m_retNum[i] = 0;
	}
}

void PCWriter::addPoint(const liblas::Point& pt) {
	if(m_returns >= maxPoints || !m_writer)
		open();

	double x = pt.GetX();
	double y = pt.GetY();
	double z = pt.GetZ();

	++m_returns;
	++m_totalReturns;
	m_retNum[pt.GetReturnNumber() - 1]++;

	if(x < m_outBounds[0]) m_outBounds[0] = x;
	if(x > m_outBounds[2]) m_outBounds[2] = x;
	if(y < m_outBounds[1]) m_outBounds[1] = y;
	if(y > m_outBounds[3]) m_outBounds[3] = y;
	if(z < m_outBounds[4]) m_outBounds[4] = z;
	if(z > m_outBounds[5]) m_outBounds[5] = z;

	m_writer->WritePoint(pt);
}

PCWriter::~PCWriter() {
	close();
	delete m_header;
	if(!m_totalReturns || m_dod) {
		for(const std::string& f : m_filenames)
			Util::rm(f);
	}
}

int even(int num) {
	if(num % 2 == 1)
		++num;
	return num;
}

geo::pc::Tile::Tile(double minx, double miny, double maxx, double maxy, double buffer) {
	m_bounds[0] = minx;
	m_bounds[1] = miny;
	m_bounds[2] = maxx;
	m_bounds[3] = maxy;
	m_bufferedBounds[0] = m_bounds[0] - buffer;
	m_bufferedBounds[1] = m_bounds[1] - buffer;
	m_bufferedBounds[2] = m_bounds[2] + buffer;
	m_bufferedBounds[3] = m_bounds[3] + buffer;
}

PCWriter* geo::pc::Tile::writer(bool release) {
	if(release) {
		return m_writer.release();
	} else {
		return m_writer.get();
	}
}

void geo::pc::Tile::writer(PCWriter* wtr) {
	m_writer.reset(wtr);
}

bool geo::pc::Tile::contains(double x, double y) {
	return x >= m_bounds[0] && x < m_bounds[2] 
		&& y >= m_bounds[1] && y < m_bounds[3];
}

bool geo::pc::Tile::containsBuffered(double x, double y) {
	return x >= m_bufferedBounds[0] && x < m_bufferedBounds[2] 
		&& y >= m_bufferedBounds[1] && y < m_bufferedBounds[3];
}


Tiler::Tiler(const std::vector<std::string> filenames) {
	for(const std::string& filename : filenames)
		files.emplace_back(filename);
}

void Tiler::tile(const std::string& outdir, double size, double buffer, int srid,
	double easting, double northing, int maxFileHandles) {

	std::cerr << std::setprecision(12);

	// Calculate the overall bounds of the file set.
	double allBounds[6] = {9999999999., 9999999999., -9999999999., -9999999999., 999999999., -999999999.};
	double fBounds[6];
	for(PCFile& f : files) {
		f.init();
		f.bounds(fBounds);
		std::cerr << "file bounds " << fBounds[0] << ", " << fBounds[1] << "; " << fBounds[2] << ", " << fBounds[3] << "\n";
		if(fBounds[0] < allBounds[0]) allBounds[0] = fBounds[0];
		if(fBounds[1] < allBounds[1]) allBounds[1] = fBounds[1];
		if(fBounds[2] > allBounds[2]) allBounds[2] = fBounds[2];
		if(fBounds[3] > allBounds[3]) allBounds[3] = fBounds[3];
		if(fBounds[4] < allBounds[4]) allBounds[4] = fBounds[4];
		if(fBounds[5] > allBounds[5]) allBounds[5] = fBounds[5];
	}

	// If the easting and northing aren't given, calculate as a
	// multiple of size.
	if(std::isnan(easting))
		easting = ((int) (allBounds[0] / size)) * size - size;
	if(std::isnan(northing))
		northing = ((int) (allBounds[1] / size)) * size - size;

	// Reset the bounds using easting, northing and multiples of size.
	allBounds[0] = easting;
	allBounds[1] = northing;
	allBounds[2] = std::ceil(allBounds[2] / size) * size + size;
	allBounds[3] = std::ceil(allBounds[3] / size) * size + size;

	// Compute the number of columns and rows of tiles.
	int cols = (int) (allBounds[2] - allBounds[0]) / size;
	int rows = (int) (allBounds[3] - allBounds[1]) / size;

	// Double the tile size until few enough writers are used.
	double size0 = size;
	int cols0 = cols;
	int rows0 = rows;
	while(cols0 * rows0 > maxFileHandles && cols0 > 1 && rows0 > 1) {
		size0 *= 2.0;
		cols0 = even(std::ceil(cols0 * 0.5));
		rows0 = even(std::ceil(rows0 * 0.5));
	}

	liblas::ReaderFactory rfact;
	std::vector<std::unique_ptr<Tile> > tiles;

	do {

		std::cerr << "cols " << cols0 << "; rows " << rows0 << "; size " << size0 << "\n";
		std::cerr << "bounds " << allBounds[0] << ", " << allBounds[1] << "; " << allBounds[2] << ", " << allBounds[3] << "\n";

		if(!tiles.empty()) {
			// This is run n>0, so we now read from the intermediate tiles.

			// To prevent the previous set of tiles from getting deleted before
			// we're done reading them.
			std::vector<std::unique_ptr<PCWriter> > tmpWriters;

			// Clear the readers list and rebuild with files from the previous
			// level of tiles.
			files.clear();
			for(std::unique_ptr<Tile>& tile : tiles) {
				PCWriter* wtr = tile->writer(true);
				// Create and emplace a LASReader.
				files.emplace_back(wtr->filenames(), wtr->x(), wtr->y());
				// Close the PCWriter.
				wtr->close();
				// Add to tmpWriters to prevent immediate destruction.
				tmpWriters.emplace_back(wtr);
			}

			// Clear the tiles list to start rebuilding it with the next level.
			tiles.clear();

			for(PCFile& file : files) {
				{
					// Get a liblas::Header from the first file to use as a template
					// for the PCWriters.
					std::ifstream istr(file.filenames()[0], std::ios::in | std::ios::binary);
					liblas::Reader irdr = rfact.CreateWithStream(istr);
					const liblas::Header& ihdr = irdr.GetHeader();

					// Create the writers for the intermediate tiles. There are four tiles
					// associated with each LASReader so only four handles at a time can be open.
					for(int r = 0; r < 2; ++r) {
						for(int c = 0; c < 2; ++c) {
							std::stringstream ss;
							double x = file.x() + c * size0;
							double y = file.y() + r * size0;
							ss << "tile_" << (int) x << "_" << (int) y << "_" << size0;
							std::string outfile = Util::pathJoin(outdir, ss.str());
							std::unique_ptr<Tile> tile(new Tile(x, y, x + size0, y + size0, buffer));
							tile->writer(new PCWriter(outfile, ihdr, x, y));
							tiles.push_back(std::move(tile));
						}
					}
				}

				// Now run over the files and reat their points into the new tiles.
				for(const std::string& filename : file.filenames()) {
					std::ifstream istr(filename, std::ios::in | std::ios::binary);
					liblas::Reader irdr = rfact.CreateWithStream(istr);
					while(irdr.ReadNextPoint()) {
						const liblas::Point& pt = irdr.GetPoint();
						for(std::unique_ptr<Tile>& tile : tiles) {
							if(tile->containsBuffered(pt.GetX(), pt.GetY()))
								tile->writer()->addPoint(pt);
						}
					}
				}

			}

			// Destroy the previous level of tiles.
			tmpWriters.clear();

		} else {
			// This is the first run, so we read all the input files into the first set of
			// intermediate files.

			{
				// Get the header from the first file to use as a template.
				std::ifstream istr(files[0].filenames()[0], std::ios::in | std::ios::binary);
				liblas::Reader irdr = rfact.CreateWithStream(istr);
				const liblas::Header& ihdr = irdr.GetHeader();

				// Create the writers for the intermediate tiles.
				for(int r = 0; r < rows0; ++r) {
					for(int c = 0; c < cols0; ++c) {
						std::stringstream ss;
						double x = allBounds[0] + c * size0;
						double y = allBounds[1] + r * size0;
						ss << "tile_" << (int) x << "_" << (int) y << "_" << size0;
						std::string outfile = Util::pathJoin(outdir, ss.str());
						std::unique_ptr<Tile> tile(new Tile(x, y, x + size0, y + size0, buffer));
						tile->writer(new PCWriter(outfile, ihdr, x, y));
						tiles.push_back(std::move(tile));
					}
				}
			}

			// Iterate over the files, filling the tiles.
			for(PCFile& file : files) {
				for(const std::string& filename : file.filenames()) {
					std::ifstream istr(filename, std::ios::in | std::ios::binary);
					liblas::Reader irdr = rfact.CreateWithStream(istr);
					while(irdr.ReadNextPoint()) {
						const liblas::Point& pt = irdr.GetPoint();
						for(std::unique_ptr<Tile>& tile : tiles) {
							if(tile->containsBuffered(pt.GetX(), pt.GetY()))
								tile->writer()->addPoint(pt);
						}
					}
				}
			}
		}

		// At the end of each loop, we halve the tile size and
		// double the col/row size to move down the pyramid.
		size0 *= 0.5;
		cols0 *= 2;
		rows0 *= 2;

	} while(size0 >= size);

	// We don't want to delete the last batch of tiles, so this is
	// set to false.
	for(std::unique_ptr<Tile>& tile : tiles)
		tile->writer()->deleteOnDestruct(false);

	// Destroy the last batch of tiles.
	tiles.clear();
}

Tiler::~Tiler() {}




geo::pc::Point::Point(const liblas::Point& pt) :
	x(pt.GetX()), y(pt.GetY()), z(pt.GetZ()),
	point(new liblas::Point(pt)){
}

geo::pc::Point::Point(double x, double y, double z) :
	x(x), y(y), z(z),
	point(nullptr) {
}

int geo::pc::Point::classId() const {
	if(point)
		return point->GetClassification().GetClass();
	return 0;
}

double geo::pc::Point::operator[](int idx) const {
	switch(idx % 2) {
	case 0:
		return x;
	case 1:
		return y;
	}
	return 0;
}

geo::pc::Point::~Point() {
}

class MeanComputer : public Computer {
public:
	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		double sum = 0;
		int count = 0;
		for(const geo::pc::Point& pt : pts) {
			sum += pt.z;
			++count;
		}
		return count > 0 ? sum / count : std::nan("");
	}
};

class PopVarianceComputer : public MeanComputer {
public:
	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		double mean = MeanComputer::compute(pts, dists, radius);
		if(std::isnan(mean))
			return mean;
		double sum = 0;
		for(const geo::pc::Point& pt : pts)
			sum += std::pow(pt.z - mean, 2.0);
		return sum / pts.size();
	}

};

class PopStdDevComputer : public PopVarianceComputer {
public:
	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		double variance = PopVarianceComputer::compute(pts, dists, radius);
		if(std::isnan(variance))
			return variance;
		return std::sqrt(variance);
	}

};

class SampVarianceComputer : public MeanComputer {
public:
	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		if(pts.size() < 2)
			return std::nan("");
		double mean = MeanComputer::compute(pts, dists, radius);
		if(std::isnan(mean))
			return mean;
		double sum = 0;
		for(const geo::pc::Point& pt : pts)
			sum += std::pow(pt.z - mean, 2.0);
		return sum / (pts.size() - 1);
	}

};

class SampStdDevComputer : public SampVarianceComputer {
public:
	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		double variance = SampVarianceComputer::compute(pts, dists, radius);
		if(std::isnan(variance))
			return variance;
		return std::sqrt(variance);
	}

};

bool pointSort(const geo::pc::Point& a, const geo::pc::Point& b) {
	return a.z < b.z;
}

class PercentileComputer : public Computer {
public:
	double percentile;

	PercentileComputer(double percentile) :
		percentile(percentile) {
		if(percentile <= 0 || percentile >= 100)
			g_runerr("Percentile must be between 0 and 100.");
	}

	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		if(pts.size() < 2)
			return std::nan("");
		std::vector<geo::pc::Point> _pts(pts.begin(), pts.end());
		std::sort(_pts.begin(), _pts.end(), pointSort);
		size_t idx = (size_t) (_pts.size() * percentile);
		return (_pts[idx + 1].z - _pts[idx].z) / 2.0;
	}

};

class MaxComputer : public Computer {
public:

	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
	double max = std::numeric_limits<double>::lowest();
		for(const geo::pc::Point& pt : pts) {
			if(pt.z > max)
				max = pt.z;
		}
		return max;
	}

};

class MinComputer : public Computer {
public:

	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		double min = std::numeric_limits<double>::max();
		for(const geo::pc::Point& pt : pts) {
			if(pt.z < min)
				min = pt.z;
		}
		return min;
	}

};

class DensityComputer : public Computer {
public:
	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		double area = radius * radius * M_PI;
		if(pts.empty())
			return std::nan("");
		return pts.size() / area;
	}

};

class CountComputer : public Computer {
public:
	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		return (double) pts.size();
	}

};

class SkewComputer : public Computer {
public:
	MeanComputer meanComp;
	SampStdDevComputer sampStdDevComp;

	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		if(pts.empty())
			return std::nan("");
		// Fisher-Pearson
		double mean = meanComp.compute(pts, dists, radius);
		if(std::isnan(mean))
			return mean;
		size_t count = pts.size();
		double sum = 0.0;
		for (const geo::pc::Point& pt : pts)
			sum += std::pow(pt.z - mean, 3.0) / count;
		double sd = sampStdDevComp.compute(pts, dists, radius);
		double skew = sum / std::pow(sd, 3.0);
		return skew;
	}

};

class KurtosisComputer : public Computer {
public:
	MeanComputer meanComp;
	SampStdDevComputer sampStdDevComp;

	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		if(pts.empty())
			return std::nan("");
		// Fisher-Pearson
		double mean = meanComp.compute(pts, dists, radius);
		if(std::isnan(mean))
			return mean;
		size_t count = pts.size();
		double sum = 0.0;
		for (const geo::pc::Point& pt : pts)
			sum += std::pow(pt.z - mean, 4.0) / count;
		double sd = sampStdDevComp.compute(pts, dists, radius);
		double kurt = sum / std::pow(sd, 4.0) - 3.0;
		return kurt;
	}

};

class CoVComputer : public Computer {
public:
	MeanComputer meanComp;
	SampStdDevComputer sampStdDevComp;

	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		if(pts.empty())
			return std::nan("");
		// Fisher-Pearson
		double mean = meanComp.compute(pts, dists, radius);
		if(std::isnan(mean))
			return mean;
		double sd = sampStdDevComp.compute(pts, dists, radius);
		double cov = mean != 0 ? sd / mean : std::nan("");
		return cov;
	}

};

class CanopyGapFractionComputer : public Computer {
public:
	MeanComputer meanComp;
	SampStdDevComputer sampStdDevComp;

	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		g_runerr("Not implemented.");
	}

};

/*
 * void fcLidarBLa(const std::list<std::unique_ptr<LiDARPoint> > &values, double *result) {
	double gnd = 0.0;
	double all = 0.0;
	for (const std::unique_ptr<LiDARPoint>& pt : values) {
		if (pt->isGround())
			gnd += pt->intensity();
		if (pt->cls() < 2) // TODO: This should perhaps be filtered by class to remove bogus points.
			all += pt->intensity();
	}
	result[0] = all != 0.0 ? 1.0 - std::sqrt(gnd / all) : -9999.0;
}

void fcLidarBLb(const std::list<std::unique_ptr<LiDARPoint> > &values, double *result) {
	double gndSingle = 0.0, gndLast = 0.0, first = 0.0, single = 0.0,
			intermediate = 0.0, last = 0.0, total = 0.0;
	for (const std::unique_ptr<LiDARPoint>& pt : values) {
		if (pt->isGround()) {
			if (pt->isSingle())
				gndSingle += pt->intensity();
			if (pt->isLast())
				gndLast += pt->intensity();
		}
		if (pt->isFirst())
			first += pt->intensity();
		if (pt->isSingle())
			single += pt->intensity();
		if (pt->isIntermediate())
			intermediate += pt->intensity();
		if (pt->isLast())
			last += pt->intensity();
		total += pt->intensity(); // TODO: This should perhaps be filtered by class to remove bogus points.
	}
	if (total == 0.0) {
		result[0] = -9999.0;
	} else {
		double denom = (first + single) / total
				+ std::sqrt((intermediate + last) / total);
		if (denom == 0.0) {
			result[0] = -9999.;
		} else {
			result[0] = (gndSingle / total + std::sqrt(gndLast / total))
					/ denom;
		}
	}
}

void fcLidarIR(const std::list<std::unique_ptr<LiDARPoint> > &values, double *result) {
	double canopy = 0.0, total = 0.0;
	for (const std::unique_ptr<LiDARPoint>& pt : values) {
		if (!pt->isGround())
			canopy += pt->intensity();
		total += pt->intensity();
	}
	result[0] = total != 0.0 ? canopy / total : -9999.0;
}

void fcLidarRR(const std::list<std::unique_ptr<LiDARPoint> > &values, double *result) {
	unsigned int canopy = 0, total = 0;
	for (const std::unique_ptr<LiDARPoint>& pt : values) {
		if (!pt->isGround())
			++canopy;
		++total;
	}
	result[0] = total != 0.0 ? (double) canopy / total : -9999.0;
}

void fcLidarFR(const std::list<std::unique_ptr<LiDARPoint> > &values, double *result) {
	unsigned int canopy = 0, total = 0;
	for (const std::unique_ptr<LiDARPoint>& pt : values) {
		if (pt->isFirst()) {
			if (!pt->isGround())
				++canopy;
			++total;
		}
	}
	result[0] = total != 0.0 ? (double) canopy / total : -9999.0;
}

void ccf(const std::list<std::unique_ptr<LiDARPoint> >& values, double *result, double threshold) {
	if (values.size() < 75) {
		result[0] = -9999.0;
	} else {
		double maxZ = -9999.0;
		for (const std::unique_ptr<LiDARPoint>& pt : values)
			maxZ = g_max(maxZ, pt->z());
		double htIncrement = (maxZ - threshold) / 20.0;
		double curHeight = threshold;
		for (int band = 0; band <= 20; ++band) {
			double count = 0;
			for (const std::unique_ptr<LiDARPoint>& pt : values) {
				if (pt->z() > curHeight)
					++count;
			}
			result[band] = (double) count / values.size();
			curHeight += htIncrement;
		}
	}
}

void gap(const std::list<std::unique_ptr<LiDARPoint> > &values, double *result, double threshold) {
	if (!values.size()) {
		result[0] = -9999.0;
	} else {
		int cnt = 0;
		for (const std::unique_ptr<LiDARPoint>& p : values) {
			if (p->z() > threshold)
				++cnt;
		}
		result[0] = 1.0 - ((double) cnt / values.size());
	}
}

CellGapFraction::CellGapFraction(unsigned char type, double threshold) :
		CellStats(), m_type(type), m_threshold(threshold) {
}

void CellGapFraction::threshold(double t) {
	m_threshold = t;
}

int CellGapFraction::bands() const {
	switch (m_type) {
	case GAP_CCF:
		return 21;
	default:
		return 1;
	}
}

void CellGapFraction::compute(double x, double y,
		const std::list<std::unique_ptr<LiDARPoint> > &values, double *result) {
	if (!values.size()) {
		result[0] = -9999.0;
	} else {
		switch (m_type) {
		case GAP_BLA:
			fcLidarBLa(values, result);
			break;
		case GAP_BLB:
			fcLidarBLb(values, result);
			break;
		case GAP_IR:
			fcLidarIR(values, result);
			break;
		case GAP_RR:
			fcLidarRR(values, result);
			break;
		case GAP_FR:
			fcLidarFR(values, result);
			break;
		case GAP_CCF:
			ccf(values, result, m_threshold);
			break;
		case GAP_GAP:
			gap(values, result, m_threshold);
			break;
		default:
			g_argerr("Unknown Gap Fraction method: " << m_type);
		}
	}
}
 *
 */

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K> Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3 Point_3;
typedef K::Plane_3 Plane_3;
typedef Delaunay::Finite_faces_iterator Finite_faces_iterator;
typedef Delaunay::Face Face;

class RugosityComputer : public Computer {
private:

	/**
	 * Compute planar area.
	 */
	double computePArea(double x1, double y1, double z1, double x2, double y2,
			double z2, double x3, double y3, double z3) {
		double side0 = std::sqrt(std::pow(x1 - x2, 2.0) + std::pow(y1 - y2, 2.0) + std::pow(z1 - z2, 2.0));
		double side1 = std::sqrt(std::pow(x2 - x3, 2.0) + std::pow(y2 - y3, 2.0) + std::pow(z2 - z3, 2.0));
		double side2 = std::sqrt(std::pow(x3 - x1, 2.0) + std::pow(y3 - y1, 2.0) + std::pow(z3 - z1, 2.0));
		double s = (side0 + side1 + side2) / 2.0;
		return std::sqrt(s * (s - side0) * (s - side1) * (s - side2));
	}

	/**
	 * Compute the area of a face.
	 */
	double computeFArea(const Face &face) {
		Point_3 p1 = face.vertex(0)->point();
		Point_3 p2 = face.vertex(1)->point();
		Point_3 p3 = face.vertex(2)->point();
		return computePArea(p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z(), p3.x(),
				p3.y(), p3.z());
	}

	double toPlane(const Point_3 &p, const Plane_3 &plane, const Point_3 &centroid) {
		return (p.x() * plane.a() + p.y() * plane.b() + plane.d()) / -plane.c();
	}

	double polyArea(const std::list<Point_3> &hull, const Plane_3 &plane,
			const Point_3 &centroid) {
		double area = 0.0;
		auto it0 = hull.begin();
		auto it1 = hull.begin();
		it1++;
		do {
			double z0 = toPlane(*it0, plane, centroid);
			double z1 = toPlane(*it1, plane, centroid);
			area += computePArea(it0->x(), it0->y(), z0, it1->x(), it1->y(), z1,
					centroid.x(), centroid.y(), centroid.z());
			it0++;
			it1++;
			if (it1 == hull.end())
				it1 = hull.begin();
		} while (it0 != hull.end());
		return area;
	}


public:

	double compute(const std::list<geo::pc::Point>& pts, const std::list<double>& dists, double radius) {
		if(pts.empty())
			return std::nan("");

		double area = radius * radius * M_PI;
		double density = pts.size() / area;

		std::list<Point_3> verts;
		for (const geo::pc::Point& pt : pts)
			verts.push_back(Point_3(pt.x, pt.y, pt.z));

		// Convex hull and POBF
		std::list<Point_3> hull;
		Plane_3 plane;
		Point_3 centroid;
		CGAL::convex_hull_2(verts.begin(), verts.end(), std::back_inserter(hull), Gt());
		CGAL::linear_least_squares_fitting_3(hull.begin(), hull.end(), plane, centroid, CGAL::Dimension_tag<0>());

		// POBF surface area.
		double parea = polyArea(hull, plane, centroid);

		// If the poly area is zero, quit.
		if(parea <= 0)
			return std::nan("");

		// Delaunay 3D surface area.
		double tarea = 0.0;
		Delaunay dt(verts.begin(), verts.end());
		for (Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it)
			tarea += computeFArea(*it);

		// TODO: This is an attempt at modelling the relationship between the ACR and
		// density. The fractal dimension is involved. Should be redone and documented.
		double densityFactor = 1.0 / (2.49127261 + 9.01659384 * std::sqrt(density * 32.65748276));

		double acr = (tarea / parea) * densityFactor;

		return acr;
	}

};

/*
 *         namespace pointstats_config {

            std::map<std::string, uint8_t> types = {
                {"Minimum", TYPE_MIN},
                {"Maximum", TYPE_MAX},
                {"Mean", TYPE_MEAN},
                {"Density", TYPE_DENSITY},
                {"Sample Variance", TYPE_VARIANCE},
                {"Sample Std. Dev.", TYPE_STDDEV},
                {"Population Variance", TYPE_PVARIANCE},
                {"Population Std. Dev.", TYPE_PSTDDEV},
                {"Count", TYPE_COUNT},
                {"Quantile", TYPE_QUANTILE},
                {"Median", TYPE_MEDIAN},
                {"Rugosity", TYPE_RUGOSITY},
                {"Kurtosis", TYPE_KURTOSIS},
                {"Skewness", TYPE_SKEW},
                {"Gap Fraction", TYPE_GAP_FRACTION},
                {"CoV", TYPE_COV}
            };
            std::map<std::string, uint8_t> attributes = {
                {"Height", ATT_HEIGHT},
                {"Intensity", ATT_INTENSITY}
            };
            std::map<std::string, uint8_t> gapFractionTypes = {
                {"IR", GAP_IR},
                {"BLa", GAP_BLA},
                {"BLb", GAP_BLB},
                {"RR", GAP_RR},
                {"FR", GAP_FR},
                {"CCF", GAP_CCF},
                {"GAP", GAP_GAP}
            };

            std::map<std::string, uint8_t> snapModes = {
                {"None" , SNAP_NONE},
                {"Grid" , SNAP_GRID},
                {"Origin" , SNAP_ORIGIN}
            };

            std::map<std::string, uint8_t> areaModes = {
                {"Full Cell", AREA_CELL},
                {"Radius", AREA_RADIUS}
            };

        } // config
 *
 */

Rasterizer::Rasterizer(const std::vector<std::string> filenames) :
	m_tree(nullptr) {
	for(const std::string& filename : filenames)
		m_files.emplace_back(filename);
	addComputer("min", new MinComputer());
	addComputer("min", new MaxComputer());
	addComputer("percentile-5", new PercentileComputer(0.05));
	addComputer("decile-1", new PercentileComputer(0.1));
	addComputer("decile-2", new PercentileComputer(0.2));
	addComputer("quartile-1", new PercentileComputer(0.25));
	addComputer("decile-3", new PercentileComputer(0.3));
	addComputer("decile-4", new PercentileComputer(0.4));
	addComputer("decile-5", new PercentileComputer(0.5));
	addComputer("quantile-2", new PercentileComputer(0.5));
	addComputer("median", new PercentileComputer(0.5));
	addComputer("decile-6", new PercentileComputer(0.6));
	addComputer("decile-7", new PercentileComputer(0.7));
	addComputer("quantile-3", new PercentileComputer(0.75));
	addComputer("decile-8", new PercentileComputer(0.8));
	addComputer("decile-9", new PercentileComputer(0.9));
	addComputer("percentile-95", new PercentileComputer(0.95));
	addComputer("mean", new MeanComputer());
	addComputer("sample-variance", new SampVarianceComputer());
	addComputer("sample-std-dev", new SampStdDevComputer());
	addComputer("population-variance", new PopVarianceComputer());
	addComputer("population-std-dev", new PopStdDevComputer());
	addComputer("rugosity-acr", new RugosityComputer());
}

bool Rasterizer::filter(const geo::pc::Point& pt) const {
	return pt.classId() == 2; // TODO: Make configurable.
}

void Rasterizer::updateTree(double x, double y, double radius) {
	std::unordered_set<int> requiredFiles;
	double tbounds[4] = {x - radius, y - radius, x + radius, y + radius};
	double fBounds[6];
	for(size_t i = 0; i < m_files.size(); ++i) {
		m_files[i].bounds(fBounds);
		if(!(tbounds[2] < fBounds[0] || tbounds[0] > fBounds[2] ||
			tbounds[3] < fBounds[1] || tbounds[1] > fBounds[3])) {
			requiredFiles.insert(i);
		}
	}
	bool changed = false;
	if(m_currentFiles.empty()) {
		m_currentFiles.insert(requiredFiles.begin(), requiredFiles.end());
		changed = true;
	} else {
		for(int a : m_currentFiles) {
			if(requiredFiles.find(a) == requiredFiles.end()) {
				changed = true;
				break;
			}
		}
		for(int a : requiredFiles) {
			if(m_currentFiles.find(a) == m_currentFiles.end()) {
				changed = true;
				break;
			}
		}
		if(changed) {
			m_currentFiles.clear();
			m_currentFiles.insert(requiredFiles.begin(), requiredFiles.end());
		}
	}
	if(changed) {
		if(m_tree) {
			delete m_tree;
			m_tree = nullptr;
		}
		if(!m_currentFiles.empty()) {
			liblas::ReaderFactory fact;
			m_tree = new geo::ds::KDTree<geo::pc::Point>(2);
			for(int a : m_currentFiles) {
				std::ifstream str(m_files[a].filenames()[0]);
				liblas::Reader rdr = fact.CreateWithStream(str);
				while(rdr.ReadNextPoint()) {
					const liblas::Point& pt = rdr.GetPoint();
					geo::pc::Point lpt(pt);
					if(filter(lpt))
						m_tree->add(lpt);
				}
			}
			try {
				m_tree->build();
			} catch(const std::exception& ex) {
				g_warn("Failed to build tree.");
			}
		}
	}
}

int Rasterizer::getPoints(double x, double y, double radius, int count, 
		std::list<geo::pc::Point>& pts, std::list<double>& dists) {
	updateTree(x, y, radius);
	geo::pc::Point pt(x, y, 0);
	int ret = 0;
	if(m_tree)
		ret = m_tree->radSearch(pt, radius, count, std::back_inserter(pts), std::back_inserter(dists));
	return ret;
}

void Rasterizer::addComputer(const std::string& name, Computer* computer) {
	if(m_computers.find(name) != m_computers.end())
		delete m_computers[name];
	m_computers[name] = computer;
}

void Rasterizer::rasterize(const std::string& filename, const std::string& type, double res, 
	double easting, double northing, double radius, int srid, int threads, double ext) {

	Computer* comp = m_computers[type];
	if(!comp)
		g_runerr("No computer for type " << type);

	double bounds[4] = {9999999999, 9999999999, -9999999999, -9999999999};
	double fBounds[6];

	for(const PCFile& f : m_files) {
		f.bounds(fBounds);
		if(fBounds[0] < bounds[0]) bounds[0] = fBounds[0];
		if(fBounds[1] < bounds[1]) bounds[1] = fBounds[1];
		if(fBounds[2] > bounds[2]) bounds[2] = fBounds[2];
		if(fBounds[3] > bounds[3]) bounds[3] = fBounds[3];
	}

	if(easting <= 0)
		easting = ((int) (bounds[0] / res)) * res;
	if(northing <= 0)
		northing = ((int) (bounds[3] / res)) * res;

	int cols = (int) ((bounds[2] - bounds[0]) / res) + 1;
	int rows = (int) ((bounds[3] - bounds[1]) / res) + 1;

	GridProps props;
	props.setTrans(easting, res, northing, -res);
	props.setSize(cols, rows);
	props.setNoData(-9999.0);
	props.setDataType(DataType::Float32);
	props.setSrid(srid);
	props.setWritable(true);
	Raster rast(filename, props);

	for(int r = 0; r < rows; ++r) {
		for(int c = 0; c < cols; ++c) {
			int count = 100;
			double x = props.toX(c) + props.resolutionX() / 2;
			double y = props.toY(r) + props.resolutionY() / 2;
			std::list<Point> pts;
			std::list<double> dists;
			int ret;
			while((ret = getPoints(x, y, radius, count, pts, dists)) >= count) {
				pts.clear();
				dists.clear();
				count *= 2;
				std::cerr << "count " << count << "; ret " << ret << "\n";
			}
			if(ret) {
				rast.setFloat(c, r, comp->compute(pts, dists, radius));
			} else {
				rast.setFloat(c, r, -9999.0);
			}
		}
	}
}

Rasterizer::~Rasterizer() {
	if(m_tree)
		delete m_tree;
	for(auto& c : m_computers)
		delete c.second;
}
