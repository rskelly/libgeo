#include <queue>
#include <string>
#include <fstream>
#include <atomic>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <condition_variable>
#include <chrono>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/filesystem.hpp>

#include <gdal.h>
#include <gdal_alg.h>
#include <ogr_feature.h>
#include <ogrsf_frmts.h>
#include <cpl_string.h>
#include <cpl_port.h>
#include <cpl_string.h>

#include <geos/geom/PrecisionModel.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/CoordinateSequenceFactory.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/Point.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/operation/union/CascadedPolygonUnion.h>

#include "omp.h"

#include "geo.hpp"
#include "util.hpp"
#include "raster.hpp"
#include "db.hpp"

using namespace geo::util;
using namespace geo::raster;
using namespace geos::geom;
using namespace geos::operation::geounion;

// Dummy cancel variable for when a cancel flag
// isn't passed.
static bool s_cancel = false;

namespace geo {

	namespace raster {

		namespace util {

			int getTypeSize(DataType type) {
				switch(type) {
				case DataType::Byte: return sizeof(uint8_t);
				case DataType::Float32: return sizeof(float);
				case DataType::Float64: return sizeof(double);
				case DataType::Int16: return sizeof(int16_t);
				case DataType::Int32: return sizeof(int32_t);
				case DataType::UInt16: return sizeof(uint16_t);
				case DataType::UInt32: return sizeof(uint32_t);
				default:
					g_runerr("No size for type: " << type);
				}
			}

			GDALDataType dataType2GDT(DataType type) {
				switch(type) {
				case DataType::Byte:  	return GDT_Byte;
				case DataType::UInt16: 	return GDT_UInt16;
				case DataType::UInt32:	return GDT_UInt32;
				case DataType::Int16:	return GDT_Int16;
				case DataType::Int32:	return GDT_Int32;
				case DataType::Float64:	return GDT_Float64;
				case DataType::Float32:	return GDT_Float32;
				case DataType::None:
				default:
					break;
				}
				return GDT_Unknown;
			}

			DataType gdt2DataType(GDALDataType type) {
				switch(type) {
				case GDT_Byte:	  	return DataType::Byte;
				case GDT_UInt16: 	return DataType::UInt16;
				case GDT_UInt32:	return DataType::UInt32;
				case GDT_Int16:		return DataType::Int16;
				case GDT_Int32:		return DataType::Int32;
				case GDT_Float64:	return DataType::Float64;
				case GDT_Float32:	return DataType::Float32;
				case GDT_Unknown:
				case GDT_CInt16:
				case GDT_CInt32:
				case GDT_CFloat32:
				case GDT_CFloat64:
				case GDT_TypeCount:
				default:
					break;
				}
				return DataType::None;
			}

			template <class T>
			inline void writeToBlock(void *block, GDALDataType type, T value, int idx) {
				switch (type) {
				case GDT_Float32:
					*(((float *)block) + idx) = (float)value;
					break;
				case GDT_Float64:
					*(((double *)block) + idx) = (double)value;
					break;
				case GDT_UInt32:
					*(((uint32_t *)block) + idx) = (uint32_t)value;
					break;
				case GDT_UInt16:
					*(((uint16_t *)block) + idx) = (uint16_t)value;
					break;
				case GDT_Int32:
					*(((int32_t *)block) + idx) = (int32_t)value;
					break;
				case GDT_Int16:
					*(((int16_t *)block) + idx) = (int16_t)value;
					break;
				case GDT_Byte:
					*(((uint8_t *)block) + idx) = (uint8_t)value;
					break;
				default:
					g_runerr("Data type not implemented: " << type);
					break;
				}
			}

			template <class T>
			inline void readFromBlock(void* block, GDALDataType type, T* value, int idx) {
				switch (type) {
				case GDT_Float32:
					*value = (double) *(((float *)block) + idx);
					break;
				case GDT_Float64:
					*value = (double) *(((double *)block) + idx);
					break;
				case GDT_UInt32:
					*value = (double) *(((uint32_t *)block) + idx);
					break;
				case GDT_UInt16:
					*value = (double) *(((uint16_t *)block) + idx);
					break;
				case GDT_Int32:
					*value = (double) *(((int32_t *)block) + idx);
					break;
				case GDT_Int16:
					*value = (double) *(((int16_t *)block) + idx);
					break;
				case GDT_Byte:
					*value = (double) *(((uint8_t *)block) + idx);
					break;
				default:
					g_runerr("Data type not implemented: " << type);
					break;
				}
			}

		} //util
	} // raster
} // geo

GridProps::GridProps() :
		m_cols(0), m_rows(0),
		m_vsrid(0), m_hsrid(0),		// Vertical and horizontal srid
		m_bands(1),           		// The number of bands
		m_writable(false),			// True if the grid is writable
		m_nodata(0),
		m_nodataSet(false),
		m_type(DataType::None) {	// The data type.
}

std::string GridProps::driver() const {
	return m_driver;
}

void GridProps::setDriver(const std::string &name) {
	m_driver = name;
}

bool GridProps::isInt() const {
	switch(m_type) {
	case DataType::Byte:
	case DataType::Int16:
	case DataType::Int32:
	case DataType::UInt16:
	case DataType::UInt32:
		return true;
	case DataType::Float32:
	case DataType::Float64:
		return false;
	default:
		g_runerr("Can't decide whether float or int: " << m_type);
	}
}

bool GridProps::isFloat() const {
	return !isInt();
}

uint64_t GridProps::size() const {
	return (uint64_t) cols() * rows();
}

double GridProps::nodata() const {
	return m_nodata;
}

void GridProps::unsetNodata() {
	m_nodataSet = false;
}

void GridProps::bounds(double* bounds) const {
	double x0 = m_trans[0];
	double y0 = m_trans[3];
	double x1 = x0 + m_trans[1] * m_cols;
	double y1 = y0 + m_trans[5] * m_rows;
	bounds[0] = g_min(x0, x1);
	bounds[1] = g_min(y0, y1);
	bounds[2] = g_max(x0, x1);
	bounds[3] = g_max(y0, y1);
}

Bounds GridProps::bounds() const {
	double x0 = m_trans[0];
	double y0 = m_trans[3];
	double x1 = x0 + m_trans[1] * m_cols;
	double y1 = y0 + m_trans[5] * m_rows;
	return Bounds(g_min(x0, x1), g_min(y0, y1), g_max(x0, x1), g_max(y0, y1));
}

void GridProps::setNoData(double nodata) {
	m_nodata = nodata;
	m_nodataSet = true;
}

bool GridProps::nodataSet() const {
	return m_nodataSet;
}

bool GridProps::hasCell(int col, int row) const {
	return !(col < 0 || row < 0 || row >= m_rows || col >= m_cols);
}

bool GridProps::hasCell(double x, double y) const {
	return hasCell(toCol(x), toRow(y));
}

int GridProps::toRow(double y) const {
	return (int) ((y - m_trans[3]) / m_trans[5]);
}

int GridProps::toCol(double x) const {
	return (int) ((x - m_trans[0]) / m_trans[1]);
}

double GridProps::toX(int col) const {
	return m_trans[0] + col * m_trans[1];
}


double GridProps::toY(int row) const {
	return m_trans[3] + row * m_trans[5];
}

// Returns the x-coordinate for the cell centroid of a given column.
double GridProps::toCentroidX(int col) const {
	return toX(col) + m_trans[1] * 0.5;
}

// Returns the y-coordinate for the cell centorid of a given row.
double GridProps::toCentroidY(int row) const {
	return toY(row) + m_trans[5] * 0.5;
}

void GridProps::setResolution(double resolutionX, double resolutionY) {
	m_trans[1] = resolutionX;
	m_trans[5] = resolutionY;
}

double GridProps::resolutionX() const {
	return m_trans[1];
}

double GridProps::resolutionY() const {
	return m_trans[5];
}

double GridProps::tlx() const {
	return m_trans[0];
}

double GridProps::tly() const {
	return m_trans[3];
}

void GridProps::setDataType(DataType type) {
	m_type = type;
}

DataType GridProps::dataType() const {
	return m_type;
}

void GridProps::setSize(int cols, int rows) {
	m_cols = cols;
	m_rows = rows;
}

int GridProps::cols() const {
	return m_cols;
}

int GridProps::rows() const {
	return m_rows;
}

void GridProps::setSrid(int hsrid, int vsrid) {
	m_vsrid = vsrid;
	m_hsrid = hsrid;
}

int GridProps::vsrid() const {
	return m_vsrid;
}

int GridProps::hsrid() const {
	return m_vsrid;
}

void GridProps::setProjection(const std::string &proj) {
	m_projection = proj;
}

std::string GridProps::projection() const {
	if(m_projection.empty() && m_hsrid > 0) {
		OGRSpatialReference* ref = new OGRSpatialReference();
		ref->importFromEPSG(m_hsrid);
		char *proj;
		ref->exportToWkt(&proj);
		ref->Release();
		return std::string(proj);
	} else {
		return m_projection;
	}
}

void GridProps::setTrans(double trans[6]) {
	for(int i = 0; i < 6; ++i)
		m_trans[i] = trans[i];
	setResolution(m_trans[1], m_trans[5]);
}

void GridProps::setTrans(double tlx, double resX, double tly, double resY) {
	double t[] = {tlx, resX, 0, tly, 0, resY};
	setTrans(t);
}

void GridProps::trans(double trans[6]) const {
	for(int i = 0; i < 6; ++i)
		trans[i] = m_trans[i];
}

void GridProps::setBands(int bands) {
	m_bands = bands;
}

int GridProps::bands() const {
	return m_bands;
}

void GridProps::setWritable(bool writable) {
	m_writable = writable;
}

bool GridProps::writable() const {
	return m_writable;
}


// Implementations for Cell

Cell::Cell(int col, int row) :
	col(col), row(row) {
}


// Implementation for TileIterator

Tile::Tile(Grid* tile, Grid* source, int cols, int rows,
		int col, int row, int buffer,
		int srcCol, int srcRow,
		int dstCol, int dstRow, int band, bool writeOnFlush) :
		m_tile(tile), m_source(source),
		m_cols(cols), m_rows(rows),
		m_col(col), m_row(row), m_buffer(buffer),
		m_srcCol(srcCol), m_srcRow(srcRow),
		m_dstCol(dstCol), m_dstRow(dstRow),
		m_band(band),
		m_writeOnFlush(writeOnFlush) {
}

int Tile::srcCol() const {
	return m_srcCol;
}

int Tile::srcRow() const {
	return m_srcRow;
}

int Tile::dstCol() const {
	return m_dstCol;
}

int Tile::dstRow() const {
	return m_dstRow;
}

int Tile::cols() const {
	return m_cols;
}

int Tile::rows() const {
	return m_rows;
}

int Tile::col() const {
	return m_col;
}

int Tile::row() const {
	return m_row;
}

Grid& Tile::grid() {
	return *m_tile;
}

void Tile::flush() {
	if(m_writeOnFlush && m_tile)
		writeTo(*m_source);
}

void Tile::writeTo(Grid& dest) {
	m_tile->writeTo(dest, m_cols, m_rows, m_buffer, m_buffer, m_col, m_row);
}

Tile::~Tile() {
	flush();
	if(m_tile)
		delete m_tile;
}

TileIterator::TileIterator(Grid& source, int cols, int rows, int buffer, int band) :
	m_source(source),
	m_cols(cols), m_rows(rows),
	m_buffer(buffer),
	m_curCol(0), m_curRow(0),
	m_band(band) {

	if(m_cols < m_buffer || m_rows < m_buffer)
		g_runerr("The column and row size must be larger than the buffer.");
}

TileIterator::TileIterator(const TileIterator& iter) :
		TileIterator(iter.m_source, iter.m_cols, iter.m_rows, iter.m_buffer, iter.m_band) {
}

bool TileIterator::hasNext() {
	std::lock_guard<std::mutex> lk(m_mtx);
	const GridProps& p = m_source.props();
	return m_curCol < p.cols() && m_curRow < p.rows();
}

int TileIterator::count() const {
	const GridProps& p = m_source.props();
	return g_max(1, (int) std::ceil((float) p.cols() / m_cols) * (int) std::ceil((float) p.rows() / m_rows));
}

Tile TileIterator::next() {

	const GridProps& props = m_source.props();

	int col, row, srcCol, srcRow, dstCol, dstRow, cols, rows;

	{
		std::lock_guard<std::mutex> lk(m_mtx);
		if(!(m_curCol < props.cols() && m_curRow < props.rows()))
			g_runerr("No more tiles.");

		col = m_curCol;
		row = m_curRow;
		srcCol = m_curCol > 0 ? m_curCol - m_buffer : m_curCol;
		srcRow = m_curRow > 0 ? m_curRow - m_buffer : m_curRow;
		dstCol = m_curCol > 0 ? 0 : m_buffer;
		dstRow = m_curRow > 0 ? 0 : m_buffer;
		cols = g_min(m_cols + m_buffer * 2, props.cols() - srcCol);
		rows = g_min(m_rows + m_buffer * 2, props.rows() - srcRow);

		m_curCol += m_cols;
		if(m_curCol > props.cols()) {
			m_curCol = 0;
			m_curRow += m_rows;
		}
	}

	GridProps p(props);
	p.setSize(m_cols + m_buffer * 2, m_rows + m_buffer * 2);
	MemRaster* tile = new MemRaster(p);

	m_source.writeTo(*tile, cols, rows, srcCol, srcRow, dstCol, dstRow, m_band, 1);

	return Tile(tile, &m_source, m_cols, m_rows, col, row, m_buffer, srcCol, srcRow, dstCol, dstRow, m_band, props.writable());
}

Tile TileIterator::create(Tile &tpl) {
	const GridProps& props = m_source.props();
	GridProps p(props);
	p.setSize(m_cols + m_buffer * 2, m_rows + m_buffer * 2);
	MemRaster* tile = new MemRaster(p);
	m_source.writeTo(*tile, tpl.m_cols, tpl.m_rows, tpl.m_srcCol, tpl.m_srcRow, tpl.m_dstCol, tpl.m_dstRow, m_band, 1);
	return Tile(tile, &m_source, m_cols, m_rows, tpl.m_col, tpl.m_row, m_buffer, tpl.m_srcCol, tpl.m_srcRow, tpl.m_dstCol, tpl.m_dstRow, m_band, props.writable());
}

TileIterator::~TileIterator() {
}


// Implementations for the Grid class

Grid::Grid() {
}

TileIterator Grid::iterator(int cols, int rows, int buffer, int band) {
	return TileIterator(*this, cols, rows, buffer, band);
}

void Grid::gaussianWeights(double *weights, int size, double sigma) {
	// If size is an even number, bump it up.
	if (size % 2 == 0) {
		++size;
		g_warn("Gaussian kernel size must be an odd number >=3. Bumping up to " << size);
	}
	for (int r = 0; r < size; ++r) {
		for (int c = 0; c < size; ++c) {
			int x = size / 2 - c;
			int y = size / 2 - r;
			weights[r * size + c] = (1 / (2 * G_PI * sigma * sigma))
					* pow(G_E, -((x * x + y * y) / (2.0 * sigma * sigma)));
		}
	}
}

GridStats Grid::stats(int band) {
	GridStats st;
	uint64_t i;
	const GridProps& gp = props();
	double nodata = gp.nodata();
	double v, m = 0, s = 0;
	int k = 1;
	st.sum = 0;
	st.count = 0;
	st.min = G_DBL_MAX_POS;
	st.max = G_DBL_MAX_NEG;
	// Welford's method for variance.
	// i has the index of the first dpata element.
	for (i = 0; i < gp.size(); ++i) {
		if ((v = getFloat(i, band)) != nodata) {
			double oldm = m;
			m = m + (v - m) / k;
			s = s + (v - m) * (v - oldm);
			st.sum += v;
			if(v < st.min) st.min = v;
			if(v > st.max) st.max = v;
			++st.count;
			++k;
		}
	}
	st.mean = st.sum / st.count;
	st.variance = s / st.count;
	st.stdDev = std::sqrt(st.variance);
	return st;
}

void Grid::normalize(int band) {
	GridStats st = stats(1);
	const GridProps& gp = props();
	double v, nodata = gp.nodata();
	double mean = st.mean;
	double stdDev = st.stdDev;
	for (uint64_t i = 0; i < gp.size(); ++i) {
		if ((v = getFloat(i, band)) != nodata && !std::isnan(v) && v < G_DBL_MAX_POS) {
			setFloat(i, ((v - mean) / stdDev), band);
		} else {
			setFloat(i, nodata, band);
		}
	}
}

void Grid::logNormalize(int band) {
	GridStats st = stats(1);
	const GridProps& gp = props();
	double n = st.min;
	double x = st.max;
	double e = std::exp(1.0) - 1.0;
	for(uint64_t i = 0; i < gp.size(); ++i)
		setFloat(i, std::log(1.0 + e * (getFloat(i) - n) / (x - n)));
}

void Grid::convert(Grid &g, int srcBand, int dstBand) {
	const GridProps& gp = props();
	if(g.props().isInt()) {
		for (uint64_t i = 0; i < gp.size(); ++i)
			g.setInt(i, getInt(i, srcBand), dstBand);
	} else {
		for (uint64_t i = 0; i < gp.size(); ++i)
			g.setFloat(i, getFloat(i, srcBand), dstBand);
	}
}

void Grid::voidFillIDW(const std::string& filename, double radius, int count, double exp, int band,
		const std::string& mask, int maskBand) {

	if(!props().isFloat())
		g_runerr("IDW fill only implemented for float rasters.");

	if (radius <= 0.0)
		throw std::invalid_argument("Radius must be larger than 0.");

	if (count <= 0)
		throw std::invalid_argument("Count must be larger than 0.");

	if (exp <= 0.0)
		throw std::invalid_argument("Exponent must be larger than 0.");

	GridProps iprops(props());
	iprops.setBands(1);
	iprops.setWritable(true);
	MemRaster input(iprops);
	writeTo(input, iprops.cols(), iprops.rows(), 0, 0, 0, 0, band);

	GridProps oprops(props());
	oprops.setBands(1);
	oprops.setWritable(true);
	MemRaster output(oprops);

	double maxDist = 100;
	bool holesOnly = true;

	double nodata = props().nodata();
	double v, d;
	int rows = props().rows();
	int cols = props().cols();

    TargetFillOperator<double, double> op1(&input, nodata, 99999);
    TargetFillOperator<double, double> op2(&input, &output, 99999, nodata);
    TargetFillOperator<double, double> op3(&input, 99999, 99998);
    int outminc, outminr, outmaxc, outmaxr;

    for (int r = 0; r < rows; ++r) {
    	if(r % 100 == 0) {
    		std::cerr << "Row " << r << " of " << rows << "\n";
    	}
		for (int c = 0; c < cols; ++c) {

			v = input.getFloat(c, r);
			if(v == 99998) {

				//output.setFloat(c, r, nodata);

			} else if (v != nodata) {

				output.setFloat(c, r, v);

			} else if(!holesOnly) {

				double dp, a = 0, b = 0;
				int cnt = 0;
				for(int r0 = g_max(0, r - maxDist); r0 < g_min(rows, r + maxDist + 1); ++r0) {
					for(int c0 = g_max(0, c - maxDist); c0 < g_min(cols, c + maxDist + 1); ++c0) {
						if((c0 == c && r0 == r) || (d = g_sq(c0 - c) + g_sq(r0 - r)) > maxDist || (v = input.getFloat(c0, r0)) == nodata)
							continue;
						dp = 1.0 / std::pow(d, exp);
						a += dp * v;
						b += dp;
						++cnt;
					}
				}
				output.setFloat(c, r, cnt ? (a / b) : nodata);

			} else {

				// Fill the hole with a unique value.
	            input.floodFill(c, r, op1, false, &outminc, &outminr, &outmaxc, &outmaxr);

	            // If it touches the edges, re-fill with nodata and continue.
	            if(outminc == 0 || outmaxc == cols - 1 || outminr == 0 || outmaxr == rows - 1) {
	            	output.floodFill(c, r, op2, false);
					input.floodFill(c, r, op3, false);
	            	continue;
	            }

	            // Find all the pixels which were filled
	            std::vector<std::tuple<int, int, double> > vpx;
	            std::vector<std::tuple<int, int> > npx;
				for(int r0 = g_max(0, outminr - 1); r0 < g_min(rows, outmaxr + 2); ++r0) {
					for(int c0 = g_max(0, outminc - 1); c0 < g_min(cols, outmaxc + 2); ++c0) {
						v = input.getFloat(c0, r0);
						if(v == 99999) {
							npx.push_back(std::make_tuple(c0, r0));
						} else if(v != nodata && v != 99998) {
							vpx.push_back(std::make_tuple(c0, r0, v));
						}
					}
				}

	            // Fill voids using the surrounding pixel values.
				int pc, pr, nc, nr, cnt;
				double dp, pv, a, b;
				for(auto& np : npx) {
					nc = std::get<0>(np);
					nr = std::get<1>(np);
					cnt = 0;
					a = 0;
					b = 0;
					for(auto& vp : vpx) {
						pc = std::get<0>(vp);
						pr = std::get<1>(vp);
						pv = std::get<2>(vp);
						d = g_sq(pc - nc) + g_sq(pr - nr);
						dp = 1.0 / std::pow(d, exp);
						a += dp * pv;
						b += dp;
						++cnt;
					}
					output.setFloat(nc, nr, cnt ? (a / b) : nodata);
				}

				// Fill again with a different value so it will be ignored.
				input.floodFill(c, r, op3, false);
			}
		}
	}

	Raster routput(filename, oprops);
	output.writeTo(routput);

}

using namespace geo::raster::util;

void Grid::smooth(Grid &smoothed, double sigma, int size, int band,
		Status* status, bool *cancel) {

	const GridProps& gp = props();

	if (!cancel)
		cancel = &s_cancel;
	if (status)
		status->update(0.01f);
	if (sigma <= 0)
		g_argerr("Sigma must be > 0.");
	if (size < 3)
		g_argerr("Kernel size must be 3 or larger.");
	if (size % 2 == 0) {
		g_warn("Kernel size must be odd. Rounding up.");
		size++;
	}

	Buffer weightsBuf(size * size * getTypeSize(DataType::Float64));
	double* weights = (double*) weightsBuf.buf;
	Grid::gaussianWeights(weights, size, sigma);

	double nd = gp.nodata();
	std::atomic<int> curTile(0);

	if (status)
		status->update(0.02f);

	TileIterator iter = iterator(512, 512, size, band);
	int tiles = iter.count();

	#pragma omp parallel for
	for (int i = 0; i < tiles; ++i) {
		if (*cancel) continue;

		Tile tile = iter.next();

		Grid& grid = tile.grid();
		const GridProps& props = grid.props();
		MemRaster buf(props);
		grid.writeTo(buf);

		// Process the entire block, even the buffer parts.
		for (int r = 0; r < props.rows() - size; ++r) {
			for (int c = 0; c < props.cols() - size; ++c) {
				double v, t = 0.0;
				bool foundNodata = false;
				for (int gr = 0; gr < size; ++gr) {
					for (int gc = 0; gc < size; ++gc) {
						v = buf.getFloat(c + gc, r + gr);
						if (v == nd) {
							foundNodata = true;
							break;
						} else {
							t += weights[gr * size + gc] * v;
						}
					}
					if(foundNodata) break;
				}
				if (!foundNodata)
					grid.setFloat(c + size / 2, r + size / 2, t);
			}
		}

		tile.writeTo(smoothed);

		if (status)
			status->update(g_min(0.99f, 0.2f + (float) curTile++ / tiles * 0.97f));
	}

	if (status)
		status->update(1.0);
}

Grid::~Grid() {
}

// Implementations for MemRaster

bool MemRaster::mmapped() const {
	return m_mmapped;
}

void MemRaster::checkInit() const {
	if (m_grid == nullptr)
		g_runerr("This instance has not been initialized.");
}

MemRaster::MemRaster() :
	m_grid(nullptr),
	m_mmapped(false) {
}

MemRaster::MemRaster(const GridProps &props, bool mapped) :
	m_grid(nullptr),
	m_mmapped(mapped) {
	init(props, mapped);
}

MemRaster::~MemRaster() {
	freeMem();
}

const GridProps& MemRaster::props() const {
	return m_props;
}

void* MemRaster::grid() {
	return m_grid;
}

void MemRaster::freeMem() {
	if (m_grid) {
		if (m_mmapped) {
			delete m_mappedFile.release();
			m_grid = nullptr;
		} else {
			free(m_grid);
			m_grid = nullptr;
		}
	}
}

void MemRaster::init(const GridProps &pr, bool mapped) {
	std::lock_guard<std::mutex> lk(m_mtx);
	m_grid = nullptr;
	m_mmapped = false;
	if (pr.cols() != m_props.cols() || pr.rows() != m_props.rows()) {
		freeMem();
		m_props = GridProps(pr);
		if(m_props.isInt()) {
			m_props.setDataType(DataType::Int32);
		} else {
			m_props.setDataType(DataType::Float64);
		}
		m_mmapped = mapped;
		m_grid = nullptr;
		size_t typeSize = getTypeSize(m_props.dataType());
		size_t size = MappedFile::fixSize(typeSize * m_props.cols() * m_props.rows());
		if (mapped) {
			m_mappedFile.reset(new MappedFile(size));
			m_grid = m_mappedFile->data();
		} else {
			m_grid = malloc(size);
		}
		if (!m_grid)
			g_runerr("Failed to allocate memory for MemRaster.");
	}
}

void MemRaster::fillFloat(double value, int band) {
	checkInit();
	if(m_props.isInt()) {
		fillInt((int) value, band);
	} else {
		size_t chunk = MappedFile::pageSize();
		size_t size = m_props.size() * sizeof(double);
		Buffer buffer(chunk);
		double* buf = (double*) buffer.buf;
		for (size_t i = 0; i < chunk / sizeof(double); ++i)
			*(buf + i) = value;
		char* grid = (char*)m_grid;
		{
			std::lock_guard<std::mutex> lk(m_mtx);
			for (uint64_t i = 0; i < size; i += chunk) {
				std::memcpy(grid, buffer.buf, chunk);
				grid += chunk;
			}
		}
	}
}

void MemRaster::fillInt(int value, int band) {
	checkInit();
	if(m_props.isFloat()) {
		fillFloat((double) value, band);
	} else {
		size_t chunk = MappedFile::pageSize();
		size_t size = m_props.size() * sizeof(int);
		Buffer buffer(chunk);
		int* buf = (int*) buffer.buf;
		for(size_t i = 0; i < chunk / sizeof(int); ++i)
			*(buf + i) = value;
		char* grid = (char*) m_grid;
		{
			std::lock_guard<std::mutex> lk(m_mtx);
			for (uint64_t i = 0; i < size; i += chunk) {
				std::memcpy(grid, buffer.buf, chunk);
				grid += chunk;
			}
		}
	}
}

double MemRaster::getFloat(uint64_t idx, int band) {
	checkInit();
	if (idx < 0 || idx >= m_props.size())
		g_argerr("Index out of bounds: " << idx << "; size: " << m_props.size());
	if(m_props.isInt()) {
		return (double) getInt(idx, band);
	} else {
		return *(((double*) m_grid) + idx);
	}
}

double MemRaster::getFloat(int col, int row, int band) {
	uint64_t idx = (uint64_t) row * m_props.cols() + col;
	return getFloat(idx, band);
}

int MemRaster::getIntRow(int row, int band, int* buf) {
	checkInit();
	if(row < 0 || row >= m_props.rows())
		g_argerr("Row index out of bounds: " << row << "; rows: " << m_props.rows());
	if(m_props.isInt()) {
		int cols = m_props.cols();
		std::memcpy(buf, ((int*) m_grid) + (row * cols), cols * sizeof(int));
		return cols;
	} else {
		g_runerr("Not an integer raster.");
	}
	return 0;
}

int MemRaster::getFloatRow(int row, int band, float* buf) {
	checkInit();
	if(row < 0 || row >= m_props.rows())
		g_argerr("Row index out of bounds: " << row << "; rows: " << m_props.rows());
	if(m_props.isFloat()) {
		int cols = m_props.cols();
		std::memcpy(buf, ((float*) m_grid) + (row * cols), cols * sizeof(double));
		return cols;
	} else {
		g_runerr("Not a float raster.");
	}
	return 0;
}

int MemRaster::getInt(uint64_t idx, int band) {
	checkInit();
	if (idx < 0 || idx >= m_props.size())
		g_argerr("Index out of bounds: " << idx << "; size: " << m_props.size());
	if(m_props.isInt()) {
		return *(((int*) m_grid) + idx);
	} else {
		return (int) getFloat(idx, band);
	}
}

int MemRaster::getInt(int col, int row, int band) {
	uint64_t idx = (uint64_t) row * m_props.cols() + col;
	return getInt(idx, band);
}

void MemRaster::setFloat(int col, int row, double value, int band) {
	uint64_t idx = (uint64_t) row * m_props.cols() + col;
	setFloat(idx, value, band);
}

void MemRaster::setFloat(uint64_t idx, double value, int band) {
	checkInit();
	uint64_t ps = m_props.size();
	if (idx >= ps)
		g_argerr("Index out of bounds: " << idx << "; size: " << m_props.size()
						<< "; value: " << value << "; col: " << (idx % m_props.cols())
						<< "; row: " << (idx / m_props.cols()));
	if(m_props.isInt()) {
		setInt(idx, (int) value, band);
	} else {
		std::lock_guard<std::mutex> lk(m_mtx);
		*(((double *) m_grid) + idx) = value;
	}
}

void MemRaster::setInt(int col, int row, int value, int band) {
	uint64_t idx = (uint64_t) row * m_props.cols() + col;
	setInt(idx, value, band);
}

void MemRaster::setInt(uint64_t idx, int value, int band) {
	checkInit();
	uint64_t ps = m_props.size();
	if (idx >= ps)
		g_argerr("Index out of bounds: " << idx << "; size: " << m_props.size()
						<< "; value: " << value << "; col: " << (idx % m_props.cols())
						<< "; row: " << (idx / m_props.cols()));
	if(m_props.isInt()) {
		std::lock_guard<std::mutex> lk(m_mtx);
		*(((int *) m_grid) + idx) = value;
	} else {
		setFloat(idx, (double) value, band);
	}
}

void MemRaster::toMatrix(
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &mtx, int band) {
	int cols = m_props.cols();
	int rows = m_props.rows();
	for (int r = 1; r < rows; ++r) {
		for (int c = 0; c < cols; ++c)
			mtx(r, c) = getFloat(c, r, band);
	}
}

void MemRaster::fromMatrix(
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &mtx, int band) {
	int cols = m_props.cols();
	int rows = m_props.rows();
	for (int r = 1; r < rows; ++r) {
		for (int c = 0; c < cols; ++c)
			setFloat(c, r, (double) mtx(r, c), band);
	}
}

void MemRaster::writeToRaster(Raster &grd,
			int cols, int rows,
			int srcCol, int srcRow,
			int dstCol, int dstRow,
			int srcBand, int dstBand) {

	if(dstBand < 1 || dstBand > grd.m_props.bands())
		g_argerr("Invalid destination band: " << dstBand);

	cols = g_abs(cols);
	rows = g_abs(rows);
	if(cols == 0) cols = m_props.cols();
	if(rows == 0) rows = m_props.rows();
	cols = g_min(m_props.cols() - srcCol, cols);
	rows = g_min(m_props.rows() - srcRow, rows);
	cols = g_min(grd.m_props.cols() - dstCol, cols);
	rows = g_min(grd.m_props.rows() - dstRow, rows);

	const GridProps& gp = m_props;
	GDALDataType gtype = dataType2GDT(gp.dataType()); // TODO: The type of the buffer, not the type of the memraster.
	int typeSize = getTypeSize(gp.dataType());
	int gcols = gp.cols();

	GDALRasterBand *band = grd.m_ds->GetRasterBand(dstBand);

	char* input  = (char*) grid();

	{
		std::lock_guard<std::mutex> lk1(grd.m_mtx);
		std::lock_guard<std::mutex> lk0(m_mtx);
		for(int r = 0; r < rows; ++r) {
			//std::memcpy(output + r * cols * typeSize, input + ((srcRow + r) * gcols + srcCol) * typeSize, cols * typeSize);
			if(CPLE_None != band->RasterIO(GF_Write, dstCol, dstRow + r, cols, 1, input + ((srcRow + r) * gcols + srcCol) * typeSize,
					cols, 1, gtype, 0, 0, 0))
				g_runerr("Failed to write to: " << grd.filename());
		}
	}

}

void MemRaster::writeToMemRaster(MemRaster &grd,
		int cols, int rows,
		int srcCol, int srcRow,
		int dstCol, int dstRow,
		int srcBand, int dstBand) {

	cols = g_abs(cols);
	rows = g_abs(rows);
	if(cols == 0) cols = grd.m_props.cols();
	if(rows == 0) rows = grd.m_props.rows();
	cols = g_min(m_props.cols() - srcCol, cols);
	rows = g_min(m_props.rows() - srcRow, rows);
	cols = g_min(grd.m_props.cols() - dstCol, cols);
	rows = g_min(grd.m_props.rows() - dstRow, rows);

	if(grd.m_props.isInt()) {
		if(m_props.isInt()) {
			std::lock_guard<std::mutex> lk(m_mtx);
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c)
					grd.setInt(c + dstCol, r + dstRow, getInt(c + srcCol, r + srcRow, srcBand), dstBand);
			}
		} else {
			std::lock_guard<std::mutex> lk(m_mtx);
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c)
					grd.setInt(c + dstCol, r + dstRow, (int) getFloat(c + srcCol, r + srcRow, srcBand), dstBand);
			}
		}
	} else {
		if(m_props.isInt()) {
			std::lock_guard<std::mutex> lk(m_mtx);
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c)
					grd.setFloat(c + dstCol, r + dstRow, (double) getInt(c + srcCol, r + srcRow, srcBand), dstBand);
			}
		} else {
			std::lock_guard<std::mutex> lk(m_mtx);
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c)
					grd.setFloat(c + dstCol, r + dstRow, getFloat(c + srcCol, r + srcRow, srcBand), dstBand);
			}
		}
	}
}

void MemRaster::writeTo(Grid &grd,
		int cols, int rows,
		int srcCol, int srcRow,
		int dstCol, int dstRow,
		int srcBand, int dstBand) {
	if(dynamic_cast<MemRaster*>(&grd)) {
		writeToMemRaster(dynamic_cast<MemRaster&>(grd), cols, rows, srcCol, srcRow, dstCol, dstRow, srcBand, dstBand);
	} else if(dynamic_cast<Raster*>(&grd)) {
		writeToRaster(dynamic_cast<Raster&>(grd), cols, rows, srcCol, srcRow, dstCol, dstRow, srcBand, dstBand);
	} else {
		g_runerr("writeTo not implemented to handle this type of grid.");
	}
}

// Implementations for Raster

std::map<std::string, std::set<std::string> > Raster::extensions() {
	GDALAllRegister();
	std::map<std::string, std::set<std::string> > extensions;
	GDALDriverManager* mgr = GetGDALDriverManager();
	for(int i = 0; i < mgr->GetDriverCount(); ++i) {
		GDALDriver* drv = mgr->GetDriver(i);
		const char* cc = drv->GetMetadataItem(GDAL_DCAP_RASTER);
		if(cc != NULL && std::strncmp(cc, "YES", 3) == 0) {
			const char* desc = drv->GetDescription();
			if(desc != NULL) {
				const char* ext = drv->GetMetadataItem(GDAL_DMD_EXTENSION);
				if(ext != NULL ) {
					std::list<std::string> lst;
					Util::splitString(std::back_inserter(lst), std::string(ext));
					for(const std::string &item : lst)
						extensions[desc].insert("." + Util::lower(item));
				}

			}
		}
	}
	return extensions;
}

std::map<std::string, std::string> Raster::drivers() {
	GDALAllRegister();
	std::map<std::string, std::string> drivers;
	GDALDriverManager *mgr = GetGDALDriverManager();
	for(int i = 0; i < mgr->GetDriverCount(); ++i) {
		GDALDriver *drv = mgr->GetDriver(i);
		const char* cc = drv->GetMetadataItem(GDAL_DCAP_RASTER);
		if(cc != NULL && std::strncmp(cc, "YES", 3) == 0) {
			const char* name = drv->GetMetadataItem(GDAL_DMD_LONGNAME);
			const char* desc = drv->GetDescription();
			if(name != NULL && desc != NULL) {
				drivers[desc] = name;
			}
		}
	}
	return drivers;
}

std::string Raster::getDriverForFilename(const std::string &filename) {
	std::string ext = Util::extension(filename);
	std::map<std::string, std::set<std::string> > drivers = extensions();
	std::string result;
	for(const auto &it : drivers) {
		if(it.second.find(ext) != it.second.end())
			result = it.first;
	}
	return result;
}

Raster::Raster(const std::string &filename, const GridProps &props) :
		m_ds(nullptr),
		m_bcols(0), m_brows(0),
		m_bcol(-1), m_brow(-1),
		m_band(1),
		m_dirty(false), 
		m_bband(0),
		m_type(GDT_Unknown) {

	if (props.resolutionX() == 0 || props.resolutionY() == 0)
		g_argerr("Resolution must not be zero.");
	if (props.cols() <= 0 || props.rows() <= 0)
		g_argerr("Columns and rows must be larger than zero.");
	if (filename.empty())
		g_argerr("Filename must be given.");

	m_props = GridProps(props);
	m_filename = filename;

	// Create GDAL dataset.
	char **opts = NULL;
	// TODO: Compress option in props, for tiffs
	/*
	if(m_props.compress())
		opts = CSLSetNameValue(opts, "COMPRESS", "LZW");
	if(m_props.bigTiff())
		opts = CSLSetNameValue(opts, "BIGTIFF", "YES");
	*/
	GDALAllRegister();
	std::string drvName = m_props.driver();
	if(drvName.empty())
		drvName = getDriverForFilename(m_filename);
	GDALDriver *drv = GetGDALDriverManager()->GetDriverByName(drvName.c_str());
	const char *create = drv->GetMetadataItem(GDAL_DCAP_CREATE);
	if(create == NULL || std::strncmp(create, "YES", 3) != 0)
		g_runerr("The " << drvName << " driver does not support dataset creation. Please specify a different driver.");
	Util::rm(filename);
	m_ds = drv->Create(
			filename.c_str(), m_props.cols(), m_props.rows(), m_props.bands(),
			dataType2GDT(m_props.dataType()), opts);
	if (!m_ds)
		g_runerr("Failed to create file.");

	// Initialize geotransform.
	double trans[6];
	m_props.trans(trans);
	m_ds->SetGeoTransform(trans);

	// Set projection.
	std::string proj = m_props.projection();
	if (!proj.empty())
		m_ds->SetProjection(proj.c_str());
	if(m_props.nodataSet()) {
		for(int i = 1; i <= m_props.bands(); ++i)
			m_ds->GetRasterBand(i)->SetNoDataValue(m_props.nodata());
	}
	m_ds->GetRasterBand(1)->GetBlockSize(&m_bcols, &m_brows);
}

Raster::Raster(const std::string &filename, bool writable) :
		m_ds(nullptr),
		m_bcols(0), m_brows(0),
		m_bcol(-1), m_brow(-1),
		m_band(1),
		m_dirty(false),
		m_bband(0),
		m_type(GDT_Unknown) {

	if (filename.empty())
		g_argerr("Filename must be given.");

	m_filename = filename;

	// Attempt to open the dataset.
	GDALAllRegister();
	m_ds = (GDALDataset *) GDALOpen(filename.c_str(), writable ? GA_Update : GA_ReadOnly);
	if (m_ds == NULL)
		g_runerr("Failed to open raster.");

	GDALDriver *drv = m_ds->GetDriver();
	if(drv == NULL)
		g_runerr("Failed to retrieve driver.");
	const char *drvName = drv->GetDescription();
	if(drvName != NULL)
		m_props.setDriver(drvName);

	m_type = m_ds->GetRasterBand(1)->GetRasterDataType();
	// Save some raster properties
	double trans[6];
	m_ds->GetGeoTransform(trans);
	m_props.setTrans(trans);
	m_props.setSize(m_ds->GetRasterXSize(), m_ds->GetRasterYSize());
	m_props.setDataType(gdt2DataType(m_type));
	m_props.setBands(m_ds->GetRasterCount());
	m_props.setWritable(writable);
	m_props.setProjection(std::string(m_ds->GetProjectionRef()));
	m_props.setNoData(m_ds->GetRasterBand(1)->GetNoDataValue()); // TODO: This might not be a real nodata value.
	m_ds->GetRasterBand(1)->GetBlockSize(&m_bcols, &m_brows);
}

void* Raster::getBlock(int band) {
	if(m_blocks.find(band) == m_blocks.end()) {
		m_blocks[band] = malloc(m_bcols * m_brows * getTypeSize(m_props.dataType()));
	}
	return m_blocks[band];
}

GDALDataset* Raster::ds() const {
	return m_ds;
}

const GridProps& Raster::props() const {
	return m_props;
}

DataType Raster::getFileDataType(const std::string &filename) {
	GDALDataset *ds = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly);
	DataType type = gdt2DataType(ds->GetRasterBand(1)->GetRasterDataType());
	GDALClose(ds);
	return type;
}

std::string Raster::filename() const {
	return m_filename;
}

void writeInt(char* buf, int size, int value) {
	for(int i = size - 1; i >= 0; --i)
		*(buf + i) = (value >> (size - i - 1)) & 0xff;

}

void Raster::fillInt(int value, int band) {
	GDALRasterBand *bnd = m_ds->GetRasterBand(band);
	std::lock_guard<std::mutex> lk(m_mtx);
	bnd->Fill((int) value);
}

void Raster::fillFloat(double value, int band) {
	GDALRasterBand *bnd = m_ds->GetRasterBand(band);
	std::lock_guard<std::mutex> lk(m_mtx);
	bnd->Fill(value);
}

void Raster::writeToRaster(Raster &grd,
			int cols, int rows,
			int srcCol, int srcRow,
			int dstCol, int dstRow,
			int srcBand, int dstBand) {

	if(m_dirty) {
		GDALRasterBand* bnd = m_ds->GetRasterBand(m_bband);
		if(!bnd)
			g_runerr("Failed to find band " << m_bband);
		std::lock_guard<std::mutex> lk(m_mtx);
		if(CPLE_None != bnd->WriteBlock(m_bcol, m_brow, getBlock(m_bband)))
			g_runerr("Failed to flush to: " << filename());
		m_dirty = false;
	}

	if (&grd == this)
		g_runerr("Recursive call to readBlock.");
	if(srcBand < 1 || srcBand > m_props.bands())
		g_argerr("Invalid source band: " << srcBand);
	if(dstBand < 1 || dstBand > grd.m_props.bands())
		g_argerr("Invalid destination band: " << srcBand);

	cols = g_abs(cols);
	rows = g_abs(rows);
	if(cols == 0) cols = grd.m_props.cols();
	if(rows == 0) rows = grd.m_props.rows();
	cols = g_min(m_props.cols() - srcCol, cols);
	rows = g_min(m_props.rows() - srcRow, rows);
	cols = g_min(grd.m_props.cols() - dstCol, cols);
	rows = g_min(grd.m_props.rows() - dstRow, rows);

	// A buffer for the entire copy.
	DataType type = grd.m_props.dataType();
	GDALDataType gtype = dataType2GDT(type);
	int typeSize = getTypeSize(type);
	Buffer buf(cols * rows * typeSize);

	GDALRasterBand* srcBnd = m_ds->GetRasterBand(srcBand);
	if(!srcBnd)
		g_runerr("Failed to find source band " << srcBand);
	GDALRasterBand* dstBnd = grd.m_ds->GetRasterBand(dstBand);
	if(!dstBnd)
		g_runerr("Failed to find destination band " << dstBand);
	{
		std::lock_guard<std::mutex> lk(m_mtx);
		if(CPLE_None != srcBnd->RasterIO(GF_Read, srcCol, srcRow, cols, rows,
				buf.buf, cols, rows, gtype, 0, 0, 0))
			g_runerr("Failed to read from: " << grd.filename());
		if(CPLE_None != dstBnd->RasterIO(GF_Write, dstCol, dstRow, cols, rows,
				buf.buf, cols, rows, gtype, 0, 0, 0))
			g_runerr("Failed to write to: " << filename());
	}
}

void Raster::writeToMemRaster(MemRaster &grd,
		int cols, int rows,
		int srcCol, int srcRow,
		int dstCol, int dstRow,
		int srcBand, int dstBand) {

	if(m_dirty) {
		GDALRasterBand* bnd = m_ds->GetRasterBand(m_bband);
		if(!bnd)
			g_runerr("Failed to find band " << bnd);
		std::lock_guard<std::mutex> lk(m_mtx);
		if(CPLE_None != bnd->WriteBlock(m_bcol, m_brow, getBlock(m_bband)))
			g_runerr("Failed to flush to: " << filename());
		m_dirty = false;
	}

	if(srcBand < 1 || srcBand > m_props.bands())
		g_argerr("Invalid source band: " << srcBand);

	cols = g_abs(cols);
	rows = g_abs(rows);
	if(cols == 0) cols = grd.props().cols();
	if(rows == 0) rows = grd.props().rows();
	cols = g_min(m_props.cols() - srcCol, cols);
	rows = g_min(m_props.rows() - srcRow, rows);
	cols = g_min(grd.props().cols() - dstCol, cols);
	rows = g_min(grd.props().rows() - dstRow, rows);

	const GridProps& gp = grd.props();
	GDALDataType gtype = dataType2GDT(gp.dataType());
	uint64_t typeSize = getTypeSize(gp.dataType());
	uint64_t gcols = gp.cols();

	Buffer buf(cols * rows * typeSize);
	GDALRasterBand *band = m_ds->GetRasterBand(srcBand);
	if(!band)
		g_runerr("Failed to find band " << srcBand);
	{
		std::lock_guard<std::mutex> lk(m_mtx);
		if(CPLE_None != band->RasterIO(GF_Read, srcCol, srcRow, cols, rows, buf.buf,
				cols, rows, gtype, 0, 0, 0))
			g_runerr("Failed to read from: " << filename());
	}

	char* output = (char*) grd.grid();
	char* input  = (char*) buf.buf;
	uint64_t len = cols * typeSize;
	for(int r = 0; r < rows; ++r) {
		uint64_t doff = ((dstRow + r) * gcols + dstCol) * typeSize;
		uint64_t soff = (uint64_t) r * cols * typeSize;
		std::memcpy(output + doff , input + soff, len);
	}

}

void Raster::writeTo(Grid &grd,
		int cols, int rows,
		int srcCol, int srcRow,
		int dstCol, int dstRow,
		int srcBand, int dstBand) {
	if(dynamic_cast<MemRaster*>(&grd)) {
		writeToMemRaster(dynamic_cast<MemRaster&>(grd), cols, rows, srcCol, srcRow, dstCol, dstRow, srcBand, dstBand);
	} else if(dynamic_cast<Raster*>(&grd)) {
		writeToRaster(dynamic_cast<Raster&>(grd), cols, rows, srcCol, srcRow, dstCol, dstRow, srcBand, dstBand);
	} else {
		g_runerr("writeTo not implemented to handle this type of grid.");
	}
}

GDALDataType Raster::getGDType() const {
	return dataType2GDT(m_props.dataType());
}

double Raster::getFloat(int col, int row, int band) {
	int bcol = col / m_bcols;
	int brow = row / m_brows;
	if(bcol != m_bcol || brow != m_brow || band != m_bband) {
		if(m_dirty) {
			GDALRasterBand* bnd = m_ds->GetRasterBand(m_bband);
			if(!bnd)
				g_runerr("Failed to find band " << bnd);
			std::lock_guard<std::mutex> lk(m_mtx);
			if(CPLE_None != bnd->WriteBlock(m_bcol, m_brow, getBlock(m_bband)))
				g_runerr("Failed to flush to: " << filename());
			m_dirty = false;
		}
		GDALRasterBand *rb = m_ds->GetRasterBand(band);
		if(!rb)
			g_argerr("Failed to find band " << band);
		{
			std::lock_guard<std::mutex> lk(m_mtx);
			if(CPLE_None != rb->ReadBlock(bcol, brow, getBlock(band)))
				g_runerr("Failed to read from: " << filename());
		}
		m_bcol = bcol;
		m_brow = brow;
		m_bband = band;
	}
	int idx = (row % m_brows) * m_bcols + (col % m_bcols);
	double v = 0;
	readFromBlock(getBlock(m_bband), getGDType(), &v, idx);
	return v;
}

double Raster::getFloat(uint64_t idx, int band) {
	return getFloat((int) (idx % m_props.cols()), (int) (idx / m_props.cols()), band);
}


double Raster::getFloat(double x, double y, int band) {
	return getFloat(m_props.toCol(x), m_props.toRow(y), band);
}

int Raster::getInt(int col, int row, int band) {
	int bcol = col / m_bcols;
	int brow = row / m_brows;
	if(bcol != m_bcol || brow != m_brow || band != m_bband) { // TODO: No effective cacheing if the band changes.
		if(m_dirty) {
			GDALRasterBand* bnd = m_ds->GetRasterBand(m_bband);
			if(!bnd)
				g_runerr("Failed to find band " << m_bband);
			std::lock_guard<std::mutex> lk(m_mtx);
			if(CPLE_None != bnd->WriteBlock(m_bcol, m_brow, getBlock(m_bband)))
				g_runerr("Failed to flush to: " << filename());
			m_dirty = false;
		}
		GDALRasterBand *rb = m_ds->GetRasterBand(band);
		if(!rb)
			g_argerr("Failed to find band " << band);
		{
			std::lock_guard<std::mutex> lk(m_mtx);
			if(CPLE_None != rb->ReadBlock(bcol, brow, getBlock(band)))
				g_runerr("Failed to read from: " << filename());
		}
		m_bcol = bcol;
		m_brow = brow;
		m_bband = band;
	}
	int idx = (row % m_brows) * m_bcols + (col % m_bcols);
	int v = 0;
	readFromBlock(getBlock(m_bband), getGDType(), &v, idx);
	return v;
}

int Raster::getInt(uint64_t idx, int band) {
	return getInt((int) (idx % m_props.cols()), (int) (idx / m_props.cols()), band);
}

int Raster::getInt(double x, double y, int band) {
	return getInt(m_props.toCol(x), m_props.toRow(y), band);
}

void Raster::setInt(uint64_t idx, int v, int band) {
	setInt((int) (idx % m_props.cols()), (int) (idx / m_props.cols()), v, band);
}

void Raster::setInt(double x, double y, int v, int band) {
	setInt(m_props.toCol(x), m_props.toRow(y), v, band);
}

void Raster::setFloat(int col, int row, double v, int band) {
	if (!m_props.writable())
		g_runerr("This raster is not writable.");
	int bcol = col / m_bcols;
	int brow = row / m_brows;
	if(bcol != m_bcol || brow != m_brow || band != m_bband) {
		if(m_dirty) {
			GDALRasterBand* bnd = m_ds->GetRasterBand(m_bband);
			if(!bnd)
				g_runerr("Failed to find band " << m_bband);
			if(CPLE_None != bnd->WriteBlock(m_bcol, m_brow, getBlock(m_bband)))
				g_runerr("Failed to flush to: " << filename());
			m_dirty = false;
		}
		GDALRasterBand *rb = m_ds->GetRasterBand(band);
		if(!rb)
			g_argerr("Failed to find band " << band);
		if(CPLE_None != rb->ReadBlock(bcol, brow, getBlock(band)))
			g_runerr("Failed to read from: " << filename());
		m_bcol = bcol;
		m_brow = brow;
		m_bband = band;
	}
	int idx = (row % m_brows) * m_bcols + (col % m_bcols);
	writeToBlock(getBlock(m_bband), getGDType(), v, idx);
	m_dirty = true;
}

void Raster::setInt(int col, int row, int v, int band) {
	if (!m_props.writable())
		g_runerr("This raster is not writable.");
	int bcol = col / m_bcols;
	int brow = row / m_brows;
	if(bcol != m_bcol || brow != m_brow || band != m_bband) {
		if(m_dirty) {
			GDALRasterBand* bnd = m_ds->GetRasterBand(m_bband);
			if(!bnd)
				g_runerr("Failed to find band " << m_bband);
			if(CPLE_None != bnd->WriteBlock(m_bcol, m_brow, getBlock(m_bband)))
				g_runerr("Failed to flush to: " << filename());
			m_dirty = false;
		}
		GDALRasterBand *rb = m_ds->GetRasterBand(band);
		if(!rb)
			g_argerr("Failed to find band " << band);
		if(CPLE_None != rb->ReadBlock(bcol, brow, getBlock(band)))
			g_runerr("Failed to read from: " << filename());
		m_bcol = bcol;
		m_brow = brow;
		m_bband = band;
	}
	int idx = (row % m_brows) * m_bcols + (col % m_bcols);
	writeToBlock(getBlock(m_bband), getGDType(), v, idx);
	m_dirty = true;
}

void Raster::setFloat(uint64_t idx, double v, int band) {
	setFloat((int) (idx % m_props.cols()), (int) (idx / m_props.cols()), v, band);
}

void Raster::setFloat(double x, double y, double v, int band) {
	setFloat(m_props.toCol(x), m_props.toRow(y), v, band);
}


/**
 * Class used in polygonization
 * implementation.
 */
class Poly {
private:

	// The final geometry.
	geos::geom::Geometry* m_geom;
	// The list of consitiuent geometries.
	std::vector<geos::geom::Geometry*> m_geoms;
	// Unique geometry ID.
	uint64_t m_id;
	// The minimum and maximum row
	// index from the source raster.
	long m_minRow, m_maxRow;
	// True if this poly has already been finalized.
	bool m_finalized;
	// True if the geometries have been collapsed after
	// the most recent update.
	bool m_collapsed;

	/**
	 * Add the list of polygons to this instance with the rows
	 * covered by them.
	 * @param u The vector of Polygons.
	 * @param minRow The lowest row index occupied by this geometry.
	 * @param maxRow The highest row index occupied by this geometry.
	 */
	void update(std::vector<Polygon*>& u, int minRow, int maxRow) {
		for (Polygon* p : u)
			m_geoms.push_back(static_cast<Geometry*>(p));
		update(minRow, maxRow);
		m_collapsed = false;
	}


	/**
	 * Updates the row indices covered by this polygon.
	 * @param minRow The lowest row index occupied by this geometry.
	 * @param maxRow The highest row index occupied by this geometry.
	 */
	void update(int minRow, int maxRow) {
		if(minRow < m_minRow) m_minRow = minRow;
		if(maxRow > m_maxRow) m_maxRow = maxRow;
	}


public:

	/**
	 * Create a polygon with the given ID and initial geometry,
	 * for the given rows. The factory is shared by all
	 * instances.
	 *
	 * @param id A unique ID.
	 */
	Poly(uint64_t id) :
		m_geom(nullptr),
		m_id(id),
		m_minRow(std::numeric_limits<long>::max()),
		m_maxRow(std::numeric_limits<long>::min()),
		m_finalized(false),
		m_collapsed(false) {
	}

	/**
	 * Add the contents of the Poly to this instance.
	 * @param p The Poly to add.
	 */
	void update(Poly& p) {
		p.collapse();
		update(p.geom(), p.minRow(), p.maxRow());
	}

	/**
	 * Add the polygon to this instance with the
	 * rows covered by it.
	 * @param u The Polygon to add.
	 * @param minRow The lowest row index occupied by this geometry.
	 * @param maxRow The highest row index occupied by this geometry.
	 */
	void update(Geometry* u, int minRow, int maxRow) {
		m_geoms.push_back(u);
		update(minRow, maxRow);
		m_collapsed = false;
	}

	/**
	 * Returns true if the range of rows given by start and end was finalized
	 * within the given thread, or by a thread whose block is completed.
	 * The checked range includes one row above and one below,
	 * which is required to guarantee that a polygon is completed.
	 * @param rowFinished The last row read.
	 */
	bool isRangeFinalized(std::vector<bool> rowsFinished) const {
		for(long i = std::max(0L, m_minRow - 1); i < std::max((long) rowsFinished.size(), m_maxRow + 1); ++i) {
			if(!rowsFinished[i])
				return false;
		}
		return true;
	}

	/**
	 * Return the pointer to the unioned polygon.
	 */
	Geometry* geom() const {
		return m_geom;
	}

	const std::vector<Geometry*>& geoms() const {
		return m_geoms;
	}

	void _collapse(std::vector<Polygon*>* group, std::vector<Geometry*>* geoms) {
		if(!group->empty()) {
			Geometry* p = geos::operation::geounion::CascadedPolygonUnion::CascadedPolygonUnion::Union(group);
			switch(p->getGeometryTypeId()) {
			case GEOS_POLYGON:
				geoms->push_back(p);
				break;
			case GEOS_MULTIPOLYGON:
				const Geometry *g;
				for(size_t i = 0; i < p->getNumGeometries(); ++i) {
					g = p->getGeometryN(i);
					if(!g || g->getGeometryTypeId() != GEOS_POLYGON)
						g_runerr("Null or invalid polygon.");
					geoms->push_back(g->clone());
				}
				break;
			default:
				break;
			}
		}
	}

	void collapse() {
		if(m_collapsed || m_geoms.empty())
			return;
		size_t geomCount = m_geoms.size();
		int change = 1;
		Geometry* g0;
		while(change) {
			std::vector<Polygon*> polys;
			for(Geometry* g : m_geoms) {
				switch(g->getGeometryTypeId()){
				case GEOS_MULTIPOLYGON:
					for(size_t i = 0; i < g->getNumGeometries(); ++i) {
						g0 = g->getGeometryN(i)->clone();
						if(!g0 || g0->getGeometryTypeId() != GEOS_POLYGON)
							g_runerr("Null or invalid polygon.");
						polys.push_back(dynamic_cast<Polygon*>(g0));
					}
					break;
				case GEOS_POLYGON:
					g0 = g->clone();
					if(!g0 || g0->getGeometryTypeId() != GEOS_POLYGON)
						g_runerr("Null or invalid polygon.");
					polys.push_back(dynamic_cast<Polygon*>(g0));
					break;
				default:
					g_runerr("Illegal geometry type: " << g->getGeometryTypeId());
					break;
				}
				delete g;
			}
			m_geoms.clear();

			std::cerr << "collapsing " << polys.size() << " polys\n";
			// The max group size is 1024; min 2.
			int groupSize = g_max(2, g_min(1024, polys.size() / 16));
			// If the group size is 2, just use one group.
			if(groupSize == 2)
				groupSize = polys.size();

			// Number of threads
			// TODO: Configurable. Causes breaks in polys; need a final merge.
			size_t tc = 4;
			// The number of items in one group.
			size_t groupCount = (int) std::ceil((float) polys.size() / groupSize);
			// Container for thread instances.
			std::vector<std::thread> threads(tc);
			// The output geometry groups.
			std::vector<std::vector<Geometry*> > output(tc);
			// The input polygon groups.
			std::vector<std::vector<Polygon*> > groups(tc);

			size_t i, j;
			for(size_t g = 0; g < groupCount; g += tc) {
				for(i = 0; i < tc; ++i) {
					size_t begin = (g + i) * groupSize;
					if(begin >= polys.size())
						break;
					size_t end = std::min((g + i + 1) * groupSize, polys.size());
					groups[i].clear(); // TODO: More efficient lists.
					output[i].clear();
					groups[i].assign(polys.begin() + begin, polys.begin() + end);
					threads[i] = std::thread(&Poly::_collapse, this, &groups[i], &output[i]);
				}
				for(j = 0; j < i; ++j) {
					threads[j].join();
					m_geoms.insert(m_geoms.end(), output[j].begin(), output[j].end());
				}
			}

			change = geomCount - m_geoms.size();
			geomCount = m_geoms.size();

			std::cerr << "collapse " << change << ", " << geomCount << "\n";
		}

		m_collapsed = true;
	}

	/**
	 * Generate the unioned polygon from its parts. Remove dangles and holes if required.
	 * @param fact The GeometryFactory.
	 * @param removeHoles True, if it is desired that holes in polygons be removed.
	 * @param removeDangles True, if it is desired that single-pixel artifacts, attached at the corners, be removed.
	 */
	void finalize(const GeometryFactory::unique_ptr& fact, bool removeHoles, bool removeDangles) {

		if(m_finalized)
			return;
		m_finalized = true;

		collapse();

		std::cerr << "finalize\n";

		Geometry* geom = fact->createMultiPolygon(m_geoms);

		// If we're removing dangles, throw away all but the
		// largest single polygon. If it was originally a polygon,
		// there are no dangles.
		if(removeDangles) {
			size_t idx = 0;
			double area = 0;
			for(size_t i = 0; i < geom->getNumGeometries(); ++i) {
				const Geometry* p = geom->getGeometryN(i);
				double a = p->getArea();
				if(a > area) {
					area = a;
					idx = i;
				}
			}
			Geometry *g = geom->getGeometryN(idx)->clone(); // Force copy.
			delete geom;
			geom = g;
		}

		// If we're removing holes, extract the exterior rings
		// of all constituent polygons.
		if(removeHoles) {
			std::vector<Geometry*>* geoms0 = new std::vector<Geometry*>();
			for(size_t i = 0; i < geom->getNumGeometries(); ++i) {
				const Polygon* p = dynamic_cast<const Polygon*>(geom->getGeometryN(i));
				const LineString* l = p->getExteriorRing();
				LinearRing* r = fact->createLinearRing(l->getCoordinates());
				geoms0->push_back(fact->createPolygon(r, nullptr));
			}
			Geometry* g = fact->createMultiPolygon(geoms0); // Do not copy -- take ownership.
			delete geom;
			geom = g;
		}

		m_geom = geom;
		std::cerr << "geom " << m_geom << "\n";
	}

	/**
	 * Return the unique ID.
	 */
	uint64_t id() const {
		return m_id;
	}

	/**
	 * Return the minimum row index covered by the geometry.
	 */
	int minRow() const {
		return m_minRow;
	}

	/**
	 * Return the maximum row index covered by the geometry.
	 */
	int maxRow() const {
		return m_maxRow;
	}

	~Poly() {
		if(m_geom)
			delete m_geom;
	}
};

typedef std::unordered_map<uint64_t, std::unique_ptr<Poly> > PolyMap;
typedef std::queue<std::unique_ptr<Poly> > PolyQueue;

class PolyContext {
private:
	Raster* m_raster;
	Raster* m_maskRaster;
	Status* m_status;
	int* m_block;
	bool* m_cancel;
	GDALDataset* m_ds;
	OGRLayer *m_layer;
	GeometryFactory::unique_ptr m_geomFactory;
	GEOSContextHandle_t m_gctx;
	OGRSpatialReference m_sr;
	PolyMap m_polyMap;
	PolyQueue m_polyQueue;
	std::condition_variable m_polyMapCond;
	std::condition_variable m_polyQueueCond;
	std::mutex m_polyMapMtx;    // For the PolyMap.
	std::mutex m_polyQueueMtx;  // For the PolyQueue.
	std::mutex m_finMtx;        // For the finished list.
	std::mutex m_blockMtx;      // For the block counter.
	std::mutex m_fidMtx;        // For the featureId.
	std::mutex m_ogrMtx;        // For the OGRLayer.
	uint64_t m_featureId;
	int m_bufSize;
	int m_band;
	bool m_removeHoles;
	bool m_removeDangles;
	bool m_queueUpdate;
	bool m_readFinish;
	bool m_transferFinish;
	std::vector<bool> m_rowsFinished;

public:
	PolyContext(Raster* raster, Status* status, int* block, bool* cancel,
			int bufSize, int band, bool removeHoles, bool removeDangles, Raster* maskRaster = nullptr) :
		m_raster(raster),
		m_maskRaster(maskRaster),
		m_status(status),
		m_block(block),
		m_cancel(cancel),
		m_ds(nullptr),
		m_layer(nullptr),
		m_geomFactory(nullptr),
		m_gctx(nullptr),
		m_featureId(0),
		m_bufSize(bufSize),
		m_band(band),
		m_removeHoles(removeHoles),
		m_removeDangles(removeDangles),
		m_queueUpdate(false),
		m_readFinish(false),
		m_transferFinish(false) {

		// Size the finished array w/ or w/o mask.
		const GridProps& props = m_raster->props();
		m_rowsFinished.resize(getRows(props, m_maskRaster));
		std::fill(m_rowsFinished.begin(), m_rowsFinished.end(), false);
	}

	~PolyContext() {
		if(m_ds)
			GDALClose(m_ds);
	}

	const GeometryFactory::unique_ptr* geomFactory() const {
		return &m_geomFactory;
	}

	void readFinish() {
		m_readFinish = true;
	}

	void transferFinish() {
		m_transferFinish = true;
	}

	void initOutput(const std::string& driver, const std::string& filename, const std::string& layerName, int srid) {
		GDALAllRegister();

		// Create the GEOS context and factory.
		m_gctx = OGRGeometry::createGEOSContext();	
		m_geomFactory = GeometryFactory::create(new PrecisionModel(1.0));

		// Get the vector driver.
		GDALDriver *drv = GetGDALDriverManager()->GetDriverByName(driver.c_str());
		if(!drv)
			g_runerr("Failed to find driver for " << driver << ".");

		// Create an output dataset for the polygons.
		char **dopts = NULL;
		if(Util::lower(driver) == "sqlite")
			dopts = CSLSetNameValue(dopts, "SPATIALITE", "YES");
		m_ds = drv->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, dopts);
		CPLFree(dopts);
		if(!m_ds)
			g_runerr("Failed to create dataset " << filename << ".");

		// Create the layer.
		m_sr.importFromEPSG(srid);
		char **lopts = NULL;
		if(Util::lower(driver) == "sqlite")
			lopts = CSLSetNameValue(lopts, "FORMAT", "SPATIALITE");
		m_layer = m_ds->CreateLayer(layerName.c_str(), &m_sr, wkbMultiPolygon, lopts);
		CPLFree(lopts);
		if(!m_layer) {
			g_runerr("Failed to create layer " << layerName << ".");
		}

		// There's only one field -- an ID.
		OGRFieldDefn field( "id", OFTInteger);
		m_layer->CreateField(&field);

		if(OGRERR_NONE != m_layer->StartTransaction())
			g_runerr("Failed to start transaction.");

	}

	void commitOutput() {
		if(OGRERR_NONE != m_layer->CommitTransaction())
			g_runerr("Failed to commit transation.");
	}

	OGRLayer* layer() {
		return m_layer;
	}

	uint64_t nextFeatureId() {
		std::lock_guard<std::mutex> lk(m_fidMtx);
		return ++m_featureId;
	}

	void notifyTransfer() {
		m_polyMapCond.notify_all();
	}

	void notifyWrite() {
		m_polyQueueCond.notify_all();
	}

	int getCols(const GridProps& props, Raster* mask) {
		if(mask) {
			const GridProps& mprops = mask->props();
			const Bounds& mbounds = mprops.bounds();
			return std::abs(props.toCol(mbounds.maxx()) - props.toCol(mbounds.minx())) + 1;
		} else {
			return props.cols();
		}
	}

	int getRows(const GridProps& props, Raster* mask) {
		if(mask) {
			const GridProps& mprops = mask->props();
			const Bounds& mbounds = mprops.bounds();
			return std::abs(props.toRow(mbounds.maxy()) - props.toRow(mbounds.miny())) + 1;
		} else {
			return props.rows();
		}
	}

	int getCol(const GridProps& props, Raster* mask) {
		if(mask) {
			const GridProps& mprops = mask->props();
			const Bounds& mbounds = mprops.bounds();
			return props.toCol(mprops.resolutionX() > 0 ? mbounds.minx() : mbounds.maxx());
		} else {
			return 0;
		}
	}

	int getRow(const GridProps& props, Raster* mask) {
		if(mask) {
			const GridProps& mprops = mask->props();
			const Bounds& mbounds = mprops.bounds();
			return props.toRow(mprops.resolutionY() ? mbounds.miny() : mbounds.maxy());
		} else {
			return 0;
		}
	}

	void polyReadBlocks() {

		// Get cols, rows an nodata.
		const GridProps& props = m_raster->props();
		int col = getCol(props, m_maskRaster);
		int row = getRow(props, m_maskRaster);
		int cols = getCols(props, m_maskRaster);
		int rows = getRows(props, m_maskRaster);
		uint64_t nd = (uint64_t) props.nodata();

		// The number of blocks int he raster.
		int blocks = (int) rows / m_bufSize + 1;
		// The current block index.
		int b;

		// Buffer for reading raster.
		GridProps gp(props);
		gp.setSize(cols, m_bufSize);
		MemRaster blockBuf(gp, false);

		// Buffer to hold individual row data.
		std::vector<int> rowBuf(cols);

		while(true) {
			{
				std::lock_guard<std::mutex> lk(m_blockMtx);
				b = *m_block;
				++*m_block;
			}

			// If cancelled or no more blocks, quit.
			if(*m_cancel || b >= blocks)
				break;

			// Work out the height of the current buffer. If it's too small quit.
			int bufHeight = g_min(m_bufSize, rows - b * m_bufSize);
			if(bufHeight < 1)
				break;

			if(m_status)
				m_status->update((float) b / blocks);

			// Write into the buffer.
			m_raster->writeTo(blockBuf, cols, bufHeight, col, row + b * m_bufSize, 0, 0, m_band);

			// Read over the rows in the buffer.
			for(int rr = 0; rr < bufHeight; ++rr) {

				if(*m_cancel)
					break;

				// The current overall row index.
				int r = b * m_bufSize + rr;

				// Read into the row buffer from the block.
				blockBuf.getIntRow(rr, m_band, rowBuf.data());

				// Read over the columns in the current row.
				for(int c = 0; c < cols; ++c) {

					// Get the current ID, skip if nodata.
					uint64_t id0 = rowBuf[c];
					if(id0 == nd || id0 == 0)
						continue;

					// Get the coord of one corner of the polygon.
					double x0 = gp.toX(c + col);
					double y0 = gp.toY(r + row);

					// Scan right...
					while(++c <= cols) {

						// Get the next ID or zero if beyond the edge.
						uint64_t id1 = c < cols ? rowBuf[c] : 0;

						// If the ID changes, capture and output the polygon.
						if(id0 > 0 && id1 != id0) {

							// Coord of the other corner.
							double x1 = gp.toX(c + col);
							double y1 = gp.toY(r + row) + gp.resolutionY();

							// Build the geometry.
							CoordinateSequence* seq = m_geomFactory->getCoordinateSequenceFactory()->create(5, 2);
							seq->setAt(Coordinate(x0, y0), 0);
							seq->setAt(Coordinate(x0, y1), 1);
							seq->setAt(Coordinate(x1, y1), 2);
							seq->setAt(Coordinate(x1, y0), 3);
							seq->setAt(Coordinate(x0, y0), 4);
							LinearRing* ring = m_geomFactory->createLinearRing(seq);
							Polygon* geom = m_geomFactory->createPolygon(ring, NULL);
							// Update the polygon list with the new poly.
							{
								std::lock_guard<std::mutex> lk(m_polyMapMtx);
								if(m_polyMap.find(id0) == m_polyMap.end())
									m_polyMap[id0].reset(new Poly(id0));
								m_polyMap[id0]->update(geom, r, r);
							}

							--c; // Back up the counter by one, to start with the new ID.
							break;
						}
						id0 = id1;
					}
				}

				// Update the finished rows.
				m_rowsFinished[r] = true;

				// Notify the transfer threads of an update.
				notifyTransfer();
			}

		} // while

		std::cerr << "read block finished\n";

		notifyTransfer();
	}

	void polyQueueTransfer() {

		while(!*m_cancel && !(m_readFinish && m_polyMap.empty())) {

			std::unique_ptr<Poly> geom;
			{
				// Wait for a wake-up.
				std::unique_lock<std::mutex> lk(m_polyMapMtx);
				while(!*m_cancel && !m_readFinish && m_polyMap.empty())
					m_polyMapCond.wait(lk);

				// Process the finished items if there are any.
				for(auto it = m_polyMap.begin(); it != m_polyMap.end(); ) {
					if(it->second->isRangeFinalized(m_rowsFinished)) {
						geom.swap(it->second);
						it = m_polyMap.erase(it);
					} else {
						++it;
					}
					break;
				}
			}

			if(geom.get()) {
				geom->finalize(m_geomFactory, m_removeHoles, m_removeDangles);
				{
					std::lock_guard<std::mutex> lk(m_polyQueueMtx);
					m_polyQueue.push(std::move(geom));
				}
			}

			notifyWrite();
		}

		std::cerr << "poly transfer queue finished\n";

		notifyWrite();
	}

	std::unique_ptr<Poly> polyFromPath(int id, const std::vector<double>& path) {
		// Build the geometry.
		CoordinateSequence* seq = m_geomFactory->getCoordinateSequenceFactory()->create((size_t) 0, 2);
		std::vector<Coordinate> coords;
		for(size_t i = 0; i < path.size(); i += 2)
			coords.push_back(Coordinate(path[i], path[i+1]));
		seq->add(&coords, false);

		// Get the ring and make a polygon.
		LinearRing* ring = m_geomFactory->createLinearRing(seq);
		Polygon* geom = m_geomFactory->createPolygon(ring, NULL);
		// Create and return the poly.
		std::unique_ptr<Poly> p(new Poly(id));
		p->update(geom, 0, 0);
		return std::move(p);
	}

	void enqueuePoly(std::unique_ptr<Poly>& poly) {
		poly->finalize(m_geomFactory, m_removeHoles, m_removeDangles);
		{
			std::lock_guard<std::mutex> lk(m_polyQueueMtx);
			m_polyQueue.push(std::move(poly));
		}
		notifyWrite();
	}

	void polyWriteQueue() {

		// Use the first thread in the group for writing polys to the DB.
		// The loop runs as long as the queue isn't "finalized" and is not empty.
		while(!*m_cancel && !(m_transferFinish && m_polyQueue.empty())) {

			std::unique_ptr<Poly> p;
			{
				// Wait for a wake-up.
				std::unique_lock<std::mutex> lk(m_polyQueueMtx);
				while(!*m_cancel && !m_transferFinish && m_polyQueue.empty())
					m_polyQueueCond.wait(lk);

				// Grab the element from the queue, if there is one.
				if(!m_polyQueue.empty()) {
					p.swap(m_polyQueue.front());
					m_polyQueue.pop();
				}
			}

			// If no element, loop.
			if(!p.get())
				continue;

			std::cerr << "write\n";

			// Retrieve the unioned geometry and write it.
			const Geometry* g = p->geom();

			{
				std::vector<Polygon*> polys;
				Geometry* p;
				switch(g->getGeometryTypeId()){
				case GEOS_MULTIPOLYGON:
					for(size_t i = 0; i < g->getNumGeometries(); ++i) {
						p = g->getGeometryN(i)->clone();
						if(!p || p->getGeometryTypeId() != GEOS_POLYGON)
							g_runerr("Null or invalid polygon.");
						polys.push_back(dynamic_cast<Polygon*>(p));
					}
					break;
				case GEOS_POLYGON:
					p = g->clone();
					if(!p || p->getGeometryTypeId() != GEOS_POLYGON)
						g_runerr("Null or invalid polygon.");
					polys.push_back(dynamic_cast<Polygon*>(p));
					break;
				default:
					g_runerr("Illegal geometry type: " << g->getGeometryTypeId());
					break;
				}
				g = geos::operation::geounion::CascadedPolygonUnion::CascadedPolygonUnion::Union(&polys);
			}

			OGRGeometry* geom = OGRGeometryFactory::createFromGEOS(m_gctx, (GEOSGeom) g);
			if(!geom)
				g_runerr("Null geometry.");
			OGRFeature feat(m_layer->GetLayerDefn());
			feat.SetGeometry(geom);
			feat.SetField("id", (GIntBig) p->id());
			feat.SetFID(nextFeatureId());
			delete geom;        // The geom is copied by the feature.

			int err;
			{
				std::lock_guard<std::mutex> lk(m_ogrMtx);
				err = m_layer->CreateFeature(&feat);
			}
			if(OGRERR_NONE != err)
				g_runerr("Failed to add geometry.");

			std::cerr << "poly\n";
		}

		std::cerr << "poly write queue finished\n";
	}

};

/*
void Raster::voidfill(const std::string& filename, int band, const std::string& method, const std::string mask) {

	bool hasMask = !mask.empty();

	MemRaster rast(props(), false);
	writeTo(rast);

	std::unique_ptr<MemRaster> mask;
	if(hasMask) {
		Raster rmask(mask);
		mask.reset(new MemRaster(rmask.props(), false));
		rmask.writeTo(*mask);
	}

	if(type == "idw") {
		fillIdw(rast, mask)	
	} else if(type == "nearest") {
		fillNearest(rast, mask);
	}

	GridProps props(props);
	props.setWritable(true);
	Raster output(filename, props());
	rast.writeTo(output);
}
*/

bool computeDirection(int tl, int tr, int bl, int br, int target, int& c, int& r, int& dir) {
	int q = (tl == target && tr == target ? 1 : 0) | (bl == target && br == target ? 2 : 0)
			| (tl == target && bl == target ? 4 : 0) | (tr == target && br == target ? 8 : 0);

	int p = (tl == target ? 1 : 0) | (tr == target ? 2 : 0) | (bl == target ? 4 : 0) | (br == target ? 8 : 0);

	// dir:
	// 0 - none
	// 1 - right
	// 2 - up
	// 4 - left
	// 8 - down

	switch(q) {
	case 1:
	case 9: 
		switch(dir) {
		case 2: ++c; dir = 1; return true;
		case 8: --c; dir = 4; return true;
		}
		break;
	case 2: 
	case 6: 
		switch(dir) {
		case 1: ++r; dir = 8; return true;
		case 4: --r; dir = 2; return true;
		}
		break;
	case 4:
	case 5: ++r; dir = 8; return true;
	case 8:
	case 10: --r; dir = 2; return true;
	case 0:
		switch(p) {
		case 1: --c; dir = 4; return true;
		case 2: --r; dir = 2; return true;
		case 4: ++r; dir = 8; return true;
		case 8: ++c; dir = 1; return true;
		case 6:
		case 9: ++c; dir = 1; return true; // Default to right when crossed.
		}
	}

	/*
	std::cerr << "tl: " << tl << ", tr: " << tr << "\n;";
	std::cerr << "bl: " << bl << ", br: " << br << "\n;";
	std::cerr << "target: " << target << "\n";
	std::cerr << "c : " << c << ", r: " << r << "\n";
	g_runerr("Failed to determine a direction.");
	*/
	return false;
}

void rightCoord(int& c, int& r, int dir) {
	switch(dir) {
	case 1: break;
	case 2: --r; break;
	case 4: --r; --c; break;
	case 8: --c; break;
	}
}

/*
template <class T>
T getXatD(T x, int d) {
	switch(d % 8) {
	case 3:
	case 4:
	case 5: return x + 1;
	case 0:
	case 1:
	case 7: return x - 1;
	default: return x;
	}
}

template <class T>
T getYatD(T y, int d) {
	switch(d % 8) {
	case 1:
	case 2:
	case 3: return y + 1;
	case 5:
	case 6:
	case 7: return y - 1;
	default: return y;
	}
}

uint64_t _tPolyId = 0;
uint64_t _tArcId = 0;

class TPoly {
public:
	uint64_t id;
	std::vector<TArc*> leftArcs;
	std::vector<TArc*> rightArcs;

	TPoly() : id(++_tPolyId) {}

	~TPoly() {
		for(TArc* a : leftArcs)
			delete a;
		for(TArc* a : rightArcs)
			delete a;
	}
}

class TArc {
public:
	uint16_t id;
	int state; // 0 - ignore; 1 - virtual; 2 solid
	std::vector<int> vertices; // Alternating x, y.

	TArc() : id(++_tArcId) {}

	void vertex(int x, int y) {
		vertices.push_back(x);
		vertices.push_back(y);
	}
}

class TArm {
public:
	int pixValue;
	int col;
	TPoly* insidePoly;
	TPoly* abovePoly;
	TPoly* leftPoly;
	TArc vertArm;
	TArc horizArm;

	TArm(int v, int c) : pixValue(v), col(c),
		insidePoly(nullptr), abovePoly(nullptr), leftPoly(nullptr) {}

	~TArm() {
		if(insidePoly) delete insidePoly;
		if(abovePoly) delete abovePoly;
		if(leftPoly) delete leftPoly;
	}
}

void Raster::polygonizeTacet(const std::string& filename, const std::string& layerName,
		const std::string& driver, uint16_t srid, uint16_t band, uint16_t threads,
		bool removeHoles, bool removeDangles, Status *status, bool *cancel) {

	std::cerr << "starting polygonize\n";

	if(!cancel)
		cancel = &s_cancel;

	const GridProps& gp = props();
	int cols = gp.cols();
	int rows = gp.rows();
	int nodata = gp.nodata();

	GridProps rp(props);
	rp.setRows(2);	
	MemRaster rowBuf(rp, false);

	//std::cerr << "initializing writer\n";

	//PolyContext ctx(this, status, 0, cancel, 0, band, removeHoles, removeDangles);
	//ctx.initOutput(driver, filename, layerName, srid);

	//std::thread transT(&PolyContext::polyQueueTransfer, &ctx);
	//std::thread writeT(&PolyContext::polyWriteQueue, &ctx);

	// Organized by ID; each edge is stored as a pair of coordinates, 
	// in i and i+1.
	//std::unordered_map<int, std::vector<double> > edges;

	std::vector<TArm*> arms(cols);

	writeto(rowBuf, 0, 0, cols, 2);

	int v, v0, lastV = nodata;

	// Build the arms using the first row of the buffer. From here
	// on, we use the second row and compare it to the first.
	for(int c = 0; c < cols; ++c) {
		if(c == 0 || (v = getInt(c, 0)) != lastV) {
			TArm* arm = arms[c] = new TArm(v, c);
			arm->vertArm.state = 2;
			arm->vertArm.vertex(c, 0);
			arm->vertArm.vertex(c, 1);
			arm->horizArm.state = c < cols - 1 ? 2 : 0;
			arm->horizArm.vertex(c, 0);
			arm->horizArm.vertex(c + 1, 0);
		}
		lastV = v;
	}

	// For each seed, build a path.
	for(int r = 0; r < rows - 1; ++r) {

		lastV = nodata;
		writeTo(rowBuf, 0, r, cols, 2);

		for(int c = 0; c < cols; ++c) {
			
			v = getInt(c, 1);
			v0 = getInt(c, 0);

			if(c == 0 || (v = getInt(c, 0)) != lastV) {

				TArm* arm;

				if(arms[c]) {
					if(arms[c].pixValue != v) {
						arm = arms[c] = new TArm(v, c);
					} else {
						
					}
				} 
				
				arm->vertArm.state = 2;
				arm->vertArm.vertex(c, 0);
				arm->vertArm.vertex(c, 1);
				arm->horizArm.state = c < cols - 1 ? 2 : 0;
				arm->horizArm.vertex(c, 0);
				arm->horizArm.vertex(c + 1, 0);
			}
			lastV = v;



			int v = getInt()
			if(c == 0 || )

			for(int d = 0; d < 8; d += 2) {

				int cc = getXatD(c, d);
				int rr = getYatD(r, d);

				if(!gp.hasCell(cc, rr) || rast.getInt(cc, rr) != id) {
					 ++edges;
					 cc = c * 2 + 1;
					 rr = r * 2 + 1;
					 erast.setInt(cc, rr, id);
					 erast.setInt(getXatD(cc, d), getYatD(rr, d), -1);
					 erast.setInt(getXatD(cc, d - 1), getYatD(rr, d - 1), -1);
				}

			}
		}
	}


	ctx.readFinish();
	ctx.notifyTransfer();
	transT.join();

	ctx.transferFinish();
	ctx.notifyWrite();
	writeT.join();

	ctx.commitOutput();

}
*/
void Raster::potrace(const std::string& filename, const std::string& layerName,
		const std::string& driver, uint16_t srid, uint16_t band, uint16_t threads,
		bool removeHoles, bool removeDangles, Status *status, bool *cancel) {

	std::cerr << "starting potrace\n";

	std::cerr << "copying raster\n";
	MemRaster rast(props(), true);
	writeTo(rast);

	int cols = props().cols();
	int rows = props().rows();
	
	int ignore = 99999;

	if(!cancel)
		cancel = &s_cancel;

	std::cerr << "initializing writer\n";

	const GridProps& gp = props();

	PolyContext ctx(this, status, 0, cancel, 0, band, removeHoles, removeDangles, nullptr);
	ctx.initOutput(driver, filename, layerName, srid);

	std::thread transT(&PolyContext::polyQueueTransfer, &ctx);
	std::thread writeT(&PolyContext::polyWriteQueue, &ctx);

	std::vector<bool> visited(cols * rows);
	std::fill(visited.begin(), visited.end(), false);

	int lastId = ignore * 2;

	// For each seed, build a path.
	for(int r = 0; r < rows; ++r) {
		for(int c = 0; c < cols; ++c) {

			if(visited[r * cols + c])
				continue;

			std::vector<double> path;

			int id = rast.getInt(c, r);
			if(id == lastId)
				continue;
			lastId = id;

			int c0 = c, r0 = r;
			size_t count = 0;

			std::cerr << "building path for " << id << "\n";

			path.push_back(gp.toX(c));
			path.push_back(gp.toY(r));

			do {
				// Get the neighbouring values. If the coord is off-raster, the value will be nodata.
				int tl = c0 == 0 || r0 == 0 ? ignore : rast.getInt(c0 - 1, r0 - 1);        // Top left.
				int tr = c0 == cols - 1 || r0 == 0 ? ignore : rast.getInt(c0, r0 - 1);     // Top right.
				int bl = c0 == 0 || r0 == rows - 1 ? ignore : rast.getInt(c0 - 1, r0);     // Bottom left.
				int br = c0 == cols - 1 || r0 == rows - 1 ? ignore : rast.getInt(c0, r0);  // Bottom right.
				int dir = 0;

				computeDirection(tl, tr, bl, br, id, c0, r0, dir);

				path.push_back(gp.toX(c0));
				path.push_back(gp.toY(r0));

				// Set the pixel to the right of the direction to true,
				// this will be filled later.
				int vc = c0, vr = r0;
				rightCoord(vc, vr, dir);
				//std::cerr << "vc: " << vc << ", vr: " << vr << "\n";
				size_t vidx = vr * cols + vc;
				if(vidx >= 0 && vidx < visited.size())
					visited[vidx] = true;

				if(++count % 10000 == 0)
					std::cerr << count << "\n";

				std::cerr << "c " << c << ", " << r << ", " << c0 << ", " << r0 << "\n";
			} while(!(c0 == c && r0 == r));

			std::cerr << "last: " << c0 << ", " << r0 << "; " << props().toX(c0) << ", " << props().toY(r0) << "\n";
			std::cerr << " -- " << path[0] << ", " << path[1] << "; " << path[path.size() - 2] << ", " << path[path.size() - 1] << "\n";

			std::unique_ptr<Poly> p = ctx.polyFromPath(id, path);
			ctx.enqueuePoly(p);

			/*
			for(size_t i = 0; i < visited.size(); ++i) {
				if(visited[i]) {
					int vc = i % cols;
					int vr = i / cols;
					if(rast.getInt(vc, vr) == id)
						rast.setInt(vc, vr, ignore);
				}
			}
			*/
		}
	}

	ctx.readFinish();
	ctx.notifyTransfer();
	transT.join();

	ctx.transferFinish();
	ctx.notifyWrite();
	writeT.join();

	ctx.commitOutput();
}

void Raster::polygonize(const std::string& filename, const std::string& layerName,
		const std::string& driver, uint16_t srid, uint16_t band, uint16_t threads,
		int bufSize, bool removeHoles, bool removeDangles, Status *status, bool *cancel, const std::string& mask) {

	if(!m_props.isInt())
		g_runerr("Only integer rasters can be polygonized.");

	// There need to be at least three threads for readin from the raster(1), 
	// writing output (1), and unioning and transferring to the write queue (n - 2).
	if(threads < 3)
		threads = 3;

	// If no cancel pointer is given, use a dummy.
	if(cancel == nullptr)
		cancel = &s_cancel;

	// Counter for the current block.
	int block = 0;

	// If the bufsize given is invalid, use the entire raster height.
	// This can use a lot of memory!
	if(bufSize <= 0 || bufSize > props().rows()) {
		g_warn("Invalid buffer size; using raster height.");
		bufSize = props().rows();
	}

	// Remove the original file; some can't be overwritten directly.
	// This will not take care of any auxillary files (e.g. shapefiles)
	Util::rm(filename);

	Raster* maskRaster = nullptr;
	if(!mask.empty())
		maskRaster = new Raster(mask);

	// Set up the shared context.
	PolyContext ctx(this, status, &block, cancel, bufSize, band, removeHoles, removeDangles, maskRaster);

	// Initialize database.
	ctx.initOutput(driver, filename, layerName, srid);

	// Start the read thread.
	std::thread readT(&PolyContext::polyReadBlocks, &ctx);

	// Start the transfer threads.
	std::thread transferT(&PolyContext::polyQueueTransfer, &ctx);

	// Start the write thread.
	std::thread writeT(&PolyContext::polyWriteQueue, &ctx);

	readT.join();

	ctx.readFinish();
	ctx.notifyTransfer();

	transferT.join();

	ctx.transferFinish();
	ctx.notifyWrite();

	writeT.join();

	g_debug("Writing...");
	if(status)
		status->update(0.99f, "Writing polygons...");

	ctx.commitOutput();

	g_debug("Done");
	if(status)
		status->update(1.0f, "Done.");
}

void Raster::flush() {
	if(m_brow > -1 && m_bcol > -1 && m_props.writable()) {
		GDALRasterBand* bnd = m_ds->GetRasterBand(m_band);
		if(!bnd)
			g_runerr("Failed to find band " << m_band);
		std::lock_guard<std::mutex> lk(m_mtx);
		if(CPLE_None != bnd->WriteBlock(m_bcol, m_brow, getBlock(m_band)))
			g_warn("Failed to write block to " << filename());
	}
	m_ds->FlushCache();
}

Raster::~Raster() {
	flush();
	for(auto& item : m_blocks)
		free(item.second);
	if(m_ds)
		GDALClose(m_ds);
}
