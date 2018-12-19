#include <queue>
#include <string>
#include <fstream>
#include <atomic>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <condition_variable>
#include <chrono>
#include <memory>

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

#include "geo.hpp"
#include "util.hpp"
#include "raster.hpp"
#include "db.hpp"

using namespace geo::util;
using namespace geo::raster;

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
		m_compress(false),
		m_bigTiff(false),
		m_type(DataType::None) {	// The data type.
}

void GridProps::setCompress(bool compress) {
	m_compress = compress;
};

bool GridProps::compress() const {
	return m_compress;
}

void GridProps::setBigTiff(bool bigTiff) {
	m_bigTiff = bigTiff;
}

bool GridProps::bigTiff() const {
	return m_bigTiff;
}

std::string GridProps::driver() const {
	return m_driver;
}

void GridProps::setDriver(const std::string& name) {
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

size_t GridProps::size() const {
	return (size_t) cols() * rows();
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

void GridProps::setBounds(const Bounds& bounds) {
	m_trans[0] = m_trans[1] > 0 ? bounds.minx() : bounds.maxx();
	m_trans[3] = m_trans[5] > 0 ? bounds.miny() : bounds.maxy();
	m_cols = (int) std::ceil(bounds.width() / std::abs(m_trans[1]));
	m_rows = (int) std::ceil(bounds.height() / std::abs(m_trans[5]));
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

void GridProps::setProjection(const std::string& proj) {
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
	std::lock_guard<std::mutex> lk(dest.mutex());
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
	const GridProps& p = m_source.props();
	return m_curCol < p.cols() && m_curRow < p.rows();
}

int TileIterator::count() const {
	const GridProps& p = m_source.props();
	return g_max(1, (int) std::ceil((float) p.cols() / m_cols) * (int) std::ceil((float) p.rows() / m_rows));
}

Tile* TileIterator::next() {

	int col, row, srcCol, srcRow, dstCol, dstRow, cols, rows;
	const GridProps& props = m_source.props();

	if(!(m_curCol < props.cols() && m_curRow < props.rows()))
		return nullptr;

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

	GridProps p(props);
	p.setSize(m_cols + m_buffer * 2, m_rows + m_buffer * 2);
	std::unique_ptr<MemRaster> tile(new MemRaster(p));

	{
		std::lock_guard<std::mutex> lk(m_source.mutex());
		m_source.writeTo(*tile, cols, rows, srcCol, srcRow, dstCol, dstRow, m_band, 1);
	}

	return new Tile(tile.release(), &m_source, m_cols, m_rows, col, row,
			m_buffer, srcCol, srcRow, dstCol, dstRow, m_band, props.writable());
}

Tile* TileIterator::create(Tile &tpl) {
	GridProps p(m_source.props());
	p.setSize(m_cols + m_buffer * 2, m_rows + m_buffer * 2);
	std::unique_ptr<MemRaster> tile(new MemRaster(p));

	{
		std::lock_guard<std::mutex> lk(m_source.mutex());
		m_source.writeTo(*tile, tpl.m_cols, tpl.m_rows, tpl.m_srcCol, tpl.m_srcRow, tpl.m_dstCol, tpl.m_dstRow, m_band, 1);
	}

	return new Tile(tile.release(), &m_source, m_cols, m_rows, tpl.m_col, tpl.m_row,
			m_buffer, tpl.m_srcCol, tpl.m_srcRow, tpl.m_dstCol, tpl.m_dstRow, m_band, p.writable());
}

TileIterator::~TileIterator() {
}


// Implementations for the Grid class

TileIterator Grid::iterator(int cols, int rows, int buffer, int band) {
	return TileIterator(*this, cols, rows, buffer, band);
}

void Grid::gaussianWeights(double *weights, int size, double sigma, double mean) {
	// If size is an even number, bump it up.
	if (size % 2 == 0) {
		++size;
		g_warn("Gaussian kernel size must be an odd number >=3. Bumping up to " << size);
	}
	for (int r = 0; r < size; ++r) {
		for (int c = 0; c < size; ++c) {
			int x = size / 2 - c;
			int y = size / 2 - r;
			weights[r * size + c] = (1 / (2 * G_PI * sigma * sigma)) * std::pow(G_E, -((x * x + y * y) / (2.0 * sigma * sigma)));
		}
	}
}

GridStats Grid::stats(int band) {
	GridStats st;
	const GridProps& gp = props();
	double nodata = gp.nodata();
	double v, m = 0, s = 0;
	int k = 1;
	st.sum = 0;
	st.count = 0;
	st.min = G_DBL_MAX_POS;
	st.max = G_DBL_MAX_NEG;
	// Welford's method for variance.
	int rows = gp.rows();
	int cols = gp.cols();
	for(int row = 0; row < rows; ++row) {
		for(int col = 0; col < cols; ++col) {
			if ((v = getFloat(col, row, band)) != nodata) {
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
	int rows = gp.rows();
	int cols = gp.cols();
	for(int row = 0; row < rows; ++row) {
		for(int col = 0; col < cols; ++col) {
			if ((v = getFloat(col, row, band)) != nodata && !std::isnan(v) && v < G_DBL_MAX_POS) {
				setFloat(col, row, ((v - mean) / stdDev), band);
			} else {
				setFloat(col, row, nodata, band);
			}
		}
	}
}

void Grid::logNormalize(int band) {
	GridStats st = stats(band);
	const GridProps& gp = props();
	double n = st.min;
	double x = st.max;
	double e = std::exp(1.0) - 1.0;
	int rows = gp.rows();
	int cols = gp.cols();
	for(int row = 0; row < rows; ++row) {
		for(int col = 0; col < cols; ++col) {
				setFloat(col, row, std::log(1.0 + e * (getFloat(col, row, band) - n) / (x - n)), band);
		}
	}
}

void Grid::convert(Grid& g, int srcBand, int dstBand) {
	const GridProps& gp = props();
	int rows = gp.rows();
	int cols = gp.cols();
	if(g.props().isInt()) {
		for(int row = 0; row < rows; ++row) {
			for(int col = 0; col < cols; ++col) {
				g.setInt(col, row, getInt(col, row, srcBand), dstBand);
			}
		}
	} else {
		for(int row = 0; row < rows; ++row) {
			for(int col = 0; col < cols; ++col) {
				g.setFloat(col, row, getFloat(col, row, srcBand), dstBand);
			}
		}
	}
}

void Grid::voidFillIDW(const std::string& filename, int band, const std::string& mask,
		int maskBand, double radius, int count, double exp) {

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

    TargetFillOperator<double, double> op1(&input, band, band, nodata, 99999);
    TargetFillOperator<double, double> op2(&input, band, &output, band, 99999, nodata);
    TargetFillOperator<double, double> op3(&input, band, band, 99999, 99998);
    int outminc, outminr, outmaxc, outmaxr;

    for (int r = 0; r < rows; ++r) {
    	if(r % 100 == 0) {
    		std::cerr << "Row " << r << " of " << rows << "\n";
    	}
		for (int c = 0; c < cols; ++c) {

			v = input.getFloat(c, r, band);
			if(v == 99998) {

				//output.setFloat(c, r, nodata);

			} else if (v != nodata) {

				output.setFloat(c, r, v, 1);

			} else if(!holesOnly) {

				double dp, a = 0, b = 0;
				int cnt = 0;
				for(int r0 = g_max(0, r - maxDist); r0 < g_min(rows, r + maxDist + 1); ++r0) {
					for(int c0 = g_max(0, c - maxDist); c0 < g_min(cols, c + maxDist + 1); ++c0) {
						if((c0 == c && r0 == r) || (d = g_sq(c0 - c) + g_sq(r0 - r)) > maxDist ||
								(v = input.getFloat(c0, r0, band)) == nodata)
							continue;
						dp = 1.0 / std::pow(d, exp);
						a += dp * v;
						b += dp;
						++cnt;
					}
				}
				output.setFloat(c, r, cnt ? (a / b) : nodata, 1);

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
						v = input.getFloat(c0, r0, band);
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
					output.setFloat(nc, nr, cnt ? (a / b) : nodata, 1);
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

/**
 * Parallel function for smoothing, called by smooth().
 * @param iter A pointer to a TileIterator.
 * @param smoothed The grid to contain smoothed output.
 * @param status The status object.
 * @param size The size of the kernel.
 * @param nodata The nodata value from the original raster.
 * @param weights A list of Gaussian weights.
 * @param rmtx The read mutex.
 * @param wmtx The write mutex.
 * @param cancel If set to true during operation, cancels the operation.
 * @param ex If the function terminates with an exception, this pointer should point to it.
 */
void _smooth(TileIterator* iter, Grid* smoothed,
		int size, double nodata, double* weights,
		std::mutex* rmtx, std::mutex* wmtx, bool* cancel, Status* status, int* curTile,
		std::exception_ptr* p) {

	try {
		std::unique_ptr<Tile> tile;
		int tileCount = iter->count();

		while(!*cancel) {

			{
				std::lock_guard<std::mutex> lk(*rmtx);
				tile.reset(iter->next());
				if(!tile.get())
					return;
				++*curTile;
			}

			// Copy the tile to a buffer, process the buffer and write back to the grid.
			Grid& grid = tile->grid();
			const GridProps& props = grid.props();
			MemRaster buf(props);
			grid.writeTo(buf);

			// Process the entire block, even the buffer parts.
			for (int r = 0; !*cancel && r < props.rows() - size; ++r) {
				for (int c = 0; !*cancel && c < props.cols() - size; ++c) {
					double v, t = 0.0;
					bool foundNodata = false;
					for (int gr = 0; gr < size; ++gr) {
						for (int gc = 0; gc < size; ++gc) {
							v = buf.getFloat(c + gc, r + gr, 1);
							if (v == nodata) {
								foundNodata = true;
								break;
							} else {
								t += weights[gr * size + gc] * v;
							}
						}
						if(foundNodata) break;
					}
					if (!foundNodata)
						grid.setFloat(c + size / 2, r + size / 2, t, 1);
				}
			}

			if(*cancel)
				break;

			{
				std::lock_guard<std::mutex> lk(*wmtx);
				tile->writeTo(*smoothed);
			}

			status->update(g_min(0.99f, 0.2f + (float) *curTile / tileCount * 0.97f));
		}
	} catch(const std::exception& e) {
		*p = std::current_exception();
	}
}

void Grid::smooth(Grid& smoothed, double sigma, int size, int band,
		bool& cancel, Status& status) {

	const GridProps& gp = props();

	status.update(0.01f);

	if (sigma <= 0)
		g_argerr("Sigma must be > 0.");
	if (size < 3)
		g_argerr("Kernel size must be 3 or larger.");
	if (size % 2 == 0) {
		g_warn("Kernel size must be odd. Rounding up.");
		size++;
	}

	// Compute the weights for Gaussian smoothing.
	Buffer weightsBuf(size * size * getTypeSize(DataType::Float64));
	double* weights = (double*) weightsBuf.buf;
	Grid::gaussianWeights(weights, size, sigma);

	double nodata = gp.nodata();

	status.update(0.02f);

	TileIterator iter = iterator(512, 512, size, band);
	std::mutex wmtx; // write mutex
	std::mutex rmtx; // read mutex
	int curTile = 0;

	// Run the smoothing jobs.
	int numThreads = 8;
	std::vector<std::thread> threads;
	std::vector<std::exception_ptr> exceptions(numThreads);
	for(int i = 0; i < numThreads; ++i) {
		threads.emplace_back(_smooth, &iter, &smoothed, size, nodata, weights,
				&rmtx, &wmtx, &cancel, &status, &curTile, &exceptions[i]);
	}

	// Wait for jobs to complete.
	for(int i = 0; i < numThreads; ++i) {
		if(threads[i].joinable())
			threads[i].join();
	}

	// Check if any exceptions were trapped. Raise the first one.
	for(int i = 0; i < numThreads; ++i) {
		if(exceptions[i])
			std::rethrow_exception(exceptions[i]);
	}

	status.update(1.0);
}

// Implementations for MemRaster

std::mutex& MemRaster::mutex() {
	return m_mtx;
}

bool MemRaster::mmapped() const {
	return m_mmapped;
}

void MemRaster::checkInit() const {
	if (m_grid == nullptr)
		g_runerr("This instance has not been initialized.");
}

MemRaster::MemRaster() : 
	m_mmapped(false),
	m_grid(nullptr),
	m_mappedFile(nullptr) {
}

MemRaster::MemRaster(const MemRaster& other) :
		MemRaster() {
	init(other.props(), other.mmapped());
}

MemRaster::MemRaster(const GridProps &props, bool mapped) :
		MemRaster() {
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
	std::lock_guard<std::mutex> lk(m_freeMtx);
	if (m_grid) {
		if (m_mmapped) {
			delete m_mappedFile;
			m_mappedFile = nullptr;
			m_grid = nullptr;
		} else {
			free(m_grid);
			m_grid = nullptr;
		}
	}
}

void MemRaster::init(const GridProps& pr, bool mapped) {
	std::lock_guard<std::mutex> lk(m_initMtx);
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
			m_mappedFile = new MappedFile(size, true);
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
		for (size_t i = 0; i < size; i += chunk) {
			std::memcpy(grid, buffer.buf, chunk);
			grid += chunk;
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
		for (size_t i = 0; i < size; i += chunk) {
			std::memcpy(grid, buffer.buf, chunk);
			grid += chunk;
		}
	}
}

double MemRaster::getFloat(int col, int row, int band) {
	checkInit();
	if(m_props.isInt()) {
		return (double) getInt(col, row, band);
	} else {
		size_t idx = (size_t) row * m_props.cols() + col;
		if (idx < 0 || idx >= m_props.size())
			g_argerr("Index out of bounds: " << idx << "; size: " << m_props.size());
		return *(((double*) m_grid) + idx);
	}
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

int MemRaster::getInt(int col, int row, int band) {
	checkInit();
	if(m_props.isInt()) {
		size_t idx = (size_t) row * m_props.cols() + col;
		if (idx < 0 || idx >= m_props.size())
			g_argerr("Index out of bounds: " << idx << "; size: " << m_props.size());
		return *(((int*) m_grid) + idx);
	} else {
		return (int) getFloat(col, row, band);
	}
}

void MemRaster::setFloat(int col, int row, double value, int band) {
	checkInit();
	if(m_props.isInt()) {
		setInt(col, row, (int) value, band);
	} else {
		size_t idx = (size_t) row * m_props.cols() + col;
		size_t ps = m_props.size();
		if (idx >= ps)
			g_argerr("Index out of bounds: " << idx << "; size: " << m_props.size()
							<< "; value: " << value << "; col: " << (idx % m_props.cols())
							<< "; row: " << (idx / m_props.cols()));
		*(((double *) m_grid) + idx) = value;
	}
}

void MemRaster::setInt(int col, int row, int value, int band) {
	checkInit();
	if(m_props.isInt()) {
		size_t idx = (size_t) row * m_props.cols() + col;
		size_t ps = m_props.size();
		if (idx >= ps)
			g_argerr("Index out of bounds: " << idx << "; size: " << m_props.size()
							<< "; value: " << value << "; col: " << (idx % m_props.cols())
							<< "; row: " << (idx / m_props.cols()));
		*(((int *) m_grid) + idx) = value;
	} else {
		setFloat(col, row, (double) value, band);
	}
}

void MemRaster::toMatrix(
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mtx, int band) {
	int cols = m_props.cols();
	int rows = m_props.rows();
	for (int r = 1; r < rows; ++r) {
		for (int c = 0; c < cols; ++c)
			mtx(r, c) = getFloat(c, r, band);
	}
}

void MemRaster::fromMatrix(
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mtx, int band) {
	int cols = m_props.cols();
	int rows = m_props.rows();
	for (int r = 1; r < rows; ++r) {
		for (int c = 0; c < cols; ++c)
			setFloat(c, r, (double) mtx(r, c), band);
	}
}

void MemRaster::writeToRaster(Raster& grd,
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

	GDALRasterBand* band = grd.m_ds->GetRasterBand(dstBand);

	char* input  = (char*) grid();

	for(int r = 0; r < rows; ++r) {
		if(CPLE_None != band->RasterIO(GF_Write, dstCol, dstRow + r, cols, 1, input + ((srcRow + r) * gcols + srcCol) * typeSize,
				cols, 1, gtype, 0, 0, 0))
			g_runerr("Failed to write to: " << grd.filename());
	}

	band->FlushCache();
}

void MemRaster::writeToMemRaster(MemRaster& grd,
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
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c)
					grd.setInt(c + dstCol, r + dstRow, getInt(c + srcCol, r + srcRow, srcBand), dstBand);
			}
		} else {
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c)
					grd.setInt(c + dstCol, r + dstRow, (int) getFloat(c + srcCol, r + srcRow, srcBand), dstBand);
			}
		}
	} else {
		if(m_props.isInt()) {
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c)
					grd.setFloat(c + dstCol, r + dstRow, (double) getInt(c + srcCol, r + srcRow, srcBand), dstBand);
			}
		} else {
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c)
					grd.setFloat(c + dstCol, r + dstRow, getFloat(c + srcCol, r + srcRow, srcBand), dstBand);
			}
		}
	}
}

void MemRaster::writeTo(Grid& grd,
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

std::mutex& Raster::mutex() {
	return m_mtx;
}

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
					for(const std::string& item : lst)
						extensions[desc].insert("." + Util::lower(item));
				}

			}
		}
	}
	return extensions;
}

std::map<std::string, std::string> Raster::drivers() {
	std::vector<std::string> f;
	return drivers(f);
}

std::map<std::string, std::string> Raster::drivers(const std::vector<std::string>& filter) {
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
				bool found = true;
				if(!filter.empty()) {
					found = false;
					for(const std::string& f : filter) {
						if(f == desc) {
							found = true;
							break;
						}
					}
				}
				if(found)
					drivers[desc] = name;
			}
		}
	}
	return drivers;
}

std::string Raster::getDriverForFilename(const std::string& filename) {
	std::string ext = Util::extension(filename);
	std::map<std::string, std::set<std::string> > drivers = extensions();
	std::string result;
	for(const auto& it : drivers) {
		if(it.second.find(ext) != it.second.end())
			result = it.first;
	}
	return result;
}

Raster::Raster(const std::string& filename, const GridProps& props) :
		m_ds(nullptr),
		m_bcols(0), m_brows(0),
		m_bcol(-1), m_brow(-1),
		m_bband(0),
		m_dirty(false), 
		m_type(GDT_Unknown) {

	if (props.resolutionX() == 0 || props.resolutionY() == 0)
		g_argerr("Resolution must not be zero.");
	if (props.cols() <= 0 || props.rows() <= 0)
		g_argerr("Columns and rows must be larger than zero.");
	if (filename.empty())
		g_argerr("Filename must be given.");

	m_props = props;
	m_filename = filename;

	// Create GDAL dataset.
	char **opts = NULL;
	if(m_props.compress()) {
		opts = CSLSetNameValue(opts, "COMPRESS", "LZW");
		opts = CSLSetNameValue(opts, "PREDICTOR", "2");
	}
	if(m_props.bigTiff())
		opts = CSLSetNameValue(opts, "BIGTIFF", "YES");
	opts = CSLSetNameValue(opts, "INTERLEAVE", "BAND");

	GDALAllRegister();
	std::string drvName = m_props.driver();
	if(drvName.empty())
		drvName = getDriverForFilename(m_filename);
	if(drvName.empty())
		g_runerr("Couldn't find driver for: " << m_filename)
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

Raster::Raster(const std::string& filename, bool writable) :
		m_ds(nullptr),
		m_bcols(0), m_brows(0),
		m_bcol(-1), m_brow(-1),
		m_bband(0),
		m_dirty(false),
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
	if(m_blocks.find(band) == m_blocks.end())
		m_blocks[band] = malloc(m_bcols * m_brows * getTypeSize(m_props.dataType()));
	return m_blocks[band];
}

GDALDataset* Raster::ds() const {
	return m_ds;
}

const GridProps& Raster::props() const {
	return m_props;
}

DataType Raster::getFileDataType(const std::string& filename) {
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
	bnd->Fill((int) value);
}

void Raster::fillFloat(double value, int band) {
	GDALRasterBand *bnd = m_ds->GetRasterBand(band);
	bnd->Fill(value);
}

void Raster::writeToRaster(Raster& grd,
			int cols, int rows,
			int srcCol, int srcRow,
			int dstCol, int dstRow,
			int srcBand, int dstBand) {

	if(m_dirty) {
		GDALRasterBand* bnd = m_ds->GetRasterBand(m_bband);
		if(!bnd)
			g_runerr("Failed to find band " << m_bband);
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
	if(CPLE_None != srcBnd->RasterIO(GF_Read, srcCol, srcRow, cols, rows,
			buf.buf, cols, rows, gtype, 0, 0, 0))
		g_runerr("Failed to read from: " << grd.filename());
	if(CPLE_None != dstBnd->RasterIO(GF_Write, dstCol, dstRow, cols, rows,
			buf.buf, cols, rows, gtype, 0, 0, 0))
		g_runerr("Failed to write to: " << filename());
}

void Raster::writeToMemRaster(MemRaster& grd,
		int cols, int rows,
		int srcCol, int srcRow,
		int dstCol, int dstRow,
		int srcBand, int dstBand) {

	if(m_dirty) {
		GDALRasterBand* bnd = m_ds->GetRasterBand(m_bband);
		if(!bnd)
			g_runerr("Failed to find band " << bnd);
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
	size_t typeSize = getTypeSize(gp.dataType());
	size_t gcols = gp.cols();

	Buffer buf(cols * rows * typeSize);
	GDALRasterBand *band = m_ds->GetRasterBand(srcBand);
	if(!band)
		g_runerr("Failed to find band " << srcBand);
	if(CPLE_None != band->RasterIO(GF_Read, srcCol, srcRow, cols, rows, buf.buf,
			cols, rows, gtype, 0, 0, 0))
		g_runerr("Failed to read from: " << filename());

	char* output = (char*) grd.grid();
	char* input  = (char*) buf.buf;
	size_t len = cols * typeSize;
	for(int r = 0; r < rows; ++r) {
		size_t doff = ((dstRow + r) * gcols + dstCol) * typeSize;
		size_t soff = (size_t) r * cols * typeSize;
		std::memcpy(output + doff , input + soff, len);
	}

}

void Raster::writeTo(Grid& grd,
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
		if(m_dirty)
			flush();
		GDALRasterBand* bnd = m_ds->GetRasterBand(band);
		if(!bnd)
			g_argerr("Failed to find band " << band);
		if(CPLE_None != bnd->ReadBlock(bcol, brow, getBlock(band)))
			g_runerr("Failed to read from: " << filename());
		m_bcol = bcol;
		m_brow = brow;
		m_bband = band;
	}
	int idx = (row % m_brows) * m_bcols + (col % m_bcols);
	double v = 0;
	readFromBlock(getBlock(m_bband), getGDType(), &v, idx);
	return v;
}

/*
double Raster::getFloat(size_t idx, int band) {
	int cols = m_props.cols();
	return getFloat((int) (idx % cols), (int) (idx / cols), band);
}
*/

double Raster::getFloat(double x, double y, int band) {
	return getFloat(m_props.toCol(x), m_props.toRow(y), band);
}

int Raster::getInt(int col, int row, int band) {
	int bcol = col / m_bcols;
	int brow = row / m_brows;
	if(bcol != m_bcol || brow != m_brow || band != m_bband) { // TODO: No effective cacheing if the band changes.
		if(m_dirty)
			flush();
		GDALRasterBand* bnd = m_ds->GetRasterBand(band);
		if(!bnd)
			g_argerr("Failed to find band " << band);
		if(CPLE_None != bnd->ReadBlock(bcol, brow, getBlock(band)))
			g_runerr("Failed to read from: " << filename());
		m_bcol = bcol;
		m_brow = brow;
		m_bband = band;
	}
	int idx = (row % m_brows) * m_bcols + (col % m_bcols);
	int v = 0;
	readFromBlock(getBlock(m_bband), getGDType(), &v, idx);
	return v;
}

/*
int Raster::getInt(size_t idx, int band) {
	int cols = m_props.cols();
	return getInt((int) (idx % cols), (int) (idx / cols), band);
}
*/

int Raster::getInt(double x, double y, int band) {
	return getInt(m_props.toCol(x), m_props.toRow(y), band);
}

/*
void Raster::setInt(size_t idx, int v, int band) {
	int cols = m_props.cols();
	setInt((int) (idx % cols), (int) (idx / cols), v, band);
}
*/

void Raster::setInt(double x, double y, int v, int band) {
	setInt(m_props.toCol(x), m_props.toRow(y), v, band);
}

void Raster::setFloat(int col, int row, double v, int band) {
	if (!m_props.writable())
		g_runerr("This raster is not writable.");
	int bcol = col / m_bcols;
	int brow = row / m_brows;
	if(bcol != m_bcol || brow != m_brow || band != m_bband) {
		if(m_dirty)
			flush();
		GDALRasterBand* bnd = m_ds->GetRasterBand(band);
		if(!bnd)
			g_argerr("Failed to find band " << band);
		if(CPLE_None != bnd->ReadBlock(bcol, brow, getBlock(band)))
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
		if(m_dirty)
			flush();
		GDALRasterBand* bnd = m_ds->GetRasterBand(band);
		if(!bnd)
			g_argerr("Failed to find band " << band);
		if(CPLE_None != bnd->ReadBlock(bcol, brow, getBlock(band)))
			g_runerr("Failed to read from: " << filename());
		m_bcol = bcol;
		m_brow = brow;
		m_bband = band;
	}
	int idx = (row % m_brows) * m_bcols + (col % m_bcols);
	writeToBlock(getBlock(m_bband), getGDType(), v, idx);
	m_dirty = true;
}

/*
void Raster::setFloat(size_t idx, double v, int band) {
	int cols = m_props.cols();
	setFloat((int) (idx % cols), (int) (idx / cols), v, band);
}
*/

void Raster::setFloat(double x, double y, double v, int band) {
	setFloat(m_props.toCol(x), m_props.toRow(y), v, band);
}

void Raster::flush() {
	if(m_dirty && m_props.writable()) {
		GDALRasterBand* bnd = m_ds->GetRasterBand(m_bband);
		if(!bnd)
			g_runerr("Failed to find band " << m_bband);
		if(CPLE_None != bnd->WriteBlock(m_bcol, m_brow, getBlock(m_bband)))
			g_warn("Failed to write block to " << filename());
	}
	m_ds->FlushCache();
	m_dirty = false;
}

Raster::~Raster() {
	flush();
	for(auto& item : m_blocks)
		free(item.second);
	if(m_ds)
		GDALClose(m_ds);
}
