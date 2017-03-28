#include <queue>
#include <string>
#include <fstream>
#include <atomic>
#include <unordered_map>
#include <unordered_set>

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

#include <geos/geom/GeometryFactory.h>
#include <geos/geom/CoordinateSequenceFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/CoordinateSequence.h>

#include "omp.h"

#include "geo.hpp"
#include "util.hpp"
#include "raster.hpp"
#include "db.hpp"

using namespace geo::util;
using namespace geo::raster;
using namespace geos::geom;


bool _cancel = false;

namespace geo {

	namespace raster {

		namespace util {

			class Poly {
			public:
				int thread;
				int minRow, maxRow;
				Geometry* poly;
				const GeometryFactory* fact;

				Poly(Geometry* poly, int thread, int row, const GeometryFactory* fact) :
						thread(thread), minRow(row), maxRow(row), poly(nullptr), fact(fact) {
					this->poly = poly;
				}

				void update(Geometry* u, int minRow, int maxRow) {
					Geometry* p = this->poly;
					poly = p->Union(u);
					delete p;
					update(minRow, maxRow);
				}

				void update(int minRow, int maxRow) {
					if(minRow < this->minRow) this->minRow = minRow;
					if(maxRow > this->maxRow) this->maxRow = maxRow;
				}

				// Returns true if the range of rows given by start and end was finalized
				// within the given thread, or by a thread whose block is completed.
				// The checked range includes one row above and one below,
				// which is required to guarantee that a polygon is completed.
				bool isRangeFinalized(const std::vector<int> &finalRows) const {
					int start = g_max(minRow - 1, 0);
					int end = g_min(maxRow + 2, (int) finalRows.size());
					for (int i = start; i < end; ++i) {
						if (finalRows[i] != -2 && finalRows[i] != thread)
							return false;
					}
					return true;
				}

				Geometry* polyAsMulti() {
					if(GEOS_MULTIPOLYGON != poly->getGeometryTypeId()) {
						std::vector<geos::geom::Geometry*>* geoms = new std::vector<Geometry*>();
						geoms->push_back(poly);
						return fact->createMultiPolygon(geoms);
					} else {
						return poly;
					}
				}

				~Poly() {
					if(poly)
						delete poly;
				}

			};
			// Keeps track of polygonization progress and provides a cancellation mechanism.
			class PolyProgressData {
			public:
				Status *status;
				bool *cancel;
				bool _cancel;

				PolyProgressData(Status *status, bool *cancel) :
					status(status), cancel(cancel), _cancel(false) {
					if (cancel == nullptr)
						cancel = &_cancel;
				}
			};

			// Callback for polygonization.
			int polyProgress(double dfComplete, const char *pszMessage, void *pProgressArg) {
				PolyProgressData *data = (PolyProgressData*) pProgressArg;
				if (data) {
					if(data->status)
						data->status->update((float) dfComplete);
					return *(data->cancel) ? 0 : 1;
				}
				return 1;
			}

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
		m_bands(0),           		// The number of bands
		m_writable(false),			// True if the grid is writable
		m_nodata(0),
		m_type(DataType::None) {	// The data type.
}


void GridProps::setMappedPath(const std::string& path) {
	m_mappedPath = path;
}

std::string GridProps::mappedPath() const {
	return m_mappedPath;
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

Bounds GridProps::bounds() const {
	double x0 = m_trans[0];
	double y0 = m_trans[3];
	double x1 = x0 + m_trans[1] * m_cols;
	double y1 = y0 + m_trans[5] * m_rows;
	return Bounds(g_min(x0, x1), g_min(y0, y1), g_max(x0, x1), g_max(y0, y1));
}

void GridProps::setNoData(double nodata) {
	m_nodata = nodata;
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


Grid::Grid() {
}

// Implementations forthe Grid class

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
	GridStats st = stats();
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
	GridStats st = stats();
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

void Grid::voidFillIDW(double radius, int count, double exp, int band) {

	if(!props().isFloat())
		g_runerr("IDW fill only implemented for float rasters.");

	if (radius <= 0.0)
		throw std::invalid_argument("Radius must be larger than 0.");

	if (count <= 0)
		throw std::invalid_argument("Count must be larger than 0.");

	if (exp <= 0.0)
		throw std::invalid_argument("Exponent must be larger than 0.");

	MemRaster tmp(props());
	double nodata = props().nodata();
	int rows = props().rows();
	int cols = props().cols();
	for (int r = 0; r < rows; ++r) {
		for (int c = 0; c < cols; ++c) {
			if (getFloat(c, r) != nodata)
				continue;
			double rad = radius;
			bool found = false;
			do {
				double d = g_sq(rad);
				double a = 0.0;
				double b = 0.0;
				int cnt = 0;
				for (int r0 = (int) g_max(0, r - rad); r0 < (int) g_min(rows, r + rad + 1); ++r0) {
					for (int c0 = (int) g_max(0, c - rad); c0 < (int) g_min(cols, c + rad + 1); ++c0) {
						double d0 = g_sq((double) c0 - c) + g_sq((double) r0 - r);
						if (d0 <= d && getFloat(c0, r0, band) != nodata) {
							double dp = 1.0 / std::pow(d0, exp);
							a += dp * getFloat(c0, r0, band);
							b += dp;
							++cnt;
						}
					}
				}

				if (cnt >= count) {
					tmp.setFloat(c, r, (a / b), band);
					found = true;
					break;
				}

				rad += 1.0;

			} while (rad < g_min(cols, rows));

			if (!found)
				g_warn("Pixel not filled at " << c << "," << r << ". Consider larger radius or smaller count.");
		}
	}
	writeTo(tmp);
}

using namespace geo::raster::util;

void Grid::smooth(Grid &smoothed, double sigma, int size, int band,
		Callbacks *status, bool *cancel) {
	const GridProps& gp = props();
	if (!cancel)
		cancel = &_cancel;
	if (status)
		status->stepCallback(0.01f);
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
	int cols = gp.cols();
	int rows = gp.rows();
	std::atomic<int> curRow(0);

	if (status)
		status->stepCallback(0.02f);

	int bufSize = g_max(size, 256);

	#pragma omp parallel
	{

		GridProps pr = GridProps(gp);
		pr.setSize(cols, bufSize + size);
		pr.setNoData(nd);
		pr.setDataType(DataType::Float64);

		MemRaster buf(pr, false);
		MemRaster smooth(pr, false);

		#pragma omp for
		for (int i = 0; i < (rows - size) / bufSize + 1; ++i) {
			if (*cancel) continue;

			buf.fillFloat(nd);
			smooth.fillFloat(nd);

			int b = i * bufSize;
			int readOffset = b > 0 ? b - size / 2 : 0;  // If this is the first row, read from zero, otherwise -(size / 2)
			int writeOffset = b > 0 ? 0 : size / 2;     // If this is the first row, write to (size / 2), otherwise 0.
			#pragma omp critical(__smooth_read)
			writeTo(buf, pr.cols(), pr.rows(), 0, readOffset, 0, writeOffset, band);

			// Process the entire block, even the buffer parts.
			for (int r = 0; r < pr.rows() - size; ++r) {
				for (int c = 0; c < pr.cols() - size; ++c) {
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
						smooth.setFloat(c + size / 2, r + size / 2, t);
				}
			}

			if (status) {
				curRow += bufSize;
				status->stepCallback(0.2f + (float) curRow / rows * 0.97f);
			}

			#pragma omp critical(__smooth_write)
			smooth.writeTo(smoothed, pr.cols(), g_min(bufSize, smoothed.props().rows() - b), 0, size / 2, 0, b); // Always write to b and read from (size / 2)
		}
	}

	if (*cancel)
		return;

	if (status)
		status->stepCallback(1.0);
}

Grid::~Grid() {
}

// Implementations for MemRaster

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
			m_mappedFile.release();
			m_grid = nullptr;
		} else {
			free(m_grid);
			m_grid = nullptr;
		}
	}
}

void MemRaster::init(const GridProps &pr, bool mapped) {
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
		size_t size = typeSize * pr.cols() * pr.rows();
		m_mappedFile.release();
		if (mapped) {
			const std::string filename = Util::tmpFile(pr.mappedPath());
			m_mappedFile = Util::mapFile(filename, size);
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
		size_t chunk = m_mmapped ? m_mappedFile->pageSize() * 1024 : 2048 * 1024;
		size_t size = m_props.size() * sizeof(double);
		size_t last = size % chunk;
		Buffer buf(chunk);
		for(size_t i = 0; i < chunk / sizeof(double); ++i)
			*((double *) buf.buf + i) = value;
		for(uint64_t i = 0; i < size; i += chunk)
			std::memcpy((char *) m_grid + i, buf.buf, g_min(chunk, last));
	}
}

void MemRaster::fillInt(int value, int band) {
	checkInit();
	if(m_props.isFloat()) {
		fillFloat((double) value, band);
	} else {
		size_t chunk = m_mmapped ? m_mappedFile->pageSize() * 1024 : 2048 * 1024;
		size_t size = m_props.size() * sizeof(int);
		size_t last = size % chunk;
		Buffer buf(chunk);
		for(size_t i = 0; i < chunk / sizeof(int); ++i)
			*((int *) buf.buf + i) = value;
		for(uint64_t i = 0; i < size; i += chunk)
			std::memcpy((char *) m_grid + i, buf.buf, g_min(chunk, last));
	}
}

double MemRaster::getFloat(uint64_t idx, int band) {
	checkInit();
	if (idx < 0 || idx >= m_props.size())
		g_argerr("Index out of bounds: " << idx << "; size: " << m_props.size());
	if(m_props.isInt()) {
		return (double) getInt(idx, band);
	} else {
		return *(((double *) m_grid) + idx);
	}
}

double MemRaster::getFloat(int col, int row, int band) {
	uint64_t idx = (uint64_t) row * m_props.cols() + col;
	return getFloat(idx, band);
}

int MemRaster::getInt(uint64_t idx, int band) {
	checkInit();
	if (idx < 0 || idx >= m_props.size())
		g_argerr("Index out of bounds: " << idx << "; size: " << m_props.size());
	if(m_props.isInt()) {
		return *(((int *) m_grid) + idx);
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

	Buffer buf(cols * rows * typeSize);
	GDALRasterBand *band = grd.m_ds->GetRasterBand(dstBand);

	char* input  = (char*) grid();
	char* output = (char*) buf.buf;
	for(int r = 0; r < rows; ++r)
		std::memcpy(output + r * cols * typeSize, input + ((srcRow + r) * gcols + srcCol) * typeSize, cols * typeSize);

	if(CPLE_None != band->RasterIO(GF_Write, dstCol, dstRow, cols, rows, output,
			cols, rows, gtype, 0, 0, 0))
		g_runerr("Failed to write to: " << grd.filename());

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
						extensions[desc].insert(Util::lower(item));
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
		m_block(nullptr),
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
	for(int i = 1; i <= m_props.bands(); ++i)
		m_ds->GetRasterBand(i)->SetNoDataValue(m_props.nodata());
	m_ds->GetRasterBand(1)->GetBlockSize(&m_bcols, &m_brows);
 	m_block = malloc(m_bcols * m_brows * getTypeSize(m_props.dataType()));
}

Raster::Raster(const std::string &filename, bool writable) :
		m_ds(nullptr),
		m_bcols(0), m_brows(0),
		m_bcol(-1), m_brow(-1),
		m_band(1),
		m_block(nullptr),
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
	m_props.setNoData(m_ds->GetRasterBand(1)->GetNoDataValue());
	m_ds->GetRasterBand(1)->GetBlockSize(&m_bcols, &m_brows);
 	m_block = malloc(m_bcols * m_brows * getTypeSize(m_props.dataType()));
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
	int size = getTypeSize(m_props.dataType());
	Buffer buf(m_bcols * m_brows * size);
	for(int i = 0; i < m_bcols * m_brows; ++i)
		writeInt(((char*) buf.buf) + i, size, value);
	GDALRasterBand *bnd = m_ds->GetRasterBand(band);
	int cols = props().cols();
	int rows = props().rows();
	for (int r = 0; r < rows / m_brows; ++r) {
		for (int c = 0; c < cols / m_bcols; ++c) {
			if (bnd->WriteBlock(c, r, buf.buf) != CE_None)
				g_runerr("Fill error.");
		}
	}
}

void Raster::fillFloat(double value, int band) {
	int size = getTypeSize(m_props.dataType());
	Buffer buf(m_bcols * m_brows * size);
	if(size == 4) {
		for(int i = 0; i < m_bcols * m_brows; ++i)
			*(((float*) buf.buf) + i) = (float) value;
	} else {
		for(int i = 0; i < m_bcols * m_brows; ++i)
			*(((double*) buf.buf) + i) = value;
	}
	GDALRasterBand *bnd = m_ds->GetRasterBand(band);
	int cols = props().cols();
	int rows = props().rows();
	for (int r = 0; r < rows / m_brows; ++r) {
		for (int c = 0; c < cols / m_bcols; ++c) {
			if (bnd->WriteBlock(c, r, buf.buf) != CE_None)
				g_runerr("Fill error.");
		}
	}
}

void Raster::writeToRaster(Raster &grd,
			int cols, int rows,
			int srcCol, int srcRow,
			int dstCol, int dstRow,
			int srcBand, int dstBand) {

	if(m_dirty) {
		if(CPLE_None != m_ds->GetRasterBand(m_bband)->WriteBlock(m_bcol, m_brow, m_block))
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

	if(CPLE_None != m_ds->GetRasterBand(srcBand)->RasterIO(GF_Read, srcCol, srcRow, cols, rows,
			buf.buf, cols, rows, gtype, 0, 0, 0))
		g_runerr("Failed to read from: " << grd.filename());
	if(CPLE_None != grd.m_ds->GetRasterBand(dstBand)->RasterIO(GF_Write, dstCol, dstRow, cols, rows,
			buf.buf, cols, rows, gtype, 0, 0, 0))
		g_runerr("Failed to write to: " << filename());
}

void Raster::writeToMemRaster(MemRaster &grd,
		int cols, int rows,
		int srcCol, int srcRow,
		int dstCol, int dstRow,
		int srcBand, int dstBand) {

	if(m_dirty) {
		if(CPLE_None != m_ds->GetRasterBand(m_bband)->WriteBlock(m_bcol, m_brow, m_block))
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
	int typeSize = getTypeSize(gp.dataType());
	int gcols = gp.cols();

	Buffer buf(cols * rows * typeSize);
	GDALRasterBand *band = m_ds->GetRasterBand(srcBand);

	if(CPLE_None != band->RasterIO(GF_Read, srcCol, srcRow, cols, rows, buf.buf,
			cols, rows, gtype, 0, 0, 0))
		g_runerr("Failed to read from: " << filename());

	char* output = (char*) grd.grid();
	char* input  = (char*) buf.buf;
	for(int r = 0; r < rows; ++r)
		std::memcpy(output + ((dstRow + r) * gcols + dstCol) * typeSize , input + r * cols * typeSize, cols * typeSize);

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
	if(bcol != m_bcol || brow != m_brow) {
		if(m_dirty) {
			if(CPLE_None != m_ds->GetRasterBand(m_bband)->WriteBlock(m_bcol, m_brow, m_block))
				g_runerr("Failed to flush to: " << filename());
			m_dirty = false;
		}
		GDALRasterBand *rb = m_ds->GetRasterBand(band);
		if(!rb)
			g_argerr("No such band: " << band);
		if(CPLE_None != rb->ReadBlock(bcol, brow, m_block))
			g_runerr("Failed to read from: " << filename());
		m_bcol = bcol;
		m_brow = brow;
		m_bband = band;
	}
	int idx = (row % m_brows) * m_bcols + (col % m_bcols);
	double v = 0;
	readFromBlock(m_block, getGDType(), &v, idx);
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
	if(bcol != m_bcol || brow != m_brow) {
		if(m_dirty) {
			if(CPLE_None != m_ds->GetRasterBand(m_bband)->WriteBlock(m_bcol, m_brow, m_block))
				g_runerr("Failed to flush to: " << filename());
			m_dirty = false;
		}
		GDALRasterBand *rb = m_ds->GetRasterBand(band);
		if(!rb)
			g_argerr("No such band: " << band);
		if(CPLE_None != rb->ReadBlock(bcol, brow, m_block))
			g_runerr("Failed to read from: " << filename());
		m_bcol = bcol;
		m_brow = brow;
		m_bband = band;
	}
	int idx = (row % m_brows) * m_bcols + (col % m_bcols);
	int v = 0;
	readFromBlock(m_block, getGDType(), &v, idx);
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
	if(bcol != m_bcol || brow != m_brow || band != m_band) {
		if(m_dirty) {
			if(CPLE_None != m_ds->GetRasterBand(m_bband)->WriteBlock(m_bcol, m_brow, m_block))
				g_runerr("Failed to flush to: " << filename());
			m_dirty = false;
		}
		GDALRasterBand *rb = m_ds->GetRasterBand(band);
		if(!rb)
			g_argerr("No such band: " << band);
		if(CPLE_None != rb->ReadBlock(bcol, brow, m_block))
			g_runerr("Failed to read from: " << filename());
		m_bcol = bcol;
		m_brow = brow;
		m_bband = band;
	}
	int idx = (row % m_brows) * m_bcols + (col % m_bcols);
	writeToBlock(m_block, getGDType(), v, idx);
	m_dirty = true;
}

void Raster::setInt(int col, int row, int v, int band) {
	if (!m_props.writable())
		g_runerr("This raster is not writable.");
	int bcol = col / m_bcols;
	int brow = row / m_brows;
	if(bcol != m_bcol || brow != m_brow || band != m_band) {
		if(m_dirty) {
			if(CPLE_None != m_ds->GetRasterBand(m_bband)->WriteBlock(m_bcol, m_brow, m_block))
				g_runerr("Failed to flush to: " << filename());
			m_dirty = false;
		}
		GDALRasterBand *rb = m_ds->GetRasterBand(band);
		if(!rb)
			g_argerr("No such band: " << band);
		if(CPLE_None != rb->ReadBlock(bcol, brow, m_block))
			g_runerr("Failed to read from: " << filename());
		m_bcol = bcol;
		m_brow = brow;
		m_bband = band;
	}
	int idx = (row % m_brows) * m_bcols + (col % m_bcols);
	writeToBlock(m_block, getGDType(), v, idx);
	m_dirty = true;
}

void Raster::setFloat(uint64_t idx, double v, int band) {
	setFloat((int) (idx % m_props.cols()), (int) (idx / m_props.cols()), v, band);
}

void Raster::setFloat(double x, double y, double v, int band) {
	setFloat(m_props.toCol(x), m_props.toRow(y), v, band);
}

void Raster::polygonize(const std::string &filename, const std::string &layerName,
		const std::string &driver, uint16_t srid, uint16_t band, Status *status, bool *cancel) {

	if(cancel == nullptr)
		cancel = &_cancel;

	// Remove the original file; some can't be overwritten directly.
	Util::rm(filename);

	GDALAllRegister();

	// Get the vector driver.
	GDALDriver *drv = GetGDALDriverManager()->GetDriverByName(driver.c_str());
	if(!drv)
		g_runerr("Failed to find driver for " << driver << ".");

	// Create an output dataset for the polygons.
	char **dopts = NULL;
	dopts = CSLSetNameValue(dopts, "SPATIALITE", "YES");
	GDALDataset *ds = drv->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, dopts);
	CPLFree(dopts);
	if(!ds)
		g_runerr("Failed to create dataset " << filename << ".");

	// Create the layer.
	OGRSpatialReference sr;
	sr.importFromEPSG(srid);
	char **lopts = NULL;
	OGRLayer *layer = ds->CreateLayer(layerName.c_str(), &sr, wkbMultiPolygon, lopts);
	CPLFree(lopts);
	if(!layer) {
		GDALClose(ds);
		g_runerr("Failed to create layer " << layerName << ".");
	}

	// There's only one field -- an ID.
	OGRFieldDefn field( "id", OFTInteger);
	layer->CreateField(&field);

	int cols = m_props.cols();
	int rows = m_props.rows();
	// Generate a unique fID for each feature.
	std::atomic<uint64_t> fid(0);
	// Keeps track of extra polys leftover after each thread completes.
	std::unordered_map<uint64_t, std::unique_ptr<Poly> > extraPolys;
	// Keeps track of which rows are completed.
	std::vector<int> finalRows(rows);
	std::fill(finalRows.begin(), finalRows.end(), -1);
	// Status tracker.
	std::atomic<int> stat(0);

	// For creating GEOS objects.
	GEOSContextHandle_t gctx = OGRGeometry::createGEOSContext();
	const GeometryFactory* gf = GeometryFactory::getDefaultInstance();

	if(OGRERR_NONE != layer->StartTransaction())
		g_runerr("Failed to start transation.");

	// TODO: Perturbation to allow polygon union without multis.
	double yShift0 = props().resolutionY() > 0 ? -0.00000001 : 0.00000001;

	#pragma omp parallel
	{

		// Buffer for reading/polygonizing.
		GridProps gp(m_props);
		gp.setSize(cols, 1);
		MemRaster rowBuf(gp);
		// Keeps polygons created in thread.
		std::unordered_map<uint64_t, std::unique_ptr<Poly> > polys;
		int thread = omp_get_thread_num();
		// To track the first and last rows in this thread's block.
		int first = -1, last = 0;

		#pragma omp for
		for(int r = 0; r < rows; ++r) {

			if(*cancel)
				continue;

			if(first == -1) first = r;
			last = r;

			// Write into the buffer.
			#pragma omp critical(__read_raster)
			writeTo(rowBuf, cols, 1, 0, r, 0, 0, band);

			for(int c = 0; c < cols; ++c) {

				// Get the current ID, skip if zero.
				uint64_t id0 = rowBuf.getInt(c, 0);
				if(!id0) continue;

				// Get the coord of one corner of the polygon.
				double x0 = gp.toX(c);
				double y0 = gp.toY(r) + yShift0;

				// Scan right...
				while(++c <= cols) {

					uint64_t id1 = c < cols ? rowBuf.getInt(c, 0) : 0;

					// If the ID changes, capture and output the polygon.
					if(id0 > 0 && id1 != id0) {

						// Coord of the other corner.
						double x1 = gp.toX(c);
						double y1 = gp.toY(r) + gp.resolutionY();

						// Build the geometry.
						CoordinateSequence* seq = gf->getCoordinateSequenceFactory()->create((size_t) 0, 2);
						seq->add(Coordinate(x0, y0));
						seq->add(Coordinate(x1, y0));
						seq->add(Coordinate(x1, y1));
						seq->add(Coordinate(x0, y1));
						seq->add(Coordinate(x0, y0));
						LinearRing* ring = gf->createLinearRing(seq);
						Geometry* geom = gf->createPolygon(ring, nullptr);

						// If it's already in the list, union it, otherwise add it.
						if(polys.find(id0) == polys.end()) {
							std::unique_ptr<Poly> p(new Poly(geom, thread, r, gf));
							polys[id0] = std::move(p);
						} else {
							polys[id0]->update(geom, r, r);
							delete geom;
						}

						--c; // Back up the counter by one, to start with the new ID.
						break;
					}
					id0 = id1;
				}

			}

			// Finalize the row.
			#pragma omp critical(__finalize_rows)
			finalRows[r] = thread;

			// At the end of the row, write any finalizable polygons.
			std::unordered_set<uint64_t> remove;
			#pragma omp critical(__write_layer)
			{
				for(const auto &it : polys) {
					if(it.second->isRangeFinalized(finalRows)) {
						// Add to the removal list
						remove.insert(it.first);
						// Retrieve the unioned geometry.
						OGRGeometry* geom = OGRGeometryFactory::createFromGEOS(gctx, (GEOSGeom) it.second->polyAsMulti());
						// Create and append the feature.
						OGRFeature feat(layer->GetLayerDefn());
						feat.SetGeometry(geom);
						feat.SetField("id", (GIntBig) it.first);
						feat.SetFID(++fid);
						delete geom;
						if(OGRERR_NONE != layer->CreateFeature(&feat))
							g_runerr("Failed to add geometry.");
					}
				}
			}
			// Remove any polys that were written.
			for(const uint64_t &id : remove)
				polys.erase(id);

			if(status)
				status->update((float) ++stat / rows);
		}

		// When a block of rows is completed, set the finalization
		// thread ID to -2, so all threads know this block is
		// finalized.
		#pragma omp critical(__finalize_rows)
		{
			for(int j = first; j < last; ++j) {
				if(finalRows[j] == thread)
					finalRows[j] = -2;
			}
		}

		// If the thread completes and polygons are left over, merge into the
		// extra dict.
		#pragma omp critical(__merge_polys)
		{
			for(auto &it : polys) {
				if(extraPolys.find(it.first) != extraPolys.end()) {
					std::unique_ptr<Poly> &p = it.second;
					extraPolys[it.first]->update(p->poly, p->minRow, p->maxRow);
				} else {
					extraPolys[it.first] = std::move(it.second);
				}
			}
		}
	}

	// Write all the remaining polygons to the output.
	for(const auto &it : extraPolys) {
		if(*cancel)
			break;
		// Retrieve the unioned geometry.
		OGRGeometry* geom = OGRGeometryFactory::createFromGEOS(gctx, (GEOSGeom) it.second->polyAsMulti());
		// Create and append the feature.
		OGRFeature feat(layer->GetLayerDefn());
		feat.SetGeometry(geom);
		feat.SetField("id", (GIntBig) it.first);
		feat.SetFID(++fid);
		delete geom;
		if(OGRERR_NONE != layer->CreateFeature(&feat))
			g_runerr("Failed to add geometry.");
	}

	if(OGRERR_NONE != layer->CommitTransaction())
		g_runerr("Failed to commit transation.");

	GDALClose(ds);
}

void Raster::flush() {
	if(m_brow > -1 && m_bcol > -1 && m_props.writable()) {
		if(CPLE_None != m_ds->GetRasterBand(m_band)->WriteBlock(m_bcol, m_brow, m_block))
			g_warn("Failed to write block to " << filename());
	}
	m_ds->FlushCache();
}

Raster::~Raster() {
	flush();
	if(m_block)
		free(m_block);
	m_block = nullptr;
	if(m_ds)
		GDALClose(m_ds);
	m_ds = nullptr;
}
