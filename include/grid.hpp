/**
 * Raster.hpp
 *
 *  Created on: Jan 21, 2016
 *  Author: rob
 */

#ifndef INCLUDE_GRID_HPP_
#define INCLUDE_GRID_HPP_

#include <sys/mman.h>

//#include <queue>
//#include <stdexcept>
//#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <vector>
#include <list>
#include <cstring>
#include <string>
#include <memory>
#include <mutex>
#include <algorithm>
#include <type_traits>
#include <condition_variable>
#include <thread>

#include <geos_c.h>

#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <ogr_geometry.h>
#include <ogr_feature.h>
#include <ogrsf_frmts.h>

#include "geo.hpp"
#include "util.hpp"

using namespace geo::util;

namespace geo {
namespace grid {


/**
 * Return the interleave type from the string representation.
 *
 * \param[in] str The interleave name.
 * \return The Interleave type.
 */
Interleave interleaveFromString(const std::string& str);

} // grid
} // geo

namespace {

	using namespace geo::grid;

	size_t poly_fid(0);					///<! A feature ID for geometries.
	std::mutex poly_gmtx;  				///<! For the geoms map.
	std::mutex poly_fmtx; 				///<! For the final list.
	std::mutex poly_omtx;  				///<! For OGR writes.
	std::condition_variable poly_cv;	///<! For waiting on the queue.

	/**
	 * Get the size (in bytes) of the given DataType.
	 *
	 * \param type The DataType.
	 * \return The number of bytes required to represent the type.
	 */
	int getTypeSize(DataType type);

	/**
	 * Return the GDAL datatype corresponding to the DataType.
	 *
	 * \param type A DataType.
	 * \return The GDAL datatype corresponding to the given DataType.
	 */
	GDALDataType dataType2GDT(DataType type);

	/**
	 * Return the DataType corresponding to the given GDAL datatype.
	 *
	 * \param type A GDAL datatype.
	 * \return A DataType corresponding to the GDAL datatype.
	 */
	DataType gdt2DataType(GDALDataType type);

	/**
	 * Update the polygonization row status.
	 *
	 * \param r The current row.
	 * \param rows The total number of rows.
	 * \return A status message.
	 */
	std::string polyRowStatus(int r, int rows);

	/**
	 * Make a dataset to contain the polygons.
	 *
	 * \param filename The output filename of the dataset.
	 * \param driver The output driver.
	 * \param layerName The output layer name.
	 * \param idField The field name of the geometry ID.
	 * \param sr The spatial reference.
	 * \param gType The WKB geometry type.
	 * \param ds The GDAL dataset.
	 * \param layer The OGR layer.
	 */
	void polyMakeDataset(const std::string& filename, const std::string& driver, const std::string& layerName,
				const std::string& idField,	OGRSpatialReference* sr, OGRwkbGeometryType gType,
				GDALDataset** ds, OGRLayer** layer);

	/**
	 * Make and return a rectangular geometry with the given corners
	 * and number of dimensions.
	 *
	 * \param x0 The minimum corner x.
	 * \param y0 The minimum corner y.
	 * \param x1 The maximum corner x.
	 * \param y1 The maximum corner y.
	 * \param dims The number of dimensions.
	 * \return A new geometry.
	 */
	GEOSGeometry* polyMakeGeom(double x0, double y0, double x1, double y1, int dims);

	void polyWriteToFile(std::unordered_map<int, std::vector<GEOSGeometry*> >* geoms, std::set<int>* finalIds,
				const std::string* idField, OGRLayer* layer, GEOSContextHandle_t* gctx,
				bool removeHoles, bool removeDangles, bool* running, bool* cancel);

	/**
	 * Fix the given coordinates so that the source does not extend
	 * past the destination and vice versa, etc.
	 *
	 * \param[inout] srcCol The source column.
	 * \param[inout] srcRow The source column.
	 * \param[inout] dstCol The destination column.
	 * \param[inout] dstRow The destination column.
	 * \param[inout] cols The number of columns to operate on.
	 * \param[inout] rows The number of rows to operate on.
	 * \param srcCols The number of columns in the source.
	 * \param srcRows The number of rows in the source.
	 * \param dstCols The number of columns in the destination.
	 * \param dstRows The number of rows in the destination.
	 */
	bool fixCoords(int& srcCol, int& srcRow, int& dstCol, int& dstRow,
			int& cols, int& rows, int srcCols, int srcRows, int dstCols, int dstRows);

	/**
	 * Return the key for the minimum value in the given map.
	 *
	 * \return The key for the minimum value in the given map.
	 */
	size_t minValue(std::unordered_map<size_t, double>& m) {
		double min = std::numeric_limits<double>::max();
		size_t key = 0;
		for(const auto& it : m) {
			if(it.second < min) {
				min = it.second;
				key = it.first;
			}
		}
		return key;
	}

} // anon



namespace geo {
namespace grid {

/**
 * A class containing the properties of a raster.
 */
class G_DLL_EXPORT GridProps {
private:
	double m_trans[6];			///<! The geotransform properties.
	int m_cols, m_rows;			///<! The number of rows and columns.
	int m_vsrid, m_hsrid;		///<! The vertical and horizontal SRIDs
	int m_bands;           		///<! The number of bands.
	bool m_writable;            ///<! True if the raster is writable
	double m_nodata;			///<! The nodata value.
	bool m_nodataSet;			///<! True if nodata is set.
	bool m_compress;			///<! True if the file is a TIF and will be compressed.
	bool m_bigTiff;				///<! Use bigtiff.
	DataType m_type;			///<! The data type.
	Interleave m_interleave;	///<! The interleave method.
	std::string m_filename;		///<! The grid filename.
	std::string m_projection;	///<! The WKT representation of the projection
	std::string m_driver;		///<! The name of the GDAL driver.

public:

	/**
	 * Construct an empty GridProps.
	 */
	GridProps() :
		m_cols(0), m_rows(0),
		m_vsrid(0), m_hsrid(0),
		m_bands(1),
		m_writable(false),
		m_nodata(0),
		m_nodataSet(false),
		m_compress(false),
		m_bigTiff(false),
		m_type(DataType::None),
		m_interleave(Interleave::BIL) {
	}

	/**
	 * Return the internal index of the pixel in the mapped data region
	 * depending on the interleave method.
	 *
	 * \param col The column.
	 * \param row The row.
	 * \param band The band.
	 * \return The index.
	 */
	size_t index(int col, int row, int band) const {
		if(band < 0 || band >= m_bands)
			g_runerr("Invalid band: " << band);
		switch(m_interleave) {
		case Interleave::BIL:
			return row * m_cols * m_bands + band * m_cols + col;
		case Interleave::BSQ:
			return band * m_cols * m_rows + row * m_cols + col;
		case Interleave::BIP:
			return row * m_cols * m_bands + col * m_bands + band;
		default:
			g_runerr("Invalid interleave: "  << (int) m_interleave);
		}
	}

	/**
	 * Return the geographic bounds of the raster.
	 *
	 * \return The geographic bounds of the raster.
	 */
	Bounds bounds() const {
		double x0 = m_trans[0];
		double y0 = m_trans[3];
		double x1 = x0 + m_trans[1] * m_cols;
		double y1 = y0 + m_trans[5] * m_rows;
		return Bounds(g_min(x0, x1), g_min(y0, y1), g_max(x0, x1), g_max(y0, y1));
	}

	/**
	 * Set the geographic bounds of the raster.
	 *
	 * \param bounds The geographic bounds of the raster.
	 */
	void setBounds(const Bounds& bounds) {
		m_trans[0] = m_trans[1] > 0 ? bounds.minx() : bounds.maxx();
		m_trans[3] = m_trans[5] > 0 ? bounds.miny() : bounds.maxy();
		m_cols = (int) std::ceil(bounds.width() / std::abs(m_trans[1]));
		m_rows = (int) std::ceil(bounds.height() / std::abs(m_trans[5]));
	}

	/**
	 * Use compression for tiff files.
	 *
	 * \param compress True to use compression for tiff files.
	 */
	void setCompress(bool compress) {
		m_compress = compress;
	};

	/**
	 * Use compression for tiff files.
	 *
	 * \return True to use compression for tiff files.
	 */
	bool compress() const {
		return m_compress;
	}

	/**
	 * Use Big Tiff setting.
	 *
	 * \param bigTuff True to use Big Tiff setting.
	 */
	void setBigTiff(bool bigTiff) {
		m_bigTiff = bigTiff;
	}

	/**
	 * Use Big Tiff setting.
	 *
	 * \return True to use Big Tiff setting.
	 */
	bool bigTiff() const {
		return m_bigTiff;
	}

	/**
	 * Populate an (at least) 4-element double array with the bounding
	 * box of this object.
	 *
	 * \param bounds A four-element double array.
	 */
	void bounds(double* bounds) const;

	/**
	 * Return the interleave method.
	 *
	 * \return The interleave method.
	 */
	Interleave interleave() const {
		return m_interleave;
	}

	/**
	 * Set the interleave method.
	 *
	 * \param interleave The interleave method.
	 */
	void setInterleave(Interleave interleave) {
		m_interleave = interleave;
	}

	/**
	 * Return the name of the GDAL driver used by the raster.
	 * Only relevant for file-based rasters.
	 *
	 * \return The name of the GDAL driver.
	 */
	std::string driver() const {
		return m_driver;
	}

	/**
	 * Set the name of the GDAL driver used by the raster.
	 * Only relevant for file-based rasters.
	 *
	 * \param name The name of the driver.
	 */
	void setDriver(const std::string& name) {
		m_driver = name;
	}

	/**
	 * Returns the no data value.
	 *
	 * \return The no data value.
	 */
	double nodata() const {
		return m_nodata;
	}

	/**
	 * Set the no data value.
	 *
	 * \param nodata The no data value.
	 */
	void setNoData(double nodata) {
		m_nodata = nodata;
		m_nodataSet = true;
	}

	/**
	 * Returns true if the no data value has been set.
	 *
	 * \return True if the no data value has been set.
	 */
	bool nodataSet() const {
		return m_nodataSet;
	}

	/**
	 * Remove the no data value.
	 */
	void unsetNodata() {
		m_nodataSet = false;
	}

	/**
	 * Return the number of columns.
	 *
	 * \return The number of columns.
	 */
	int cols() const{
		return m_cols;
	}

	/*
	 * Return the number of rows.
	 *
	 * \param The number of rows.
	 */
	int rows() const{
		return m_rows;
	}

	/**
	 * Returns true if the cell is in the raster.
	 *
	 * \param col The column.
	 * \param row The row.
	 * \return True if the cell is in the raster.
	 */
	bool hasCell(int col, int row) const {
		return !(col < 0 || row < 0 || row >= m_rows || col >= m_cols);
	}

	/**
	 * Returns true if the cell is in the raster.
	 *
	 * \param x The geographic x or longitude coordinate.
	 * \param y The geographic y or latitude coordinate.
	 * \return True if the cell is in the raster.
	 */
	bool hasCell(double x, double y) const {
		return hasCell(toCol(x), toRow(y));
	}

	/**
	 * Returns the row for a given y-coordinate.
	 *
	 * \param y The geographic y or latitude coordinate.
	 * \return The row index.
	 */
	int toRow(double y) const {
		return (int) ((y - m_trans[3]) / m_trans[5]);
	}

	/**
	 * Returns the column for a given x-coordinate.
	 *
	 * \param x The geographic x or longitude coordinate.
	 * \return The column index.
	 */
	int toCol(double x) const {
		return (int) ((x - m_trans[0]) / m_trans[1]);
	}

	/**
	 * Returns the x-coordinate for the cell centroid of a given column.
	 *
	 * \param col The column.
	 * \return The x-coordinate at the centre of the column.
	 */
	double toX(int col) const {
		return m_trans[0] + col * m_trans[1] + m_trans[1] * 0.5;
	}

	/**
	 * Returns the y-coordinate for the cell centorid of a given row.
	 *
	 * \param row The row.
	 * \return The y-coordinate at the centre of the row.
	 */
	double toY(int row) const {
		return m_trans[3]  + row * m_trans[5] + m_trans[5] * 0.5;
	}

	/**
	 * Returns the number of pixels a single band in the raster.
	 *
	 * \return The number of pixels in the raster.
	 */
	size_t size() const {
		return m_cols * m_rows;
	}

	/**
	 * Set the data type of the raster.
	 *
	 * \param type The data type.
	 */
	void setDataType(DataType type) {
		m_type = type;
	}

	/**
	 * Get the data type of the raster.
	 *
	 * \return The data type.
	 */
	DataType dataType() const {
		return m_type;
	}

	/**
	 * Set the size of the raster in columns, rows.
	 *
	 * \param col The column.
	 * \param row The row.
	 */
	void setSize(int cols, int rows){
		m_cols = cols;
		m_rows = rows;
	}

	/**
	 * Set the horizontal and vertical (optional) SRID.
	 *
	 * \param hsrid The horizontal SRID.
	 * \param vsrid The vertical SRID.
	 */
	void setSrid(int hsrid, int vsrid = 0) {
		m_vsrid = vsrid;
		m_hsrid = hsrid;
	}

	/**
	 * Get the vertical SRID.
	 *
	 * \return The vertical SRID.
	 */
	int vsrid() const {
		return m_vsrid;
	}

	/**
	 * Get the horizontal SRID.
	 *
	 * \return The horizontal SRID.
	 */
	int hsrid() const {
		return m_hsrid;
	}


	/**
	 * Set the WKT projection.
	 *
	 * \param The projection string (proj or WKT format).
	 */
	void setProjection(const std::string& proj) {
		m_projection = proj;
	}


	/**
	 * Get the WKT projection (proj or WKT format).
	 *
	 * \return The WKT projection.
	 */
	std::string projection() const {
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

	/**
	 * Set the geo transform properties.
	 *
	 * \param trans The six-element transformation matrix.
	 */
	void setTrans(double trans[6]) {
		for(int i = 0; i < 6; ++i)
			m_trans[i] = trans[i];
		setResolution(m_trans[1], m_trans[5]);
	}

	/**
	 * Set the geo transform properties. The third and fifth elements are set to zero.
	 *
	 * \param tlx The top left x-coordinate.
	 * \param resX The horizontal resolution.
	 * \param tly The top left y-coordinate.
	 * \param resY The vertical resolution (negative for UTM (etc.) projections).
	 */
	void setTrans(double tlx, double resX, double tly, double resY) {
		double t[] = {tlx, resX, 0, tly, 0, resY};
		setTrans(t);
	}

	/**
	 * Gets the geo transform properties by setting them in the given array.
	 *
	 * \param trans The six-element transformation matrix.
	 */
	void trans(double trans[6]) const {
		for(int i = 0; i < 6; ++i)
			trans[i] = m_trans[i];
	}

	/**
	 * Set the vertical and horizontal resolution.
	 *
	 * \param resolutionX The horizontal resolution.
	 * \param resolutionY The vertical resolution (negative for UTM (etc.) projections).
	 */
	void setResolution(double resolutionX, double resolutionY) {
		m_trans[1] = resolutionX;
		m_trans[5] = resolutionY;
	}

	/**
	 * Get the horizontal resolution.
	 *
	 * \return The horizontal resolution.
	 */
	double resX() const {
		return m_trans[1];
	}

	/**
	 * Get the vertical resolution.
	 *
	 * \return The vertical resolution.
	 */
	double resY() const {
		return m_trans[5];
	}

	/**
	 * Return the top-left horizontal coordinate of the raster.
	 *
	 * \return The top-left horizontal coordinate of the raster.
	 */
	double tlx() const {
		return m_trans[0];
	}

	/**
	 * Return the top-left vertical coordinate of the raster.
	 *
	 * \return The top-left vertical coordinate of the raster.
	 */
	double tly() const {
		return m_trans[3];
	}

	/**
	 * Set the number of bands.
	 *
	 * \param bands The number of bands.
	 */
	void setBands(int bands) {
		m_bands = bands;
	}

	/**
	 * Get the number of bands.
	 *
	 * \return The number of bands.
	 */
	int bands() const {
		return m_bands;
	}

	/**
	 * Set the writable state of the raster. If the raster is not writable,
	 * attempting to write to it will throw an exception.
	 *
	 * \param writable True, if the raster should be writable.
	 */
	void setWritable(bool writable) {
		m_writable = writable;
	}

	/**
	 * Get the writable state of the raster.
	 *
	 * \param The writable state of the raster.
	 */
	bool writable() const  {
		return m_writable;
	}

	/**
	 * Return the GDAL datatype of the raster.
	 *
	 * \return The GDAL datatype of the raster.
	 */
	GDALDataType gdalDataType() const {
		return dataType2GDT(dataType());
	}

	/**
	 * Set the filename of the raster.
	 *
	 * \param filename The filename of the raster.
	 */
	void setFilename(const std::string& filename) {
		m_filename = filename;
	}

	/**
	 * Return the filename for this raster.
	 *
	 * \return The filename for this raster.
	 */
	std::string filename() const {
		return m_filename;
	}

};

/**
 * A class to contain statistics of the raster.
 */
class G_DLL_EXPORT GridStats {
public:
	double min;
	double max;
	double mean;
	double stdDev;
	double variance;
	double sum;
	size_t count;
};

/**
 * A simple class to represent a single grid cell.
 */
class G_DLL_EXPORT Cell {
public:
	int col;
	int row;

	Cell(int col, int row) :
		col(col), row(row) {}
};

/**
 * Used by Grid::floodFill to determine whether a pixel should be filled.
 * Subclasses will implement clever ways to detect fillable
 * pixels.
 */
template <class T, class U>
class G_DLL_EXPORT FillOperator {
public:

	/**
	 * Return the properties of the source raster.
	 *
	 * \return The properties of the source raster.
	 */
	virtual const GridProps& srcProps() const = 0;

	/**
	 * Return the properties of the destination raster.
	 *
	 * \return The properties of the destination raster.
	 */
	virtual const GridProps& dstProps() const = 0;

	/**
	 * Return true if if the current pixel should be filled.
	 *
	 * \param col The column.
	 * \param row The row.
	 * \return True if if the current pixel should be filled.
	 */
	virtual bool shouldFill(int col, int row) const = 0;

	/**
	 * Fill the current column.
	 *
	 * \param col The column.
	 * \param row The row.
	 */
	virtual void fill(int col, int row) const = 0;

	virtual ~FillOperator() {};
};

// Forward declaration.
template <class T, class U>
class G_DLL_EXPORT TargetFillOperator;

/**
 * Grid (raster).
 */
template <class T>
class G_DLL_EXPORT Grid {
private:
	GDALDataset *m_ds;          				///<! GDAL data set pointer.
	GDALDataType m_type;        				///<! GDALDataType -- limits the possible template types.
	GridProps m_props;							///<! Properties of the raster.
	std::unique_ptr<TmpFile> m_mapFile;
	bool m_mapped;
	size_t m_size;
	T* m_data;
	bool m_dirty;								///<! True if the current block has been written to and must be flushed.

	void initMapped() {
		if(m_data)
			destroy();
		m_size = props().cols() * props().rows() * props().bands() * sizeof(T);
		m_mapFile.reset(new TmpFile(m_size));
		m_data = (T*) mmap(0, m_size, PROT_READ|PROT_WRITE, MAP_SHARED, m_mapFile->fd, 0);
		if(!m_data)
			g_runerr("Failed to map " << m_size << " bytes for grid.");
		m_mapped = true;
	}

	void initMem() {
		if(m_data)
			destroy();
		m_size = props().cols() * props().rows() * props().bands() * sizeof(T);
		m_data = (T*) malloc(m_size);
		if(!m_data)
			g_runerr("Failed to allocate " << m_size << " bytes for grid.");
	}

	DataType type() const {
		if(std::is_same<T, double>::value) {
			return DataType::Float64;
		} else if(std::is_same<T, float>::value) {
			return DataType::Float32;
		} else if(std::is_same<T, int>::value) {
			return DataType::Int32;
		} else if(std::is_same<T, unsigned int>::value) {
			return DataType::UInt32;
		} else if(std::is_same<T, short>::value) {
			return DataType::Int16;
		} else if(std::is_same<T, unsigned short>::value) {
			return DataType::UInt16;
		} else if(std::is_same<T, char>::value) {
			return DataType::Byte;
		} else if(std::is_same<T, unsigned char>::value) {
			return DataType::Byte;
		} else {
			return DataType::None;
		}
	}

	GDALDataType gdalType() const {
		if(std::is_same<T, double>::value) {
			return GDT_Float64;
		} else if(std::is_same<T, float>::value) {
			return GDT_Float32;
		} else if(std::is_same<T, int>::value) {
			return GDT_Int32;
		} else if(std::is_same<T, unsigned int>::value) {
			return GDT_UInt32;
		} else if(std::is_same<T, short>::value) {
			return GDT_Int16;
		} else if(std::is_same<T, unsigned short>::value) {
			return GDT_UInt16;
		} else if(std::is_same<T, char>::value) {
			return GDT_Byte;
		} else if(std::is_same<T, unsigned char>::value) {
			return GDT_Byte;
		} else {
			return GDT_Unknown;
		}
	}

	bool isFloat() const {
		DataType t = type();
		return t == DataType::Float32 || t == DataType::Float64;
	}

public:

	Grid() :
		m_ds(nullptr),
		m_mapped(false),
		m_size(0),
		m_data(nullptr),
		m_dirty(false),
		m_type(GDT_Unknown) {}

	/**
	 * Create an anonymous grid of the given size with the given
	 * properties.
	 *
	 * \param props The properties of the grid.
	 * \param mapped If true, the grid is created in a mapped memory segment.
	 */
	Grid(const GridProps& props, bool mapped = false) : Grid() {
		init(props, mapped);
	}

	void init(const GridProps& props, bool mapped = false) {
		m_props = props;
		if(mapped) {
			initMapped();
		} else {
			initMem();
		}
	}

	/**
	 * Create a new raster with the given properties. The raster will be created.
	 *
	 * \param filename The path to the file.
	 * \param props A GridProps instance containing a descriptor for the raster.
	 */
	Grid(const std::string& filename, const GridProps& props) : Grid() {
		init(filename, props);
	}

	/**
	 * Initialize the raster by loading it and setting up the grid properties.
	 *
	 * \param filename The path to the file.
	 * \param props A GridProps instance containing a descriptor for the raster.
	 */
	void init(const std::string& filename, const GridProps& props) {

		if (props.resX() == 0 || props.resY() == 0)
			g_argerr("Resolution must not be zero.");
		if (props.cols() <= 0 || props.rows() <= 0)
			g_argerr("Columns and rows must be larger than zero.");
		if (filename.empty())
			g_argerr("Filename must be given.");

		m_props = props;
		m_props.setFilename(filename);

		// Create GDAL dataset.
		char **opts = NULL;
		if(m_props.compress()) {
			opts = CSLSetNameValue(opts, "COMPRESS", "LZW");
			opts = CSLSetNameValue(opts, "PREDICTOR", "2");
		}
		if(m_props.bigTiff())
			opts = CSLSetNameValue(opts, "BIGTIFF", "YES");
		if(m_props.interleave() == Interleave::BIL) {
			opts = CSLSetNameValue(opts, "INTERLEAVE", "BAND");
		} else if(m_props.interleave() == Interleave::BIP){
			opts = CSLSetNameValue(opts, "INTERLEAVE", "PIXEL");
		}

		GDALAllRegister();

		// Figure out and load the driver.

		std::string drvName = m_props.driver();
		if(drvName.empty())
			drvName = getDriverForFilename(filename);

		if(drvName.empty()) {
			g_runerr("Couldn't find driver for: " << filename);
		} else {
			m_props.setDriver(drvName);
		}

		GDALDriver *drv = GetGDALDriverManager()->GetDriverByName(drvName.c_str());
		if(!drv)
			g_runerr("Failed to get driver for " << drvName);

		const char *create = drv->GetMetadataItem(GDAL_DCAP_CREATE);
		if(create == NULL || std::strncmp(create, "YES", 3) != 0)
			g_runerr("The " << drvName << " driver does not support dataset creation. Please specify a different driver.");


		// Remove and create the raster.

		rem(filename);

		m_ds = drv->Create(filename.c_str(), m_props.cols(), m_props.rows(), m_props.bands(), dataType2GDT(m_props.dataType()), opts);

		if(opts)
			CSLDestroy(opts);

		if (!m_ds)
			g_runerr("Failed to create file: " << filename);

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

		// Map the raster into virtual memory.

		initMapped();
	}


	/**
	 * Open the given extant raster. Set the writable argument to true
	 * to enable writing.
	 *
	 * \param filename The path to the file.
	 * \param writable True if the file is to be writable.
	 */
	Grid(const std::string& filename, bool writable = false) : Grid() {

		init(filename, writable);

	}

	/**
	 * Initalize with an extant raster. Set the writable argument to true
	 * to enable writing.
	 *
	 * \param filename The path to the file.
	 * \param writable True if the file is to be writable.
	 */
	void init(const std::string& filename, bool writable = false) {

		if (filename.empty())
			g_argerr("Filename must be given.");

		GDALAllRegister();

		// Attempt to open the dataset.

		m_ds = (GDALDataset *) GDALOpen(filename.c_str(), writable ? GA_Update : GA_ReadOnly);
		if (m_ds == NULL)
			g_runerr("Failed to open raster.");

		GDALDriver *drv = m_ds->GetDriver();
		if(drv == NULL)
			g_runerr("Failed to retrieve driver.");

		const char *drvName = drv->GetDescription();
		if(drvName != NULL)
			m_props.setDriver(drvName);

		char** interleave = m_ds->GetMetadata("INTERLEAVE");

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
		m_props.setFilename(filename);
		m_props.setInterleave(interleaveFromString(*interleave));

	}

	/**
	 * Return a map containing the raster driver short name and extension.
	 *
	 * \return A map containing the raster driver short name and extension.
	 */
	static std::map<std::string, std::set<std::string> > extensions() {
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
						split(std::back_inserter(lst), std::string(ext));
						for(const std::string& item : lst)
							extensions[desc].insert("." + lowercase(item));
					}

				}
			}
		}
		return extensions;
	}

	/**
	 * Return a map containing the raster driver short name and long name.
	 *
	 * \return A map containing the raster driver short name and long name.
	 */
	static std::map<std::string, std::string> drivers() {
		std::vector<std::string> f;
		return drivers(f);
	}

	/**
	 * Return a map containing the raster driver short name and long name. Use filter
	 * to filter the returns on short name.
	 *
	 * \param filter A vector containing the short names of drivers to include.
	 * \return A map containing the raster driver short name and long name.
	 */
	static std::map<std::string, std::string> drivers(const std::vector<std::string>& filter) {
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


	/**
	 * Get the name of the driver that would be used to open a file
	 * with the given path.
	 *
	 * \param filename The path to an existing raster.
	 * \return The name of the griver used to open the file.
	 */
	static std::string getDriverForFilename(const std::string& filename) {
		std::string ext = extension(filename);
		std::map<std::string, std::set<std::string> > drivers = extensions();
		std::string result;
		for(const auto& it : drivers) {
			if(it.second.find(ext) != it.second.end())
				result = it.first;
		}
		return result;
	}

	/**
	 * Creates a virtual raster using the given files and writes it to a file
	 * with the given name.
	 *
	 * \param files A list of files to include in the raster.
	 * \param outfile The path to the virtual raster.
	 * \param nodata The nodata value for the virtual raster.
	 */
	static void createVirtualRaster(const std::vector<std::string>& files, const std::string& outfile, double nodata);

	/**
	 * Creates a virtual raster using the given files and writes it to a file
	 * with the given name.
	 *
	 * \param begin An iterator into a list of files to include in the raster.
	 * \param files The end iterator of the list of files to include in the raster.
	 * \param outfile The path to the virtual raster.
	 * \param nodata The nodata value for the virtual raster.
	 */
	template <class U>
	static void createVirtualRaster(U begin, U end, const std::string& outfile, double nodata) {
		std::vector<std::string> files(begin, end);
		return createVirtualRaster(files, outfile, nodata);
	}

	/**
	 * Compute the table of Gaussian weights given the size of
	 * the table and the standard deviation.
	 *
	 * \param weights The list of weights.
	 * \param size The size of the weights list.
	 * \param sigma The standard deviation.
	 * \param mean The centre of the curve.
	 */
	static void gaussianWeights(double *weights, int size, double sigma, double mean = 0) {
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

	/**
	 * Attempts to return the data type of the raster
	 * with the given filename.
	 *
	 * \param filename The path to an existing raster.
	 * \return The data type.
	 */
	static DataType dataType(const std::string& filename) {
		return DataType::None;
	}

	/**
	 * Return the GDAL data set pointer.
	 *
	 * \return The GDAL data set pointer.
	 */
	GDALDataset* gdalDataset() const {
		return m_ds;
	}

	/**
	 * Copies the image data from an entire row into the buffer
	 * which must be pre-allocated.
	 *
	 * \param row The row index.
	 * \param band The band number.
	 * \param buf A pre-allocated buffer to store the data.
	 */
	void getRow(int row, int band, T* buf) {

		for(int c = 0; c < props().cols(); ++c)
			buf[c] = get(c, row, band);

	}

	/**
	 * Copies the image data from an entire row into the buffer
	 * which must be pre-allocated.
	 *
	 * \param row The row index.
	 * \param band The band number.
	 * \param buf A pre-allocated buffer to store the data.
	 */
	void setRow(int row, int band, T* buf) {

		for(int c = 0; c < props().cols(); ++c)
			set(c, row, buf[c], band);

	}

	/**
	 * Copies the image data from an entire column into the buffer
	 * which must be pre-allocated.
	 *
	 * \param column The column index.
	 * \param band The band number.
	 * \param buf A pre-allocated buffer to store the data.
	 */
	void getColumn(int col, int band, T* buf) {

		for(int r = 0; r < props().rows(); ++r)
			buf[r] = get(col, r, band);

	}

	/**
	 * Copies the image data from an entire column into the buffer
	 * which must be pre-allocated.
	 *
	 * \param col The column index.
	 * \param band The band number.
	 * \param buf A pre-allocated buffer to store the data.
	 */
	void setColumn(int col, int band, T* buf) {

		for(int r = 0; r < props().rows(); ++r)
			set(col, r, buf[r], band);

	}

	/**
	 * Copies the image data from an entire pixel into the buffer
	 * which must be pre-allocated. The buffer represents
	 * the band stack for that pixel.
	 *
	 * \param col The column index.
	 * \param row The row index.
	 * \param buf A pre-allocated buffer to store the data.
	 */
	void getPixel(int col, int row, T* buf) {

		for(int b = 0; b < props().bands(); ++b)
			buf[b] = get(col, row, b);

	}

	/**
	 * Copies the image data from an entire column into the buffer
	 * which must be pre-allocated.
	 *
	 * \param band The band index.
	 * \param band The band number.
	 * \param buf A pre-allocated buffer to store the data.
	 */
	void setPixel(int col, int row, T* buf) {

		for(int b = 0; b < props().bands(); ++b)
			set(col, row, buf[b], b);

	}

	/**
	 * Return the properties of this Grid.
	 *
	 * \return The properties of this Grid.
	 */
	const GridProps& props() const {
		return m_props;
	}

	/**
	 * Compute and return the statistics for the band.
	 *
	 * \param band The raster band.
	 * \return A GridStats instance containing computed statistics.
	 */
	GridStats stats(int band) {
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
				if ((v = get<double>(col, row, band)) != nodata) {
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

	/**
	 * Fill the entire dataset with the given value.
	 *
	 * \param value The value to fill the raster with.
	 * \param band The band to fill.
	 */
	template <class U>
	void fill(U value, int band) {

		if(band >= 0 && band < props().bands()){

			int cols = props().cols();
			int rows = props().rows();

			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c)
					set<U>(c, r, value, band);
			}

		} else {

			g_runerr("Invalid fill band: " << band);

		}

	}

	/**
	 * Fill the entire dataset with the given value.
	 *
	 * \param value The value to fill the raster with.
	 * \param band The band to fill.
	 */
	void fill(T value, int band) {
		fill<T>(value, band);
	}

	/**
	 * Fill all bands with the given value.
	 *
	 * \param value The fill value.
	 */
	void fill(T value) {

		for(int i = 0; i < props().bands(); ++i)
			fill(value, i);

	}

	/**
	 * Return a the value held at the given position in the grid.
	 *
	 * \param col The column.
	 * \param row The row.
	 * \param band The band.
	 * \return The value held at the given index in the grid.
	 */
	template <class U>
	U get(int col, int row, int band) {

		size_t idx = props().index(col, row, band);
		return (U) m_data[idx];

	}

	template <class U>
	U get(double x, double y, int band) {

		return get<U>(props().toCol(x), props().toRow(y), band);

	}

	/**
	 * Set the value held at  the given index in the grid.
	 *
	 * \param col The column.
	 * \param row The row.
	 * \param value The value to set.
	 * \param band The band.
	 */
	template <class U>
	void set(int col, int row, U value, int band) {

		if (!m_props.writable())
			g_runerr("This raster is not writable.");

		size_t idx = props().index(col, row, band);
		m_data[idx] =(T) value;

		m_dirty = true;

	}

	template <class U>
	void set(double x, double y, U value, int band) {

		set<U>(m_props.toCol(x), m_props.toRow(y), value, band);

	}

	/**
	 * Return a the value held at the given position in the grid.
	 *
	 * \param col The column.
	 * \param row The row.
	 * \param band The band.
	 * \return The value held at the given index in the grid.
	 */
	T get(int col, int row, int band) {

		return get<T>(col, row, band);

	}

	/**
	 * Return a the value held at the given position in the grid.
	 *
	 * \param x The x coordinate.
	 * \param y The y coordinate.
	 * \param band The band.
	 * \return The value held at the given index in the grid.
	 */
	T get(double x, double y, int band) {

		return get(props().toCol(x), props().toRow(y), band);

	}

	/**
	 * Set the value held at  the given index in the grid.
	 *
	 * \param col The column.
	 * \param row The row.
	 * \param value The value to set.
	 * \param band The band.
	 */
	void set(int col, int row, T value, int band) {

		set<T>(col, row, value, band);

	}

	/**
	 * Set the value held at  the given index in the grid.
	 *
	 * \param x The x coordinate.
	 * \param y The y coordinate.
	 * \param value The value to set.
	 * \param band The band.
	 */
	void set(double x, double y, T value, int band) {

		set(m_props.toCol(x), m_props.toRow(y), value, band);

	}

	/**
	 * Write data from the current Grid instance to the given grid.
	 *
	 * \param grd The target grid.
	 * \param cols The number of columns to write.
	 * \param rows The number of rows to write.
	 * \param srcCol The source column to read from.
	 * \param srcRow The source row to read from.
	 * \param dstCol The destination column to write to.
	 * \param dstRow The destination row to write to.
	 * \param srcBand The source band.
	 * \param dstBand The destination band.
	 */
	void writeTo(Grid<T>& grd,
			int cols = 0, int rows = 0,
			int srcCol = 0, int srcRow = 0,
			int dstCol = 0, int dstRow = 0,
			int srcBand = 1, int dstBand = 1) {

		int srcCols = props().cols();
		int srcRows = props().rows();
		int dstCols = grd.props().cols();
		int dstRows = grd.props().rows();

		fixCoords(srcCol, srcRow, dstCol, dstRow, cols, rows, srcCols, srcRows, dstCols, dstRows);

		for(int r = srcRow; r < srcRow + rows; ++r) {
			for(int c = srcCol; c < srcCol + cols; ++c)
				set(c - srcCol + dstCol, r - srcRow + dstRow, get(c, r, srcBand), dstBand);
		}

	}

	/**
	 * Return the pixel values as a vector, for the given band.
	 */
	std::vector<T> asVector(int band) {

		int cols = props().cols();
		int rows = props().rows();

		std::vector<T> data(cols * rows);

		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c)
				data[r * cols + c] = get(c, r, band);
		}

		return data;

	}

	/**
	 * Normalize the grid so that one standard deviation is +-1.
	 *
	 * \param band The target band.
	 */
	void normalize(int band)  {
		GridStats st = stats(1);
		const GridProps& gp = props();
		double v, nodata = gp.nodata();
		double mean = st.mean;
		double stdDev = st.stdDev;
		int rows = gp.rows();
		int cols = gp.cols();
		for(int row = 0; row < rows; ++row) {
			for(int col = 0; col < cols; ++col) {
				if ((v = get<double>(col, row, band)) != nodata && !std::isnan(v) && v < G_DBL_MAX_POS) {
					set(col, row, ((v - mean) / stdDev), band);
				} else {
					set(col, row, nodata, band);
				}
			}
		}
	}

	/**
	 * Normalize the grid so that the max value is equal to 1, and the
	 * minimum is zero.
	 *
	 * \param band The target band.
	 */
	void logNormalize(int band) {
		GridStats st = stats(band);
		const GridProps& gp = props();
		double n = st.min;
		double x = st.max;
		double e = std::exp(1.0) - 1.0;
		int rows = gp.rows();
		int cols = gp.cols();
		for(int row = 0; row < rows; ++row) {
			for(int col = 0; col < cols; ++col) {
					set(col, row, std::log(1.0 + e * (get<double>(col, row, band) - n) / (x - n)), band);
			}
		}
	}

	/**
	 * Convert a Grid to some other type.
	 *
	 * \param g The destination Grid.
	 * \param srcBand The source band.
	 * \param dstBand The destination band.
	 */
	template <class U>
	void convert(Grid<U>& g, int srcBand, int dstBand) {
		const GridProps& gp = props();
		int rows = gp.rows();
		int cols = gp.cols();
		for(int row = 0; row < rows; ++row) {
			for(int col = 0; col < cols; ++col) {
				g.set(col, row, get<U>(col, row, srcBand), dstBand);
			}
		}
	}

	/**
	 * Fill the grid, beginning with the target cell, where any contiguous cell
	 * satisfies the given FillOperator. The other grid is actually filled,
	 * and the present grid is unchanged *unless* the present grid is passed
	 * as other.
	 *
	 * \param col   The column to start on.
	 * \param row   The row to start on.
	 * \param op    A FillOperator instance which will determine
	 *              whether a pixel should be filled.
	 * \param other The grid whose cells will actually be filled.
	 * \param fill  The value to fill cells with.
	 * \param d8    Whether to enable diagonal fills.
	 * \param out*  Pointer to variables that hold min and max rows and columns
	 *              plus the area of the fill's bounding box.
	 */
	template <class U>
	static void floodFill(int col, int row,
		FillOperator<T, U>& op,  bool d8 = false,
		int *outminc = nullptr, int *outminr = nullptr,
		int *outmaxc = nullptr, int *outmaxr = nullptr,
		int *outarea = nullptr) {

		const GridProps& gp = op.srcProps();

		int cols = gp.cols();
		int rows = gp.rows();
		size_t size = gp.size();
		int minc = cols + 1;
		int minr = rows + 1;
		int maxc = -1;
		int maxr = -1;
		int area = 0;

		std::queue<Cell> q;
		q.emplace(col, row);

		std::vector<bool> visited(size, false); // Tracks visited pixels.

		while (q.size()) {

			const Cell& cel = q.front();
			row = cel.row;
			col = cel.col;
			q.pop();

			size_t idx = (size_t) row * cols + col;

			if (!visited[idx] && op.shouldFill(col, row)) {

				minc = g_min(col, minc);
				maxc = g_max(col, maxc);
				minr = g_min(row, minr);
				maxr = g_max(row, maxr);
				++area;
				op.fill(col, row);
				visited[idx] = true;

				if (row > 0)
					q.push(Cell(col, row - 1));
				if (row < rows - 1)
					q.push(Cell(col, row + 1));

				int c;
				for (c = col - 1; c >= 0; --c) {
					idx = (size_t) row * cols + c;
					if (!visited[idx] && op.shouldFill(c, row)) {
						minc = g_min(c, minc);
						++area;
						op.fill(c, row);
						visited[idx] = true;
						if (row > 0)
							q.push(Cell(c, row - 1));
						if (row < rows - 1)
							q.push(Cell(c, row + 1));
					} else {
						break;
					}
				}
				if(d8) {
					if (row > 0)
						q.push(Cell(c, row - 1));
					if (row < rows - 1)
						q.push(Cell(c, row + 1));
				}
				for (c = col + 1; c < cols; ++c) {
					idx = (size_t) row * cols + c;
					if (!visited[idx] && op.shouldFill(c, row)) {
						maxc = g_max(c, maxc);
						++area;
						op.fill(c, row);
						visited[idx] = true;
						if (row > 0)
							q.push(Cell(c, row - 1));
						if (row < rows - 1)
							q.push(Cell(c, row + 1));
					} else {
						break;
					}
				}
				if(d8) {
					if (row > 0)
						q.push(Cell(c, row - 1));
					if (row < rows - 1)
						q.push(Cell(c, row + 1));
				}
			}
		}
		if(outminc != nullptr)
			*outminc = minc;
		if(outminr != nullptr)
			*outminr = minr;
		if(outmaxc != nullptr)
			*outmaxc = maxc;
		if(outmaxr != nullptr)
			*outmaxr = maxr;
		if(outarea != nullptr)
			*outarea = area;
	}

	/**
	 * Smooth the raster and write the smoothed version to the output raster.
	 * Callback is an optional function reference with a single float
	 * between 0 and 1, for status tracking.
	 *
	 * Sigma defaults to 0.84089642, window size to 3.
	 *
	 * \param smoothed The smoothed grid.
	 * \param sigma    The standard deviation.
	 * \param size     The window size.
	 * \param band     The target band.
	 * \param monitor  A reference to the Monitor.
	 */
	void smooth(Grid<T>& smoothed, double sigma, int size, int band, geo::Monitor* monitor = nullptr) {

		static Monitor _monitor;

		if(!monitor)
			monitor = &_monitor;

		m_dirty = true;

		const GridProps& gp = props();

		monitor->status(0.01);

		if (sigma <= 0)
			g_argerr("Sigma must be > 0.");
		if (size < 3)
			g_argerr("Kernel size must be 3 or larger.");
		if (size % 2 == 0) {
			g_warn("Kernel size must be odd. Rounding up.");
			size++;
		}

		// Compute the weights for Gaussian smoothing.
		std::vector<double> weights(size * size * getTypeSize(DataType::Float64));
		Grid::gaussianWeights(weights.data(), size, sigma);

		double nodata = gp.nodata();

		monitor->status(0.02);

		std::mutex wmtx; // write mutex
		std::mutex rmtx; // read mutex
		int row = 0;

		// Run the smoothing jobs.
		int numThreads = 1;
		std::vector<std::thread> threads;
		std::vector<std::exception_ptr> exceptions(numThreads);
		for(int i = 0; i < numThreads; ++i) {
			threads.emplace_back(smooth, this, &smoothed, size, nodata, weights.data(), &row,
					&rmtx, &wmtx, monitor, &exceptions[i]);
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

		monitor->status(1.0);
	}

	/**
	 * The radius is given with cells as the unit, but can be rational.
	 * When determining which cells to include in the calculation,
	 * any cell which partially falls in the radius will be included.
	 *
	 * \param filename 	The output filename.
	 * \param band   	The source band.
	 * \param mask		A raster to use as a mask; invalid pixels will be ignored.
	 * \param maskBand  The band to use in the mask.
	 * \param radius 	The search radius.
	 * \param count  	The number of pixels to use for calculations.
	 * \param exp    	The exponent.
	 */
	void voidFillIDW(const std::string& filename, int band, const std::string& mask,
			int maskBand, double radius, int count = 4, double exp = 2.0) {

		if(!isFloat())
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
		Grid<T> input(iprops);
		writeTo(input, iprops.cols(), iprops.rows(), 0, 0, 0, 0, band);

		GridProps oprops(props());
		oprops.setBands(1);
		oprops.setWritable(true);
		Grid<T> output(oprops);

		double maxDist = 100;
		bool holesOnly = true;

		double nodata = props().nodata();
		double v, d;
		int rows = props().rows();
		int cols = props().cols();

		TargetFillOperator<T, T> op1(&input, band, band, nodata, 99999);
		TargetFillOperator<T, T> op2(&input, band, &output, band, 99999, nodata);
		TargetFillOperator<T, T> op3(&input, band, band, 99999, 99998);
		int outminc, outminr, outmaxc, outmaxr;

		for (int r = 0; r < rows; ++r) {
			if(r % 100 == 0) {
				std::cerr << "Row " << r << " of " << rows << "\n";
			}
			for (int c = 0; c < cols; ++c) {

				v = input.get(c, r, band);
				if(v == 99998) {

					//output.setFloat(c, r, nodata);

				} else if (v != nodata) {

					output.set(c, r, v, 1);

				} else if(!holesOnly) {

					double dp, a = 0, b = 0;
					int cnt = 0;
					for(int r0 = g_max(0, r - maxDist); r0 < g_min(rows, r + maxDist + 1); ++r0) {
						for(int c0 = g_max(0, c - maxDist); c0 < g_min(cols, c + maxDist + 1); ++c0) {
							if((c0 == c && r0 == r) || (d = g_sq(c0 - c) + g_sq(r0 - r)) > maxDist ||
									(v = input.get(c0, r0, band)) == nodata)
								continue;
							dp = 1.0 / std::pow(d, exp);
							a += dp * v;
							b += dp;
							++cnt;
						}
					}
					output.set(c, r, cnt ? (a / b) : nodata, 1);

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
							v = input.get(c0, r0, band);
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
						output.set(nc, nr, cnt ? (a / b) : nodata, 1);
					}

					// Fill again with a different value so it will be ignored.
					input.floodFill(c, r, op3, false);
				}
			}
		}

		Grid<T> routput(filename, oprops);
		output.writeTo(routput);

	}
	// TODO: Document me.
	template <class V>
	void writeAStarPath(size_t start, std::unordered_map<size_t, size_t>& parents, V inserter) {
		*inserter = start;
		++inserter;
		while(parents.find(start) != parents.end()) {
			start = parents[start];
			*inserter = start;
			++inserter;
		}
	}

	/**
	 * Finds the least-cost path from the start cell to the goal cell,
	 * using the given heuristic. Populates the given iterator with
	 * the optimal path between the start cell and the goal.
	 *
	 * If the search fails for some reason, like exceeding the maxCost, returns
	 * false. Otherwise returns true.
	 *
	 * \param startCol The starting column.
	 * \param startrow The starting row.
	 * \param goalCol The column of the goal.
	 * \param goalRow The row of the goal.
	 * \param heuristic Used by the algorithm to estimate the future cost of the path.
	 * \param inserter Used to accumulate the path results.
	 * \param maxCost If the total cost exceeds this amount, just quit and return false.
	 * \return True if the search succeeded, false otherwise.
	 */
	template <class U, class V>
	bool searchAStar(int startCol, int startRow, int goalCol, int goalRow, U heuristic, V inserter, double maxCost = std::numeric_limits<double>::infinity()) {

		static double offsets[4][2] = {{0, -1}, {-1, 0}, {1, 0}, {0, 1}};

		size_t goal = ((size_t) goalCol << 32) | goalRow;
		size_t start = ((size_t) startCol << 32) | startRow;

		std::unordered_map<size_t, size_t> parents;
		std::unordered_map<size_t, double> gscore;
		std::unordered_map<size_t, double> fscore;

		std::unordered_set<size_t> openSet;
		std::unordered_set<size_t> closedSet;

		openSet.insert(start);
		gscore[start] = 0; 						// Distance from start to neighbour
		fscore[start] = heuristic(start, goal); // Distance from neighbour to goal.

		int cols = props().cols();
		int rows = props().rows();

		while(!openSet.empty()) {

			if(openSet.size() % 10000 == 0) {
				std::cerr << openSet.size() << "\n";
			}

			size_t top = minValue(fscore);

			if(top == goal) {
				writeAStarPath(top, parents, inserter);
				return true;
			}

			double gscore0 = gscore[top];

			fscore.erase(top);
			gscore.erase(top);

			openSet.erase(top);
			closedSet.insert(top);

			int qcol = (top >> 32) & 0xffffffff;
			int qrow = top & 0xffffffff;

			for(int i = 0; i < 4; ++i) {
				int col = qcol + offsets[i][0];
				int row = qrow + offsets[i][1];

				if(col < 0 || row < 0 || col >= cols || row >= rows)
					continue;

				size_t n = ((size_t) col << 32) | row;

				if(closedSet.find(n) != closedSet.end())
					continue;

				double tgscore = gscore0 + heuristic(top, n);

				if(tgscore > maxCost)
					return false;

				if(openSet.find(n) == openSet.end()) {
					openSet.insert(n);
				} else if(tgscore >= gscore[n]) {
					continue;
				}

				parents[n] = top;
				gscore[n] = tgscore;
				fscore[n] = tgscore + heuristic(n, goal);
			}
		}

		return true;
	}

	void flush(Monitor* monitor = nullptr) {

		if(!m_ds || !props().writable())
			return;

		int cols = props().cols();
		int rows = props().rows();
		int bands = props().bands();

		std::vector<T> buf(cols);

		for(int b = 0; b < bands; ++b) {
			GDALRasterBand* band = m_ds->GetRasterBand(b + 1);
			for(int r = 0; r < rows; ++r) {
				getRow(r, b, buf.data());
				if(CE_None != band->RasterIO(GF_Write, 0, r, cols, 1, buf.data(), cols, 1, gdalType(), 0, 0, 0)) {
					if(monitor)
						monitor->error("Failed to write to raster.");
				}

				if(monitor)
					monitor->status((float) r / rows, "Flushing grid to file.");
			}
		}
	}

	void destroy()  {

		if(m_ds) {

			GDALClose(m_ds);
			m_ds = nullptr;

		}

		if(m_mapped) {

			munmap(m_data, m_size);
			m_mapped = false;
			m_data = nullptr;

		} else if(m_data) {

			free(m_data);
			m_data = nullptr;

		}

		m_size = 0;
	}

	bool mapped() const {
		return m_mapped;
	}

	void remap(Interleave interleave) {

		if(props().interleave() == interleave)
			return;

		int cols = props().cols();
		int rows = props().rows();
		int bands = props().bands();

		T* mem = (T*) mmap(0, m_size, PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_SHARED, 0, 0);
		if(!mem)
			g_runerr("Failed to reserve space for remap.")

		GridProps nprops(props());
		nprops.setInterleave(interleave);

		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c) {
				for(int b = 0; b < bands; ++b)
					mem[nprops.index(c, r, b)] = m_data[props().index(c, r, b)];
			}
		}

		std::memcpy(m_data, mem, m_size);

		m_props.setInterleave(interleave);

	}

	/**
	 * Vectorize the raster.
	 *
	 * \param filename The filename of the output vector.
	 * \param layerName The name of the output layer.
	 * \param idField A field name for the ID.
	 * \param driver The name of the output driver. Any of the GDAL options.
	 * \param srid The SRID of the projection for the database.
	 * \param band The band to vectorize.
	 * \param removeHoles Remove holes from the polygons.
	 * \param removeDangles Remove small polygons attached to larger ones diagonally.
	 * \param mask The name of a raster file that will be used to set the bounds for vectorization.
	 * \param maskBand The band from the mask raster.
	 * \param threads The number of threads to use.
	 * \param d3 Set to true for 3D geometries; 2D otherwise.
	 * \param fields A list of fields to add to the dataset.
	 * \param status A Status object to receive progress updates.
	 * \param cancel A boolean that will be set to true if the algorithm should quit.
	 */
	void polygonize(const std::string& filename, const std::string& layerName, const std::string& idField,
			const std::string& driver, int srid, int band,
			bool removeHoles = false, bool removeDangles = false,
			const std::string& mask = "", int maskBand = 0, int threads = 1, bool d3 = false,
			const std::vector<std::pair<std::string, OGRFieldType> >& fields = {},
			Monitor* monitor = nullptr) {

		static Monitor _monitor;

		if(!monitor)
			monitor = &_monitor;

		OGRSpatialReference sr;
		sr.importFromEPSG(srid);
		std::vector<char> buf(2048);
		char* bufc = buf.data();
		sr.exportToWkt(&bufc);
		std::string projection(buf.data());

		polygonize(filename, layerName, idField, driver, projection, band, removeHoles, removeDangles,
				mask, maskBand, threads, d3, fields, monitor);
	}

	/**
	 * Vectorize the raster.
	 *
	 * \param filename The filename of the output vector.
	 * \param layerName The name of the output layer.
	 * \param idField A field name for the ID.
	 * \param driver The name of the output driver. Any of the GDAL options.
	 * \param projection The WKT projection for the database.
	 * \param band The band to vectorize.
	 * \param removeHoles Remove holes from the polygons.
	 * \param removeDangles Remove small polygons attached to larger ones diagonally.
	 * \param mask The name of a raster file that will be used to set the bounds for vectorization.
	 * \param maskBand The band from the mask raster.
	 * \param threads The number of threads to use.
	 * \param d3 Set to true for 3D geometries; 2D otherwise.
	 * \param fields A list of fields to add to the dataset.
	 * \param status A Status object to receive progress updates.
	 * \param cancel A boolean that will be set to true if the algorithm should quit.
	 */
	void polygonize(const std::string& filename, const std::string& layerName, const std::string& idField,
			const std::string& driver, const std::string& projection, int band,
			bool removeHoles = false, bool removeDangles = false,
			const std::string& mask = "", int maskBand = 0, int threads = 1, bool d3 = false,
			const std::vector<std::pair<std::string, OGRFieldType> >& fields = {},
			Monitor* monitor = nullptr) {

		static Monitor _monitor;

		if(!monitor) monitor = &_monitor;

		if(isFloat())
			g_runerr("Only int rasters can be polygonized.");

		// It's faster to work on a band-interleaved raster.
		if(props().interleave() == Interleave::BIP)
			remap(Interleave::BIL);

		if(threads < 1)
			threads = 1;

		initGEOS(0, 0);

		// Extract some grid properties.
		int cols = props().cols();
		int rows = props().rows();
		double resX = props().resX();
		double resY = props().resY();
		const Bounds& bounds = props().bounds();

		// Create a buffer for the row.
		std::vector<T> buf(cols);

		// The starting corner coordinates.
		double startX = resX > 0 ? bounds.minx() : bounds.maxx();
		double startY = resY > 0 ? bounds.miny() : bounds.maxy();

		// Create the output dataset
		GEOSContextHandle_t gctx = OGRGeometry::createGEOSContext();

		OGRSpatialReference* sr = nullptr;
		if(!projection.empty())
			sr = new OGRSpatialReference(projection.c_str());

		GDALDataset* ds;
		OGRLayer* layer;
		polyMakeDataset(filename, driver, layerName, idField, sr, d3 ? wkbMultiPolygon25D : wkbMultiPolygon, &ds, &layer);

		if(sr)
			sr->Release();

		if(!fields.empty()) {
			for(const std::pair<std::string, OGRFieldType>& field : fields) {
				if(idField == field.first) {
					g_warn("The ID field matches one of the additional field names. Ignored.")
					continue;
				}
				OGRFieldDefn dfn(field.first.c_str(), field.second);
				layer->CreateField(&dfn);
			}
		}

		double xeps = 0.001 * (resX > 0 ? 1 : -1);
		double yeps = 0.001 * (resY > 0 ? 1 : -1);

		// Data containers.
		std::unordered_map<int, std::vector<GEOSGeometry*> > geoms;
		std::set<int> activeIds;
		std::set<int> finalIds;

		// Thread control features.
		bool running = true;
		std::list<std::thread> ths;

		// Start output threads.
		for(int i = 0; i < threads; ++i)
			ths.emplace_back(&polyWriteToFile, &geoms, &finalIds, &idField, layer, &gctx, removeHoles, removeDangles, monitor);

		// Process raster.
		for(int r = 0; r < rows; ++r) {

			if(monitor->cancel()) break;

			monitor->status((float) r / rows, polyRowStatus(r, rows));

			// Load the row buffer.
			getRow(r, band, buf);

			// Initialize the corner coordinates.
			double x0 = startX;
			double y0 = startY + r * resY;
			double x1 = x0;
			double y1 = y0 + resY;

			// For tracking cell values. TODO: An unsigned int is possible here: overflow.
			int v0 = buf[0];
			int v1 = -1;

			// Reset the list of IDs extant in the current row.
			activeIds.clear();

			for(int c = 1; c < cols; ++c) {

				if(monitor->cancel()) break;

				// If the current cell value differs from the previous one...
				if((v1 = buf[c]) != v0) {
					// Update the right x coordinate.
					x1 = startX + c * resX;
					// If the value is a valid ID, create and the geometry and save it for writing.
					if(v0 > 0) {
						geoms[v0].push_back(polyMakeGeom(x0, y0, x1 + xeps, y1 + yeps, d3 ? 3 : 2));
						activeIds.insert(v0);
					}
					// Update values for next loop.
					v0 = v1;
					x0 = x1;
				}
			}

			// IDs that are in the geoms array and not in the current row are ready to be finalized.
			if(!monitor->cancel()) {
				std::lock_guard<std::mutex> lk0(poly_gmtx);
				std::lock_guard<std::mutex> lk1(poly_fmtx);
				for(const auto& it : geoms) {
					if(activeIds.find(it.first) == activeIds.end())
						finalIds.insert(it.first);
				}
			}
			poly_cv.notify_all();

			while(!monitor->cancel() && !finalIds.empty())
				std::this_thread::yield();
		}

		// Finalize all remaining geometries.
		if(!monitor->cancel()) {
			std::lock_guard<std::mutex> lk0(poly_gmtx);
			std::lock_guard<std::mutex> lk1(poly_fmtx);
			for(auto& it : geoms)
				finalIds.insert(it.first);
		}

		// Let the threads shut down when they run out of geometries.
		running = false;
		poly_cv.notify_all();

		for(std::thread& th : ths) {
			if(th.joinable())
				th.join();
		}

		OGRGeometry::freeGEOSContext(gctx);

		// Release the layer -- GDAL will take care of it. But close the dataset so that can happen.
		layer->Dereference();
		GDALClose(ds);

		monitor->status(1.0, "Finished polygonization");

		finishGEOS();

	}

	/**
	 * Destroy the grid.
	 */
	~Grid() {

		flush();
		destroy();

	}

};

/**
 * Used by flood fill to determine whether a pixel should be filled.
 * Identifies pixels that match a given value.
 */
template <class T, class U>
class G_DLL_EXPORT TargetFillOperator : public FillOperator<T, U> {
private:
	Grid<T>* m_src;
	int m_srcBand;
	Grid<U>* m_dst;
	int m_dstBand;
	T m_target;
	U m_fill;
public:

	/**
	 * Construct a TargetFillOperator to fill a different raster.
	 *
	 * \param src The source raster.
	 * \param dst The destination raster.
	 * \param target The target value.
	 * \param fill The fill value.
	 * \param band The band to fill.
	 */
	TargetFillOperator(Grid<T>* src, int srcBand, Grid<U>* dst, int dstBand, T target, U fill) :
		m_src(src), m_srcBand(srcBand),
		m_dst(dst), m_dstBand(dstBand),
		m_target(target), m_fill(fill) {
	}

	/**
	 * Construct a TargetFillOperator. To fill the same raster.
	 *
	 * \param grd The source and destination raster.
	 * \param srcBand The source band.
	 * \param dstBand The destination band.
	 * \param target The target value.
	 * \param fill The fill value.
	 */
	TargetFillOperator(Grid<T>* grd, int srcBand, int dstBand, T target, U fill) :
		m_src(grd), m_srcBand(srcBand),
		m_dst(grd), m_dstBand(dstBand),
		m_target(target), m_fill(fill) {
	}

	/**
	 * Return the source grid properties.
	 *
	 * \return The source grid properties.
	 */
	const GridProps& srcProps() const {
		return m_src->props();
	}

	/**
	 * Return the destination grid properties.
	 *
	 * \return The destination grid properties.
	 */
	const GridProps& dstProps() const {
		return m_dst->props();
	}

	/**
	 * Return true if the given cell should be filled.
	 *
	 * \param col A column index.
	 * \param row A row index.
	 * \return True if the given cell should be filled.
	 */
	bool shouldFill(int col, int row) const {
		return m_src->get(col, row, m_srcBand) == m_target;
	}

	/**
	 * Fill the given cell.
	 *
	 * \param col A column index.
	 * \param row A row index.
	 */
	void fill(int col, int row) const {
		m_dst->set(col, row, (U) m_fill, m_dstBand);
	}

	~TargetFillOperator() {}
};


} // grid
} // geo

namespace {

	// Write to the file from the map of geometry lists.
	// On each loop, extracts a single finalized poly ID and loads those polys for unioning.
	void polyWriteToFile(std::unordered_map<int, std::vector<GEOSGeometry*> >* geoms, std::set<int>* finalIds,
			const std::string* idField, OGRLayer* layer, GEOSContextHandle_t* gctx,
			bool removeHoles, bool removeDangles, bool* running, bool* cancel) {

		std::vector<GEOSGeometry*> polys;
		GEOSGeometry* geom = nullptr;
		int id;

		while(!*cancel && (*running || !finalIds->empty())) {

			// Get an ID and the list of polys from the queue.
			{
				std::unique_lock<std::mutex> lk(poly_fmtx);
				// Wait for a notification if the queue is empty.
				while(!*cancel && *running && finalIds->empty())
					poly_cv.wait(lk);
				// If the wakeup is spurious, skip.
				if(finalIds->empty())
					continue;
				// Get the ID.
				id = *(finalIds->begin());
				finalIds->erase(id);
			}
			// Get the list of polys and remove the list from the map.
			{
				std::lock_guard<std::mutex> lk(poly_gmtx);
				// If the geoms list is still here, grab it.
				if(geoms->find(id) == geoms->end())
					continue;
				polys = geoms->at(id);
				geoms->erase(id);
			}

			if(polys.empty())
				continue;

			if(*cancel) {
				for(GEOSGeometry* p : polys)
					GEOSGeom_destroy(p);
				continue;
			}

			// Union the polys.
			geom = polys.front();
			for(size_t i = 1; i < polys.size(); ++i) {
				GEOSGeometry* tmp = GEOSUnion_r(*gctx, geom, polys[i]); //
				GEOSGeom_destroy(geom);
				GEOSGeom_destroy(polys[i]);
				geom = tmp;
			}
			polys.clear();

			if(*cancel) {
				GEOSGeom_destroy(geom);
				continue;
			}

			// If we're removing dangles, throw away all but the
			// largest single polygon. If it was originally a polygon, there are no dangles.
			size_t numGeoms;
			if(!*cancel && removeDangles && (numGeoms = GEOSGetNumGeometries(geom)) > 1) {
				size_t idx = 0;
				double a, area = 0;
				for(size_t i = 0; i < numGeoms; ++i) {
					const GEOSGeometry* p = GEOSGetGeometryN(geom, i);
					GEOSArea(p, &a);
					if(a > area) {
						area = a;
						idx = i;
					}
				}
				GEOSGeometry *g = GEOSGeom_clone(GEOSGetGeometryN(geom, idx)); // Force copy.
				GEOSGeom_destroy(geom);
				geom = g;
			}

			// If we're removing holes, extract the exterior rings of all constituent polygons.
			if(!*cancel && removeHoles) {
				std::vector<GEOSGeometry*> geoms0;
				for(int i = 0; i < GEOSGetNumGeometries(geom); ++i) {
					const GEOSGeometry* p = GEOSGetGeometryN(geom, i);
					const GEOSGeometry* l = GEOSGetExteriorRing(p);
					const GEOSCoordSequence* seq = GEOSGeom_getCoordSeq(l);
					GEOSGeometry* r = GEOSGeom_createLinearRing(GEOSCoordSeq_clone(seq));
					GEOSGeometry* npoly = GEOSGeom_createPolygon(r, 0, 0);
					geoms0.push_back(npoly);
				}
				GEOSGeometry* g = GEOSGeom_createCollection(GEOSGeomTypes::GEOS_MULTIPOLYGON, geoms0.data(), geoms0.size()); // Do not copy -- take ownership.
				GEOSGeom_destroy(geom);
				geom = g;
			}

			// If the result is not a multi, make it one.
			if(!*cancel && GEOSGeomTypeId(geom) != GEOSGeomTypes::GEOS_MULTIPOLYGON) {
				std::vector<GEOSGeometry*> gs;
				gs.push_back(geom);
				geom = GEOSGeom_createCollection(GEOSGeomTypes::GEOS_MULTIPOLYGON, gs.data(), gs.size());
			}

			if(!geom)
				g_runerr("Null geometry.");

			if(*cancel) {
				GEOSGeom_destroy(geom);
				continue;
			}

			// Create and write the OGR geometry.
			OGRGeometry* ogeom = OGRGeometryFactory::createFromGEOS(*gctx, (GEOSGeom) geom);
			GEOSGeom_destroy(geom);

			// Create and configure the feature. Feature owns the OGRGeometry.
			OGRFeatureDefn* fdef = layer->GetLayerDefn();
			OGRFeature* feat = OGRFeature::CreateFeature(fdef); // Creates on the OGR heap.
			feat->SetGeometryDirectly(ogeom);
			feat->SetField(idField->c_str(), (GIntBig) id);
			feat->SetFID(++poly_fid);

			// Write to the output file.
			int err;
			{
				std::lock_guard<std::mutex> lk(poly_omtx);
				err = layer->CreateFeature(feat);
			}

			OGRFeature::DestroyFeature(feat);

			if(OGRERR_NONE != err)
				g_runerr("Failed to add geometry.");
		}
	}

	// Produce a rectangular polygon from the four corners.
	GEOSGeometry* polyMakeGeom(double x0, double y0, double x1, double y1, int dims) {

		// Build the geometry.
		GEOSCoordSequence* seq = GEOSCoordSeq_create(5, dims);
		GEOSCoordSeq_setX(seq, 0, x0);
		GEOSCoordSeq_setY(seq, 0, y0);
		GEOSCoordSeq_setX(seq, 1, x0);
		GEOSCoordSeq_setY(seq, 1, y1);
		GEOSCoordSeq_setX(seq, 2, x1);
		GEOSCoordSeq_setY(seq, 2, y1);
		GEOSCoordSeq_setX(seq, 3, x1);
		GEOSCoordSeq_setY(seq, 3, y0);
		GEOSCoordSeq_setX(seq, 4, x0);
		GEOSCoordSeq_setY(seq, 4, y0);
		GEOSGeometry* ring = GEOSGeom_createLinearRing(seq);
		GEOSGeometry* poly = GEOSGeom_createPolygon(ring, 0, 0);
		return poly;
	}

	// Make or open an OGR database to write the polygons to.
	void polyMakeDataset(const std::string& filename, const std::string& driver, const std::string& layerName, const std::string& idField,
			OGRSpatialReference* sr, OGRwkbGeometryType gType,
			GDALDataset** ds, OGRLayer** layer) {

		*ds = (GDALDataset*) GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE, nullptr, nullptr, nullptr);

		std::string drvl = lowercase(driver);

		if(!*ds) {
			// Get the vector driver.
			GDALDriver *drv = GetGDALDriverManager()->GetDriverByName(driver.c_str());
			if(!drv)
				g_runerr("Failed to find driver for " << driver << ".");

			// Create an output dataset for the polygons.
			char** dopts = NULL;
			if(drvl == "sqlite")
				dopts = CSLSetNameValue(dopts, "SPATIALITE", "YES");

			*ds = drv->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, dopts);
			CPLFree(dopts);
			if(!*ds)
				g_runerr("Failed to create dataset " << filename << ".");

		}

		// Create the layer.
		char** lopts = NULL;
		if(drvl == "sqlite") {
			lopts = CSLSetNameValue(lopts, "FORMAT", "SPATIALITE");
		} else if(drvl == "esri shapefile") {
			lopts = CSLSetNameValue(lopts, "2GB_LIMIT", "YES");
		}
		*layer = (*ds)->CreateLayer(layerName.c_str(), sr, gType, lopts);
		CPLFree(lopts);

		if(!*layer)
			g_runerr("Failed to create layer " << layerName << ".");

		// There's only one field -- an ID.
		OGRFieldDefn field(idField.c_str(), OFTInteger);
		(*layer)->CreateField(&field, TRUE);

	}

	// Produce a status message from the current and total row counds.
	std::string polyRowStatus(int r, int rows) {
		std::stringstream ss;
		ss << "Polygonizing row " << (r + 1) << " of " << rows;
		return ss.str();
	}

}



namespace {

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
	template <class T>
	void smooth(Grid<T>* source, Grid<T>* smoothed,
			int size, double nodata, double* weights, int* row,
			std::mutex* rmtx, std::mutex* wmtx, geo::Monitor* monitor = nullptr) {

		/*
		try {
			std::unique_ptr<Tile<T>> tile;
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
				Grid<T>& grid = tile->grid();
				const GridProps& props = grid.props();
				Grid<T> buf(props);
				grid.writeTo(buf);

				// Process the entire block, even the buffer parts.
				double v, t;
				for (int r = 0; !*cancel && r < props.rows(); ++r) {
					for (int c = 0; !*cancel && c < props.cols(); ++c) {
						t = 0.0;
						for (int gr = 0; gr < size; ++gr) {
							for (int gc = 0; gc < size; ++gc) {
								int cc = c - size / 2 + gc;
								int rr = r - size / 2 + gr;
								if(props.hasCell(cc, rr) && (v = buf.get(cc, rr, 1)) != nodata)
									t += weights[gr * size + gc] * v;
							}
						}
						grid.set(c, r, t, 1);
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
		*/
	}

	using namespace geo::grid;

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
			g_runerr("No size for type: " << (int) type);
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

	void fixWriteBounds(int& cols, int& rows, int& srcCol, int& srcRow, int& dstCol, int& dstRow, int rcols, int rrows, int gcols, int grows) {

		if(cols <= 0)
			cols = gcols;
		if(rows <= 0)
			rows = grows;

		if(srcCol < 0) {
			dstCol -= srcCol;
			cols += srcCol;
			srcCol = 0;
		}
		if(srcRow < 0) {
			dstRow -= srcRow;
			rows += srcRow;
			srcRow = 0;
		}
		if(dstCol < 0) {
			cols -= dstCol;
			dstCol = 0;
		}
		if(dstRow < 0) {
			rows -= dstRow;
			dstRow = 0;
		}

		if(srcCol + cols >= rcols)
			cols = rcols - srcCol;
		if(srcRow + rows >= rrows)
			rows = rrows - srcRow;

		if(dstCol + cols > gcols)
			cols = gcols - dstCol;
		if(dstRow + rows > grows)
			rows = grows - dstRow;

	}

	bool fixCoords(int& srcCol, int& srcRow, int& dstCol, int& dstRow, int& cols, int& rows, int srcCols, int srcRows, int dstCols, int dstRows) {

		if(cols <= 0) cols = srcCols;
		if(rows <= 0) rows = srcRows;

		if(srcCol >= srcCols || srcRow >= srcRows || srcCol + cols < 0 || srcRow + rows < 0) {
			g_warn("Col/row out of range." << srcCol << ", " << srcRow << ", " << srcCols << ", " << srcRows << ", " << cols << ", " << rows);
			return false;
		}
		if(srcCol < 0) {
			cols += srcCol;
			srcCol = 0;
		}
		if(srcRow < 0) {
			rows += srcRow;
			srcRow = 0;
		}
		if(srcCol + cols > srcCols) {
			cols = srcCols - srcCol;
		}
		if(srcRow + rows > srcRows) {
			rows = srcRows - srcRow;
		}
		if(dstCol < 0) {
			cols += dstCol;
			dstCol = 0;
		}
		if(dstRow < 0) {
			rows += dstRow;
			dstRow = 0;
		}
		if(dstCol + cols > dstCols) {
			cols = dstCols - dstCol;
		}
		if(dstRow + rows > dstRows) {
			rows = dstRows - dstRow;
		}
		if(cols <= 0 || rows <= 0) {
			g_warn("Zero copy area.")
			return false;
		}
		return true;
	}


} // anonymous



#endif
