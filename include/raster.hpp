/**
 * Raster.hpp
 *
 *  Created on: Jan 21, 2016
 *  Author: rob
 */

#ifndef INCLUDE_RASTER_HPP_
#define INCLUDE_RASTER_HPP_

#include <queue>
#include <stdexcept>
#include <map>
#include <unordered_map>
#include <vector>
#include <cstring>
#include <string>
#include <memory>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <eigen3/Eigen/Core>

#include "geotools.hpp"
#include "util.hpp"

using namespace geotools::util;

namespace geotools {

    namespace raster {

		enum DataType {
			Float64 = 7, Float32 = 6, UInt32 = 5, UInt16 = 4, Byte = 3, Int32 = 2, Int16 = 1, None = 0
		};


    	class G_DLL_EXPORT GridProps {
    	private:
    		double m_trans[6];			// The geotransform properties.
    		int m_cols, m_rows;			// The number of rows and columns.
    		int m_vsrid, m_hsrid;		// The vertical and horizontal SRIDs
            int m_bands;           		// The number of bands.
            bool m_writable;            // True if the raster is writable
            double m_nodata;			// The nodata value.
    		DataType m_type;
    		std::string m_projection;	// The WKT representation of the projection
    		std::string m_driver;

    	public:

    		GridProps();

    		Bounds bounds() const;

    		std::string driver() const;
    		void setDriver(const std::string &name);

    		bool isInt() const;
    		bool isFloat() const;

    		double nodata() const;

    		void setNoData(double nodata);

    		// Return the number of columns.
            int cols() const;

            // Return the number of rows.
            int rows() const;

            // Returns true if the cell is in the raster.
            bool hasCell(int col, int row) const;
            bool hasCell(double x, double y) const;

            // Returns the row for a given y-coordinate.
            int toRow(double y) const;

            // Returns the column for a given x-coordinate.
            int toCol(double x) const;

            // Returns the x-coordinate for a given column.
            double toX(int col) const;

            // Returns the y-coordinate for a given row.
            double toY(int row) const;

            // Returns the x-coordinate for the cell centroid of a given column.
            double toCentroidX(int col) const;

            // Returns the y-coordinate for the cell centorid of a given row.
            double toCentroidY(int row) const;

            // Returns the number of pixels.
            long size() const;

            // Set the data type of the raster.
    		void setDataType(DataType type);

    		// Get the data type of the raster.
    		DataType dataType() const;

    		// Set the size of the raster in columns, rows.
    		void setSize(int cols, int rows);

    		// Set the horizontal and vertical (optional) SRID.
    		void setSrid(int hsrid, int vsrid = 0);

    		// Get the vertical SRID.
    		int vsrid() const;

    		// Get the horizontal SRID.
    		int hsrid() const;

    		// Set the WKT projection.
    		void setProjection(const std::string &proj);

    		// Get the WKT projection.
    		std::string projection() const;

    		// Set the geo transform properties.
    		void setTrans(double m_trans[6]);

    		// Get the geo transform properties.
    		void trans(double m_trans[6]) const;

    		// Set the vertical and horizontal resolution.
    		void setResolution(double resolutionX, double resolutionY);

    		// Get the horizontal resolution.
    		double resolutionX() const;

    		// Get the vertical resolution.
    		double resolutionY() const;

    		double tlx() const;

    		double tly() const;

    		// Set the number of bands.
    		void setBands(int bands);

    		// Get the number of bands.
    		int bands() const;

    		// Set the writable state of the raster.
    		void setWritable(bool writable);

    		// Get the writable state of the raster.
    		bool writable() const;
    	};

    	class G_DLL_EXPORT GridStats {
    	public:
            double min;
            double max;
            double mean;
            double stdDev;
            double variance;
            double sum;
            long count;
    	};

        // Simple class to represent a single grid cell.
        class G_DLL_EXPORT Cell {
        public:
            int col;
            int row;
            Cell(int col, int row);
        };

        // Used by Grid::floodFill to determine whether
        // a pixel should be filled.
        class G_DLL_EXPORT FillOperator {
        public:
            virtual bool fill(int value) const = 0;
            virtual ~FillOperator() = 0;
        };

        class G_DLL_EXPORT TargetOperator : public FillOperator {
        private:
            int m_match;
        public:
            TargetOperator(int match);

            bool fill(int value) const;

            ~TargetOperator();
        };

        // Abstract class for grids (rasters).
        class G_DLL_EXPORT Grid {
        public:
            Grid();

            virtual ~Grid() = 0;
            
            // Compute the table of Gaussian weights given the size of the table
            // and the std. deviation.
            static void gaussianWeights(double *weights, int size, double sigma);

            GridStats stats();

            // Returns the grid properties
            virtual const GridProps &props() const = 0;

            // Fill the entire dataset with the given value.
            virtual void fillFloat(double value, int band = 1) = 0;
            virtual void fillInt(int value, int band = 1) = 0;

            // Return a the value held at the given index in the grid.
            virtual int getInt(long idx, int band = 1) = 0;
            virtual int getInt(int col, int row, int band = 1) = 0;
            virtual double getFloat(long idx, int band = 1) = 0;
            virtual double getFloat(int col, int row, int band = 1) = 0;

            // Set the value held at  the given index in the grid.
            virtual void setInt(long idx, int value, int band = 1) = 0;
            virtual void setInt(int col, int row, int value, int band = 1) = 0;
            virtual void setFloat(long idx, double value, int band = 1) = 0;
            virtual void setFloat(int col, int row, double value, int band = 1) = 0;

            // Write data from Grid instance.
            virtual void write(Grid &grd,
            		int cols = 0, int rows = 0,
            		int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
					int srcBand = 1, int dstBand = 1) = 0;

            // Normalize the grid so that one standard deviation is +-1.
            void normalize(int band = 1);

            // Normalize the grid so that the max value is equal to 1, and
            // the minimum is zero.
            void logNormalize(int band = 1);

            // Convert a Grid to some other type.
            void convert(Grid &g, int srcBand = 1, int dstBand = 1);

            // Fill the grid, beginning with the target cell, where any contiguous cell
            // satisfies the given FillOperator. The other grid is actually filled,
            // and the present grid is unchanged *unless* the present grid is passed
            // as other.
            // col, row -- The column and row to start on.
            // op       -- A FillOperator instance which will determine
            //             whether a pixel should be filled.
            // other    -- The grid whose cells will actually be filled.
            // fill     -- The value to fill cells with.
            // d8       -- Whether to enable diagonal fills.
            // out*     -- Pointer to variables that hold min and max rows and columns
            //             plus the area of the fill's bounding box.
            void floodFill(int col, int row,
                FillOperator &op, Grid &other, int fill, bool d8 = false,
				int *outminc = nullptr, int *outminr = nullptr,
				int *outmaxc = nullptr, int *outmaxr = nullptr,
				int *outarea = nullptr);

            // Begin flood fill at the given cell; fill cells equal to the target value.
            void floodFill(int col, int row, int target, int fill, bool d8 = false,
    				int *outminc = nullptr, int *outminr = nullptr,
    				int *outmaxc = nullptr, int *outmaxr = nullptr,
    				int *outarea = nullptr);

            // Begin flood fill at the given cell; fill cells that satisfy the operator.
            void floodFill(int col, int row, FillOperator &op, int fill, bool d8 = false,
    				int *outminc = nullptr, int *outminr = nullptr,
    				int *outmaxc = nullptr, int *outmaxr = nullptr,
    				int *outarea = nullptr);

            // Smooth the raster and write the smoothed version to the output raster.
            // Callback is an optional function reference with a single float
            // between 0 and 1, for status tracking.
            void smooth(Grid &smoothed, double sigma, int size, int band = 1,
                geotools::util::Callbacks *status = nullptr, 
                bool *cancel = nullptr);

            // The radius is given with cells as the unit, but
            // can be rational. When determining which cells to
            // include in the calculation, any cell which partially
            // falls in the radius will be included.
            void voidFillIDW(double radius, int count = 4, double exp = 2.0, int band = 1);

        };

        class Raster;

        // A convenience class for managing a grid of values.
        // Handles allocation and deallocation of memory.
        class G_DLL_EXPORT MemRaster : public Grid {
        private:
            void *m_grid;
            bool m_mmapped;
            GridProps m_props;
            std::unique_ptr<geotools::util::MappedFile> m_mappedFile;
            std::unique_ptr<boost::interprocess::mapped_region> m_region;
            std::unique_ptr<boost::interprocess::file_mapping> m_mapping;

            // Checks if the grid has been initialized. Throws exception otherwise.
            void checkInit() const;

            void freeMem();

        public:
            MemRaster();

            MemRaster(const GridProps &props, bool mapped = false);

            ~MemRaster();

            // Return a pointer to the allocated memory.
            void *grid();

            const GridProps& props() const;

            // Initialize with the given number of cols and rows.
            // (Re)allocates memory for the internal grid.
            void init(const GridProps &props, bool mapped = false);

            // Fill the entire dataset with the given value.
            void fillFloat(double value, int band = 1);
            void fillInt(int value, int band = 1);

            // Return a the value held at the given index in the grid.
            int getInt(long idx, int band = 1);
            int getInt(int col, int row, int band = 1);
            double getFloat(long idx, int band = 1);
            double getFloat(int col, int row, int band = 1);

            // Set the value held at  the given index in the grid.
            void setInt(long idx, int value, int band = 1);
            void setInt(int col, int row, int value, int band = 1);
            void setFloat(long idx, double value, int band = 1);
            void setFloat(int col, int row, double value, int band = 1);

			// Write data from Grid instance.
			void write(Grid &grd,
					int cols = 0, int rows = 0,
					int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
					int srcBand = 1, int dstBand = 1);
            void writeMemRaster(MemRaster &grd,
            		int cols = 0, int rows = 0,
            		int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
            		int srcBand = 1, int dstBand = 1);
            void writeRaster(Raster &grd,
            		int cols = 0, int rows = 0,
            		int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
					int srcBand = 1, int dstBand = 1);

            // Convert the grid to matrix.
            void toMatrix(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &mtx, int band = 1);

            // Initialize the grid from a matrix.
            void fromMatrix(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &mtx, int band = 1);

        };

        class G_DLL_EXPORT Raster : public Grid {
        	friend class MemRaster;
        private:
            GDALDataset *m_ds;          // GDAL dataset
            int m_bcols, m_brows;		// The size of the GDAL block.
            int m_bcol, m_brow;			// The current loaded block position.
            int m_band;
    		void *m_block;
            std::string m_filename;     // Raster filename
            GridProps m_props;
            GDALDataType m_type;        // GDALDataType -- limits the possible template types.

            GDALDataType getGDType() const;

        protected:
            GDALDataset* ds() const;

        public:

            // Create a new raster for writing with a template.
            Raster(const std::string &filename, const GridProps &props);

            // Open the given raster. Set the writable argument to true
            // to enable writing.
            Raster(const std::string &filename, bool writable = false);

            // Return the grid properties object.
            const GridProps& props() const;

            // Attempts to return the datatype of the raster
            // with the given filename.
            static DataType getFileDataType(const std::string &filename);

            // Return a map containing the raster driver short name and extension.
            static std::map<std::string, std::set<std::string> > extensions();

            // Return a map containing the raster driver short name and long name.
            static std::map<std::string, std::string> drivers();

            static std::string getDriverForFilename(const std::string &filename);

            // Return the filename for this raster.
            std::string filename() const;

            // Fill the given band with the given value.
            void fillInt(int value, int band = 1);
            void fillFloat(double value, int band = 1);

			// Write data from Grid instance.
			void write(Grid &grd,
					int cols = 0, int rows = 0,
					int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
					int srcBand = 1, int dstBand = 1);
            void writeMemRaster(MemRaster &grd,
            		int cols = 0, int rows = 0,
            		int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
            		int srcBand = 1, int dstBand = 1);
            void writeRaster(Raster &grd,
            		int cols = 0, int rows = 0,
            		int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
					int srcBand = 1, int dstBand = 1);

            // Returns a pixel value.
            int getInt(double x, double y, int band = 1);
            int getInt(int col, int row, int band = 1);
            int getInt(long idx, int band = 1);

            double getFloat(double x, double y, int band = 1);
            double getFloat(int col, int row, int band = 1);
            double getFloat(long idx, int band = 1);

            // Set an pixel value.
            void setInt(double x, double y, int v, int band = 1);
            void setInt(int col, int row, int v, int band = 1);
            void setInt(long idx, int v, int band = 1);

            void setFloat(double x, double y, double v, int band = 1);
            void setFloat(int col, int row, double v, int band = 1);
            void setFloat(long idx, double v, int band = 1);

            // Returns true if the raster is square.
            bool isSquare() const;

            // Flush the current block to the dataset.
            void flush();

            // Vectorize the raster.
            void polygonize(const std::string &filename, const std::string &layerName, 
                const std::string &driver, uint16_t srid = 0, uint16_t band = 1,
				geotools::util::Status *status = nullptr, bool *cancel = nullptr);

            ~Raster();

        };

    } // raster

} // geotools


#endif
