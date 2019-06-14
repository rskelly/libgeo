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
#include <unordered_set>
#include <queue>
#include <vector>
#include <cstring>
#include <string>
#include <memory>
#include <mutex>
#include <set>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <ogr_geometry.h>

#include <Eigen/Core>

#include "util.hpp"

using namespace geo::util;

namespace geo {

    namespace raster {

    	/**
    	 * The allowable types for a raster.
    	 */
		enum DataType {
			Float64 = 7, Float32 = 6, UInt32 = 5, UInt16 = 4, Byte = 3, Int32 = 2, Int16 = 1, None = 0
		};

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
    		std::string m_projection;	///<! The WKT representation of the projection
    		std::string m_driver;		///<! The name of the GDAL driver.

		public:

    		/**
    		 * Construct an empty GridProps.
    		 */
    		GridProps();

    		/**
    		 * Return the geographic bounds of the raster.
    		 *
    		 * @return The geographic bounds of the raster.
    		 */
    		Bounds bounds() const;

    		/**
    		 * Set the geographic bounds of the raster.
    		 *
    		 * @param bounds The geographic bounds of the raster.
    		 */
    		void setBounds(const Bounds& bounds);

    		/**
    		 * Use compression for tiff files.
    		 *
    		 * @param compress True to use compression for tiff files.
    		 */
    		void setCompress(bool compress);

    		/**
    		 * Use compression for tiff files.
    		 *
    		 * @return True to use compression for tiff files.
    		 */
    		bool compress() const;

    		/**
    		 * Use Big Tiff setting.
    		 *
    		 * @param bigTuff True to use Big Tiff setting.
    		 */
    		void setBigTiff(bool bigTiff);

    		/**
    		 * Use Big Tiff setting.
    		 *
    		 * @return True to use Big Tiff setting.
    		 */
    		bool bigTiff() const;

    		/**
    		 * Populate an (at least) 4-element double array with the bounding
    		 * box of this object.
    		 *
    		 * @param bounds A four-element double array.
    		 */
    		void bounds(double* bounds) const;

    		/**
    		 * Return the name of the GDAL driver used by the raster.
    		 * Only relevant for file-based rasters.
    		 *
    		 * @return The name of the GDAL driver.
    		 */
    		std::string driver() const;

    		/**
    		 * Set the name of the GDAL driver used by the raster.
    		 * Only relevant for file-based rasters.
    		 *
    		 * @param name The name of the driver.
    		 */
    		void setDriver(const std::string &name);

    		/**
    		 * Returns true if the raster contains integer values.
    		 */
    		bool isInt() const;

    		/**
    		 * Returns true if the raster contains double values.
    		 */
    		bool isFloat() const;

    		/**
    		 * Returns the no data value.
    		 *
    		 * @return The no data value.
    		 */
    		double nodata() const;

    		/**
    		 * Set the no data value.
    		 *
    		 * @param nodata The no data value.
    		 */
    		void setNoData(double nodata);

    		/**
    		 * Returns true if the no data value has been set.
    		 *
    		 * @return True if the no data value has been set.
    		 */
    		bool nodataSet() const;

    		/**
    		 * Remove the no data value.
    		 */
    		void unsetNodata();

    		/**
    		 * Return the number of columns.
    		 *
    		 * @return The number of columns.
    		 */
            int cols() const;

            /*
             * Return the number of rows.
             *
             * @param The number of rows.
             */
            int rows() const;

            /**
             * Returns true if the cell is in the raster.
             *
             * @param col The column.
             * @param row The row.
             * @return True if the cell is in the raster.
             */
            bool hasCell(int col, int row) const;

            /**
             * Returns true if the cell is in the raster.
             *
             * @param x The geographic x or longitude coordinate.
             * @param y The geographic y or latitude coordinate.
             * @return True if the cell is in the raster.
             */
            bool hasCell(double x, double y) const;

            /**
             * Returns the row for a given y-coordinate.
             *
             * @param y The geographic y or latitude coordinate.
             * @return The row index.
             */
            int toRow(double y) const;

            /**
             * Returns the column for a given x-coordinate.
             *
             * @param x The geographic x or longitude coordinate.
             * @return The column index.
             */
            int toCol(double x) const;

            /**
             * Returns the x-coordinate at the minimum corner of a given column.
             *
             * @param col The column.
             * @return The x-coordinate of the column corner.
             */
            double toX(int col) const;

            /**
             * Returns the y-coordinate at the minimum corner of a given row.
             *
             * @param row The row.
             * @return The y-coordinate of the row corner.
             */
            double toY(int row) const;

            /**
             * Returns the x-coordinate for the cell centroid of a given column.
             *
             * @param col The column.
             * @return The x-coordinate at the centre of the column.
             */
            double toCentroidX(int col) const;

            /**
             * Returns the y-coordinate for the cell centorid of a given row.
             *
             * @param row The row.
             * @return The y-coordinate at the centre of the row.
             */
            double toCentroidY(int row) const;

            /**
             * Returns the number of pixels in the raster.
             *
             * @return The number of pixels in the raster.
             */
            size_t size() const;

            /**
             * Set the data type of the raster.
             *
             * @param type The data type.
             */
    		void setDataType(DataType type);

    		/**
    		 * Get the data type of the raster.
    		 *
    		 * @return The data type.
    		 */
    		DataType dataType() const;

    		/**
    		 * Set the size of the raster in columns, rows.
    		 *
             * @param col The column.
             * @param row The row.
    		 */
    		void setSize(int cols, int rows);

    		/**
    		 * Set the horizontal and vertical (optional) SRID.
    		 *
    		 * @param hsrid The horizontal SRID.
    		 * @param vsrid The vertical SRID.
    		 */
    		void setSrid(int hsrid, int vsrid = 0);

    		/**
    		 * Get the vertical SRID.
    		 *
    		 * @return The vertical SRID.
    		 */
    		int vsrid() const;

    		/**
    		 * Get the horizontal SRID.
    		 *
    		 * @return The horizontal SRID.
    		 */
    		int hsrid() const;

    		/**
    		 * Set the WKT projection.
    		 *
    		 * @param The projection string (proj or WKT format).
    		 */
    		void setProjection(const std::string &proj);

    		/**
    		 * Get the WKT projection (proj or WKT format).
    		 *
    		 * @return The WKT projection.
    		 */
    		std::string projection() const;

    		/**
    		 * Set the geo transform properties.
    		 *
    		 * @param trans The six-element transformation matrix.
    		 */
    		void setTrans(double trans[6]);

    		/**
    		 * Set the geo transform properties. The third and fifth elements are set to zero.
    		 *
    		 * @param tlx The top left x-coordinate.
    		 * @param resX The horizontal resolution.
    		 * @param tly The top left y-coordinate.
    		 * @param resY The vertical resolution (negative for UTM (etc.) projections).
    		 */
    		void setTrans(double tlx, double resX, double tly, double resY);

    		/**
    		 * Gets the geo transform properties by setting them in the given array.
    		 *
    		 * @param trans The six-element transformation matrix.
    		 */
    		void trans(double trans[6]) const;

    		/**
    		 * Set the vertical and horizontal resolution.
    		 *
    		 * @param resolutionX The horizontal resolution.
    		 * @param resolutionY The vertical resolution (negative for UTM (etc.) projections).
    		 */
    		void setResolution(double resolutionX, double resolutionY);

    		/**
    		 * Get the horizontal resolution.
    		 *
    		 * @return The horizontal resolution.
    		 */
    		double resolutionX() const;

    		/**
    		 * Get the vertical resolution.
    		 *
    		 * @return The vertical resolution.
    		 */
    		double resolutionY() const;

    		/**
    		 * Return the top-left horizontal coordinate of the raster.
    		 *
    		 * @return The top-left horizontal coordinate of the raster.
    		 */
    		double tlx() const;

    		/**
    		 * Return the top-left vertical coordinate of the raster.
    		 *
    		 * @return The top-left vertical coordinate of the raster.
    		 */
    		double tly() const;

    		/**
    		 * Set the number of bands.
    		 *
    		 * @param bands The number of bands.
    		 */
    		void setBands(int bands);

    		/**
    		 * Get the number of bands.
    		 *
    		 * @return The number of bands.
    		 */
    		int bands() const;

    		/**
    		 * Set the writable state of the raster. If the raster is not writable,
    		 * attempting to write to it will throw an exception.
    		 *
    		 * @param writable True, if the raster should be writable.
    		 */
    		void setWritable(bool writable);

    		/**
    		 * Get the writable state of the raster.
    		 *
    		 * @param The writable state of the raster.
    		 */
    		bool writable() const;


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
            Cell(int col, int row);
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
        	 * @return The properties of the source raster.
        	 */
            virtual const GridProps& srcProps() const = 0;

            /**
             * Return the properties of the destination raster.
             *
             * @return The properties of the destination raster.
             */
            virtual const GridProps& dstProps() const = 0;

            /**
             * Return true if if the current pixel should be filled.
             *
             * @param col The column.
             * @param row The row.
             * @return True if if the current pixel should be filled.
             */
            virtual bool shouldFill(int col, int row) const = 0;

            /**
             * Fill the current column.
             *
             * @param col The column.
             * @param row The row.
             */
            virtual void fill(int col, int row) const = 0;

            virtual ~FillOperator() {};
        };

        // Forward declarations; see below.
        class Grid;
        class TileIterator;

        /**
         * A tile represents a discrete region of a raster.
         */
        class G_DLL_EXPORT Tile {
        	friend class TileIterator;
        private:
        	Grid* m_tile;			///<! Contains the tile grid.
        	Grid* m_source;			///<! Contains the source grid.
        	int m_cols;				///<! The number of columns in the tile.
        	int m_rows;				///<! The number of rows in the tile.
        	int m_col; 				///<! The column index of the ROI in the source.
        	int m_row;				///<! The row index of the ROI in the source.
        	int m_buffer;			///<! The number of buffer pixels.
        	int m_srcCol;			///<! The buffered column index in the source.
        	int m_srcRow;			///<! The buffered row index in the source.
        	int m_dstCol;			///<! The column index written to in the tile.
        	int m_dstRow;			///<! The row index written to in the tile.
        	int m_band;				///<! The raster band in the source raster.
        	bool m_writeOnFlush;	///<! If true, writes back to the raster on destruction.

        protected:

        	/**
        	 * Create the tile.
        	 *
        	 * @param tile The tile data.
        	 * @param source The source data.
        	 * @param cols The number of columns in the tile.
        	 * @param rows The number of rows in the tile.
        	 * @param col The source column.
        	 * @param row The row column.
        	 * @param buffer The size of the buffer in pixels.
        	 * @param srcCol The position of the tile in the source, including buffer.
        	 * @param srcRow The position of the tile in the source, including buffer.
        	 * @param dstCol The position of the grid information with respect to the tile.
        	 * @param dstRow The position of the grid information with respect to the tile.
        	 * @param band The raster band in the source.
        	 * @param writeOnFlush If true, writes the tile data to the source on destruction.
        	 */
        	Tile(Grid* tile, Grid* source, int cols, int rows,
        			int col, int row, int buffer,
        			int srcCol, int srcRow,
					int dstCol, int dstRow, int band, bool writeOnFlush);

        public:

        	/**
        	 * Return the source column.
        	 *
        	 * @return The source column.
        	 */
        	int srcCol() const;

        	/**
        	 * Return the destination column.
        	 *
        	 * @return The destination column.
        	 */
        	int dstCol() const;

        	/**
        	 * Return the source row.
        	 *
        	 * @return The source row.
        	 */
        	int srcRow() const;

        	/**
        	 * Return the destination row.
        	 *
        	 * @return The destination row.
        	 */
        	int dstRow() const;

        	/**
        	 * Return the number of columns.
        	 *
        	 * @return The number of columns.
        	 */
        	int cols() const;

        	/**
        	 * Return the number of rows.
        	 *
        	 * @return The number of rows.
        	 */
        	int rows() const;

        	/**
        	 * Return the column index of the tile in the source raster.
        	 *
        	 * @return The column index of the tile in the source raster.
        	 */
        	int col() const;

        	/**
        	 * Return the row index of the tile in the source raster.
        	 *
        	 * @return The row index of the tile in the source raster.
        	 */
        	int row() const;

        	/**
        	 * Return a mutable reference to the grid containing tile data.
        	 *
        	 * @return A mutable reference to the grid containing tile data.
        	 */
        	Grid& grid();

        	/**
        	 * Write the tile's data to the destination grid.
        	 *
        	 * @param dest The destination grid.
        	 */
        	void writeTo(Grid& dest);

        	/**
        	 * Flush tile contents to the source. Called implicitly on destruction.
        	 */
        	void flush();

        	/**
        	 * Destroy the Tile. Implicitly flushes any changes.
        	 */
        	~Tile();
        };

        /**
         * Provides a way of iterating over a raster tile by tile.
         */
        class G_DLL_EXPORT TileIterator {
        	friend class Grid;
        private:
        	Grid& m_source;				///<! The data source.
        	int m_cols;					///<! The number of columns in the tile.
        	int m_rows;					///<! The number of rows in the tile.
        	int m_buffer;				///<! The size of the buffer region around the tile.
        	int m_curCol;				///<! The current column index.
        	int m_curRow;				///<! The current row index.
        	int m_band;					///<! The band index in the data source.

        protected:

        	/**
        	 * Create a TileIterator of the given size. If a buffer is given,
        	 * The tile size in increased, and pixels are read from the source
        	 * to fill the buffer, however only pixels within the unbuffered
        	 * window are written back from a Tile.
        	 *
        	 * @param source The source raster.
        	 * @param cols The number of columns in each tile.
        	 * @param rows The number of rows in each tile.
        	 * @param buffer The size of the buffer around the tile.
        	 * @param band The band in the source raster.
        	 */
        	TileIterator(Grid& source, int cols, int rows, int buffer, int band);

        	/**
        	 * Construct a TileIterator.
        	 */
        	TileIterator();

        public:

        	/**
        	 * Copy constructor.
        	 */
        	TileIterator(const TileIterator& iter);

        	/**
        	 * Returns true if there's another tile to be retrieved.
        	 *
        	 * @return True if there's another tile to be retrieved.
        	 */
        	bool hasNext();

        	/**
        	 * Returns the number of tiles.
        	 *
        	 * @return The number of tiles.
        	 */
        	int count() const;

        	/**
        	 * Returns the next tile or nullptr if there isn't one. The
        	 * caller is responsible for disposing of the returned pointer.
        	 *
        	 * @return The next tile or nullptr if there isn't one.
        	 */
        	Tile* next();

        	/**
        	 * Create a Tile using the given tile as a Template. Does not
        	 * interfere with the iterator. The caller is responsible for
        	 * disposing of the returned pointer.
        	 *
        	 * @param tpl A tile template.
        	 */
        	Tile* create(Tile &tpl);

        	/**
        	 * Destroy the TileIterator.
        	 */
        	~TileIterator();
        };

        /**
         * Abstract class for grids (rasters).
         */
        class G_DLL_EXPORT Grid {
        public:

            /**
             * Destroy the grid.
             */
            virtual ~Grid() {}
            
            /**
             * Return a mutex that can be used to protect the resource.
             *
             * @return A mutex.
             */
            virtual std::mutex& mutex() = 0;

            /**
             * Return a TileIterator.
             *
             * @param cols The number of columns in each tile.
             * @param rows The number of rows in each tile.
             * @param buffer The buffer around the tile.
             * @param band The band in the source raster.
             * @return A TileIterator.
             */
            TileIterator iterator(int cols, int rows, int buffer = 0, int band = 1);

            /**
             * Compute the table of Gaussian weights given the size of
             * the table and the standard deviation.
             *
             * @param weights The list of weights.
             * @param size The size of the weights list.
             * @param sigma The standard deviation.
             * @param mean The centre of the curve.
             */
            static void gaussianWeights(double *weights, int size, double sigma, double mean = 0);

            /**
             * Compute and return the statistics for the band.
             *
             * @param band The raster band.
             * @return A GridStats instance containing computed statistics.
             */
            GridStats stats(int band);

            /**
             * Returns the grid properties.
             *
             * @return A reference to this raster's grid properties.
             */
            virtual const GridProps& props() const = 0;

            /**
             * Fill the entire dataset with the given value.
             *
             * @param value The value to fill the raster with.
             * @param band The band to fill.
             */
            virtual void fillFloat(double value, int band) = 0;

            /**
             * Fill the entire dataset with the given value.
             *
             * @param value The value to fill the raster with.
             * @param band The band to fill.
             */
            virtual void fillInt(int value, int band) = 0;

            /**
             * Return a the value held at the given index in the grid.
             *
             * @param idx The index in the raster; left to right, top to bottom.
             * @param band The band.
             * @return The value held at the given index in the grid.
             */
            //virtual int getInt(size_t idx, int band) = 0;

            /**
             * Return a the value held at the given position in the grid.
             *
             * @param col The column.
             * @param row The row.
             * @param band The band.
             * @return The value held at the given index in the grid.
             */
            virtual int getInt(int col, int row, int band) = 0;

            /**
             * Return a the value held at the given index in the grid.
             *
             * @param idx The index in the raster; left to right, top to bottom.
             * @param band The band.
             * @return The value held at the given index in the grid.
             */
            //virtual double getFloat(size_t idx, int band) = 0;

            /**
             * Return a the value held at the given position in the grid.
             *
             * @param col The column.
             * @param row The row.
             * @param band The band.
             * @param The value held at the given index in the grid.
             */
            virtual double getFloat(int col, int row, int band) = 0;

            /**
             * Set the value held at  the given index in the grid.
             *
             * @param idx The index in the raster; left to right, top to bottom.
             * @param value The value to set.
             * @param band The band.
             */
            //virtual void setInt(size_t idx, int value, int band) = 0;

            /**
             * Set the value held at  the given index in the grid.
             *
             * @param col The column.
             * @param row The row.
             * @param value The value to set.
             * @param band The band.
             */
            virtual void setInt(int col, int row, int value, int band) = 0;

            /**
             * Set the value held at  the given index in the grid.
             *
             * @param idx The index in the raster; left to right, top to bottom.
             * @param value The value to set.
             * @param band The band.
             */
            //virtual void setFloat(size_t idx, double value, int band) = 0;

            /**
             * Set the value held at  the given index in the grid.
             *
             * @param col The column.
             * @param row The row.
             * @param value The value to set.
             * @param band The band.
             */
            virtual void setFloat(int col, int row, double value, int band) = 0;

            /**
             * Write data from the current Grid instance to the given grid.
             *
             * @param grd The target grid.
             * @param cols The number of columns to write.
             * @param rows The number of rows to write.
             * @param srcCol The source column to read from.
             * @param srcRow The source row to read from.
             * @param dstCol The destination column to write to.
             * @param dstRow The destination row to write to.
             * @param srcBand The source band.
             * @param dstBand The destination band.
             */
            virtual void writeTo(Grid &grd,
            		int cols = 0, int rows = 0,
            		int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
					int srcBand = 1, int dstBand = 1) = 0;

            /**
             * Normalize the grid so that one standard deviation is +-1.
             *
             * @param band The target band.
             */
            void normalize(int band);

            /**
             * Normalize the grid so that the max value is equal to 1, and the
             * minimum is zero.
             *
             * @param band The target band.
             */
            void logNormalize(int band);

            /**
             * Convert a Grid to some other type.
             *
             * @param g The destination Grid.
             * @param srcBand The source band.
             * @param dstBand The destination band.
             */
            void convert(Grid &g, int srcBand, int dstBand);

            /**
             * Fill the grid, beginning with the target cell, where any contiguous cell
             * satisfies the given FillOperator. The other grid is actually filled,
             * and the present grid is unchanged *unless* the present grid is passed
             * as other.
             *
             * @param col   The column to start on.
             * @param row   The row to start on.
             * @param op    A FillOperator instance which will determine
             *              whether a pixel should be filled.
             * @param other The grid whose cells will actually be filled.
             * @param fill  The value to fill cells with.
             * @param d8    Whether to enable diagonal fills.
             * @param out*  Pointer to variables that hold min and max rows and columns
             *              plus the area of the fill's bounding box.
             */
            template <class T, class U>
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
             * @param smoothed The smoothed grid.
             * @param sigma    The standard deviation.
             * @param size     The window size.
             * @param band     The target band.
             * @param cancel   A reference to a variable that will be true if the
             *                 process should be cancelled.
             * @param status   A reference to a Status object which will be updated with
             *                 the current status of the process.
			 */
            void smooth(Grid &smoothed, double sigma, int size, int band,
            		bool& cancel, geo::util::Status& status);

            void smooth(Grid &smoothed, double sigma, int size, int band);

            /**
             * Smooth the raster and write the smoothed version to the output raster.
             * Callback is an optional function reference with a single float
             * between 0 and 1, for status tracking.
             *
             * Sigma defaults to 0.84089642, window size to 3.
             *
             * @param smoothed The smoothed grid.
             * @param band     The target band.
             * @param cancel   A reference to a variable that will be true if the
             *                 process should be cancelled.
             * @param status   A reference to a Status object which will be updated with
             *                 the current status of the process.
			 */
            void smooth(Grid &smoothed, int band, bool& cancel, geo::util::Status& status);

            /**
             * The radius is given with cells as the unit, but can be rational.
             * When determining which cells to include in the calculation,
             * any cell which partially falls in the radius will be included.
             *
             * @param filename 	The output filename.
             * @param band   	The source band.
             * @param mask		A raster to use as a mask; invalid pixels will be ignored.
             * @param maskBand  The band to use in the mask.
             * @param radius 	The search radius.
             * @param count  	The number of pixels to use for calculations.
             * @param exp    	The exponent.
			 */
            void voidFillIDW(const std::string& filename, int band, const std::string& mask, int maskBand, double radius, int count = 4, double exp = 2.0);

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
        	 * Return the key for the minimum value in the given map.
        	 *
        	 * @return The key for the minimum value in the given map.
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

            /**
             * Finds the least-cost path from the start cell to the goal cell,
             * using the given heuristic. Populates the given iterator with
             * the optimal path between the start cell and the goal.
             *
             * If the search fails for some reason, like exceeding the maxCost, returns
             * false. Otherwise returns true.
             *
             * @param startCol The starting column.
             * @param startrow The starting row.
             * @param goalCol The column of the goal.
             * @param goalRow The row of the goal.
             * @param heuristic Used by the algorithm to estimate the future cost of the path.
             * @param inserter Used to accumulate the path results.
             * @param maxCost If the total cost exceeds this amount, just quit and return false.
             * @return True if the search succeeded, false otherwise.
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

            /**
             * Vectorize the raster.
             *
             * @param filename The filename of the output vector.
             * @param layerName The name of the output layer.
             * @param idField A field name for the ID.
             * @param driver The name of the output driver. Any of the GDAL options.
             * @param projection The WKT projection for the database.
             * @param band The band to vectorize.
             * @param removeHoles Remove holes from the polygons.
             * @param removeDangles Remove small polygons attached to larger ones diagonally.
             * @param mask The name of a raster file that will be used to set the bounds for vectorization.
             * @param maskBand The band from the mask raster.
             * @param threads The number of threads to use.
             * @param d3 Set to true for 3D geometries; 2D otherwise.
             * @param status A Status object to receive progress updates.
             * @param cancel A boolean that will be set to true if the algorithm should quit.
             */
            void polygonize(const std::string &filename, const std::string &layerName, const std::string& idField,
                const std::string &driver, const std::string& projection, int band, bool removeHoles, bool removeDangles,
				const std::string& mask, int maskBand, int threads, bool d3,
				bool& cancel, geo::util::Status& status);

            /**
             * Vectorize the raster.
             *
             * @param filename The filename of the output vector.
             * @param layerName The name of the output layer.
             * @param idField A field name for the ID.
             * @param driver The name of the output driver. Any of the GDAL options.
             * @param projection The WKT projection for the database.
             * @param band The band to vectorize.
             * @param removeHoles Remove holes from the polygons.
             * @param removeDangles Remove small polygons attached to larger ones diagonally.
             * @param mask The name of a raster file that will be used to set the bounds for vectorization.
             * @param maskBand The band from the mask raster.
             * @param threads The number of threads to use.
             * @param d3 Set to true for 3D geometries; 2D otherwise.
             * @param fields A list of pairs of names of columns to add. Will not be populated.
             * @param status A Status object to receive progress updates.
             * @param cancel A boolean that will be set to true if the algorithm should quit.
             */
            void polygonize(const std::string &filename, const std::string &layerName, const std::string& idField,
                const std::string &driver, const std::string& projection, int band, bool removeHoles, bool removeDangles,
				const std::string& mask, int maskBand, int threads, bool d3, const std::vector<std::pair<std::string, OGRFieldType> >& fields,
				bool& cancel, geo::util::Status& status);
        };

        /**
         * Used by flood fill to determine whether a pixel should be filled.
         * Identifies pixels that match a given value.
         */
        template <class T, class U>
        class G_DLL_EXPORT TargetFillOperator : public FillOperator<T, U> {
        private:
            Grid* m_src;
            int m_srcBand;
            Grid* m_dst;
            int m_dstBand;
            bool m_sint;
            bool m_dint;
            T m_target;
            U m_fill;
        public:

            /**
             * Construct a TargetFillOperator to fill a different raster.
             *
             * @param src The source raster.
             * @param dst The destination raster.
             * @param target The target value.
             * @param fill The fill value.
             * @param band The band to fill.
             */
            TargetFillOperator(Grid* src, int srcBand, Grid* dst, int dstBand, T target, U fill) :
                m_src(src), m_srcBand(srcBand),
				m_dst(dst), m_dstBand(dstBand),
				m_target(target), m_fill(fill) {
                    m_sint = m_src->props().isInt();
                    m_dint = m_dst->props().isInt();
            }

            /**
             * Construct a TargetFillOperator. To fill the same raster.
             *
             * @param src The source raster.
             * @param dst The destination raster.
             * @param target The target value.
             * @param fill The fill value.
             * @param band The band to fill.
             */
            TargetFillOperator(Grid* grd, int srcBand, int dstBand, T target, U fill) :
                m_src(grd), m_srcBand(srcBand),
				m_dst(grd), m_dstBand(dstBand),
				m_target(target), m_fill(fill) {
                    m_sint = m_src->props().isInt();
                    m_dint = m_dst->props().isInt();
            }

            /**
             * Return the source grid properties.
             *
             * @return The source grid properties.
             */
            const GridProps& srcProps() const {
                return m_src->props();
            }

            /**
             * Return the destination grid properties.
             *
             * @return The destination grid properties.
             */
            const GridProps& dstProps() const {
                return m_dst->props();
            }

            /**
             * Return true if the given cell should be filled.
             *
             * @param col A column index.
             * @param row A row index.
             * @return True if the given cell should be filled.
             */
            bool shouldFill(int col, int row) const {
                if(m_sint) {
                    return m_src->getInt(col, row, m_srcBand) == m_target;
                } else {
                    return m_src->getFloat(col, row, m_srcBand) == m_target;
                }
            }

            /**
             * Fill the given cell.
             *
             * @param col A column index.
             * @param row A row index.
             */
            void fill(int col, int row) const {
                if(m_dint) {
                    m_dst->setInt(col, row, (int) m_fill, m_dstBand);
                } else {
                    m_dst->setFloat(col, row, (double) m_fill, m_dstBand);
                }
            }

            ~TargetFillOperator() {}
        };



        class Raster;

        /**
         * A convenience class for managing a grid of values.
         * Handles allocation and deallocation of memory.
         */
        class G_DLL_EXPORT MemRaster : public Grid {
        private:
            bool m_mmapped;							///<! True if the memory is mapped to a file.
            void *m_grid;							///<! The allocated memory pointer.
            geo::util::MappedFile* m_mappedFile;	///<! The mapped file container.
            GridProps m_props;						///<! Grid properties.
            std::mutex m_freeMtx;					///<! Mutex to protect memory resources.
            std::mutex m_initMtx;					///<! Mutex to protect memory resources.
            std::mutex m_mtx;						///<! Mutex for use by callers.

            /**
             * Checks if the grid has been initialized. Throws exception otherwise.
             */
            void checkInit() const;

            /**
             * Free reserved memory.
             */
            void freeMem();

        public:

            MemRaster();

            /**
             * Copy another MemRaster.
             *
             * @param other Another MemRaster.
             */
            MemRaster(const MemRaster& other);
            
            /**
             * Create a MemRaster with the given properties.
             *
             * @param props The grid properties.
             * @param mapped True if the raster should be mapped to a file.
             */
            MemRaster(const GridProps &props, bool mapped = false);

            ~MemRaster();

            /**
             * Return a pointer to the allocated memory.
             *
             * @return A pointer to the allocated memory.
             */
            void* grid();

            /**
             * Return true if the memory is mapped to a file.
             *
             * @param True if the memory is mapped to a file.
             */
            bool mmapped() const;

            /**
             * Initialize with the given number of cols and rows.
             * (Re)allocates memory for the internal grid.
             *
             * @param props The grid properties.
             * @param mapped True if the memory should be mapped to a file.
             */
            void init(const GridProps& props, bool mapped = false);

            /**
             * Copies the image data from an entire row into the buffer
             * which must be pre-allocated.
             * @param row The row index.
             * @param band The band number.
             * @param buf A pre-allocated buffer to store the data.
             */
            int getIntRow(int row, int band, int* buf);

            /**
             * Copies the image data from an entire row into the buffer
             * which must be pre-allocated.
             * @param row The row index.
             * @param band The band number.
             * @param buf A pre-allocated buffer to store the data.
             */
            int getFloatRow(int row, int band, double* buf);

            GridProps readIntoVector(std::vector<double>& data);

            void writeFromVector(std::vector<double>& data);


            /**
             * Convert the grid to a matrix.
             *
             * @param mtx The matrix to copy to.
             * @param band The band to copy from.
             */
            void toMatrix(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mtx, int band);

            /**
             * Initialize the grid from a matrix.
             *
             * @param mtx The matrix to copy from.
             * @param band The band to copy to.
             */
            void fromMatrix(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mtx, int band);

            std::mutex& mutex();

            const GridProps& props() const;

            void fillFloat(double value, int band);

            void fillInt(int value, int band);

            //int getInt(uint64_t idx, int band);

            int getInt(int col, int row, int band);

            //double getFloat(uint64_t idx, int band);

            double getFloat(int col, int row, int band);

            //void setInt(uint64_t idx, int value, int band);

            void setInt(int col, int row, int value, int band);

            //void setFloat(uint64_t idx, double value, int band);

            void setFloat(int col, int row, double value, int band);

			void writeTo(Grid& grd,
					int cols = 0, int rows = 0,
					int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
					int srcBand = 1, int dstBand = 1);

            void writeToMemRaster(MemRaster& grd,
            		int cols = 0, int rows = 0,
            		int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
            		int srcBand = 1, int dstBand = 1);

            void writeToRaster(Raster& grd,
            		int cols = 0, int rows = 0,
            		int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
					int srcBand = 1, int dstBand = 1);

            bool writeToVector(std::vector<double>& grd, int col, int row, int cols, int rows, int band = 1, double invalid = nan(""));

            bool writeFromVector(std::vector<double>& grd, int col, int row, int cols, int rows, int band = 1);
        };

        /**
         * Represents a file-backed raster.
         */
        class G_DLL_EXPORT Raster : public Grid {
        	friend class MemRaster;
        private:
            GDALDataset *m_ds;          				///<! GDAL data set pointer.
            int m_bcols, m_brows;						///<! The size of the GDAL block.
            int m_bcol, m_brow;							///<! The current loaded block position.
            int m_bband;								///<! The band corresponding to the current block.
            bool m_dirty;								///<! True if the current block has been written to and must be flushed.
            std::string m_filename;     				///<! Raster filename
            GridProps m_props;							///<! Properties of the raster.
            GDALDataType m_type;        				///<! GDALDataType -- limits the possible template types.
            std::unordered_map<int, void*> m_blocks;	///<! The block cache.
            std::mutex m_mtx;							///<! Mutex for use by callers.

            /**
             * Returns the GDAL data type.
             *
             * @return The GDAL data type.
             */
            GDALDataType getGDType() const;

            /**
             * Get an allocated block of memory from the map of cached
             * blocks corresponding to the given band. If it has not
             * already been allocated, allocated it.
             *
             * @param band The raster band.
             * @return A pointer to the allocated block.
             */
            void* getBlock(int band);

        protected:

            /**
             * Return the GDAL data set pointer.
             *
             * @return The GDAL data set pointer.
             */
            GDALDataset* ds() const;

            /**
             * Write the raster into another grid instance.
             *
             * @param cols The number of columns to write.
             * @param rows The number of rows to write.
             * @param srcCol The column to begin reading from.
             * @param srcRow The row to begin reading from.
             * @param dstCol The column to begin writing to.
             * @param dstRow The row to begin writing to.
             * @param srcBand The band to read from.
             * @param dstBand The band to write to.
             */
			void writeToMemRaster(MemRaster &grd,
            		int cols = 0, int rows = 0,
            		int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
            		int srcBand = 1, int dstBand = 1);

			void flushDirtyBlock();

            /**
             * Write the raster into another grid instance.
             *
             * @param cols The number of columns to write.
             * @param rows The number of rows to write.
             * @param srcCol The column to begin reading from.
             * @param srcRow The row to begin reading from.
             * @param dstCol The column to begin writing to.
             * @param dstRow The row to begin writing to.
             * @param srcBand The band to read from.
             * @param dstBand The band to write to.
             */
			void writeToRaster(Raster &grd,
            		int cols = 0, int rows = 0,
            		int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
					int srcBand = 1, int dstBand = 1);

        public:

            /**
             * Create a new raster with a template. The raster will be created.
             *
             * @param filename The path to the file.
             * @param props A GridProps instance containing a descriptor for the raster.
             */
            Raster(const std::string &filename, const GridProps &props);

            /**
             * Open the given extant raster. Set the writable argument to true
             * to enable writing.
             *
             * @param filename The path to the file.
             * @param writable True if the file is to be writable.
             */
            Raster(const std::string &filename, bool writable = false);

            /**
             * Attempts to return the data type of the raster
             * with the given filename.
             *
             * @param filename The path to an existing raster.
             * @return The data type.
             */
            static DataType getFileDataType(const std::string &filename);

            /**
             * Return a map containing the raster driver short name and extension.
             *
             * @return A map containing the raster driver short name and extension.
             */
            static std::map<std::string, std::set<std::string> > extensions();

            /**
             * Return a map containing the raster driver short name and long name.
             *
             * @return A map containing the raster driver short name and long name.
             */
            static std::map<std::string, std::string> drivers();

            /**
             * Return a map containing the raster driver short name and long name. Use filter
             * to filter the returns on short name.
             *
             * @param filter A vector containing the short names of drivers to include.
             * @return A map containing the raster driver short name and long name.
             */
            static std::map<std::string, std::string> drivers(const std::vector<std::string>& filter);

            /**
             * Get the name of the driver that would be used to open a file
             * with the given path.
             *
             * @param filename The path to an existing raster.
             * @return The name of the griver used to open the file.
             */
            static std::string getDriverForFilename(const std::string &filename);

            /**
             * Creates a virtual raster using the given files and writes it to a file
             * with the given name.
             *
             * @param files A list of files to include in the raster.
             * @param outfile The path to the virtual raster.
             * @param nodata The nodata value for the virtual raster.
             */
            static void createVirtualRaster(const std::vector<std::string>& files, const std::string& outfile, double nodata);

            /**
             * Creates a virtual raster using the given files and writes it to a file
             * with the given name.
             *
             * @param begin An iterator into a list of files to include in the raster.
             * @param files The end iterator of the list of files to include in the raster.
             * @param outfile The path to the virtual raster.
             * @param nodata The nodata value for the virtual raster.
             */
            template <class T>
            static void createVirtualRaster(T begin, T end, const std::string& outfile, double nodata) {
            	std::vector<std::string> files(begin, end);
            	return createVirtualRaster(files, outfile, nodata);
            }

            static GridProps readIntoVector(const std::string& filename, int band, std::vector<double>& data);

            /**
             * Return the filename for this raster.
             * @return The filename for this raster.
             */
            std::string filename() const;

            /**
             * Copies the image data from an entire row into the buffer
             * which must be pre-allocated.
             * @param row The row index.
             * @param band The band number.
             * @param buf A pre-allocated buffer to store the data.
             */
            int getIntRow(int row, int band, int* buf);

            /**
             * Copies the image data from an entire row into the buffer
             * which must be pre-allocated.
             * @param row The row index.
             * @param band The band number.
             * @param buf A pre-allocated buffer to store the data.
             * @return The number of columns read.
             */
            int getFloatRow(int row, int band, double* buf);

            /**
             * Returns true if the raster is square.
             * @return True if the raster is square.
             */
            bool isSquare() const;

            /**
             * Flush the current block to the dataset.
             */
            void flush();

            /**
             * Flush a dirty read/write block to the dataset.
             */
            void flushDirty();

            std::mutex& mutex();

            const GridProps& props() const;

            void fillInt(int value, int band = 1);

            void fillFloat(double value, int band = 1);

			void writeTo(Grid &grd,
					int cols = 0, int rows = 0,
					int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
					int srcBand = 1, int dstBand = 1);

            int getInt(double x, double y, int band);

            int getInt(int col, int row, int band);

            //int getInt(uint64_t idx, int band);

            double getFloat(double x, double y, int band);

            double getFloat(int col, int row, int band);

            //double getFloat(uint64_t idx, int band);

            void setInt(double x, double y, int v, int band);

            void setInt(int col, int row, int v, int band);

            //void setInt(uint64_t idx, int v, int band);

            void setFloat(double x, double y, double v, int band);

            void setFloat(int col, int row, double v, int band);

            //void setFloat(uint64_t idx, double v, int band);

            ~Raster();

        };

    } // raster

} // geo


#endif
