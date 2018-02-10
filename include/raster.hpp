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

#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>

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
    		double m_trans[6];			// The geotransform properties.
    		int m_cols, m_rows;			// The number of rows and columns.
    		int m_vsrid, m_hsrid;		// The vertical and horizontal SRIDs
            int m_bands;           		// The number of bands.
            bool m_writable;            // True if the raster is writable
            double m_nodata;			// The nodata value.
            bool m_nodataSet;			// True if nodata is set.
    		DataType m_type;			// The data type.
    		std::string m_projection;	// The WKT representation of the projection
    		std::string m_driver;		// The name of the GDAL driver.

		public:

    		/**
    		 * Construct an empty GridProps.
    		 */
    		GridProps();

    		/**
    		 * Return the geographic bounds of the raster.
    		 */
    		Bounds bounds() const;

    		/**
    		 * Return the name of the GDAL driver used by the raster.
    		 * Only relevant for file-based rasters.
    		 */
    		std::string driver() const;

    		/**
    		 * Set the name of the GDAL driver used by the raster.
    		 * Only relevant for file-based rasters.
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
    		 */
    		double nodata() const;

    		/**
    		 * Set the no data value.
    		 * @param nodata The no data value.
    		 */
    		void setNoData(double nodata);

    		/**
    		 * Returns true if the no data value has been set.
    		 */
    		bool nodataSet() const;

    		/**
    		 * Remove the no data value.
    		 */
    		void unsetNodata();

    		/**
    		 * Return the number of columns.
    		 */
            int cols() const;

            /*
             * Return the number of rows.
             */
            int rows() const;

            /**
             * Returns true if the cell is in the raster.
             * @param col The column.
             * @param row The row.
             */
            bool hasCell(int col, int row) const;

            /**
             * Returns true if the cell is in the raster.
             * @param x The geographic x or longitude coordinate.
             * @param y The geographic y or latitude coordinate.
             */
            bool hasCell(double x, double y) const;

            /**
             * Returns the row for a given y-coordinate.
             * @param y The geographic y or latitude coordinate.
             */
            int toRow(double y) const;

            /**
             * Returns the column for a given x-coordinate.
             * @param x The geographic x or longitude coordinate.
             */
            int toCol(double x) const;

            /**
             * Returns the x-coordinate for a given column.
             * @param col The column.
             */
            double toX(int col) const;

            /**
             * Returns the y-coordinate for a given row.
             * @param row The row.
             */
            double toY(int row) const;

            /**
             * Returns the x-coordinate for the cell centroid of a given column.
             * @param col The column.
             */
            double toCentroidX(int col) const;

            /**
             * Returns the y-coordinate for the cell centorid of a given row.
             * @param row The row.
             */
            double toCentroidY(int row) const;

            /**
             * Returns the number of pixels.
             */
            uint64_t size() const;

            /**
             * Set the data type of the raster.
             * @param type The data type.
             */
    		void setDataType(DataType type);

    		/**
    		 * Get the data type of the raster.
    		 */
    		DataType dataType() const;

    		/**
    		 * Set the size of the raster in columns, rows.
             * @param col The column.
             * @param row The row.
    		 */
    		void setSize(int cols, int rows);

    		/**
    		 * Set the horizontal and vertical (optional) SRID.
    		 * @param hsrid The horizontal SRID.
    		 * @param vsrid The vertical SRID.
    		 */
    		void setSrid(int hsrid, int vsrid = 0);

    		/**
    		 * Get the vertical SRID.
    		 */
    		int vsrid() const;

    		/**
    		 * Get the horizontal SRID.
    		 */
    		int hsrid() const;

    		/**
    		 * Set the WKT projection.
    		 * @param The projection string (proj or WKT format).
    		 */
    		void setProjection(const std::string &proj);

    		/**
    		 * Get the WKT projection (proj or WKT format).
    		 */
    		std::string projection() const;

    		/**
    		 * Set the geo transform properties.
    		 * @param trans The six-element transformation matrix.
    		 */
    		void setTrans(double trans[6]);

    		/**
    		 * Set the geo transform properties. The third and fourth
    		 * elements are set to zero.
    		 * @param tlx The top left x-coordinate.
    		 * @param resX The horizontal resolution.
    		 * @param tly The top left y-coordinate.
    		 * @param resY The vertical resolution (negative for UTM (etc.) projections).
    		 */
    		void setTrans(double tlx, double resX, double tly, double resY);

    		/**
    		 * Gets the geo transform properties by setting them
    		 * in the given array.
    		 * @param trans The six-element transformation matrix.
    		 */
    		void trans(double trans[6]) const;

    		/**
    		 * Set the vertical and horizontal resolution.
    		 * @param resolutionX The horizontal resolution.
    		 * @param resolutionY The vertical resolution (negative for UTM (etc.) projections).
    		 */
    		void setResolution(double resolutionX, double resolutionY);

    		/**
    		 * Get the horizontal resolution.
    		 */
    		double resolutionX() const;

    		/**
    		 * Get the vertical resolution.
    		 */
    		double resolutionY() const;

    		/**
    		 * Return the top-left horizontal coordinate of the raster.
    		 */
    		double tlx() const;

    		/**
    		 * Return the top-left vertical coordinate of the raster.
    		 */
    		double tly() const;

    		/**
    		 * Set the number of bands.
    		 * @param bands The number of bands.
    		 */
    		void setBands(int bands);

    		/**
    		 * Get the number of bands.
    		 */
    		int bands() const;

    		/**
    		 * Set the writable state of the raster.
    		 * If the raster is not writable, attempting to write
    		 * to it will throw an exception.
    		 * @param writable True, if the raster should be writable.
    		 */
    		void setWritable(bool writable);

    		/**
    		 * Get the writable state of the raster.
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
            uint64_t count;
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
        	 */
            virtual const GridProps& srcProps() const = 0;

            /**
             * Return the properties of the destination raster.
             */
            virtual const GridProps& dstProps() const = 0;

            /**
             * Return true if if the current pixel should be filled.
             * @param col The column.
             * @param row The row.
             */
            virtual bool shouldFill(int col, int row) const = 0;

            /**
             * Fill the current column.
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
        	Grid* m_tile;			// Contains the tile grid.
        	Grid* m_source;			// Contains the source grid.
        	int m_cols;				// The number of columns in the tile.
        	int m_rows;				// The number of rows in the tile.
        	int m_col; 				// The position of the ROI in the source.
        	int m_row;
        	int m_buffer;			// The number of buffer pixels.
        	int m_srcCol;			// The buffered position in the source.
        	int m_srcRow;
        	int m_dstCol;			// The position written to in the tile.
        	int m_dstRow;
        	int m_band;				// The raster band in the source raster.
        	bool m_writeOnFlush;	// If true, writes back to the raster on destruction.

        protected:

        	/**
        	 * Create the tile.
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
        	 */
        	int srcCol() const;

        	/**
        	 * Return the destination column.
        	 */
        	int dstCol() const;

        	/**
        	 * Return the source row.
        	 */
        	int srcRow() const;

        	/**
        	 * Return the destination row.
        	 */
        	int dstRow() const;

        	/**
        	 * Return the number of columns.
        	 */
        	int cols() const;

        	/**
        	 * Return the number of rows.
        	 */
        	int rows() const;

        	/**
        	 * Return the position of the tile in the source raster.
        	 */
        	int col() const;

        	/**
        	 * Return the position of the tile in the source raster.
        	 */
        	int row() const;

        	/**
        	 * Return the grid containing tile data.
        	 */
        	Grid& grid();

        	/**
        	 * Write the tile's data to the destination grid.
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
        	Grid& m_source;
        	int m_cols;
        	int m_rows;
        	int m_buffer;
        	int m_curCol;
        	int m_curRow;
        	int m_band;
        	mutable std::mutex m_mtx;

        protected:

        	/**
        	 * Create a TileIterator of the given size. If a buffer is given,
        	 * The tile size in increased, and pixels are read from the source
        	 * to fill the buffer, however only pixels within the unbuffered
        	 * window are written back.
        	 * @param source The source raster.
        	 * @param cols The number of columns in each tile.
        	 * @param rows The number of rows in each tile.
        	 * @param buffer The size of the buffer around the tile.
        	 * @param band The band in the source raster.
        	 */
        	TileIterator(Grid& source, int cols, int rows, int buffer, int band);

        	/**
        	 * Copy constructor.
        	 */
        	TileIterator(const TileIterator& iter);

        	/**
        	 * Construct a TileIterator.
        	 */
        	TileIterator();

        public:

        	/**
        	 * Returns true if there's another tile to be retrieved.
        	 */
        	bool hasNext();

        	/**
        	 * Returns the number of tiles.
        	 */
        	int count() const;

        	/**
        	 * Returns the next tile. Writes the previous tile to source
        	 * (if there is one) and reads the next one.
        	 */
        	Tile next();

        	/**
        	 * Create a Tile using the given tile as a Template. Does not
        	 * interfere with the iterator.
        	 * @param tpl A tile template.
        	 */
        	Tile create(Tile &tpl);

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
        	 * Construct the grid.
        	 */
            Grid();

            /**
             * Destroy the grid.
             */
            virtual ~Grid() = 0;
            
            /**
             * Return a tile iterator.
             * @param cols The number of columns in each tile.
             * @param rows The number of rows in each tile.
             * @param buffer The buffer around the tile.
             * @param band The band in the source raster.
             */
            TileIterator iterator(int cols, int rows, int buffer = 0, int band = 1);

            /**
             * Compute the table of Gaussian weights given the size of
             * the table and the standard deviation.
             * @param weights The list of weights.
             * @param size The size of the weights list.
             * @param sigma The standard deviation.
             */
            static void gaussianWeights(double *weights, int size, double sigma);

            /**
             * Compute and return the statistics for the band.
             * @param band The raster band.
             */
            GridStats stats(int band);

            /**
             * Returns the grid properties.
             */
            virtual const GridProps &props() const = 0;

            /**
             * Fill the entire dataset with the given value.
             * @param value The value to fill the raster with.
             * @param band The band to fill.
             */
            virtual void fillFloat(double value, int band = 1) = 0;

            /**
             * Fill the entire dataset with the given value.
             * @param value The value to fill the raster with.
             * @param band The band to fill.
             */
            virtual void fillInt(int value, int band = 1) = 0;

            /**
             * Return a the value held at the given index in the grid.
             * @param idx The index in the raster; left to right, top to bottom.
             * @param band The band.
             */
            virtual int getInt(uint64_t idx, int band = 1) = 0;

            /**
             * Return a the value held at the given position in the grid.
             * @param col The column.
             * @param row The row.
             * @param band The band.
             */
            virtual int getInt(int col, int row, int band = 1) = 0;

            /**
             * Return a the value held at the given index in the grid.
             * @param idx The index in the raster; left to right, top to bottom.
             * @param band The band.
             */
            virtual double getFloat(uint64_t idx, int band = 1) = 0;

            /**
             * Return a the value held at the given position in the grid.
             * @param col The column.
             * @param row The row.
             * @param band The band.
             */
            virtual double getFloat(int col, int row, int band = 1) = 0;

            /**
             * Set the value held at  the given index in the grid.
             * @param idx The index in the raster; left to right, top to bottom.
             * @param value The value to set.
             * @param band The band.
             */
            virtual void setInt(uint64_t idx, int value, int band = 1) = 0;

            /**
             * Set the value held at  the given index in the grid.
             * @param col The column.
             * @param row The row.
             * @param value The value to set.
             * @param band The band.
             */
            virtual void setInt(int col, int row, int value, int band = 1) = 0;

            /**
             * Set the value held at  the given index in the grid.
             * @param idx The index in the raster; left to right, top to bottom.
             * @param value The value to set.
             * @param band The band.
             */
            virtual void setFloat(uint64_t idx, double value, int band = 1) = 0;

            /**
             * Set the value held at  the given index in the grid.
             * @param col The column.
             * @param row The row.
             * @param value The value to set.
             * @param band The band.
             */
            virtual void setFloat(int col, int row, double value, int band = 1) = 0;

            /**
             * Write data from the current Grid instance to the given grid.
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
             * @param band The target band.
             */
            void normalize(int band = 1);

            /**
             * Normalize the grid so that the max value is equal to 1, and the
             * minimum is zero.
             * @param band The target band.
             */
            void logNormalize(int band = 1);

            /**
             * Convert a Grid to some other type.
             */
            void convert(Grid &g, int srcBand = 1, int dstBand = 1);

            /**
             * Fill the grid, beginning with the target cell, where any contiguous cell
             * satisfies the given FillOperator. The other grid is actually filled,
             * and the present grid is unchanged *unless* the present grid is passed
             * as other.
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
                FillOperator<T, U> &op,  bool d8 = false,
				int *outminc = nullptr, int *outminr = nullptr,
				int *outmaxc = nullptr, int *outmaxr = nullptr,
				int *outarea = nullptr) {

                const GridProps& gp = op.srcProps();

                int cols = gp.cols();
                int rows = gp.rows();
                uint64_t size = gp.size();
                int minc = cols + 1;
                int minr = rows + 1;
                int maxc = -1;
                int maxr = -1;
                int area = 0;

                std::queue<Cell> q;
                q.push(Cell(col, row));

                std::vector<bool> visited(size, false); // Tracks visited pixels.

                while (q.size()) {

                    const Cell& cel = q.front();
                    row = cel.row;
                    col = cel.col;
                    q.pop();

                    uint64_t idx = (uint64_t) row * cols + col;

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
                            idx = (uint64_t) row * cols + c;
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
                            idx = (uint64_t) row * cols + c;
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
             * @param smoothed The smoothed grid.
             * @param sigma    The standard deviation.
             * @param size     The window size.
             * @param band     The target band.
             * @param status   A pointer to a Status object which will be updated with
             *                 the current status of the process.
             * @param cancel   A pointer to a variable that will be true if the
             *                 process should be cancelled.
			 */
            void smooth(Grid &smoothed, double sigma = 0.84089642, int size = 3, int band = 1,
                Status* status = nullptr,
                bool *cancel = nullptr);

            /**
             * The radius is given with cells as the unit, but
             * can be rational. When determining which cells to
             * include in the calculation, any cell which partially
             * falls in the radius will be included.
             * @param filename 	The output filename.
             * @param radius 	The search radius.
             * @param count  	The number of pixels to use for calculations.
             * @param exp    	The exponent.
             * @param band   	The source band.
             * @param mask		A raster to use as a mask; invalid pixels will be ignored.
             * @param maskBand  The band to use in the mask.
			 */
            void voidFillIDW(const std::string& filename, double radius, int count = 4, double exp = 2.0, int band = 1,
            		const std::string& mask = "", int maskBand = 1);


            // TODO: Document me.
        	template <class V>
        	void writeAStarPath(uint64_t start, std::unordered_map<uint64_t, uint64_t>& parents, V inserter) {
        		*inserter = start;
        		++inserter;
        		while(parents.find(start) != parents.end()) {
        			start = parents[start];
        			*inserter = start;
        			++inserter;
        		}
            }

        	uint64_t minValue(std::unordered_map<uint64_t, double>& m) {
        		double min = G_DBL_MAX_POS;
        		uint64_t key = 0;
        		for(const auto& it : m) {
        			if(it.second < min) {
        				min = it.second;
        				key = it.first;
        			}
        		}
        		return key;
        	}

            // Finds the least-cost path from the start cell to the goal cell,
            // using the given heuristic. Returns the optimal path between the
            // start cell and the goal.
            template <class U, class V>
            void searchAStar(int startCol, int startRow, int goalCol, int goalRow, U heuristic, V inserter) {

            	uint64_t goal = ((uint64_t) goalCol << 32) | goalRow;

            	std::vector<std::pair<int, int> > offsets = Util::squareKernel(3, false);

            	std::unordered_map<uint64_t, uint64_t> parents;
            	std::unordered_map<uint64_t, double> gscore;
            	std::unordered_map<uint64_t, double> fscore;

            	std::unordered_set<uint64_t> openSet;
            	std::unordered_set<uint64_t> closedSet;

            	uint64_t start = ((uint64_t) startCol << 32) | startRow;

            	openSet.insert(start);
            	gscore[start] = 0; 						// Distance from start to neighbour
            	fscore[start] = heuristic(start, goal); // Distance from neighbour to goal.

        		int cols = props().cols();
        		int rows = props().rows();

            	while(!openSet.empty()) {

            		uint64_t top = minValue(fscore);

            		if(top == goal) {
            			writeAStarPath(top, parents, inserter);
            			break;
            		}

            		double gscore0 = gscore[top];

            		fscore.erase(top);
            		gscore.erase(top);

            		openSet.erase(top);
            		closedSet.insert(top);

            		int qcol = (top >> 32) & 0xffffffff;
            		int qrow = top & 0xffffffff;

            		for(const auto& it : offsets) {

            			if(qcol + it.first < 0 || qrow + it.second < 0 || qcol + it.first >= cols || qrow + it.second >= rows)
            				continue;

            			uint64_t n = ((uint64_t) (qcol + it.first) << 32) | (qrow + it.second);

            			if(closedSet.find(n) != closedSet.end())
            				continue;

            			double tgscore = gscore0 + heuristic(top, n);

            			if(openSet.find(n) == openSet.end()) {
            				openSet.insert(n);
            			} else if(tgscore >= gscore[n]) {
            				continue;
            			}

            			parents[n] = top;
            			gscore[n] = tgscore;
            			double h = heuristic(n, goal);
            			fscore[n] = tgscore + h;
            		}
            	}
            }

        };

       
        template <class T, class U>
        class G_DLL_EXPORT TargetFillOperator : public FillOperator<T, U> {
        private:
            Grid* m_src;
            Grid* m_dst;
            bool m_sint, m_dint;
            T m_target;
            U m_fill;
            int m_band;
        public:
            TargetFillOperator(Grid* src, Grid* dst, T target, U fill, int band = 1) :
                m_src(src), m_dst(dst), m_target(target), m_fill(fill), m_band(band) {
                    m_sint = m_src->props().isInt();
                    m_dint = m_dst->props().isInt();
            }
            TargetFillOperator(Grid* grd, T target, U fill, int band = 1) :
                m_src(grd), m_dst(grd), m_target(target), m_fill(fill), m_band(band) {
                    m_sint = m_src->props().isInt();
                    m_dint = m_dst->props().isInt();
            }
            const GridProps& srcProps() const {
                return m_src->props();
            }
            const GridProps& dstProps() const {
                return m_dst->props();
            }
            bool shouldFill(int col, int row) const {
                if(m_sint) {
                    return m_src->getInt(col, row, m_band) == m_target;
                } else {
                    return m_src->getFloat(col, row, m_band) == m_target;
                }
            }
            void fill(int col, int row) const {
                if(m_dint) {
                    m_dst->setInt(col, row, (int) m_fill, m_band);
                } else {
                    m_dst->setFloat(col, row, (double) m_fill, m_band);
                }
            }
            ~TargetFillOperator() {}
        };



        class Raster;

        // A convenience class for managing a grid of values.
        // Handles allocation and deallocation of memory.
        class G_DLL_EXPORT MemRaster : public Grid {
        private:
            void *m_grid;
            bool m_mmapped;
            GridProps m_props;
            std::unique_ptr<geo::util::MappedFile> m_mappedFile;
            std::mutex m_mtx;

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

            bool mmapped() const;

            // Initialize with the given number of cols and rows.
            // (Re)allocates memory for the internal grid.
            void init(const GridProps &props, bool mapped = false);

            // Fill the entire dataset with the given value.
            void fillFloat(double value, int band = 1);
            void fillInt(int value, int band = 1);

            // Return a the value held at the given index in the grid.
            int getInt(uint64_t idx, int band = 1);
            int getInt(int col, int row, int band = 1);
            double getFloat(uint64_t idx, int band = 1);
            double getFloat(int col, int row, int band = 1);

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
            int getFloatRow(int row, int band, float* buf);

            // Set the value held at  the given index in the grid.
            void setInt(uint64_t idx, int value, int band = 1);
            void setInt(int col, int row, int value, int band = 1);
            void setFloat(uint64_t idx, double value, int band = 1);
            void setFloat(int col, int row, double value, int band = 1);

			void writeTo(Grid &grd,
					int cols = 0, int rows = 0,
					int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
					int srcBand = 1, int dstBand = 1);
            void writeToMemRaster(MemRaster &grd,
            		int cols = 0, int rows = 0,
            		int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
            		int srcBand = 1, int dstBand = 1);
            void writeToRaster(Raster &grd,
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
    		//void *m_block;
            bool m_dirty;
            int m_bband;
            std::string m_filename;     // Raster filename
            GridProps m_props;
            GDALDataType m_type;        // GDALDataType -- limits the possible template types.
            std::mutex m_mtx;
            std::unordered_map<int, void*> m_blocks;

            GDALDataType getGDType() const;

            /**
             * Get an allocated block of memory from the map of cached
             * blocks corresponding to the given band. If it has not
             * already been allocated, allocated it.
             * @param band The raster band.
             * @return A pointer to the allocated block.
             */
            void* getBlock(int band);

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

            static void createVirtualRaster(const std::vector<std::string>& files, const std::string& outfile, double nodata);

            template <class T>
            static void createVirtualRaster(T begin, T end, const std::string& outfile, double nodata) {
            	std::vector<std::string> files(begin, end);
            	return createVirtualRaster(files, outfile, nodata);
            }

            // Return the filename for this raster.
            std::string filename() const;

            // Fill the given band with the given value.
            void fillInt(int value, int band = 1);
            void fillFloat(double value, int band = 1);

			void writeTo(Grid &grd,
					int cols = 0, int rows = 0,
					int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
					int srcBand = 1, int dstBand = 1);
            void writeToMemRaster(MemRaster &grd,
            		int cols = 0, int rows = 0,
            		int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
            		int srcBand = 1, int dstBand = 1);
            void writeToRaster(Raster &grd,
            		int cols = 0, int rows = 0,
            		int srcCol = 0, int srcRow = 0,
					int dstCol = 0, int dstRow = 0,
					int srcBand = 1, int dstBand = 1);

            // Returns a pixel value.
            int getInt(double x, double y, int band = 1);
            int getInt(int col, int row, int band = 1);
            int getInt(uint64_t idx, int band = 1);

            /**
             * Copies the image data from an entire row into the buffer
             * which must be pre-allocated.
             * @param row The row index.
             * @param band The band number.
             * @param buf A pre-allocated buffer to store the data.
             */
            int getIntRow(int row, int band, int* buf);

            double getFloat(double x, double y, int band = 1);
            double getFloat(int col, int row, int band = 1);
            double getFloat(uint64_t idx, int band = 1);

            /**
             * Copies the image data from an entire row into the buffer
             * which must be pre-allocated.
             * @param row The row index.
             * @param band The band number.
             * @param buf A pre-allocated buffer to store the data.
             */
            int getFloatRow(int row, int band, double* buf);

            // Set an pixel value.
            void setInt(double x, double y, int v, int band = 1);
            void setInt(int col, int row, int v, int band = 1);
            void setInt(uint64_t idx, int v, int band = 1);

            void setFloat(double x, double y, double v, int band = 1);
            void setFloat(int col, int row, double v, int band = 1);
            void setFloat(uint64_t idx, double v, int band = 1);

            // Returns true if the raster is square.
            bool isSquare() const;

            // Flush the current block to the dataset.
            void flush();

            // Flush a dirty read/write block to the dataset.
            void flushDirty();

            void potrace(const std::string& filename, const std::string& layerName,
                    const std::string& driver, uint16_t srid, uint16_t band = 1, uint16_t threads = 1,
                    bool removeHoles = false, bool removeDangles = false,
                    geo::util::Status *status = nullptr, bool *cancel = nullptr);

            void polygonize2(const std::string& filename, const std::string& layerName,
                    const std::string& driver, uint16_t srid, uint16_t band = 1, uint16_t threads = 1,
                    bool removeHoles = false, bool removeDangles = false,
                    geo::util::Status *status = nullptr, bool *cancel = nullptr);

            // Vectorize the raster.
            void polygonize(const std::string &filename, const std::string &layerName, 
                const std::string &driver, uint16_t srid = 0, uint16_t band = 1, uint16_t threads = 1,
				int bufSize = 0, bool removeHoles = false, bool removeDangles = false,
				geo::util::Status *status = nullptr, bool *cancel = nullptr);

            ~Raster();

        };

    } // raster

} // geo


#endif
