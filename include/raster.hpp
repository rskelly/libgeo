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

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <eigen3/Eigen/Core>

#include "geo.hpp"
#include "util.hpp"

using namespace geo::util;

namespace geo {

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
            bool m_nodataSet;			// True if nodata is set.
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
    		bool nodataSet() const;
    		void unsetNodata();

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
            uint64_t size() const;

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
            uint64_t count;
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
        template <class T, class U>
        class G_DLL_EXPORT FillOperator {
        public:
            virtual const GridProps& srcProps() const = 0;
            virtual const GridProps& dstProps() const = 0;
            virtual bool shouldFill(int col, int row) const = 0;
            virtual void fill(int col, int row) const = 0;
            virtual ~FillOperator() {};
        };

        class Grid;
        class TileIterator;

        class G_DLL_EXPORT Tile {
        	friend class TileIterator;
        private:
        	Grid* m_tile;
        	Grid* m_source;
        	int m_cols;
        	int m_rows;
        	int m_col; 		// The position of the ROI in the source.
        	int m_row;
        	int m_buffer;
        	int m_srcCol;	// The buffered position in the source.
        	int m_srcRow;
        	int m_dstCol;	// The position written to in the tile.
        	int m_dstRow;
        	int m_band;
        	bool m_writeOnFlush;

        protected:

        	// Create the tile.
        	Tile(Grid* tile, Grid* source, int cols, int rows,
        			int col, int row, int buffer,
        			int srcCol, int srcRow,
					int dstCol, int dstRow, int band, bool writeOnFlush);

        public:

        	int srcCol() const;
        	int dstCol() const;
        	int srcRow() const;
        	int dstRow() const;
        	int cols() const;
        	int rows() const;
        	int col() const;
        	int row() const;

        	// Return the grid containing tile data.
        	Grid& grid();

        	// Write the tile's data to the destination grid.
        	void writeTo(Grid& dest);

        	// Flush tile contents to the source. Called implicitly
        	// on destruction.
        	void flush();

        	// Destroy the Tile. Implicitly flushes any changes.
        	~Tile();
        };

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

        	// Create a TileIterator of the given size. If a buffer is given,
        	// The tile size in increased, and pixels are read from the source
        	// to fill the buffer, however only pixels within the unbuffered
        	// window are written back.
        	TileIterator(Grid& source, int cols, int rows, int buffer = 0, int band = 1);

        public:

        	// Returns true if there's another tile to be
        	// retrieved
        	bool hasNext();

        	// Returns the number of tiles.
        	int count() const;

        	// Returns the next tile. Writes the previous tile to
        	// source (if there is one) and reads the next one.
        	Tile next();

        	// Creates a Tile using the given template. This Tile
        	// is independent of the iterator.
        	Tile create(Tile &tpl);

        	// Destroy the TileIterator.
        	~TileIterator();
        };

        // Abstract class for grids (rasters).
        class G_DLL_EXPORT Grid {
        public:
            Grid();

            virtual ~Grid() = 0;
            
            // Return a tile iterator;
            std::unique_ptr<TileIterator> iterator(int cols, int rows, int buffer = 0, int band = 1);

            // Compute the table of Gaussian weights given the size of the table
            // and the std. deviation.
            static void gaussianWeights(double *weights, int size, double sigma);

            GridStats stats(int band);

            // Returns the grid properties
            virtual const GridProps &props() const = 0;

            // Fill the entire dataset with the given value.
            virtual void fillFloat(double value, int band = 1) = 0;
            virtual void fillInt(int value, int band = 1) = 0;

            // Return a the value held at the given index in the grid.
            virtual int getInt(uint64_t idx, int band = 1) = 0;
            virtual int getInt(int col, int row, int band = 1) = 0;
            virtual double getFloat(uint64_t idx, int band = 1) = 0;
            virtual double getFloat(int col, int row, int band = 1) = 0;

            // Set the value held at  the given index in the grid.
            virtual void setInt(uint64_t idx, int value, int band = 1) = 0;
            virtual void setInt(int col, int row, int value, int band = 1) = 0;
            virtual void setFloat(uint64_t idx, double value, int band = 1) = 0;
            virtual void setFloat(int col, int row, double value, int band = 1) = 0;

            // Write data from the current Grid instance to the given grid.
            virtual void writeTo(Grid &grd,
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
            template <class T, class U>
            static void floodFill(int col, int row,
                FillOperator<T, U> &op,  bool d8 = false,
				int *outminc = nullptr, int *outminr = nullptr,
				int *outmaxc = nullptr, int *outmaxr = nullptr,
				int *outarea = nullptr) {

                const GridProps& gp = op.srcProps();

                int cols = gp.cols();
                int rows = gp.rows();
                int size = gp.size();
                int minc = cols + 1;
                int minr = rows + 1;
                int maxc = -1;
                int maxr = -1;
                int area = 0;
                std::queue<std::unique_ptr<Cell> > q;
                q.push(std::unique_ptr<Cell>(new Cell(col, row)));

                std::vector<bool> visited(size, false); // Tracks visited pixels.

                while (q.size()) {

                    std::unique_ptr<Cell> cel = std::move(q.front());
                    row = cel->row;
                    col = cel->col;
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
                            q.push(std::unique_ptr<Cell>(new Cell(col, row - 1)));
                        if (row < rows - 1)
                            q.push(std::unique_ptr<Cell>(new Cell(col, row + 1)));

                        int c;
                        for (c = col - 1; c >= 0; --c) {
                            idx = (uint64_t) row * cols + c;
                            if (!visited[idx] && op.shouldFill(c, row)) {
                                minc = g_min(c, minc);
                                ++area;
                                op.fill(c, row);
                                visited[idx] = true;
                                if (row > 0)
                                    q.push(std::unique_ptr<Cell>(new Cell(c, row - 1)));
                                if (row < rows - 1)
                                    q.push(std::unique_ptr<Cell>(new Cell(c, row + 1)));
                            } else {
                                break;
                            }
                        }
                        if(d8) {
                            if (row > 0)
                                q.push(std::unique_ptr<Cell>(new Cell(c, row - 1)));
                            if (row < rows - 1)
                                q.push(std::unique_ptr<Cell>(new Cell(c, row + 1)));
                        }
                        for (c = col + 1; c < cols; ++c) {
                            idx = (uint64_t) row * cols + c;
                            if (!visited[idx] && op.shouldFill(c, row)) {
                                maxc = g_max(c, maxc);
                                ++area;
                                op.fill(c, row);
                                visited[idx] = true;
                                if (row > 0)
                                    q.push(std::unique_ptr<Cell>(new Cell(c, row - 1)));
                                if (row < rows - 1)
                                    q.push(std::unique_ptr<Cell>(new Cell(c, row + 1)));
                            } else {
                                break;
                            }
                        }
                        if(d8) {
                            if (row > 0)
                                q.push(std::unique_ptr<Cell>(new Cell(c, row - 1)));
                            if (row < rows - 1)
                                q.push(std::unique_ptr<Cell>(new Cell(c, row + 1)));
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

            // Smooth the raster and write the smoothed version to the output raster.
            // Callback is an optional function reference with a single float
            // between 0 and 1, for status tracking.
            void smooth(Grid &smoothed, double sigma, int size, int band = 1,
                Status* status = nullptr,
                bool *cancel = nullptr);

            // The radius is given with cells as the unit, but
            // can be rational. When determining which cells to
            // include in the calculation, any cell which partially
            // falls in the radius will be included.
            void voidFillIDW(double radius, int count = 4, double exp = 2.0, int band = 1);

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
        public:
            TargetFillOperator(Grid* src, Grid* dst, T target, U fill) :
                m_src(src), m_dst(dst), m_target(target), m_fill(fill) {
                    m_sint = m_src->props().isInt();
                    m_dint = m_dst->props().isInt();
            }
            TargetFillOperator(Grid* grd, T target, U fill) :
                m_src(grd), m_dst(grd), m_target(target), m_fill(fill) {
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
                    return m_src->getInt(col, row) == m_target;
                } else {
                    return m_src->getFloat(col, row) == m_target;
                }
            }
            void fill(int col, int row) const {
                if(m_dint) {
                    m_dst->setInt(col, row, (int) m_fill);
                } else {
                    m_dst->setFloat(col, row, (double) m_fill);
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
            std::unique_ptr<boost::interprocess::mapped_region> m_region;
            std::unique_ptr<boost::interprocess::file_mapping> m_mapping;
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
    		void *m_block;
            bool m_dirty;
            int m_bband;
            std::string m_filename;     // Raster filename
            GridProps m_props;
            GDALDataType m_type;        // GDALDataType -- limits the possible template types.
            std::mutex m_mtx;

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

            double getFloat(double x, double y, int band = 1);
            double getFloat(int col, int row, int band = 1);
            double getFloat(uint64_t idx, int band = 1);

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

            // Vectorize the raster.
            void polygonize(const std::string &filename, const std::string &layerName, 
                const std::string &driver, uint16_t srid = 0, uint16_t band = 1,
				bool removeHoles = false, bool removeDangles = false,
				geo::util::Status *status = nullptr, bool *cancel = nullptr);

            ~Raster();

        };

    } // raster

} // geo


#endif
