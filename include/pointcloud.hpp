#ifndef __PC_HPP__
#define __PC_HPP__

/**
 * Classes for working with LiDAR point clouds, specifically LAS files.

 *  Created on: Apr 13, 2017
 *  Author: rob
 */

#include <grid.hpp>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <memory>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdint>

#include <liblas/liblas.hpp>

#include "util.hpp"
#include "ds/mkdtree.hpp"

#define NODATA -9999.0
#define D_MAX std::numeric_limits<double>::max()
#define D_MIN std::numeric_limits<double>::lowest()

using namespace geo::util;
using namespace geo::ds;

namespace geo {
namespace pc {

class Point;

/**
 * A class representing a source point cloud file. Maintains
 * the position of the tile, by the minimum x and y coordinates, 
 * the tiles bounding box and a list of filenames of the tile's
 * constituent files.
 */
class PCFile {
private:
	double m_x;								///< Minimum corner x-coordinate.
	double m_y;								///< Minimum corner y-coordinate.
	double m_fileBounds[6];					///< The bounding box of the actual point cloud (may differ from tile bounds.)
	double m_bounds[4];						///< The nominal bounding box.
	double m_bufferedBounds[4];				///< The buffered bounding box.
	size_t m_pointCount;					///< The number of points in the file set.
	bool m_inited;                          ///< True if the PCFile has already been initialized.

	size_t m_index;
	std::ifstream* m_instr;
	liblas::Reader* m_reader;

	std::vector<std::string> m_filenames;	///< The list of filenames of constituent files.

	bool openReader();

	void closeReader();

	bool isReaderOpen() const;

public:

	/**
	 * Construct a PCFile object using the given filename and corner coordinate.
	 * \param filename 	The filename of a consituent file.
	 * \param x 		The minimum x-coordinate.
	 * \param y 		The minimum y-coordinate.
	 */
	PCFile(const std::string& filename, double x = 0, double y = 0, double size = 0, double buffer = 0);

	/**
	 * Construct a PCFile object using the given filenames and corner coordinate.
	 * \param filenames	A list of filenames of consituent files.
	 * \param x 		The minimum x-coordinate.
	 * \param y 		The minimum y-coordinate.
	 */
	PCFile(const std::vector<std::string>& filenames, double x = 0, double y = 0, double size = 0, double buffer = 0);

	void resize(double x, double y, double size, double buffer);

	/**
	 * Get the minimum x-coordinate.
	 * \return The minimum x-coordinate.
	 */
	double x() const;

	/**
	 * Get the minimum y-coordinate.
	 * \return The minimum y-coordinate.
	 */
	double y() const;

	/**
	 * Populates the given array with the bounds of this instance.
	 * There are six elements.
	 * \param bounds A six-element array of doubles, which is overwritten.
	 */
	void fileBounds(double* bounds) const;

	/**
	 * Populates the given array with the nominal bounds of this instance.
	 * There are four elements.
	 * \param bounds A four-element array of doubles, which is overwritten.
	 */
	void bounds(double* bounds) const;

	/**
	 * Populates the given array with the buffered nominal bounds of this instance.
	 * There are four elements.
	 * \param bounds A four-element array of doubles, which is overwritten.
	 */
	void bufferedBounds(double* bounds) const;

	/**
	 * Returns a reference to the filenames list.
	 */
	const std::vector<std::string>& filenames() const;

	/**
	 * Returns true if the point is within the unbuffered bounds of the tile.
	 */
	bool contains(double x, double y) const;

	/**
	 * Returns true if the file bounds intersect the given boundary.
	 * \param bounds An array of coordinates: (minx, miny, maxx, maxy).
	 * \return True if the bounds intersect.
	 */
	bool intersects(double* bounds) const;

	/**
	 * Returns true if the point is within the buffered bounds of the tile.
	 */
	bool containsBuffered(double x, double y) const;

	/**
	 * Returns the point count.
	 * \return The point count.
	 */
	size_t pointCount() const;

	/**
	 * Initialize the last file. Currently just computes its bounds. If
	 * the {useHeader} parameter is false, computes the actual bounds
	 * from the point cloud, otherwise uses the header bounds.
	 * \param useHeader Use the header to determine bounds (fast), rather than the point cloud (slow).
	 */
	void init(bool useHeader = true);

	/**
	 * Populates the point with the next point in the source file and returns
	 * true. If there is no next point, returns false and the point's state
	 * is undefined.
	 * \param pt A point to populate with data.
	 * \return True if there was a point to be read.
	 */
	bool next(geo::pc::Point& pt);

	/**
	 * Destroy the PCFile.
	 */
	~PCFile();

};

/**
 * A class that is responsible for writing new tiles.
 */
class PCWriter {
private:
	int m_fileIdx;							///< The current file index.
	int m_returns;							///< The number of returns (points) in the current file.
	int m_retNum[5];						///< The number of points for each return in the current file.
	long m_totalReturns;					///< The total count of returns across all files.
	double m_bounds[4];						///< The bounding box of the point cloud (may differ from the tile's bounds).
	double m_bufferedBounds[4];				///< The bounding box of the point cloud (may differ from the tile's bounds).
	double m_outBounds[6];					///< The bounding box of the point cloud (may differ from the tile's bounds).
	double m_x;								///< The minimum x-coordinate of the tile.
	double m_y;								///< The minimum y-coordinate of the tile.
	std::string m_filename;					///< The output filename template.
	std::vector<std::string> m_filenames;	///< The list of output filenames.
	liblas::Writer* m_writer;				///< The current {liblas::Writer} instance.
	liblas::Header* m_header;				///< The liblas {liblas::Header} instance.
	std::ofstream m_str;					///< The current output stream.
	bool m_dod;								///< Set to true if constituent files should be deleted on destruction.
	double m_buffer;						///< The size of the buffer on all sides.
	double m_size;

	/**
	 * Creates and returns the next filename in the series for this tile.
	 * \return The new filename.
	 */
	std::string nextFile();

public:
	/**
	 * Construct an instance of a PCWriter using the given output filename, {liblas::Header} and corner
	 * coordinate.
	 * \param filename 	The output filename template. This is the file path without an index part or extension.
	 * \param hdr 		A {liblas::Header} to be used as a template for the new file.
	 * \param x 		The minimum x-coordinate of the tile.
	 * \param y 		The minimum y-coordinate of the tile.
	 */
	PCWriter(const std::string& filename, const liblas::Header& hdr, double x, double y, double size, double buffer);

	PCWriter(PCWriter&& other);

	/**
	 * Get the minimum x-coordinate.
	 * \return The minimum x-coordinate.
	 */
	double x() const;

	/**
	 * Get the minimum y-coordinate.
	 * \return The minimum y-coordinate.
	 */
	double y() const;

	/**
	 * Get the tile size for this writer.
	 * \return The tile size.
	 */
	double size() const;

	/**
	 * The number of points written.
	 * \return The number of points written.
	 */
	size_t count() const;

	/**
	 * Populates the given array with the bounds of this instance.
	 * There are six elements.
	 * \param bounds A six-element array of doubles, which is overwritten.
	 */
	void outBounds(double* bounds) const;

	/**
	 * Populates the given array with the buffered bounds of this instance.
	 * There are six elements.
	 * \param bounds A four-element array of doubles, which is overwritten.
	 */
	void bufferedBounds(double* bounds) const;

	/**
	 * Populates the given array with the nominal bounds of this instance.
	 * There are six elements.
	 * \param bounds A four-element array of doubles, which is overwritten.
	 */
	void bounds(double* bounds) const;

	/**
	 * Returns a reference to the filenames list.
	 */
	const std::vector<std::string>& filenames() const;

	/**
	 * Open a new output stream with a new file. If an output
	 * stream is currently open, it is closed.
	 */
	void open();

	/**
	 * If set to true, the files associated with this tile will be
	 * deleted on destruction of the instnce.
	 */
	void deleteOnDestruct(bool dod);

	/**
	 * Close the current output stream and destroy the writer. Prepare
	 * the instance for the opening of a new file.
	 */
	void close();

	/**
	 * Add a point to the tile represented by this instance.
	 * The point will be added to the currently-open file.
	 * \param pt A point.
	 */
	void addPoint(const liblas::Point& pt);

	/**
	 * Returns true if the given coordinate is within the bounds of this tile's
	 * geographic extent. "Within" means greater than or equal to the minimum
	 * bound and less than the maximum.
	 * \param x The x-coordinate.
	 * \param y The y-coordinate.
	 */
	bool contains(double x, double y) const;

	/**
	 * Returns true if the given coordinate is within the bounds of this tile's
	 * buffered extent. "Within" means greater than or equal to the minimum
	 * bound and less than the maximum.
	 * \param x The x-coordinate.
	 * \param y The y-coordinate.
	 */
	bool containsBuffered(double x, double y) const;

	/**
	 * Destroy the writer. Will close any open files, and, if {dod} is set to
	 * true, will delete constituent files.
	 */
	~PCWriter();

};

/**
 * Represents a single tile, its geographic extent and its buffered extent.
 * Contains a PCWriter instance for writing to the file(s) associated
 * with the tile.
 */
class Tile {
private:
	double m_bounds[4];						///< The geographic bounds of the tile.
	double m_bufferedBounds[4];				///< The buffered extent of the tile.
	std::unique_ptr<PCWriter> m_writer;		///< The PCWriter used for writing points.

public:
	/**
	 * Construct a tile with the given bounds and buffer.
	 * \param minx 		The minimum x-coordinate.
	 * \param miny 		The minimum y-coordinate.
	 * \param maxx 	 	The maximum x-coordinate.
	 * \param maxy 		The maximum y-coordinate.
	 * \param buffer 	The buffer distance.
	 */
	Tile(double minx, double miny, double maxx, double maxy, double buffer = 0);

	/**
	 * Return a pointer to the PCWriter owned by this tile.
	 * \param release If true, the tile will relinquish ownership of the
	 * PCWriter instance.
	 * \return A pointer to the \ling PCWriter owned by this tile.
	 */
	PCWriter* writer(bool release = false);

	/**
	 * Set a PCWriter instance on this tile. This tile is
	 * the unique owner of the instance. Failure will
	 * result if the caller destroys the PCWriter.
	 * \param wtr A PCWriter pointer.
	 */
	void writer(PCWriter* wtr);

	/**
	 * Returns true if the given coordinate is within the bounds of this tile's
	 * geographic extent. "Within" means greater than or equal to the minimum 
	 * bound and less than the maximum.
	 * \param x The x-coordinate.
	 * \param y The y-coordinate.
	 */
	bool contains(double x, double y);

	/**
	 * Returns true if the given coordinate is within the bounds of this tile's
	 * buffered extent. "Within" means greater than or equal to the minimum 
	 * bound and less than the maximum.
	 * \param x The x-coordinate.
	 * \param y The y-coordinate.
	 */
	bool containsBuffered(double x, double y);

};

/**
 * Performs the tiling of a set of point cloud files.
 */
class Tiler {
public:
	std::vector<PCFile> files;	///< The list of PCFile instances.

	/**
	 * Construct a tiler with the given list of initial filenames.
	 * \param filenames A list of filenames of point cloud files to process.
	 */
	Tiler(const std::vector<std::string> filenames);

	/**
	 * Tile the list of point cloud files; write the output to the given output directory.
	 * \param outdir			The output directory.
	 * \param size 				The length of one side of the tile.
	 * \param buffer 			The size of the buffer around each tile. Default 0.
	 * \param srid 				The spatial reference ID. Default 0.
	 * \param easting   		The minimum x-coordinate of the tile grid. If NaN is given, will round
	 *							the minimum cordinate of the data extent down to a multiple of the tile size.
	 *							Default NaN.
	 * \param northing  		The minimum y-coordinate of the tile grid. If nan is given, will round
	 *							the minimum cordinate of the data extent down to a multiple of the tile size.
	 *							Default NaN.
	 * \param maxFileHandles 	The maximum number of output files that can be open at once. This determines
	 *                          the number of intermediate tiles that will be produced. Default 64.
	 */
	void tile(const std::string& outdir, double size, double buffer = 0, int srid = 0, 
		double easting = std::nan(""), double northing = std::nan(""), int maxFileHandles = 64);

	/**
	 * Destroy the Tiler.
	 */
	~Tiler();

};

/**
 * Represents a 3D point that can be used in a tree structure.
 * can be instantiated with raw coordinates or with a
 * liblas::Point instance.
 */
class Point {
private:
	double m_x; 			///< The x-coordinate.
	double m_y; 			///< The y-coordinate.
	double m_z; 			///< The z-coordinate.
	short m_intensity;
	char m_angle;
	char m_cls;
	char m_edge;
	char m_numReturns;
	char m_returnNum;

public:

	/**
	 * Construct a Point using a liblas::Point. Will read the 3D
	 * coordinates and save a pointer to a copy of the object.
	 * The copy is deleted on destruction.
	 * \param pt A liblas::Point object.
	 */
	Point(const liblas::Point& pt);

	/**
	 * Construct a point using the three raw coordinates.
	 * \param x The x-coordinate.
	 * \param y The y-coordinate.
	 * \param z The z-coordinate.
	 */
	Point(double x, double y, double z);

	/**
	 * Construct a point using the three raw coordinates.
	 * \param x The x-coordinate.
	 * \param y The y-coordinate.
	 * \param z The z-coordinate.
	 */
	Point(double x, double y, double z, double intensity, double angle,
			int cls, int returnNum, int numReturns, bool isEdge);

	Point(const geo::pc::Point& pt);

	/**
	 * Construct an empty point.
	 */
	Point();

	/**
	 * Populate this point with values from a liblas point.
	 */
	void setPoint(const liblas::Point& pt);

	/**
	 * Returns true if the point is a last return.
	 * \return True if the point is a last return.
	 */
	bool isLast() const;

	/**
	 * Returns true if the point is a first return.
	 * \return True if the point is a first return.
	 */
	bool isFirst() const;


	/**
	 * Returns the class ID of this point. If no liblas::Point
	 * or other was provided, this returns zero.
	 */
	int classId() const;

	/**
	 * Returns the 2D coordinate for the given index.
	 * The modulus of the index is taken; if index % 2 == 0,
	 * then x is returned, otherwise y.
	 * \param idx The index.
	 */
	double operator[](int idx) const;

	/**
	 * Return the x-coordinate.
	 * \return The x-coordinate.
	 */
	double x() const;

	/**
	 * Return the y-coordinate.
	 * \return The y-coordinate.
	 */
	double y() const;

	/**
	 * Return the z-coordinate.
	 * \return The z-coordinate.
	 */
	double z() const;

	/**
	 * Return the intensity.
	 * \return The intensity.
	 */
	double intensity() const;

	/**
	 * Return the scan angle.
	 * \return The scane angle.
	 */
	double scanAngle() const;

	/**
	 * Return the return number.
	 * \return The return number.
	 */
	int returnNum() const;

	int numReturns() const;

	/**
	 * Return true if the point is marked as a flight line edge.
	 * \return True if the point is marked as a flight line edge.
	 */
	bool isEdge() const;

	/**
	 * This is the getter for values used in competition.
	 * The source of the value may be changed.
	 */
	double value() const;

	/**
	 * Destroy the point.
	 */
	~Point();

};

class PCPointFilter;

class Rasterizer;

/**
 * Used by Rasterizer to compute cell values from the points found
 * in the neighbourhood of a cell.
 */
class Computer {
private:
	geo::pc::Rasterizer* m_rasterizer;		///<! Pointer to instance of Rasterizer.

public:

	/**
	 * Set a pointer to the Rasterizer, so the Computer can glean information
	 * about its processes.
	 *
	 * \param rasterizer A pointer to the Rasterizer.
	 */
	void setRasterizer(geo::pc::Rasterizer* rasterizer);

	/**
	 * Get a pointer to the Rasterizer.
	 *
	 * \return A pointer to the Rasterizer.
	 */
	geo::pc::Rasterizer* rasterizer() const;

	/**
	 * Compute and return the statistic for the points within the neighbourhood
	 * defined by the lists of points.
	 * \param type 	 The statistic to compute.
	 * \param pts 	 The list of Points.
	 * \param radius The neighbourhood to calculate within.
	 * \param out    The output list.
	 * \return The the number of bands computed.
	 */
	virtual int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr) = 0;

	virtual int compute(double x, double y, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) = 0;

	/**
	 * Return the number of bands that should be returned by the computer.
	 * \return The number of bands that should be returned by the computer.
	 */
	virtual int bandCount() const = 0;

	/**
	 * Destroy the computer.
	 */
	virtual ~Computer() {}
};

class PointFilter {
public:
	virtual bool keep(const geo::pc::Point& pt) const = 0;
	virtual void print() const = 0;
	virtual ~PointFilter() {}
};

/**
 * A configurable filter for point clouds.
 */
class PCPointFilter {
public:
	std::vector<PointFilter*> filters;

	/**
	 * Construct a PCPointFilter.
	 */
	PCPointFilter();

	~PCPointFilter();

	/**
	 * Print the command-line parameters that can be used to
	 * configure a PCPointFilter.
	 * \param str An output stream.
	 */
	static void printHelp(std::ostream& str);

	/**
	 * Configure this object by reading the arguments from the command line.
	 * Use PCPointFilter::printHelp to print the expected arguments.
	 * If the argument at the current index is relevant to this item,
	 * return true to notify the caller. Updates the index if
	 * more than one item is read. The caller is responsible for
	 * incrementing the index before the next item.
	 * \param idx  The index of the current argument.
	 * \param argv The array of arguments.
	 * \return True, if the current argument was consumed.
	 */
	bool parseArgs(int& idx, char** argv);

	/**
	 * Returns true if the point satisfied the filter conditions.
	 * \param pt A Point instance.
	 * \return True if the point should be kept.
	 */
	bool keep(const geo::pc::Point& pt) const;

	/**
	 * Filters a container of geo::pc::Point, adding them to the given iterator.
	 * \param begin An iterator.
	 * \param end The end iterator.
	 * \param iter A back_inserter or other iterator.
	 * \return The number of points added to iter.
	 */
	template <class T, class U>
	int filter(T begin, T end, U iter) const {
		int i = 0;
		for(;begin != end; ++begin) {
			if(keep(*begin)) {
				iter = *begin;
				++i;
			}
		}
		return i;
	}

	void addClassFilter(int cls);

	void addClassFilter(const std::vector<int>& cls);

	void addIntensityFilter(double min, double max);

	void addZRangeFilter(double min, double max);

	double minZRange() const;

	double maxZRange() const;

	void addScanAngleFilter(double min, double max);

	void addKeepLastFilter();

	void addKeepFirstFilter();

	void addRejectEdgeFilter();

	void print() const;

};


/**
 * Given a terrain model, normalizes a point cloud so that all point
 * elevations are relative to zero.
 */
class Normalizer {
private:
	std::vector<std::string> m_filenames;	///< A list of input point cloud files.
	PCPointFilter* m_filter;				///< A point cloud filter instance.

public:

	/**
	 * Create a Normalizer with a list of input files.
	 * \param filenames A list of input files.
	 */
	Normalizer(const std::vector<std::string> filenames);

	/**
	 * Set the PCPointFilter instance. The Normalizer copies the instance.
	 * \param filter A PCPointFilter instance.
	 */
	void setFilter(const PCPointFilter& filter);

	/**
	 * Normalize the point cloud.
	 * \param dtm 		The filename of a digital terrain model.
	 * \param outdir 	The output directory for normalized point cloud files.
	 * \param band		The DTM raster band to use for normalization.
	 * \param force     Force overwrite of existing files.
	 */
	void normalize(const std::string& dtm, const std::string& outdir, int band = 1, bool force = false);

	/**
	 * Destroy the normalizer.
	 */
	~Normalizer();

};

/**
 * Turns a set of point cloud files into a raster.
 */
class Rasterizer {
private:
	std::vector<PCFile> m_files;		///< A list of PCFile instances.
	PCPointFilter* m_filter;
	int m_thin;

	// Used by finalizer.
	std::vector<std::unique_ptr<Computer> > m_computers;
	std::vector<geo::grid::Grid<double>> m_rasters;
	std::vector<geo::pc::Point> m_filtered;
	std::vector<double> m_out;

	void finalize(int row, double radius,
			std::unordered_map<size_t, std::vector<geo::pc::Point> >& cells,
			geo::grid::Grid<double>& outrast);

public:

	/**
	 * Construct a Rasterizer using the given point cloud file names.
	 * \param filenames A vector containing file names.
	 */
	Rasterizer(const std::vector<std::string> filenames);

	/**
	 * Returns a mapping of computer names and brief descriptions.
	 */
	static const std::unordered_map<std::string, std::string>& availableComputers();

	/**
	 * Rasterize the point cloud using the given output filename, statistic type, resolution
	 * and bounds. The radius gives the size of the neighbourhood around each pixel
	 * centre.
	 * \param filename 	The output (raster) filename.
	 * \param types		The list of statistics to compute.
	 * \param resX 		The raster x resolution.
	 * \param resY 		The raster y resolution.
	 * \param easting 	The minimum corner coordinate of the raster.
	 * \param northing 	The minimum corner coordinate of the raster.
	 * \param radius 	The size of the neighbourhood around each cell centre.
	 * \param srid 		The spatial reference ID of the output.
	 * \param memory    The number of bytes to consume in RAM before changing to disk-backed storage.
	 * \param useHeader True to trust the LAS file headers for things like bounds. Otherwise, read the points.

	 */
	void rasterize(const std::string& filename, const std::vector<std::string>& types, double resX, double resY,
		double easting, double northing, double radius, int srid, int memory, bool useHeader);

	/**
	 * Esitmate the point density (per cell) given the source files, resolution and search radius.
	 * \param The grid resolution.
	 * \param The search radius. Set to zero to use the cell bounds.
	 * \return The estimated point density per cell.
	 */
	double density(double resolution, double radius);

	/**
	 * If given and larger than zero, any cell with more than this number
	 * of points will be randomly thinned. Any cell with fewer will be zeroed
	 * out.
	 *
	 * \param thin The number of points in each cell.
	 */
	void setThin(int thin);

	/**
	 * Set a point filter to use for filtering points. Removes the old filter.
	 * \param filter A PointFilter.
	 */
	void setFilter(PCPointFilter* filter);

	/**
	 * Return a pointer to the PointFilter.
	 *
	 * \return A pointer to the PointFilter.
	 */
	PCPointFilter* filter() const;

	/**
	 * Destroy the Rasterizer.
	 */
	~Rasterizer();

};

/**
 * Provides a means of iterating over a point cloud by means of successive
 * kd-tree instances, covering geographic tiles with a given buffer zone.
 * This enables the program to seamlessly scan over a set of point cloud
 * files without edge effects and without exhausting the machine's memory.
 */
class PCTreeIterator {
private:
	int m_idx;
	int m_cols;
	int m_rows;
	double m_minX;
	double m_minY;
	double m_size;
	double m_buffer;
	std::vector<geo::pc::PCFile> m_files;

public:
	/**
	 * Build a PCTreeIterator with the given source files, size and buffer. The buffer
	 * should be sized to allow operations that use a neighbourhood without the
	 * neighbourhood extending past the edge of the tile, where possible. The
	 * size should be chosen to make efficient use of the machine's physical memory
	 * without exhausting it.
	 *
	 * The main body of the tile, within the square defined by the size, will not overlap
	 * any other tile. The buffers will overlap.
	 *
	 * \param files   The filenames of the point cloud files; preferrably non-overlapping tiles.
	 * \param size    The size of the "tile" represented by each tree, not including buffer.
	 * \param buffer  The size of the buffer around the tile.
	 */
	PCTreeIterator(const std::vector<std::string>& files, double size, double buffer);

	void init();

	void reset();

	bool next(geo::ds::mKDTree<geo::pc::Point>& tree);

};

bool pointSort(const geo::pc::Point& a, const geo::pc::Point& b);

template <class T, class U>
int pointFilter(U begin, U end, T insert, geo::pc::PCPointFilter* filter) {
	int count = 0;
	for(auto& it = begin; it < end; ++it) {
		if(filter && !filter->keep(*it)) continue;
		insert = *it;
		++count;
	}
	return count;
}

} // pc
} // geo

#endif
