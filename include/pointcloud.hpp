#ifndef __PC_HPP__
#define __PC_HPP__

/**
 * Classes for working with LiDAR point clouds, specifically LAS files.

 *  Created on: Apr 13, 2017
 *  Author: rob
 */

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
#include "raster.hpp"
#include "ds/kdtree.hpp"

#define DBL_MAX std::numeric_limits<double>::max()
#define DBL_MIN std::numeric_limits<double>::lowest()

using namespace geo::util;
using namespace geo::ds;

namespace geo {
namespace pc {

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
	std::vector<std::string> m_filenames;	///< The list of filenames of consituent files.

public:

	/**
	 * Construct a PCFile object using the given filename and corner coordinate.
	 * @param filename 	The filename of a consituent file.
	 * @param x 		The minimum x-coordinate.
	 * @param y 		The minimum y-coordinate.
	 */
	PCFile(const std::string& filename, double x = 0, double y = 0, double size = 0, double buffer = 0);

	/**
	 * Construct a PCFile object using the given filenames and corner coordinate.
	 * @param filenames	A list of filenames of consituent files.
	 * @param x 		The minimum x-coordinate.
	 * @param y 		The minimum y-coordinate.
	 */
	PCFile(const std::vector<std::string>& filenames, double x = 0, double y = 0, double size = 0, double buffer = 0);

	void resize(double x, double y, double size, double buffer);

	/**
	 * Get the minimum x-coordinate.
	 * @return The minimum x-coordinate.
	 */
	double x() const;

	/**
	 * Get the minimum y-coordinate.
	 * @return The minimum y-coordinate.
	 */
	double y() const;

	/**
	 * Populates the given array with the bounds of this instance.
	 * There are six elements.
	 * @param bounds A six-element array of doubles, which is overwritten.
	 */
	void fileBounds(double* bounds) const;

	/**
	 * Populates the given array with the nominal bounds of this instance.
	 * There are four elements.
	 * @param bounds A four-element array of doubles, which is overwritten.
	 */
	void bounds(double* bounds) const;

	/**
	 * Populates the given array with the buffered nominal bounds of this instance.
	 * There are four elements.
	 * @param bounds A four-element array of doubles, which is overwritten.
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
	 * Returns true if the point is within the buffered bounds of the tile.
	 */
	bool containsBuffered(double x, double y) const;

	/**
	 * Returns the point count.
	 * @return The point count.
	 */
	size_t pointCount() const;

	/**
	 * Initialize the last file. Currently just computes its bounds. If
	 * the {useHeader} parameter is false, computes the actual bounds
	 * from the point cloud, otherwise uses the header bounds.
	 * @param useHeader Use the header to determine bounds (fast), rather than the point cloud (slow).
	 */
	void init(bool useHeader = true);

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
	 * @return The new filename.
	 */
	std::string nextFile();

public:
	/**
	 * Construct an instance of a PCWriter using the given output filename, {liblas::Header} and corner
	 * coordinate.
	 * @param filename 	The output filename template. This is the file path without an index part or extension.
	 * @param hdr 		A {liblas::Header} to be used as a template for the new file.
	 * @param x 		The minimum x-coordinate of the tile.
	 * @param y 		The minimum y-coordinate of the tile.
	 */
	PCWriter(const std::string& filename, const liblas::Header& hdr, double x, double y, double size, double buffer);

	PCWriter(PCWriter&& other);

	/**
	 * Get the minimum x-coordinate.
	 * @return The minimum x-coordinate.
	 */
	double x() const;

	/**
	 * Get the minimum y-coordinate.
	 * @return The minimum y-coordinate.
	 */
	double y() const;

	/**
	 * Get the tile size for this writer.
	 * @return The tile size.
	 */
	double size() const;

	/**
	 * The number of points written.
	 * @return The number of points written.
	 */
	size_t count() const;

	/**
	 * Populates the given array with the bounds of this instance.
	 * There are six elements.
	 * @param bounds A six-element array of doubles, which is overwritten.
	 */
	void outBounds(double* bounds) const;

	/**
	 * Populates the given array with the buffered bounds of this instance.
	 * There are six elements.
	 * @param bounds A four-element array of doubles, which is overwritten.
	 */
	void bufferedBounds(double* bounds) const;

	/**
	 * Populates the given array with the nominal bounds of this instance.
	 * There are six elements.
	 * @param bounds A four-element array of doubles, which is overwritten.
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
	 * @param pt A point.
	 */
	void addPoint(const liblas::Point& pt);

	/**
	 * Returns true if the given coordinate is within the bounds of this tile's
	 * geographic extent. "Within" means greater than or equal to the minimum
	 * bound and less than the maximum.
	 * @param x The x-coordinate.
	 * @param y The y-coordinate.
	 */
	bool contains(double x, double y) const;

	/**
	 * Returns true if the given coordinate is within the bounds of this tile's
	 * buffered extent. "Within" means greater than or equal to the minimum
	 * bound and less than the maximum.
	 * @param x The x-coordinate.
	 * @param y The y-coordinate.
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
	 * @param minx 		The minimum x-coordinate.
	 * @param miny 		The minimum y-coordinate.
	 * @param maxx 	 	The maximum x-coordinate.
	 * @param maxy 		The maximum y-coordinate.
	 * @param buffer 	The buffer distance.
	 */
	Tile(double minx, double miny, double maxx, double maxy, double buffer = 0);

	/**
	 * Return a pointer to the PCWriter owned by this tile.
	 * @param release If true, the tile will relinquish ownership of the 
	 * PCWriter instance.
	 * @return A pointer to the \ling PCWriter owned by this tile.
	 */
	PCWriter* writer(bool release = false);

	/**
	 * Set a PCWriter instance on this tile. This tile is
	 * the unique owner of the instance. Failure will
	 * result if the caller destroys the PCWriter.
	 * @param wtr A PCWriter pointer.
	 */
	void writer(PCWriter* wtr);

	/**
	 * Returns true if the given coordinate is within the bounds of this tile's
	 * geographic extent. "Within" means greater than or equal to the minimum 
	 * bound and less than the maximum.
	 * @param x The x-coordinate.
	 * @param y The y-coordinate.
	 */
	bool contains(double x, double y);

	/**
	 * Returns true if the given coordinate is within the bounds of this tile's
	 * buffered extent. "Within" means greater than or equal to the minimum 
	 * bound and less than the maximum.
	 * @param x The x-coordinate.
	 * @param y The y-coordinate.
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
	 * @param filenames A list of filenames of point cloud files to process.
	 */
	Tiler(const std::vector<std::string> filenames);

	/**
	 * Tile the list of point cloud files; write the output to the given output directory.
	 * @param outdir			The output directory.
	 * @param size 				The length of one side of the tile.
	 * @param buffer 			The size of the buffer around each tile. Default 0.
	 * @param srid 				The spatial reference ID. Default 0.
	 * @param easting   		The minimum x-coordinate of the tile grid. If NaN is given, will round
	 *							the minimum cordinate of the data extent down to a multiple of the tile size.
	 *							Default NaN.
	 * @param northing  		The minimum y-coordinate of the tile grid. If nan is given, will round
	 *							the minimum cordinate of the data extent down to a multiple of the tile size.
	 *							Default NaN.
	 * @param maxFileHandles 	The maximum number of output files that can be open at once. This determines
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
	double m_intensity;
	double m_angle;
	int m_cls;
	int m_returnNum;
	int m_numReturns;
	bool m_isEdge;
	liblas::Point* m_point;	///< The (optional) liblas::Point instance. Deleted on destruction.

public:

	/**
	 * Construct a Point using a liblas::Point. Will read the 3D
	 * coordinates and save a pointer to a copy of the object.
	 * The copy is deleted on destruction.
	 * @param pt A liblas::Point object.
	 */
	Point(const liblas::Point& pt);

	/**
	 * Construct a point using the three raw coordinates.
	 * @param x The x-coordinate.
	 * @param y The y-coordinate.
	 * @param z The z-coordinate.
	 */
	Point(double x, double y, double z);

	/**
	 * Construct a point using the three raw coordinates.
	 * @param x The x-coordinate.
	 * @param y The y-coordinate.
	 * @param z The z-coordinate.
	 */
	Point(double x, double y, double z, double intensity, double angle, int cls, int returnNum, int numReturns, bool isEdge);

	/**
	 * Construct an empty point.
	 */
	Point();

	/**
	 * Returns true if the point is a last return.
	 * @return True if the point is a last return.
	 */
	bool isLast() const;

	/**
	 * Returns true if the point is a first return.
	 * @return True if the point is a first return.
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
	 * @param idx The index.
	 */
	double operator[](int idx) const;

	/**
	 * Return the x-coordinate.
	 * @return The x-coordinate.
	 */
	double x() const;

	/**
	 * Return the y-coordinate.
	 * @return The y-coordinate.
	 */
	double y() const;

	/**
	 * Return the z-coordinate.
	 * @return The z-coordinate.
	 */
	double z() const;

	/**
	 * Return the intensity.
	 * @return The intensity.
	 */
	double intensity() const;

	/**
	 * Return the scan angle.
	 * @return The scane angle.
	 */
	double scanAngle() const;

	/**
	 * Return the return number.
	 * @return The return number.
	 */
	int returnNum() const;

	int numReturns() const;

	/**
	 * Return true if the point is marked as a flight line edge.
	 * @return True if the point is marked as a flight line edge.
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

/**
 * Used by Rasterizer to compute cell values from the points found
 * in the neighbourhood of a cell.
 */
class Computer {
public:

	/**
		 * Compute and return the statistic for the points within the neighbourhood
		 * defined by the lists of points.
		 * @param type 	 The statistic to compute.
		 * @param pts 	 The list of Points.
		 * @param radius The neighbourhood to calculate within.
		 * @param out    The output list.
		 * @return The the number of bands computed.
		 */
		virtual int compute(double x, double y, const std::vector<geo::pc::Point>& pts, double radius, std::vector<double>& out, geo::pc::PCPointFilter* filter = nullptr) = 0;

		/**
		 * Return the number of bands that should be returned by the computer.
		 * @return The number of bands that should be returned by the computer.
		 */
		virtual int bandCount() const = 0;

		/**
		 * Destroy the computer.
		 */
		virtual ~Computer() {}
};

/**
 * A configurable filter for point clouds.
 */
class PCPointFilter {
public:
	std::vector<int> classes;	///< A list of classes to keep. If this list is empty, all classes are kept. Default empty.
	double minScanAngle;		///< The minimum scan angle. Default -90.
	double maxScanAngle;		///< The maximum scan angle. Default 90.
	bool keepEdges;				///< Keep the points marked as flight line edges. Default false.
	double minZ;				///< The maximum elevation to keep. Default minimum double value.
	double maxZ;				///< The minimum elevation to keep. Default maximum double value.
	double minIntensity;		///< The minimum intensity to keep. Default minimum double value.
	double maxIntensity;		///< The maximum intensity to keep.	Default maximum double value.
	bool lastOnly;				///< Only keep the last returns. Default false.
	bool firstOnly;				///< Only keep the first returns. Default false.

	/**
	 * Construct a PCPointFilter.
	 */
	PCPointFilter();

	/**
	 * Print the command-line parameters that can be used to
	 * configure a PCPointFilter.
	 * @param str An output stream.
	 */
	static void printHelp(std::ostream& str);

	/**
	 * Configure this object by reading the arguments from the command line.
	 * Use PCPointFilter::printHelp to print the expected arguments.
	 * If the argument at the current index is relevant to this item,
	 * return true to notify the caller. Updates the index if
	 * more than one item is read. The caller is responsible for
	 * incrementing the index before the next item.
	 * @param idx  The index of the current argument.
	 * @param argv The array of arguments.
	 * @return True, if the current argument was consumed.
	 */
	bool parseArgs(int& idx, char** argv);

	/**
	 * Returns true if the point satisfied the filter conditions.
	 * @param pt A Point instance.
	 * @return True if the point should be kept.
	 */
	bool keep(const geo::pc::Point& pt) const;

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
	 * @param filenames A list of input files.
	 */
	Normalizer(const std::vector<std::string> filenames);

	/**
	 * Set the PCPointFilter instance. The Normalizer copies the instance.
	 * @param filter A PCPointFilter instance.
	 */
	void setFilter(const PCPointFilter& filter);

	/**
	 * Normalize the point cloud.
	 * @param dtm 		The filename of a digital terrain model.
	 * @param outdir 	The output directory for normalized point cloud files.
	 * @param band		The DTM raster band to use for normalization.
	 */
	void normalize(const std::string& dtm, const std::string& outdir, int band = 1);

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

public:

	/**
	 * Construct a Rasterizer using the given point cloud file names.
	 * @param filenames A vector containing file names.
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
	 * @param filename 	The output (raster) filename.
	 * @param types		The list of statistics to compute.
	 * @param res 		The raster resolution.
	 * @param easting 	The minimum corner coordinate of the raster.
	 * @param northing 	The minimum corner coordinate of the raster.
	 * @param radius 	The size of the neighbourhood around each cell centre.
	 * @param srid 		The spatial reference ID of the output.
	 * @param density   The estimated number of points per cell; used for reserving mapped memory. Default 32.
	 * @param ext 		An extra filter value. For example, for percentiles, this is the cutoff value.
	 * @param mapFile   An existing point-map file that can be loaded. When this is given, and
	 *                  it exists, the point file loading is skipped and the file is not deleted
	 *                  at completion. If it does not exist, it is created and used for loading points
	 *                  and not deleted. If it is not given, a temporary file is mapped and deleted
	 *                  on completion.
	 */
	void rasterize(const std::string& filename, const std::vector<std::string>& types, double res,
		double easting, double northing, double radius, int srid, int density = 32, double ext = 0,
		const std::string& mapFile = "");

	/**
	 * Esitmate the point density (per cell) given the source files, resolution and search radius.
	 * @param The grid resolution.
	 * @param The search radius. Set to zero to use the cell bounds.
	 * @return The estimated point density per cell.
	 */
	double density(double resolution, double radius);

	/**
	 * Set a point filter to use for filtering points. Removes the old filter.
	 * @param filter A PointFilter.
	 */
	void setFilter(const PCPointFilter& filter);

	/**
	 * Destroy the Rasterizer.
	 */
	~Rasterizer();

};

} // las
} // geo

#endif
