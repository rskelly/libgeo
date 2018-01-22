#ifndef __LAS_HPP__
#define __LAS_HPP__

/**
 * Classes for working with LiDAR point clouds, specifically LAS files.

 *  Created on: Apr 13, 2017
 *  Author: rob
 */

#include <fstream>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <memory>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdint>

#include <liblas/liblas.hpp>

#include "util.hpp"

using namespace geo::util;

namespace geo {
namespace las {

/**
 * A class representing a source LAS file. Maintains
 * the position of the tile, by the minimum x and y coordinates, 
 * the tiles bounding box and a list of filenames of the tile's
 * constituent files.
 */
class LASFile {
private:
	double m_x;								///< Minimum corner x-coordinate.
	double m_y;								///< Minimum corner y-coordinate.
	double m_bounds[6];						///< The bounding box of the actual point cloud (may differ from tile bounds.)
	std::vector<std::string> m_filenames;	///< The list of filenames of consituent files.

public:

	/**
	 * Construct a LASFile object using the given filename and corner coordinate.
	 * @param filename 	The filename of a consituent file.
	 * @param x 		The minimum x-coordinate.
	 * @param y 		The minimum y-coordinate.
	 */
	LASFile(const std::string& filename, double x = 0, double y = 0);

	/**
	 * Construct a LASFile object using the given filenames and corner coordinate.
	 * @param filenames	A list of filenames of consituent files.
	 * @param x 		The minimum x-coordinate.
	 * @param y 		The minimum y-coordinate.
	 */
	LASFile(const std::vector<std::string>& filenames, double x = 0, double y = 0);

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
	void bounds(double* bounds) const;

	/**
	 * Returns a reference to the filenames list.
	 */
	const std::vector<std::string>& filenames() const;

	/**
	 * Initialize the last file. Currently just computes its bounds. If
	 * the {useHeader} parameter is false, computes the actual bounds
	 * from the point cloud, otherwise uses the header bounds.
	 * @param useHeader Use the header to determine bounds (fast), rather than the point cloud (slow).
	 */
	void init(bool useHeader = true);

	~LASFile();

};

/**
 * A class that is responsible for writing new tiles.
 */
class LASWriter {
private:
	int m_fileIdx;							///< The current file index.
	int m_returns;							///< The number of returns (points) in the current file.
	int m_retNum[5];						///< The number of points for each return in the current file.
	long m_totalReturns;					///< The total count of returns across all files.
	double m_outBounds[6];					///< The bounding box of the point cloud (may differ from the tile's bounds).
	double m_x;								///< The minimum x-coordinate of the tile.
	double m_y;								///< The minimum y-coordinate of the tile.
	std::string m_filename;					///< The output filename template.
	std::vector<std::string> m_filenames;	///< The list of output filenames.
	liblas::Writer* m_writer;				///< The current {liblas::Writer} instance.
	liblas::Header* m_header;				///< The liblas {liblas::Header} instance.
	std::ofstream m_str;					///< The current output stream.
	bool m_dod;								///< Set to true if constituent files should be deleted on destruction.

	/**
	 * Creates and returns the next filename in the series for this tile.
	 * @return The new filename.
	 */
	std::string nextFile();

public:
	/**
	 * Construct an instance of a LASWriter using the given output filename, {liblas::Header} and corner
	 * coordinate.
	 * @param filename 	The output filename template. This is the file path without an index part or extension.
	 * @param hdr 		A {liblas::Header} to be used as a template for the new file.
	 * @param x 		The minimum x-coordinate of the tile.
	 * @param y 		The minimum y-coordinate of the tile.
	 */
	LASWriter(const std::string& filename, const liblas::Header& hdr, double x = 0, double y = 0);

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
	 * Destroy the writer. Will close any open files, and, if {dod} is set to
	 * true, will delete constituent files.
	 */
	~LASWriter();

};

/**
 * Represents a single tile, its geographic extent and its buffered extent.
 * Contains a LASWriter instance for writing to the file(s) associated
 * with the tile.
 */
class Tile {
private:
	double m_bounds[4];						///< The geographic bounds of the tile.
	double m_bufferedBounds[4];				///< The buffered extent of the tile.
	std::unique_ptr<LASWriter> m_writer;	///< The LASWriter used for writing points.

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
	 * Return a pointer to the LASWriter owned by this tile.
	 * @param release If true, the tile will relinquish ownership of the 
	 * LASWriter instance.
	 * @return A pointer to the \ling LASWriter owned by this tile.
	 */
	LASWriter* writer(bool release = false);

	/**
	 * Set a LASWriter instance on this tile. This tile is
	 * the unique owner of the instance. Failure will
	 * result if the caller destroys the LASWriter.
	 * @param wtr A LASWriter pointer.
	 */
	void writer(LASWriter* wtr);

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
 * Performs the tiling of a set of LAS files.
 */
class Tiler {
public:
	std::vector<LASFile> files;	///< The list of LASFile instances.

	/**
	 * Construct a tiler with the given list of initial filenames.
	 * @param filenames A list of filenames of LAS files to process.
	 */
	Tiler(const std::vector<std::string> filenames);

	/**
	 * Tile the list of LAS files; write the output to the given output directory.
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
	 * @param maxFileHandles 	The maximum number of output files that can be open at once. 
	 *							Default 64.
	 */
	void tile(const std::string& outdir, double size, double buffer = 0, int srid = 0, 
		double easting = std::nan(""), double northing = std::nan(""), int maxFileHandles = 64);

	~Tiler();

};

} // las
} // geo

#endif
