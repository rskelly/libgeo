/*
 * pc_memgrid.hpp
 *
 *  Created on: Feb 25, 2018
 *      Author: rob
 */

#ifndef INCLUDE_PC_MEMGRID_HPP_
#define INCLUDE_PC_MEMGRID_HPP_

#include <util_old.hpp>
#include "pointcloud.hpp"

namespace geo {
namespace pc {

/**
 * Struct used to store the line headers in the point cloud cache file.
 */
typedef struct {
	size_t idx;			// The index of the cell: row * cols + col
	size_t nextLine;	// The address in the mapped memory of the next row containing data for this cell.
	size_t count;		// The number of points in the cell.
} MappedLine;

/**
 * Struct used to store point in the point cache file.
 */
/*
typedef struct {
	float x;
	float y;
	float z;
	float intensity;
	float angle;
	int cls;
	int retNum;
	int numRets;
	int edge;
} MappedPoint;
*/
class MemGrid {
private:
	char* m_mem;
	geo::util::MappedFile m_mapped;
	size_t m_lineCount;
	size_t m_cellCount;
	size_t m_lineLength;
	size_t m_totalLength;
	size_t m_currentLine;						// The next available line index.
	size_t m_pointCount;
	size_t m_memLimit;
	std::string m_mapFile;
	std::list<size_t> m_finalized; 				// The line offsets of finalized cells available for re-use.
	std::unordered_map<size_t, size_t> m_map;   // Maps cell indices to current line indices.
	std::unordered_map<size_t, size_t> m_countMap; // The number of points in the current row.

	/**
	 * Resize the mapped file.
	 * @param size The size of the mapped file in bytes.
	 */
	void resize(size_t size);

	/**
	 * Write a point.
	 * @param idx The cell index to write the point to.
	 * @param x The x coordinate.
	 * @param y The y coordinate.
	 * @param z The z coordinate.
	 * @param intensity The intensity.
	 * @param angle The scan angle.
	 * @param cls The classification ID.
	 * @param retNum The return number.
	 * @param numRets The number of returns.
	 * @param edge Whether the point is classified as an edge (1) or not (0).
	 */
	void writePoint(size_t idx, const geo::pc::Point& pt); //float x, float y, float z, float intensity, float angle, int cls, int retNum, int numRets, int edge);

	/**
	 * Finalize the cell with the given index. The mapped slot is marked for re-use.
	 * @param idx The cell index.
	 */
	void finalize(size_t idx);

	/**
	 * Read the points for the given cell index into the vector. Finalize
	 * if needed.
	 * @param idx The cell index.
	 * @param pts The point vector.
	 * @param final True if the cell should be finalized.
	 * @return The number of points read.
	 */
	size_t readPoints(size_t idx, std::vector<geo::pc::Point>& pts, bool final);

	/**
	 * Return the data pointer which may be mapped or RAM.
	 */
	char* data();

public:

	/**
	 * Construct a MemGrid with defaults.
	 */
	MemGrid();

	/**
	 * Construct a MemGrid with the given number of cells.
	 * @param cellCount The number of cells to start with.
	 * @param memLimit The amount of memory (bytes) to reserve upfront.
	 */
	MemGrid(size_t cellCount, size_t memLimit);

	/**
	 * Initialize a MemGrid with the given number of cells.
	 * @param cellCount The number of cells to start with.
	 * @param memLimit The amount of memory (bytes) to reserve upfront.
	 */
	void init(size_t cellCount, size_t memLimit);

	/**
	/**
	 * Initialize a MemGrid with the given number of cells.
	 * @param mapFile The filename of the map file, if there is one. Empty string otherwise.
	 * @param cellCount The number of cells to start with.
	 * @param memLimit The amount of memory (bytes) to reserve upfront.
	 */
	void init(const std::string& mapFile, size_t cellCount, size_t memLimit);

	/**
	 * Return the number of points stored in the grid.
	 * @return The number of points stored in the grid.
	 */
	size_t pointCount() const;

	/**
	 * Return true if there are unread points for the given index.
	 * @param idx The cell index.
	 * @return True if there are unread points for the given index.
	 */
	bool hasUnread(size_t idx) const;

	/**
	 * Return a reference to the index map.
	 * @return A reference to the index map.
	 */
	const std::unordered_map<size_t, size_t>& indexMap() const;

	/**
	 * Add a point.
	 * @param idx The cell index to write the point to.  This should help with spatial locality, like z-order, etc.
	 * @param x The x coordinate.
	 * @param y The y coordinate.
	 * @param z The z coordinate.
	 * @param intensity The intensity.
	 * @param angle The scan angle.
	 * @param cls The classification ID.
	 * @param retNum The return number.
	 * @param numRets The number of returns.
	 * @param edge Whether the point is classified as an edge (1) or not (0).
	 */
	//void add(size_t idx, double x, double y, double z, double intensity, double angle, int cls, int retNum, int numRets, int edge);

	void add(size_t idx, const geo::pc::Point& pt);

	//void add(const std::vector<geo::pc::Point>& pts);


	/**
	 * Get the points for the given cell index into the vector. Finalize
	 * if needed.
	 * @param idx The cell index.
	 * @param pts The point vector.
	 * @param final True if the cell should be finalized.
	 * @return The number of points read.
	 */
	size_t get(size_t idx, std::vector<geo::pc::Point>& out, bool final = true);

	/**
	 * Flush the information to disk.
	 */
	void flush();

	/**
	 * Destroy the MemGrid.
	 */
	~MemGrid();
};


} // pc
} // geo


#endif /* INCLUDE_PC_MEMGRID_HPP_ */
