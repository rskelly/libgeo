/*
 * pc_memgrid.cpp
 *
 *  Created on: Feb 25, 2018
 *      Author: rob
 */

#include "pc_memgrid.hpp"

using namespace geo::pc;

char* MemGrid::data() {
	// Return data pointer either from the mapped set or the in-memory set.
	return m_mem == nullptr ? (char*) m_mapped.data() : m_mem;
}

void MemGrid::resize(size_t size) {
	if(size > m_memLimit) {
		g_debug("MemGrid resize (mapped): " << size << "; limit: " << m_memLimit);
		if(!m_mapped.data()) {
			if(m_mapFile.empty()) {
				m_mapped.init(size, true);
			} else {
				m_mapped.init(m_mapFile, size, true);
			}
			if(m_mem) {
				std::memcpy(m_mapped.data(), m_mem, m_totalLength);
				free(m_mem);
				m_mem = nullptr;
			}
		} else {
			m_mapped.reset(size);
		}
	} else if(m_mem) {
		g_debug("MemGrid resize (realloc): " << size << "; limit: " << m_memLimit);
		m_mem = (char*) realloc((void*) m_mem, size);
	} else {
		g_debug("MemGrid resize (alloc): " << size << "; limit: " << m_memLimit);
		m_mem = (char*) malloc(size);
	}
	m_totalLength = size;
}

void MemGrid::writePoint(size_t idx, const geo::pc::Point& pt) { //float x, float y, float z, float intensity, float angle, int cls, int retNum, int numRets, int edge) {

	size_t offset;

	if(m_map.find(idx) != m_map.end()) {
		// If the index is in the map, find the offset.
		offset = m_map[idx];
	} else {
		// The offset is equal to the index * the line length.
		offset = m_map[idx] = idx;
		m_countMap[idx] = 0;
	}

	// Repurpose the offset as a pointer into memory for the current line.
	offset *= m_lineLength;

	// If the new offset is too far, resize the data.
	if(offset + m_lineLength > m_totalLength)
		resize((offset + m_lineLength) * 2);

	size_t count = m_countMap[idx];
	MappedLine ml = {idx, 0, count};
	//MappedPoint mp = {x, y, z, intensity, angle, cls, retNum, numRets, edge};

	// If the line is full, update the nextLine and grab a finalized one or start a new one.
	if(count == m_lineCount) {

		// Update the nextline field and copy the line header to the previous line
		ml.nextLine = m_currentLine++;
		std::memcpy(data() + offset, &ml, sizeof(MappedLine));

		// If the new offset is too far, resize the memory.
		if(ml.nextLine * m_lineLength + m_lineLength > m_totalLength)
			resize((ml.nextLine * m_lineLength + m_lineLength) * 2);

		// Create the new line.
		m_map[idx] = ml.nextLine;
		offset = ml.nextLine * m_lineLength;
		ml.nextLine = 0;
		ml.count = 1;
		std::memcpy(data() + offset, &ml, sizeof(MappedLine));
		std::memcpy(data() + offset + sizeof(MappedLine), &pt, sizeof(geo::pc::Point));
		m_countMap[idx] = 1;
		++m_pointCount;

	} else {

		// Add the point to the line.
		ml.count++;
		std::memcpy(data() + offset, &ml, sizeof(MappedLine));
		offset += sizeof(MappedLine) + count * sizeof(geo::pc::Point);
		std::memcpy(data() + offset, &pt, sizeof(geo::pc::Point));
		++m_countMap[idx];
		++m_pointCount;

	}

}

void MemGrid::finalize(size_t idx) {
	m_map.erase(idx);
	m_countMap.erase(idx);
}

size_t MemGrid::readPoints(size_t idx, std::vector<geo::pc::Point>& pts, bool final) {
	if(!hasUnread(idx))
		return 0; //g_runerr("No unread pixel for that index.");
	size_t offset = idx * m_lineLength;
	size_t count = 0;
	MappedLine ml;
	std::vector<geo::pc::Point> ptbuf(m_lineCount);
	std::vector<char> buf(m_lineLength);
	do {
		std::memcpy(buf.data(), data() + offset, m_lineLength);
		std::memcpy(&ml, buf.data(), sizeof(MappedLine));
		std::memcpy(ptbuf.data(), buf.data() + sizeof(MappedLine), sizeof(geo::pc::Point) * m_lineCount);
		pts.insert(pts.end(), ptbuf.begin(), ptbuf.begin() + ml.count);
		offset = ml.nextLine * m_lineLength;
	} while(ml.nextLine);
	if(final) {
		finalize(idx);
		m_pointCount -= count;
	}
	return count;
}

MemGrid::MemGrid() :
	m_mem(nullptr),
	m_lineCount(0),
	m_cellCount(0),
	m_lineLength(0),
	m_totalLength(0),
	m_currentLine(0),
	m_pointCount(0),
	m_memLimit(0) {
}

MemGrid::MemGrid(size_t cellCount, size_t memLimit) :
	MemGrid() {
	init(cellCount, memLimit);
}


void MemGrid::init(size_t cellCount, size_t memLimit) {
	init("", cellCount, memLimit);
}

/**
 * Initialize the MemGrid.
 * @param mapFile The filename of the map file, if there is one. Empty string otherwise.
 * @param cellCount An initial estimate of the number of rows.
 */
void MemGrid::init(const std::string& mapFile, size_t cellCount, size_t memLimit) {
	m_lineLength = sysconf(_SC_PAGESIZE); //sizeof(MappedLine) + sizeof(MappedPoint) * m_lineCount;
	m_lineCount = (m_lineLength - sizeof(MappedLine)) / sizeof(geo::pc::Point); // TODO: Compute expected point count for each index to minimize jumping.
	m_cellCount = cellCount;
	m_totalLength = cellCount * m_lineLength;
	m_mapFile = mapFile;
	m_memLimit = std::max((size_t) 0, memLimit);
	m_currentLine = m_cellCount + 1;
	g_trace("MemGrid.init: cellCount: " << cellCount << "; lineCount: " << m_lineCount << "; lineLength: " << m_lineLength);

	resize(m_totalLength);
}

size_t MemGrid::pointCount() const {
	return m_pointCount;
}

bool MemGrid::hasUnread(size_t idx) const {
	return m_map.find(idx) != m_map.end();
}

const std::unordered_map<size_t, size_t>& MemGrid::indexMap() const {
	return m_map;
}

/*
void MemGrid::add(size_t idx, double x, double y, double z, double intensity, double angle, int cls, int retNum, int numRets, int edge) {
	writePoint(idx, x, y, z, intensity, angle, cls, retNum, numRets, edge);
}
*/

void MemGrid::add(size_t idx, const geo::pc::Point& pt) {
	writePoint(idx, pt);//pt.x(), pt.y(), pt.z(), pt.intensity(), pt.scanAngle(), pt.classId(), pt.returnNum(), pt.numReturns(), pt.isEdge());
}

size_t MemGrid::get(size_t idx, std::vector<geo::pc::Point>& out, bool final) {
	return readPoints(idx, out, final);
}

void MemGrid::flush() {
	m_mapped.flush();
}

MemGrid::~MemGrid() {
	if(m_mem)
		free(m_mem);
	if(m_mapFile.empty())
		Util::rm(m_mapped.name());
}

