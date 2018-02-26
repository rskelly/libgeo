/*
 * pc_memgrid.cpp
 *
 *  Created on: Feb 25, 2018
 *      Author: rob
 */

#include "pc_memgrid.hpp"

using namespace geo::pc;

void MemGrid::resize(size_t size) {
	g_debug("MemGrid resize: " << size);
	m_mapped.reset(size);
	m_totalLength = size;
}

void MemGrid::writePoint(size_t idx, float x, float y, float z, float intensity, float angle, int cls, int retNum, int numRets, int edge) {

	char* data = (char*) m_mapped.data();
	size_t offset;
	bool newLine;

	if(m_map.find(idx) != m_map.end()) {
		offset = m_map[idx];
		newLine = false;
	} else {
		if(!m_finalized.empty()) {
			offset = m_finalized.front();
			m_finalized.pop_front();
		} else {
			offset = m_currentLine++;
		}
		m_map[idx] = offset;
		newLine = true;
	}
	offset *= m_lineLength;

	if(offset + m_lineLength >= m_totalLength) {
		resize(m_totalLength * 2);
		data = (char*) m_mapped.data();
	}

	MappedPoint mp = {x, y, z, intensity, angle, cls, retNum, numRets, edge};
	MappedLine ml = {0, 0, 0};

	if(newLine) {

		// There is no record in the points file; start one.
		ml.idx = idx;
		ml.count = 1;
		ml.nextLine = 0;
		std::memcpy(data + offset, &ml, sizeof(MappedLine));
		std::memcpy(data + offset + sizeof(MappedLine), &mp, sizeof(MappedPoint));
		++m_pointCount;

	} else {

		// Retrieve the line header at the offset.
		std::memcpy(&ml, data + offset, sizeof(MappedLine));

		// Find the last line for this cell.
		while(ml.nextLine) {
			offset = ml.nextLine * m_lineLength;
			std::memcpy(&ml, data + offset, sizeof(MappedLine));
		}

		// If the line is full, update the nextLine and grab a finalized one or start a new one.
		if(ml.count == m_lineCount) {

			if(!m_finalized.empty()) {
				ml.nextLine = m_finalized.front();
				m_finalized.pop_front();
			} else {
				ml.nextLine = m_currentLine++;
			}

			// If the new offset is too far, resize the memory.
			if(ml.nextLine * m_lineLength + m_lineLength >= m_totalLength) {
				resize(m_totalLength * 2);
				data = (char*) m_mapped.data();
			}

			// Copy the line header to the previous line
			std::memcpy(data + offset, &ml, sizeof(MappedLine));

			// Create the new line.
			offset = ml.nextLine * m_lineLength;
			ml.nextLine = 0;
			ml.count = 1;
			std::memcpy(data + offset, &ml, sizeof(MappedLine));
			std::memcpy(data + offset + sizeof(MappedLine), &mp, sizeof(MappedPoint));
			++m_pointCount;

		} else {

			// Add the point toe the line.
			size_t count = ml.count++;
			std::memcpy(data + offset, &ml, sizeof(MappedLine));
			offset += sizeof(MappedLine) + count * sizeof(MappedPoint);
			std::memcpy(data + offset, &mp, sizeof(MappedPoint));
			++m_pointCount;

		}
	}
}

void MemGrid::finalize(size_t idx) {
	size_t offset = m_map[idx];
	m_finalized.push_back(offset);
	char* data = (char*) m_mapped.data();
	MappedLine ml = {0, 0, 0};
	std::memcpy(&ml, data + offset * m_lineLength, sizeof(MappedLine));
	while(ml.nextLine) {
		m_finalized.push_back(ml.nextLine);
		std::memcpy(&ml, data + ml.nextLine * m_lineLength, sizeof(MappedLine));
	}
	m_map.erase(idx);
}

size_t MemGrid::readPoints(size_t idx, std::vector<geo::pc::Point>& pts, bool final) {
	if(!hasUnread(idx))
		return 0; //g_runerr("No unread pixel for that index.");
	char* data = (char*) m_mapped.data();
	size_t offset = m_map[idx] * m_lineLength;
	size_t count = 0;
	MappedLine ml;
	MappedPoint mp;
	do {
		char* buf = data + offset;
		std::memcpy(&ml, buf, sizeof(MappedLine));
		buf += sizeof(MappedLine);
		for(size_t i = 0; i < ml.count; ++i) {
			std::memcpy(&mp, buf, sizeof(MappedPoint));
			pts.emplace_back(mp.x, mp.y, mp.z, mp.intensity, mp.angle, mp.cls, mp.retNum, mp.numRets, mp.edge);
			buf += sizeof(MappedPoint);
			++count;
		}
		offset = ml.nextLine * m_lineLength;
	} while(ml.nextLine);
	if(final)
		finalize(idx);
	m_pointCount -= count;
	return count;
}

MemGrid::MemGrid() :
	m_lineCount(0),
	m_cellCount(0),
	m_lineLength(0),
	m_totalLength(0),
	m_currentLine(0),
	m_pointCount(0) {
}

MemGrid::MemGrid(size_t cellCount) :
	MemGrid() {
	init(cellCount);
}


void MemGrid::init(size_t cellCount) {
	init("", cellCount);
}

/**
 * Initialize the MemGrid.
 * @param mapFile The filename of the map file, if there is one. Empty string otherwise.
 * @param cellCount An initial estimate of the number of rows.
 */
void MemGrid::init(const std::string& mapFile, size_t cellCount) {
	m_lineCount = (sysconf(_SC_PAGESIZE) - sizeof(MappedLine)) / sizeof(MappedPoint);
	m_cellCount = cellCount;
	m_lineLength = sizeof(MappedLine) + sizeof(MappedPoint) * m_lineCount;
	m_totalLength = cellCount * m_lineLength;
	m_mapFile = mapFile;

	g_trace("MemGrid.init: cellCount: " << cellCount << "; lineCount: " << m_lineCount << "; lineLength: " << m_lineLength);

	if(mapFile.empty()) {
		m_mapped.init(m_totalLength, true);
	} else {
		m_mapped.init(mapFile, m_totalLength, true);
	}
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

void MemGrid::add(size_t idx, double x, double y, double z, double intensity, double angle, int cls, int retNum, int numRets, int edge) {
	writePoint(idx, x, y, z, intensity, angle, cls, retNum, numRets, edge);
}

size_t MemGrid::get(size_t idx, std::vector<geo::pc::Point>& out, bool final) {
	return readPoints(idx, out, final);
}

void MemGrid::flush() {
	m_mapped.flush();
}

MemGrid::~MemGrid() {
	if(m_mapFile.empty())
		Util::rm(m_mapped.name());
}

