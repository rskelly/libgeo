/*
 * externalmergesort.hpp
 *
 *  Created on: Apr 8, 2018
 *      Author: rob
 */

#ifndef INCLUDE_EXTERNALMERGESORT_HPP_
#define INCLUDE_EXTERNALMERGESORT_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <cstdio>
#include <thread>

#include <liblas/liblas.hpp>

#include "pointcloud.hpp"
#include "util.hpp"
#include "geo.hpp"
#include "raster.hpp"

using namespace geo::util;
using namespace geo::pc;
using namespace geo::raster;

/**
 * Interleave the coordinates' bits to produce an index appropriate for z-ordering.
 *
 * @param x The x-coordinate. A pre-scaled integer.
 * @param y The y-coordinate. A pre-scaled integer.
 * @return The index.
 */
size_t interleave(size_t x, size_t y) {
	size_t out = 0;
	for(size_t i = 0; i < sizeof(int) * 8; ++i) {
		out |= ((x >> i) & 1) << (i * 2);
		out |= ((y >> i) & 1) << (i * 2 + 1);
	}
	return out;
}

/**
 * Reverse the interleaving process, converting an index into
 * scaled x and y coordinates.
 *
 * @param index The Morton index.
 * @param x A reference to the x coordinate.
 * @param y A reference to the y coordinate.
 */
void uninterleave(size_t index, size_t& x, size_t& y) {
	x = 0; y = 0;
	for(size_t i = 0; i < sizeof(size_t) * 8; i += 2) {
		x |= ((index >> i) & 1) << (i >> 1);
		y |= ((index >> (i + 1)) & 1) << (i >> 1);
	}
}

size_t toIdx(int col, int row, int cols) {
	return row * cols + col;
}

void fromIdx(size_t idx, int cols, int& col, int& row) {
	col = idx % cols;
	row = idx / cols;
}

/**
 * Sort function for Points. Uses the Morton index for sorting, using
 * the normalized and scaled coordinates.
 *
 * @param a A Point.
 * @param b A Point.
 * @param scale The normalized coordinates are divided by this amount.
 * @param xmin The minimum x coordinate of the bounding box. Subtracted from x.
 * @param ymin The minimum y coordinate of the bounding box. Subtracted from y.
 * @return True if a is "less" than b according to the Morton ordering.
 */
bool pointSort(const geo::pc::Point& a, const geo::pc::Point& b, const GridProps& props) {
	return toIdx(props.toCol(a.x()), props.toRow(a.y()), props.cols()) < toIdx(props.toCol(b.x()), props.toRow(b.y()), props.cols());
}

/**
 * Sorts the vector of points beginning and ending at the given indices.
 * Points are normalized and scaled and sorted according to the Morton
 * ordering.
 *
 * @param pts A vector of Points.
 * @param begin The start index.
 * @param end The end index.
 * @param scale The normalized coordinates are divided by this amount.
 * @param xmin The minimum x coordinate of the bounding box. Subtracted from x.
 * @param ymin The minimum y coordinate of the bounding box. Subtracted from y.
 */
void sortPoints(std::vector<geo::pc::Point>* pts, size_t begin, size_t end, const GridProps& props) {
	std::sort(pts->begin() + begin, pts->begin() + end,
			std::bind(pointSort, std::placeholders::_1, std::placeholders::_2, props));
}

/**
 * Abstract class for loading points and returning them as
 * vectors.
 */
class PointLoader {
public:
	virtual int getPoints(std::vector<geo::pc::Point>& pts, int count) = 0;
	virtual void computeBounds(double& xmin, double& ymin, bool useHeader) = 0;
	virtual ~PointLoader() {}
};

/**
 * Implementation of PointLoader to load LAS files.
 * TODO: Duplicate of PFile.
 */
class LASPointLoader : public PointLoader {
private:
	liblas::Reader* m_rdr;
	std::ifstream* m_instr;

public:
	LASPointLoader(const std::string& filename) {
		liblas::ReaderFactory rf;
		m_instr = new std::ifstream(filename, std::ios::binary|std::ios::in);
		m_rdr = new liblas::Reader(rf.CreateWithStream(*m_instr));

	}

	void computeBounds(double& xmin, double& ymin, bool useHeader) {
		if(useHeader) {
			const liblas::Header& hdr = m_rdr->GetHeader();
			xmin = hdr.GetMinX();
			ymin = hdr.GetMinY();
		} else {
			xmin = G_DBL_MAX_POS;
			ymin = G_DBL_MAX_POS;
			while(m_rdr->ReadNextPoint()) {
				const liblas::Point& pt = m_rdr->GetPoint();
				xmin = std::min(xmin, pt.GetX());
				ymin = std::min(ymin, pt.GetY());
			}
		}
	}

	int getPoints(std::vector<geo::pc::Point>& pts, int count) {
		int i = 0;
		while(--count >= 0 && m_rdr->ReadNextPoint()) {
			pts.emplace_back(m_rdr->GetPoint());
			++i;
		}
		return i;
	}

	~LASPointLoader() {
		delete m_rdr;
		delete m_instr;
	}
};

/**
 * Loads a cache file and returns its points one by one.
 */
class PointStream {
private:
	std::vector<geo::pc::Point> m_pts;
	std::string m_filename;
	std::ifstream* m_instr;
	size_t m_size;
	size_t m_current;
	size_t m_chunkSize;
	size_t m_ptr;

	bool loadNext() {
		if(!m_instr) {
			// If no stream is active, start one.
			m_instr = new std::ifstream(m_filename, std::ios::binary|std::ios::in);
			m_instr->read((char*) &m_size, sizeof(size_t));
			// Move the pointer ahead past the count element.
			m_current += sizeof(size_t);
		}
		// If past the end of the stream, quit.
		if(m_current >= m_size)
			return false;
		// Get chunk size or number of remaining elements and read into pts.
		size_t size = m_current + m_chunkSize > m_size ? m_size - m_current : m_chunkSize;
		m_pts.resize(size);
		m_instr->read((char*) m_pts.data(), size * sizeof(geo::pc::Point));
		// Advance the pointer.
		m_current += size;
		// Reset the current pointer.
		m_ptr = 0;
		return true;
	}

public:

	PointStream(const std::string& filename, size_t chunkSize) :
		m_filename(filename),
		m_instr(nullptr),
		m_size(0),
		m_current(0),
		m_chunkSize(chunkSize),
		m_ptr(0) {
	}

	bool empty() {
		return m_ptr >= m_size;
	}

	void reset() {
		m_ptr = 0;
	}

	~PointStream() {
		if(m_instr) {
			m_instr->close();
			delete m_instr;
		}
	}

	geo::pc::Point* current() {
		if((m_pts.empty() || m_ptr >= m_pts.size())  && !loadNext())
			return nullptr;
		return &(m_pts[m_ptr]);
	}

	void pop() {
		++m_ptr;
	}

};

/**
 * Sorts point clouds and provides methods for returning the sorted points.
 */
class ExternalMergeSort {
private:
	std::vector<std::string> m_inputFiles;
	std::vector<std::string> m_chunkFiles;
	std::string m_outputFile;
	std::string m_tmpDir;
	int m_chunkSize;						///< The number of points in a single chunk.
	int m_numChunks;						///< The number of chunks to process simultaneously.
	int m_chunk; 							///< The ID of the current chunk file.

	std::vector<geo::pc::Point> m_buffer;
	std::ifstream* m_instr;
	size_t m_count;
	size_t m_position;
	size_t m_index;

	GridProps m_props;

	/**
	 * Return the next available name for a chunk file.
	 */
	std::string nextChunkFile() {
		return Util::pathJoin(m_tmpDir, "chunk_" + std::to_string(++m_chunk) + ".tmp");
	}

	/**
	 * Merge the chunks given by the files into a single chunk. The points
	 * will be merged in order according to their Morton index.
	 *
	 * @param chunkFiles A list of point files to merge.
	 * @return The name of the new merged chunk.
	 */
	std::string mergeChunks(const std::vector<std::string>& chunkFiles) {
		// Get a file name for the output and open a stream.
		std::string outFile = nextChunkFile();
		std::ofstream outStream(outFile, std::ios::binary|std::ios::out);
		// Skip ahead to make room for the count.
		outStream.seekp(sizeof(size_t), std::ios::beg);

		// Create a point stream for each input chunk file.
		std::list<PointStream> inStreams;
		for(const std::string& file : chunkFiles)
			inStreams.emplace_back(file, m_chunkSize);

		size_t minIdx;
		size_t count = 0;
		geo::pc::Point* minPt;
		std::vector<geo::pc::Point> buffer;
		PointStream* source;
		buffer.reserve(m_chunkSize);
		while(true) {
			minPt = nullptr;
			for(PointStream& ps : inStreams) {
				geo::pc::Point* a = ps.current();
				if(a && !minPt) {
					minPt = a;
					source = &ps;
				} else if(a && minPt && pointSort(*a, *minPt, m_props)) {
					minPt = a;
					source = &ps;
				}
			}
			if(minPt) {
				source->pop();
				// A point was found at the index. Remove it.
				buffer.push_back(*minPt);
				if(buffer.size() == m_chunkSize) {
					// Write the points to the output.
					g_debug("merge write");
					outStream.write((char*) buffer.data(), buffer.size() * sizeof(geo::pc::Point));
					count += buffer.size();
					buffer.clear();
				}
			} else {
				break;
			}
		}

		if(!buffer.empty()) {
			// Write the points to the output.
			g_debug("merge write final")
			outStream.write((char*) buffer.data(), buffer.size() * sizeof(geo::pc::Point));
			count += buffer.size();
			buffer.clear();
		}

		// Write the final point count.
		outStream.seekp(0, std::ios::beg);
		outStream.write((char*) &count, sizeof(size_t));

		return outFile;
	}

	/**
	 * Read chunks of each file each file, sort them
	 * and write them to chunk files.
	 */
	void phase1() {
		int threadCount = 8;
		int size;
		// Iterate over the list of source files.
		for(const std::string& inputFile : m_inputFiles) {
			LASPointLoader ldr(inputFile);
			while(true) {
				std::list<std::thread> threads;
				std::vector<geo::pc::Point> pts;
				pts.reserve(threadCount * m_chunkSize);
				if((size = ldr.getPoints(pts, threadCount * m_chunkSize))) {
					for(int i = 0; i < size / m_chunkSize; ++i) {
						size_t begin = i * m_chunkSize;
						size_t end = std::min(size, (i + 1) * m_chunkSize);
						threads.emplace_back(sortPoints, &pts, begin, end, m_props);
					}
				} else {
					break;
				}
				for(std::thread& t : threads)
					t.join();
				for(int i = 0; i <= size / m_chunkSize; ++i) {
					size_t begin = i * m_chunkSize;
					size_t end = std::min(size, (i + 1) * m_chunkSize);
					// Create and open an output stream.
					std::string chunkFile = nextChunkFile();
					std::ofstream ostr(chunkFile, std::ios::binary|std::ios::out);
					std::vector<geo::pc::Point> data(pts.begin() + begin, pts.begin() + end);
					// Write the points to the output stream.
					size_t len = end - begin;
					ostr.write((char*) &len, sizeof(size_t));
					ostr.write((char*) data.data(), len * sizeof(geo::pc::Point));
					// Save the chunk file name.
					m_chunkFiles.push_back(chunkFile);
				}
			}
		}
	}

	/**
	 * Merge the chunk files. If there are more of these than numChunks,
	 * this must be done in multiple passes.
	 */
	void phase2() {
		// At the end there will be only one file.
		while(m_chunkFiles.size() > 1) {
			// Track the new merged chunks.
			std::vector<std::string> newChunks;
			std::vector<std::string> tmp;
			for(size_t i = 0; i < m_chunkFiles.size(); i += m_numChunks) {
				// Get a slice of the chunk file list.
				if(i + m_numChunks > m_chunkFiles.size()) {
					tmp.assign(m_chunkFiles.begin() + i, m_chunkFiles.end());
				} else {
					tmp.assign(m_chunkFiles.begin() + i, m_chunkFiles.begin() + i + m_numChunks);
				}
				// Merge the chunks into a new chunk.
				g_debug("merge " << i << "-" << i + m_numChunks << " of " << m_chunkFiles.size())
				newChunks.push_back(mergeChunks(tmp));
			}
			for(const std::string& file : m_chunkFiles)
				Util::rm(file);
			// Replace the chunk files list with the new chunks list.
			m_chunkFiles.assign(newChunks.begin(), newChunks.end());
		}
	}

	/**
	 * Load the next sorted point file and prepare for reading.
	 */
	bool loadNext() {
		if(!m_instr)
			reset();
		if(m_position >= m_count)
			return false;
		size_t count = m_position + m_chunkSize >= m_count ? m_count - m_position : m_chunkSize;
		m_buffer.resize(count);
		m_instr->read((char*) m_buffer.data(), count * sizeof(geo::pc::Point));
		m_index = 0;
		m_position += count;
		return true;
	}

	/**
	 * Close and destroy the current stream.
	 */
	void unload() {
		if(m_instr) {
			m_instr->close();
			delete m_instr;
		}
	}

public:

	/**
	 * Create a sort object for point clouds.
	 *
	 * @param outputFile The final output file.
	 * @param inputFiles A list of input point cloud files.
	 * @param tmpDir A directory to store temporary files.
	 * @param scale A factor by which to scale normalized points. To generate indices,
	 *              coordinates are normalized to the bounding box, then scaled (divided)
	 *              and truncated.
	 * @param xmin The minimum x coordinate for normalizing points.
	 * @param ymin The minimum y coordinate for normalizing points.
	 * @param chunkSize The number of points in a chunk. (RAM)
	 * @param numChunks The number of chunks that can be processed simultaneously. (Processors)
	 */
	ExternalMergeSort(const std::string& outputFile, const std::vector<std::string>& inputFiles,
			const std::string tmpDir, const GridProps& props,
			int chunkSize = 4 * 1024 * 1024, int numChunks = 4) :
		m_outputFile(outputFile),
		m_inputFiles(inputFiles),
		m_tmpDir(tmpDir),
		m_chunkSize(chunkSize),
		m_numChunks(numChunks),
		m_chunk(0),

		m_instr(nullptr),
		m_count(0),
		m_position(0),
		m_index(0),

		m_props(props) {
	}

	/**
	 * Run over the files and figure out the minimum bounds.
	 * @param useHeader True if the header can be used to determine bounds. False to
	 *                  force reading through every file.
	 */
	/*
	void computeMinimums(bool useHeader = false) {
		m_xmin = G_DBL_MAX_POS;
		m_ymin = G_DBL_MAX_POS;
		double x, y;
		for(const std::string& f : m_inputFiles) {
			LASPointLoader ldr(f);
			ldr.computeBounds(x, y, useHeader);
			if(x < m_xmin) m_xmin = x;
			if(y < m_ymin) m_ymin = y;
		}
	}
	*/

	/**
	 * Return the Morton index of the Point. The coordinate will be
	 * normalized and scaled.
	 *
	 * @param pt A Point.
	 * @return The Morton index.
	 */
	size_t index(const geo::pc::Point& pt) const {
		return index(pt.x(), pt.y());
	}

	/**
	 * Return the Morton index of the coordinate. The coordinate will be
	 * normalized and scaled.
	 *
	 * @param x The x coordinate.
	 * @param y The y coordinate.
	 * @return The Morton index.
	 */
	size_t index(double x, double y) const {
		return toIdx(m_props.toCol(x), m_props.toRow(y), m_props.cols());
	}

	/**
	 * Return the scale value.
	 * @return The scale value.
	 */
	/*
	double xscale() const {
		return m_scale;
	}
	*/

	/**
	 * Restores the un-scaled and un-normalized coordinate that
	 * was used to create the Morton index. May lose some precision.
	 *
	 * @param index The index.
	 * @param x The original x-coordinate (out).
	 * @param y The original y-coordinate (out).
	 */
	/*
	void position(size_t index, double& x, double& y) const {
		size_t ax, ay;
		uninterleave(index, ax, ay);
		x = ax / m_scale;
		y = ay / m_scale;
	}
	*/

	/**
	 * Restores the scaled and normalized coordinate that
	 * was used to create the Morton index. Essentially a column
	 * and row.
	 *
	 * @param index The index.
	 * @param col The column (out).
	 * @param row The row (out).
	 */
	void position(size_t index, size_t& col, size_t& row) const {
		uninterleave(index, col, row);
	}

	/**
	 * Sort the point cloud.
	 *
	 * @param force If true, forces the entire process, ignoring
	 * any intermediate state.
	 */
	void sort(bool force) {
		if(!Util::exists(m_outputFile) || force) {
			phase1();
			phase2();
		}
	}

	/**
	 * Reset the point reader.
	 */
	void reset() {
		unload();
		m_instr = new std::ifstream(m_chunkFiles[0], std::ios::binary|std::ios::in);
		m_instr->read((char*) &m_count, sizeof(size_t));
		m_position = 0;
	}

	/**
	 * Populate the given point with the data from the next point. Return
	 * true on success, false otherwise (e.g. if there are no more points.)
	 *
	 * @param pt A point (will be modified.)
	 * @return True if there was a point available and it was successfuly copied.
	 */
	bool next(geo::pc::Point& pt) {
		if(m_index >= m_buffer.size() && !loadNext()) {
			return false;
		} else {
			std::memcpy(&pt, &m_buffer[m_index], sizeof(geo::pc::Point));
			++m_index;
			return true;
		}
	}

	~ExternalMergeSort() {
		unload();
	}
};



#endif /* INCLUDE_EXTERNALMERGESORT_HPP_ */
