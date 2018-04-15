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

using namespace geo::util;
using namespace geo::pc;


namespace geo {
namespace pc {
namespace sort {

/**
 * Holds offset and scale information for transforming
 * points to a sortable representation.
 */
class PointTransformer {
private:
	double m_x;
	double m_y;
	double m_xscale;
	double m_yscale;
public:
	PointTransformer() :
		PointTransformer(0, 0, 1, 1) {}
	PointTransformer(double x, double y, double xscale, double yscale) :
		m_x(x), m_y(y),
		m_xscale(xscale), m_yscale(yscale) {}
	void transform(const geo::pc::Point& pt, double& x, double& y) const {
		x = (pt.x() - m_x) / m_xscale;
		y = (pt.y() - m_y) / m_yscale;
	}
};


/**
 * Abstract class for loading points and returning them as
 * vectors.
 */
class PointLoader {
public:
	virtual int getPoints(std::vector<geo::pc::Point>& pts, int count) = 0;
	virtual void load(const std::string& filename) = 0;
	virtual ~PointLoader() {}
};


/**
 * Implementation of PointLoader to load LAS files.
 */
class LASPointLoader : public PointLoader {
private:

	liblas::Reader* m_rdr;
	std::ifstream* m_instr;

	double m_xscale;
	double m_yscale;
	double m_zscale;

public:

	LASPointLoader() :
		m_rdr(nullptr),
		m_instr(nullptr),
		m_xscale(1),
		m_yscale(1),
		m_zscale(1) {
	}

	void load(const std::string& filename) {
		// Set up the reader.
		liblas::ReaderFactory rf;
		m_instr = new std::ifstream(filename, std::ios::binary|std::ios::in);
		m_rdr = new liblas::Reader(rf.CreateWithStream(*m_instr));
		const liblas::Header& hdr = m_rdr->GetHeader();
		m_xscale = hdr.GetScaleX();
		m_yscale = hdr.GetScaleY();
		m_zscale = hdr.GetScaleZ();
	}

	double xscale() const {
		return m_xscale;
	}

	double yscale() const {
		return m_yscale;
	}

	double zscale() const {
		return m_zscale;
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

	/**
	 * Returns true if there are no (more) points in the list.
	 * @return True if there are no (more) points in the list.
	 */
	bool empty() {
		return m_ptr >= m_size;
	}

	/**
	 * Reset the pointer to the beginning of the list.
	 */
	void reset() {
		m_ptr = 0;
	}

	~PointStream() {
		if(m_instr) {
			m_instr->close();
			delete m_instr;
		}
	}

	/**
	 * Return the point currently at the front of the list.
	 * @return nullptr if there are no more points.
	 */
	geo::pc::Point* current() {
		if((m_pts.empty() || m_ptr >= m_pts.size())  && !loadNext())
			return nullptr;
		return &(m_pts[m_ptr]);
	}

	/**
	 * Remove the point at the front of the list.
	 */
	void pop() {
		++m_ptr;
	}

};


/**
 * Sort function for Points. Points will be transformed before sorting.
 *
 * @param a A Point.
 * @param b A Point.
 * @param trans A PointTransformer.
 * @return True if a is "less" than b.
 */
bool pointSort(const geo::pc::Point& a, const geo::pc::Point& b, PointTransformer* trans) {
	double ax, ay, bx, by;
	trans->transform(a, ax, ay);
	trans->transform(b, bx, by);
	return ay < by || (ay == by && ax < bx);
}

/**
 * Sorts the vector of points beginning and ending at the given indices.
 *
 * @param pts A vector of Points.
 * @param begin The start index.
 * @param end The end index.
 * @param trans a PointTransformer.
 */
void sortPoints(std::vector<geo::pc::Point>* pts, size_t begin, size_t end, PointTransformer* trans) {
	std::sort(pts->begin() + begin, pts->begin() + end,
			std::bind(pointSort, std::placeholders::_1, std::placeholders::_2, trans));
}



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

	PointLoader* m_loader;
	PointTransformer* m_trans;

	std::vector<geo::pc::Point> m_buffer;
	std::ifstream* m_instr;
	size_t m_count;
	size_t m_position;
	size_t m_index;

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

		size_t count = 0;
		PointStream* source; // The source that produced minPt.
		geo::pc::Point* minPt;
		std::vector<geo::pc::Point> buffer;
		buffer.reserve(m_chunkSize);

		while(true) {
			minPt = nullptr;
			for(PointStream& ps : inStreams) {
				geo::pc::Point* a = ps.current();
				if(a && (!minPt || pointSort(*a, *minPt, m_trans))) {
					minPt = a;
					source = &ps;
				}
			}
			if(minPt) {
				// A point was found. Remove it from the source and add to the buffer.
				source->pop();
				buffer.push_back(*minPt);
				if(buffer.size() == m_chunkSize) {
					// Write the points to the output.
					outStream.write((char*) buffer.data(), buffer.size() * sizeof(geo::pc::Point));
					count += buffer.size();
					buffer.clear();
				}
			} else {
				// All the sources are empty.
				break;
			}
		}

		if(!buffer.empty()) {
			// Write the points to the output.
			outStream.write((char*) buffer.data(), buffer.size() * sizeof(geo::pc::Point));
			count += buffer.size();
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
		int size;
		// Iterate over the list of source files.
		for(const std::string& inputFile : m_inputFiles) {
			m_loader->load(inputFile);
			while(true) {
				std::list<std::thread> threads;
				std::vector<geo::pc::Point> pts;
				pts.reserve(m_numChunks * m_chunkSize);
				if((size = m_loader->getPoints(pts, m_numChunks * m_chunkSize))) {
					for(int i = 0; i < (int) std::ceil((double) size / m_chunkSize); ++i) {
						size_t begin = i * m_chunkSize;
						size_t end = std::min(size, (i + 1) * m_chunkSize);
						threads.emplace_back(sortPoints, &pts, begin, end, m_trans);
					}
				} else {
					break;
				}
				for(std::thread& t : threads)
					t.join();
				for(int i = 0; i < (int) std::ceil((double) size / m_chunkSize); ++i) {
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
			const std::string tmpDir, int chunkSize = 1024 * 1024, int numChunks = 4) :
		m_outputFile(outputFile),
		m_inputFiles(inputFiles),
		m_tmpDir(tmpDir),
		m_chunkSize(chunkSize),
		m_numChunks(numChunks),
		m_chunk(0),
		m_loader(nullptr),
		m_trans(new PointTransformer()),
		m_instr(nullptr),
		m_count(0),
		m_position(0),
		m_index(0) {
	}

	/**
	 * Set the PointLoader. This class takes ownership and
	 * will destroy the loader when finished.
	 * @param loader A PointLoader.
	 */
	void setLoader(PointLoader* loader) {
		if(m_loader && m_loader != loader)
			delete m_loader;
		m_loader = loader;
	}

	void setTransformer(PointTransformer* trans) {
		if(m_trans && m_trans != trans)
			delete m_trans;
		m_trans = trans;
	}

	/**
	 * Sort the point cloud.
	 *
	 * @param force If true, forces the entire process, ignoring
	 * any intermediate state.
	 */
	void sort(bool force) {
		if(!Util::exists(m_outputFile) || force) {
			if(!m_loader)
				g_runerr("No point loader has been configured.")
			if(!m_trans)
				g_runerr("No point transformer has been configured.")
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
		if(m_loader)
			delete m_loader;
		if(m_trans)
			delete m_trans;
	}

};

} // sort
} // pc
} // geo

#endif /* INCLUDE_EXTERNALMERGESORT_HPP_ */
