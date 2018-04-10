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

size_t interleave(const geo::pc::Point& a) {
	int ax = (int) (a.x() / 100.0);
	int ay = (int) (a.y() / 100.0);
	size_t out = 0;
	for(size_t i = 0; i < sizeof(int) * 8; ++i) {
		out |= ((ax >> i) & 1) << (i * 2);
		out |= ((ay >> i) & 1) << (i * 2 + 1);
	}
	return out;
}

bool pointSort(const geo::pc::Point& a, const geo::pc::Point& b) {
	return interleave(a) < interleave(b);
}

void sortPoints(std::vector<geo::pc::Point>* pts, size_t begin, size_t end) {
	std::sort(pts->begin() + begin, pts->begin() + end, pointSort);
}

/**
 * Abstract class for loading points and returning them as
 * vectors.
 */
class PointLoader {
public:
	virtual int getPoints(std::vector<geo::pc::Point>& pts, int count) = 0;
	virtual ~PointLoader() {}
};

/**
 * Implementation of PointLoader to load LAS files.
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

	std::string nextChunkFile() {
		return Util::pathJoin(m_tmpDir, "chunk_" + std::to_string(++m_chunk) + ".tmp");
	}

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
				} else if(a && minPt && pointSort(*a, *minPt)) {
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
						threads.emplace_back(sortPoints, &pts, begin, end);
					}
				} else {
					break;
				}
				for(std::thread& t : threads)
					t.join();
				for(int i = 0; i < size / m_chunkSize; ++i) {
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

public:
	ExternalMergeSort(const std::string& outputFile, const std::vector<std::string>& inputFiles,
			const std::string tmpDir, int chunkSize, int numChunks) :
		m_outputFile(outputFile),
		m_inputFiles(inputFiles),
		m_tmpDir(tmpDir),
		m_chunkSize(chunkSize),
		m_numChunks(numChunks),
		m_chunk(0),

		m_instr(nullptr),
		m_count(0),
		m_position(0),
		m_index(0) {
	}

	void sort(bool force) {
		if(!Util::exists(m_outputFile) || force) {
			phase1();
			phase2();
		}
	}

	void reset() {
		unload();
		m_instr = new std::ifstream(m_chunkFiles[0], std::ios::binary|std::ios::in);
		m_instr->read((char*) &m_count, sizeof(size_t));
		m_position = 0;
	}

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

	bool next(geo::pc::Point& pt) {
		if(m_index >= m_buffer.size() && !loadNext()) {
			return false;
		} else {
			std::memcpy(&pt, &m_buffer[m_index], sizeof(geo::pc::Point));
			++m_index;
			return true;
		}
	}

	void unload() {
		if(m_instr) {
			m_instr->close();
			delete m_instr;
		}
	}

	~ExternalMergeSort() {
		unload();
	}
};



#endif /* INCLUDE_EXTERNALMERGESORT_HPP_ */
