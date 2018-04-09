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

class PointLoader {
private:
	liblas::Reader* m_rdr;
	std::ifstream* m_instr;

public:
	PointLoader(const std::string& filename) {
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

	~PointLoader() {
		delete m_rdr;
		delete m_instr;
	}
};

class PointStream {
private:
	std::list<geo::pc::Point> m_pts;
	std::vector<geo::pc::Point> m_tmp;
	std::string m_filename;
	std::ifstream* m_instr;
	size_t m_size;
	size_t m_current;
	size_t m_chunkSize;

	bool loadNext() {
		if(!m_instr) {
			m_instr = new std::ifstream(m_filename, std::ios::binary|std::ios::in);
			m_instr->read((char*) &m_size, sizeof(size_t));
			m_current += sizeof(size_t);
		}
		if(m_current >= m_size)
			return false;
		size_t size = m_current + m_chunkSize > m_size ? m_size - m_current : m_chunkSize;
		m_tmp.resize(size);
		m_instr->read((char*) m_tmp.data(), size * sizeof(geo::pc::Point));
		m_pts.assign(m_tmp.begin(), m_tmp.end());
		m_current += size;
		return true;
	}

public:

	PointStream(const std::string& filename) :
		m_filename(filename),
		m_instr(nullptr),
		m_size(0),
		m_current(0),
		m_chunkSize(1024 * 1024) {
	}

	~PointStream() {
		delete m_instr;
	}

	geo::pc::Point* current() {
		if(m_pts.empty() && !loadNext())
			return nullptr;
		return &(m_pts.front());
	}

	void pop() {
		if(!m_pts.empty())
			m_pts.pop_front();
	}

};

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
		std::vector<PointStream> inStreams;
		for(const std::string& file : chunkFiles)
			inStreams.emplace_back(file);

		size_t minIdx;
		size_t count = 0;
		geo::pc::Point* minPt;
		std::vector<geo::pc::Point> buffer;
		buffer.reserve(m_chunkSize);
		while(true) {
			minPt = nullptr;
			for(size_t i = 0; i < inStreams.size(); ++i) {
				geo::pc::Point* a = inStreams[i].current();
				if(a && !minPt) {
					minPt = a;
					inStreams[i].pop();
				} else if(a && minPt && pointSort(*a, *minPt)) {
					minPt = a;
					inStreams[i].pop();
				}
			}
			if(minPt) {
				// A point was found at the index. Remove it.
				buffer.push_back(*minPt);
				if(buffer.size() == m_chunkSize) {
					// Write the points to the output.
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
		size_t size;
		std::vector<geo::pc::Point> pts;
		// Iterate over the list of source files.
		for(const std::string& inputFile : m_inputFiles) {
			PointLoader ldr(inputFile);
			// Load a section of the file into the vector.
			while((size = ldr.getPoints(pts, m_chunkSize))) {
				// Sort the points.
				//std::sort(pts.begin(), pts.end(), pointSort);
				// Create and open an output stream.
				std::string chunkFile = nextChunkFile();
				std::ofstream ostr(chunkFile, std::ios::binary|std::ios::out);
				// Write the points to the output stream.
				ostr.write((char*) &size, sizeof(size_t));
				ostr.write((char*) pts.data(), pts.size() * sizeof(geo::pc::Point));
				// Save the chunk file name.
				m_chunkFiles.push_back(chunkFile);
				pts.clear();
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
