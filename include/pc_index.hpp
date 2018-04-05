/*
 * pc_index.hpp
 *
 * Any point in the point cloud may be the last point in
 * one or more cells. This index works by triggering finalization
 * calls for each point that is the last point for one or more cells.
 *
 *  Created on: Apr 4, 2018
 *      Author: rob
 */

#ifndef INCLUDE_PC_INDEX_HPP_
#define INCLUDE_PC_INDEX_HPP_

#include <unordered_map>
#include <list>
#include <vector>
#include <fstream>

#include "liblas/liblas.hpp"

#include "geo.hpp"
#include "pointcloud.hpp"

class PointCloud {
private:
	// The key is the index of the point in the file.
	// The value is a the cell index to finalize after that point.
	// The actual number of finalized cells must be determined by
	// the caller using the computation radius.
	// TODO: THIS IS WRONG. Must be determined here, using radius.
	std::unordered_map<size_t, size_t> m_indices;	///< The file->cell indices for finalization.
	std::vector<std::string> m_files;				///< The list of point files.

	liblas::Reader* m_reader;						///< The current point file reader.
	std::ifstream* m_input;							///< The current input file stream.
	size_t m_readIdx;								///< The index of the current point.
	size_t m_fileIdx;								///< The index of the current file.
	size_t m_final;									///< The index of the last finalized cell.

	size_t m_rows;									///< The number of cells in the finalization grid.
	size_t m_cols;									///< The number of rows in the finalization grid.
	double m_height;								///< The height of the point cloud (y).
	double m_width;									///< The width of the point cloud (x).
	double m_cellSize;								///< The size of a finalization cell.
	double m_bounds[4]; // minx, miny, maxx, maxy	///< The bounds of the point cloud.

	/**
	 * Delete and close the reader and file.
	 */
	void unload() {
		if(m_reader) {
			delete m_reader;
			m_reader = nullptr;
		}
		if(m_input) {
			m_input->close();
			delete m_input;
			m_input = nullptr;
		}
	}

	/**
	 * Unload the reader and load the next one. Return
	 * true on success and if there's a point waiting to be read.
	 */
	bool loadNext() {
		unload();
		liblas::ReaderFactory rf;
		m_input = new std::ifstream(m_files[++m_fileIdx], std::ios::binary|std::ios::in);
		m_reader = new liblas::Reader(rf.CreateWithStream(*m_input));
		return m_reader->ReadNextPoint();
	}

public:

	PointCloud(std::vector<std::string>& files, double cellSize = 10.0) :
		m_files(files),
		m_reader(nullptr),
		m_input(nullptr),
		m_readIdx(0), m_fileIdx(0),
		m_final(0),
		m_rows(0), m_cols(0),
		m_height(0), m_width(0),
		m_cellSize(cellSize) {
	}

	size_t cols() const {
		return m_cols;
	}

	size_t rows() const {
		return m_rows;
	}

	void init(bool useHeader = true) {

		liblas::ReaderFactory rf;

		m_bounds = {G_DBL_MAX_POS, G_DBL_MAX_POS, G_DBL_MAX_NEG, G_DBL_MAX_NEG};

		for(const std::string& filename : m_files) {
			std::ifstream str(filename, std::ios::binary|std::ios::in);
			liblas::Reader rdr = rf.CreateWithStream(str);
			if(useHeader) {
				const liblas::Header& hdr = rdr.GetHeader();
				m_bounds[0] = std::min(hdr.GetMinX(), m_bounds[0]);
				m_bounds[1] = std::min(hdr.GetMinY(), m_bounds[1]);
				m_bounds[2] = std::max(hdr.GetMaxX(), m_bounds[2]);
				m_bounds[3] = std::max(hdr.GetMaxX(), m_bounds[3]);
			} else {
				while(rdr.ReadNextPoint()) {
					const liblas::Point& pt = rdr.GetPoint();
					m_bounds[0] = std::min(pt.GetX(), m_bounds[0]);
					m_bounds[1] = std::min(pt.GetY(), m_bounds[1]);
					m_bounds[2] = std::max(pt.GetX(), m_bounds[2]);
					m_bounds[3] = std::max(pt.GetX(), m_bounds[3]);
				}
			}
		}

		m_width = m_bounds[2] - m_bounds[0];
		m_height = m_bounds[3] - m_bounds[1];
		m_cols = (size_t) std::ceil(m_width / m_cellSize);
		m_rows = (size_t) std::ceil(m_height / m_cellSize);

		size_t lidx = 0;
		std::unordered_map<size_t, size_t> indices;
		for(const std::string& filename : m_files) {
			std::ifstream str(filename, std::ios::binary|std::ios::in);
			liblas::Reader rdr = rf.CreateWithStream(str);
			while(rdr.ReadNextPoint()) {
				const liblas::Point& pt = rdr.GetPoint();
				double px = pt.GetX();
				double py = pt.GetY();
				size_t col = (size_t) (px - m_bounds[0]) / m_width;
				size_t row = (size_t) (py - m_bounds[1]) / m_height;
				size_t idx = (row << 32) | col;
				indices[idx] = lidx;
				++lidx;
			}
		}

		// Flip the index map; we have cell idx -> file idx.
		// We want file idx -> cell idx.
		for(const auto& p : indices)
			m_indices[p.second] = p.first;

	}

	void reset() {
		m_readIdx = 0;
		m_fileIdx = 0;
		unload();
	}

	/**
	 * Populate the point with the next point's data and return true. If there is no
	 * point return false; the point's state is indeterminate.
	 * If the hasFinal value is true, a cell has been finalized with the reading of
	 * this point. Use getFinal() to get its index.
	 */
	bool next(geo::pc::Point& pt, bool& hasFinal) {
		if(!m_reader && !m_reader->ReadNextPoint() && !loadNext())
			return false;
		if(m_indices.find(m_readIdx) != m_indices.end()) {
			m_final = m_indices[m_readIdx];
			hasFinal = true;
		}
		pt.setPoint(m_reader->GetPoint());
		++m_readIdx;
		return true;
	}

	/**
	 * If the hasFinal value from next() is true, use getFinal
	 * to return the index of the cell that is finalized.
	 */
	size_t getFinal() const {
		return m_final;
	}

	~PointCloud() {
		unload();
	}

};



#endif /* INCLUDE_PC_INDEX_HPP_ */
