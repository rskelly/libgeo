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
	// The key is the index of the point in the file. The values are a the cell
	// indices to finalize after that point.The actual number of finalized cells 
	// must be determined by by calculating a set of cells that cover the target
	// cell plus any given radius.
	std::unordered_map<size_t, std::vector<size_t> > m_indices;	///< The file->cell indices for finalization.
	std::vector<std::string> m_files;				///< The list of point files.

	liblas::Reader* m_reader;						///< The current point file reader.
	std::ifstream* m_input;							///< The current input file stream.
	size_t m_readIdx;								///< The index of the current point.
	size_t m_fileIdx;								///< The index of the current file.

	int m_rows;									///< The number of cells in the finalization grid.
	int m_cols;									///< The number of rows in the finalization grid.
	double m_height;								///< The height of the point cloud (y).
	double m_width;									///< The width of the point cloud (x).
	double m_resX;									///< The size of a finalization cell.
	double m_resY;
	double m_bounds[4]; // minx, miny, maxx, maxy	///< The bounds of the point cloud.
	int m_offset;

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
		m_input = new std::ifstream(m_files[m_fileIdx], std::ios::binary|std::ios::in);
		m_reader = new liblas::Reader(rf.CreateWithStream(*m_input));
		++m_fileIdx;
		return m_reader->ReadNextPoint();
	}

	void fixBounds(double resX, double resY, double* easting, double* northing) {

		double aresX = std::abs(resX);
		double aresY = std::abs(resY);

		{
			int a = 0, b = 2;
			if(resX < 0)
				a = 2, b = 0;
			m_bounds[a] = std::floor(m_bounds[a] / resX) * resX;
			m_bounds[b] = std::ceil(m_bounds[b] / resX) * resX;
			a = 1, b = 3;
			if(resY < 0)
				a = 3, b = 1;
			m_bounds[a] = std::floor(m_bounds[a] / resY) * resY;
			m_bounds[b] = std::ceil(m_bounds[b] / resY) * resY;
		}

		if(!std::isnan(*easting)) {
			if((resX > 0 && *easting < m_bounds[0]) || (resX < 0 && *easting > m_bounds[2]))
				g_argerr("The easting is within the data boundary.");
			double w = m_bounds[2] - m_bounds[0];
			if(resX > 0) {
				while(*easting + w < m_bounds[2])
					w += resX;
				m_bounds[0] = *easting;
				m_bounds[2] = *easting + w;
			} else {
				while(*easting - w > m_bounds[1])
					w += aresX;
				m_bounds[2] = *easting;
				m_bounds[0] = *easting - w;
			}
		} else {
			*easting = m_bounds[resX > 0 ? 0 : 2];
		}

		if(!std::isnan(*northing)) {
			if((resY > 0 && *northing < m_bounds[1]) || (resY < 0 && *northing > m_bounds[3]))
				g_argerr("The *northing is within the data boundary.");
			double h = m_bounds[3] - m_bounds[1];
			if(resY > 0) {
				while(*northing + h < m_bounds[3])
					h += resY;
				m_bounds[1] = *northing;
				m_bounds[3] = *northing + h;
			} else {
				while(*northing - h > m_bounds[1])
					h += aresY;
				m_bounds[3] = *northing;
				m_bounds[1] = *northing - h;
			}
		} else {
			*northing = m_bounds[resY > 0 ? 1 : 3];
		}
	}

public:

	PointCloud(const std::vector<std::string>& files) :
		m_files(files),
		m_reader(nullptr),
		m_input(nullptr),
		m_readIdx(0), m_fileIdx(0),
		m_rows(0), m_cols(0),
		m_height(0), m_width(0),
		m_resX(0), m_resY(0),
		m_offset(0) {
	}

	int cols() const {
		return m_cols;
	}

	int rows() const {
		return m_rows;
	}

	void bounds(double* b) const {
		for(int i = 0; i < 4; ++i)
			b[i] = m_bounds[i];
	}

	void init(double radius, double resX, double resY, double* easting, double* northing, bool useHeader = true) {

		liblas::ReaderFactory rf;

		m_resX = resX;
		m_resY = resY;
		m_bounds[0] = m_bounds[1] = geo::maxvalue<double>();
		m_bounds[2] = m_bounds[3] = geo::minvalue<double>();

		for(const std::string& filename : m_files) {
			std::ifstream str(filename, std::ios::binary|std::ios::in);
			liblas::Reader rdr = rf.CreateWithStream(str);
			if(useHeader) {
				const liblas::Header& hdr = rdr.GetHeader();
				m_bounds[0] = std::min(hdr.GetMinX(), m_bounds[0]);
				m_bounds[1] = std::min(hdr.GetMinY(), m_bounds[1]);
				m_bounds[2] = std::max(hdr.GetMaxX(), m_bounds[2]);
				m_bounds[3] = std::max(hdr.GetMaxY(), m_bounds[3]);
			} else {
				while(rdr.ReadNextPoint()) {
					const liblas::Point& pt = rdr.GetPoint();
					m_bounds[0] = std::min(pt.GetX(), m_bounds[0]);
					m_bounds[1] = std::min(pt.GetY(), m_bounds[1]);
					m_bounds[2] = std::max(pt.GetX(), m_bounds[2]);
					m_bounds[3] = std::max(pt.GetY(), m_bounds[3]);
				}
			}
		}

		fixBounds(resX, resY, easting, northing);

		m_width = m_bounds[2] - m_bounds[0];
		m_height = m_bounds[3] - m_bounds[1];
		m_cols = (int) std::ceil(m_width / std::abs(m_resX));
		m_rows = (int) std::ceil(m_height / std::abs(m_resY));

		m_offset = (int) std::ceil(radius / std::max(std::abs(m_resX), std::abs(m_resY)));

		size_t lidx = 0;
		std::unordered_map<size_t, size_t> indices;
		for(const std::string& filename : m_files) {
			std::ifstream str(filename, std::ios::binary|std::ios::in);
			liblas::Reader rdr = rf.CreateWithStream(str);
			while(rdr.ReadNextPoint()) {
				const liblas::Point& pt = rdr.GetPoint();
				int col = (int) ((pt.GetX() - m_bounds[0]) / m_width * m_cols);
				int row = (int) ((pt.GetY() - m_bounds[1]) / m_height * m_rows);
				for(int r = std::max(0, row - m_offset); r < std::min((int) m_rows, row + m_offset + 1); ++r) {
					for(int c = std::max(0, col - m_offset); c < std::min((int) m_cols, col + m_offset + 1); ++c)
						indices[toIndex(col, row)] = lidx;
				}
				++lidx;
			}
		}

		// Flip the index map; we have cell idx -> file idx.
		// We want file idx -> cell idx. The mapping is not 1:1 so
		// we use a vector.
		for(const auto& p : indices)
			m_indices[p.second].push_back(p.first);

	}

	void fromIndex(size_t idx, int& col, int& row) {
		col = idx & 0xffffffff;
		row = (idx >> 32) & 0xffffffff;
	}

	void fromIndex(size_t idx, double& x, double& y) {
		int col = idx & 0xffffffff;
		int row = (idx >> 32) & 0xffffffff;
		x = m_bounds[0] + m_resX * col + m_resX * 0.5;
		y = m_bounds[1] + m_resY * row + m_resY * 0.5;
	}

	size_t toIndex(int col, int row) const {
		return (((size_t) row) << 32)|col;
	}

	size_t toIndex(const geo::pc::Point& pt) const {
		return  toIndex(toCol(pt), toRow(pt));
	}

	int toCol(const geo::pc::Point& pt) const {
		return (int) (pt.x() - m_bounds[0]) / m_width;
	}

	int toRow(const geo::pc::Point& pt) const {
		return (int) (pt.y() - m_bounds[1]) / m_height;
	}

	int getIndices(const geo::pc::Point& pt, std::vector<size_t>& indices) {
		double px = pt.x();
		double py = pt.y();
		int col = (int) ((px - m_bounds[0]) / m_width * m_cols);
		int row = (int) ((py - m_bounds[1]) / m_height * m_rows);
		size_t count = 0;
		for(int r = std::max(0, row - m_offset); r < std::min((int) m_rows, row + m_offset + 1); ++r) {
			for(int c = std::max(0, col - m_offset); c < std::min((int) m_cols, col + m_offset + 1); ++c) {
				size_t idx = toIndex(c, r);
				indices.push_back(idx);
				++count;
			}
		}
		return count;
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
	bool next(geo::pc::Point& pt, std::vector<size_t>& final) {
		if((!m_reader || !m_reader->ReadNextPoint()) && !loadNext())
			return false;
		final.clear();
		if(m_indices.find(m_readIdx) != m_indices.end())
			final.assign(m_indices[m_readIdx].begin(), m_indices[m_readIdx].end());
		pt.setPoint(m_reader->GetPoint());
		++m_readIdx;
		return true;
	}

	~PointCloud() {
		unload();
	}

};



#endif /* INCLUDE_PC_INDEX_HPP_ */
