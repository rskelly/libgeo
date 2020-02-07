/*
 * pointcloud_reader.hpp
 *
 *  Created on: Feb 5, 2020
 *      Author: rob
 */

#ifndef LIBGEO_INCLUDE_POINTCLOUD_READER_HPP_
#define LIBGEO_INCLUDE_POINTCLOUD_READER_HPP_

#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/io/LasReader.hpp>
#include <pdal/io/BufferReader.hpp>
#include <pdal/io/LasWriter.hpp>
#include <pdal/io/LasHeader.hpp>
#include <pdal/Options.hpp>

#include "geo.hpp"

size_t __idx;

class pointcloud_reader {
private:
	std::string m_projection;
	double m_bounds[6];

	pdal::PipelineManager m_mgr;
	pdal::PointTable* m_table;
	pdal::Stage* m_rdr;
	pdal::LasHeader m_hdr;
	pdal::PointViewSet m_viewSet;
	pdal::PointViewPtr m_view;
	pdal::Dimension::IdList m_dims;

	pdal::PointId m_idx;
	bool m_hasBounds;
	bool m_open;
	size_t m_size;

public:

	pointcloud_reader() :
		m_rdr(nullptr),
		m_table(nullptr),
		m_idx(0),
		m_hasBounds(false),
		m_open(false), m_size(0) {

	}

	size_t size() const {
		return m_size;
	}

	void open(const std::string& filename, bool computeBounds) {

		close();

		using namespace pdal::Dimension;

		m_rdr = &m_mgr.makeReader(filename, "readers.las");
		m_table = new pdal::PointTable();
		m_rdr->prepare(*m_table);
		m_viewSet = m_rdr->execute(*m_table);
		m_view = *m_viewSet.begin();
		m_dims = m_view->dims();

		m_projection = m_hdr.srs().getWKT();

		m_idx = 0;
		m_size = m_view->size();
		m_hasBounds = false;

		if(computeBounds) {
			for(int i = 0; i < 3; ++i) {
				m_bounds[i] = G_DBL_MAX_POS;
				m_bounds[i + 3] = G_DBL_MAX_NEG;
			}
			for (pdal::PointId idx = 0; idx < m_view->size(); ++idx) {
				double x = m_view->getFieldAs<double>(Id::X, idx);
				double y = m_view->getFieldAs<double>(Id::Y, idx);
				double z = m_view->getFieldAs<double>(Id::Z, idx);
				if(x < m_bounds[0]) m_bounds[0] = x;
				if(y < m_bounds[1]) m_bounds[1] = y;
				if(z < m_bounds[2]) m_bounds[2] = z;
				if(x > m_bounds[3]) m_bounds[3] = x;
				if(y > m_bounds[4]) m_bounds[4] = y;
				if(z > m_bounds[5]) m_bounds[5] = z;
			}
			m_hasBounds = true;
		}

		m_open = true;
	}

	void close() {
		if(m_open) {
			m_mgr.destroyStage(m_rdr);
			m_viewSet.clear();
			delete m_table;
			m_table = nullptr;
			m_hasBounds = false;
			m_open = false;
		}
	}

	bool bounds(double* bounds) const {
		if(!m_hasBounds)
			return false;
		for(int i = 0; i < 6; ++i)
			bounds[i] = m_bounds[i];
		return true;
	}

	bool extendBounds(double* bounds) const {
		if(!m_hasBounds)
			return false;
		for(int i = 0; i < 3; ++i) {
			if(m_bounds[i] < bounds[i]) bounds[i] = m_bounds[i];
			if(m_bounds[i + 3] > bounds[i + 3]) bounds[i + 3] = m_bounds[i + 3];
		}
		return true;
	}

	bool next(double& x, double& y, double& z, int& cls) {
		return next(x, y, z, cls, __idx);
	}

	bool next(double& x, double& y, double& z, int& cls, size_t& idx) {
		if(!m_open || m_idx >= m_view->size())
			return false;

		using namespace pdal::Dimension;

		x = m_view->getFieldAs<double>(Id::X, m_idx);
		y = m_view->getFieldAs<double>(Id::Y, m_idx);
		z = m_view->getFieldAs<double>(pdal::Dimension::Id::Z, m_idx);
		if(z > 1000)
			std::cerr << z;
		cls = m_view->getFieldAs<int>(Id::Classification, m_idx);
		idx = m_idx;
		++m_idx;
		return true;
	}

	bool next(double& x, double& y, double& z) {
		return next(x, y, z, __idx);
	}
	bool next(double& x, double& y, double& z, size_t& idx) {
		if(!m_open || m_idx >= m_view->size())
			return false;

		using namespace pdal::Dimension;

		x = m_view->getFieldAs<double>(Id::X, m_idx);
		y = m_view->getFieldAs<double>(Id::Y, m_idx);
		z = m_view->getFieldAs<double>(pdal::Dimension::Id::Z, m_idx);
		idx = m_idx;
		++m_idx;
		return true;
	}

	template <class V>
	bool update(pdal::Dimension::Id dimension, size_t idx, V value) {
		if(!m_open || idx >= m_view->size())
			return false;
		m_view->setField(dimension, idx, value);
		return true;
	}

	bool setX(size_t idx, double x) {
		return update(pdal::Dimension::Id::X, idx, x);
	}

	bool setY(size_t idx, double y) {
		return update(pdal::Dimension::Id::Y, idx, y);
	}

	bool setZ(size_t idx, double z) {
		return update(pdal::Dimension::Id::Z, idx, z);
	}

	bool setX(double x) {
		return update(pdal::Dimension::Id::X, m_idx, x);
	}

	bool setY(double y) {
		return update(pdal::Dimension::Id::Y, m_idx, y);
	}

	bool setZ(double z) {
		return update(pdal::Dimension::Id::Z, m_idx, z);
	}

	double x(size_t idx) {
		if(!m_open || idx >= m_view->size())
			return std::nan("");
		return m_view->getFieldAs<double>(pdal::Dimension::Id::X, idx);
	}

	double y(size_t idx) {
		if(!m_open || idx >= m_view->size())
			return std::nan("");
		return m_view->getFieldAs<double>(pdal::Dimension::Id::Y, idx);
	}

	double z(size_t idx) {
		if(!m_open || idx >= m_view->size())
			return std::nan("");
		return m_view->getFieldAs<double>(pdal::Dimension::Id::Z, idx);
	}

	double x() {
		return x(m_idx);
	}

	double y() {
		return y(m_idx);
	}

	double z() {
		return z(m_idx);
	}

	bool save(const std::string& outfile) {
		if(!m_open)
			return false;

		pdal::PointTable table;
		for(pdal::Dimension::Id dim : m_dims)
			table.layout()->registerDim(dim);

		pdal::PointViewPtr view(new pdal::PointView(*m_table));
		view->append(*m_view);

		pdal::BufferReader reader;
		reader.addView(view);

		pdal::Stage& wtr = m_mgr.makeWriter(outfile, "writers.las");
		wtr.setInput(reader);
		wtr.prepare(*m_table);
		wtr.execute(*m_table);

		m_mgr.destroyStage(&wtr);

		return true;
	}

	~pointcloud_reader() {
		close();
	}
};


#endif /* LIBGEO_INCLUDE_POINTCLOUD_READER_HPP_ */
