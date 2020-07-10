/*
 * pc.cpp
 *
 *  Created on: Jul 9, 2020
 *      Author: rob
 */

#include <pdal/Options.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/SpatialReference.hpp>
#include <pdal/io/LasReader.hpp>
#include <pdal/io/LasHeader.hpp>

#include "geo.hpp"
#include "pc.hpp"
#include "util.hpp"

using namespace geo::util;
using namespace geo::pc;

namespace {

	// TODO: Don't load the entire cloud into memory.
	class LasFile : public PointFile {
	private:
		size_t m_idx;
		pdal::LasReader* m_reader;
		pdal::PointTable* m_table;
		pdal::PointViewSet* m_viewSet;
		pdal::PointViewPtr m_view;
		pdal::LasHeader* m_header;

	protected:

		void scan() {
			geo::pc::Point pt;
			reset();
			m_bounds.collapse(3);
			while(next(pt))
				m_bounds.extend(pt.x(), pt.y(), pt.z());
			reset();
		}

	public:
		LasFile(const std::string& file) :
			m_idx(0) {
			pdal::Option opt("filename", file);
			pdal::Options opts;
			opts.add(opt);
			m_table = new pdal::PointTable;
			m_reader = new pdal::LasReader;
			m_reader->setOptions(opts);
			m_reader->prepare(*m_table);
			m_viewSet = m_reader->execute(*m_table);
			m_view = m_viewSet->begin();
			m_header = m_reader->header();
			m_bounds.collapse(3);
			m_bounds.extend(m_header->minX(), m_header->minY(), m_header->minZ());
			m_bounds.extend(m_header->maxX(), m_header->maxY(), m_header->maxZ());
		}

		void reset() {
			m_idx = 0;
		}

		size_t size() const {
			return m_header->pointCount();
		}

		const geo::util::Bounds& bounds() const {
			return m_bounds;
		}

		std::string projection() const {
			return m_table->spatialReference().getWKT();
		}

		bool next(geo::pc::Point& pt) {
			if(m_idx < m_view->size()) {
				pt.m_x = m_view->getFieldAs<double>(pdal::Dimension::Id::X, m_idx);
				pt.m_y = m_view->getFieldAs<double>(pdal::Dimension::Id::Y, m_idx);
				pt.m_z = m_view->getFieldAs<double>(pdal::Dimension::Id::Z, m_idx);
				pt.m_cls = m_view->getFieldAs<int>(pdal::Dimension::Id::Classification, m_idx);
				++m_idx;
				return true;
			}
			return false;
		}

		~LasFile() {
			delete m_viewSet;
			delete m_reader;
		}
};

} // Anon.

Reader::Reader() :
		m_size(0),
		m_fileIdx(0),
		m_trustHeader(true) {}

Reader::Reader(const std::string& file, bool trustHeader) :
		Reader() {
	init({file}, trustHeader);
}

Reader::Reader(const std::vector<std::string>& files, bool trustHeader) :
		Reader() {
	init(files, trustHeader);
}

// TODO: Implement more file formats.
// TODO: Implement full scan (don't trust headers).
void Reader::init(const std::vector<std::string>& files, bool trustHeader) {
	m_bounds.collapse(3);
	m_size = 0;
	m_files = files;
	bool first = true;
	for(const std::string& file : files) {
		std::unique_ptr<PointFile> pf(new LasFile(file));
		// If the header isn't to be trusted, scan the point cloud.
		if(!trustHeader)
			pf->scan();
		m_bounds.extend(pf->bounds());
		m_size += pf->size();
		if(first) {
			// On the first one, get the projection and preserve
			// the file for first iteration via next().
			m_projection = pf->projection();
			m_currentFile.swap(pf);
			m_currentFile->reset();
		}
	}
}

size_t Reader::size() const {
	return m_size;
}

const geo::util::Bounds& Reader::bounds() const {
	return m_bounds;
}

const std::string& Reader::projection() const {
	return m_projection;
}

bool Reader::next(geo::pc::Point& pt) {
	while(true) {
		// Try to get the next point.
		if(m_currentFile->next(pt))
			return true;

		// If there's no point and there's no other point, we're done.
		if(m_fileIdx >= m_files.size())
			return false;

		// Load the next file and go around.
		m_currentFile.reset(new LasFile(m_files[m_fileIdx]));
		++m_fileIdx;
	}
}

void Reader::reset() {
	if(m_fileIdx > 0)
		m_currentFile.reset(new LasFile(m_files[0]));
	m_currentFile->reset();
}
