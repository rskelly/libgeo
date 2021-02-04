/*
 * pc_normalizer.cpp
 *
 *  Created on: Mar 17, 2018
 *      Author: rob
 */

#include <vector>
#include <string>

#include "pointcloud.hpp"
#include "grid.hpp"

using namespace geo::grid;
using namespace geo::pc;

Normalizer::Normalizer(const std::vector<std::string> filenames) :
		m_filenames(filenames),
		m_filter(nullptr) {
}

void Normalizer::setFilter(const PCPointFilter& filter) {
	m_filter = new PCPointFilter(filter);
}

Normalizer::~Normalizer() {
	if(m_filter)
		delete m_filter;
}

void Normalizer::normalize(const std::string& dtmpath, const std::string& outdir, int band, bool force) {

	Band<float> dtm(dtmpath, band - 1, false, true);
	const GridProps& props = dtm.props();
	double nodata = props.nodata();

	for(const std::string& filename : m_filenames) {

		std::string outfile = join(outdir, basename(filename) + ".las");
		if(!force && isfile(outfile)) {
			g_trace("File " << outfile << " exists. Skipping.")
			continue;
		}

		geo::pc::PCFile rdr(filename);

		geo::pc::PCWriter wtr(outfile, filename);

		double minZ = geo::maxvalue<double>();
		double maxZ = -geo::maxvalue<double>();
		double minX = geo::maxvalue<double>();
		double maxX = -geo::maxvalue<double>();
		double minY = geo::maxvalue<double>();
		double maxY = -geo::maxvalue<double>();

		geo::pc::Point pt;
		while(rdr.next(pt)) {

			if(m_filter && !m_filter->keep(pt))
				continue;

			double x = pt.x();
			double y = pt.y();

			if(!props.hasCell(x, y)) continue;

			double t = dtm.get(x, y);

			if(t == nodata)
				continue;

			double z = pt.x();
			double z0 = z - t;

			if(z0 < 0) z0 = 0;

			if(z0 < minZ) minZ = z0;
			if(z0 > maxZ) maxZ = z0;
			if(x < minX) minX = x;
			if(x > maxX) maxX = x;
			if(y < minY) minY = y;
			if(y > maxY) maxY = y;

			pt.z(z0);
			wtr.addPoint(pt);
		}

	}

}

