/*
 * pc_normalizer.cpp
 *
 *  Created on: Mar 17, 2018
 *      Author: rob
 */

#include <vector>
#include <string>

#include "raster.hpp"
#include "pointcloud.hpp"

using namespace geo::raster;
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

	Raster dtm(dtmpath);
	const GridProps& props = dtm.props();
	double nodata = props.nodata();

	liblas::ReaderFactory fact;

	for(const std::string& filename : m_filenames) {

		std::string outfile = Util::pathJoin(outdir, Util::basename(filename) + ".las");
		if(!force && Util::exists(outfile)) {
			g_trace("File " << outfile << " exists. Skipping.")
			continue;
		}

		std::ifstream str(filename);
		liblas::Reader rdr = fact.CreateWithStream(str);

		std::ofstream ostr(outfile, std::ios::binary | std::ios::trunc | std::ios::out);
		liblas::Header hdr(rdr.GetHeader());
		liblas::Writer wtr(ostr, hdr);

		double minZ = DBL_MAX;
		double maxZ = DBL_MIN;
		double minX = DBL_MAX;
		double maxX = DBL_MIN;
		double minY = DBL_MAX;
		double maxY = DBL_MIN;

		while(rdr.ReadNextPoint()) {

			const liblas::Point& pt = rdr.GetPoint();

			if(m_filter && !m_filter->keep(pt))
				continue;

			double x = pt.GetX();
			double y = pt.GetY();

			int col = props.toCol(x);
			int row = props.toRow(y);

			if(col < 0 || col >= props.cols() || row < 0 || row >= props.rows())
				continue;

			double t = dtm.getFloat(props.toCol(x), props.toRow(y), band);

			if(t == nodata)
				continue;

			double z = pt.GetZ();
			double z0 = z - t;

			if(z0 < 0) z0 = 0;

			if(z0 < minZ) minZ = z0;
			if(z0 > maxZ) maxZ = z0;
			if(x < minX) minX = x;
			if(x > maxX) maxX = x;
			if(y < minY) minY = y;
			if(y > maxY) maxY = y;

			liblas::Point npt(pt);
			npt.SetZ(z0);
			wtr.WritePoint(npt);
		}

		hdr.SetMin(minX, minY, minZ);
		hdr.SetMax(maxX, maxY, maxZ);
		wtr.SetHeader(hdr);

	}

}

