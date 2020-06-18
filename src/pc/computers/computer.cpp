#include <vector>

#include "pointcloud.hpp"
#include "pc_computer.hpp"

using namespace geo::pc;
using namespace geo::pc::compute;


void Computer::setRasterizer(geo::pc::Rasterizer* rasterizer) {
	m_rasterizer = rasterizer;
}

Rasterizer* Computer::rasterizer() const {
	return m_rasterizer;
}

void Computer::setFilters(const std::vector<PointFilter*>& filters) {
	m_filters.assign(filters.begin(), filters.end());
}

const std::vector<PointFilter*> Computer::filters() const {
	return m_filters;
}

size_t Computer::filter(const std::vector<geo::pc::Point>& pts, std::vector<geo::pc::Point>& out) const {
	size_t count = 0;
	if(m_filters.empty()) {
		count = pts.size();
		out.assign(pts.begin(), pts.end());
	} else {
		bool keep;
		for(const geo::pc::Point& pt : pts) {
			keep = true;
			for(const PointFilter* f : m_filters) {
				if(!f->keep(pt)) {
					keep = false;
					break;
				}
			}
			if(keep) {
				out.push_back(pt);
				++count;
			}
		}
	}
	return count;
}

int DensityComputer::compute(double, double, const std::vector<geo::pc::Point>& pts, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {
	if(!filtered.empty()) {
		double area = radius * radius * M_PI;
		out.push_back(pts.size() / area);
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int DensityComputer::bandCount() const {
	return 1;
}

std::vector<std::pair<std::string,std::string>> DensityComputer::bandMeta() const {
	return {{"none", "None"}};
}

