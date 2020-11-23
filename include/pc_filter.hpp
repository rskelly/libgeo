/*
 * pc_filter.hpp
 *
 *  Created on: Feb 10, 2020
 *      Author: rob
 */

#ifndef LIBGEO_INCLUDE_PC_FILTER_HPP_
#define LIBGEO_INCLUDE_PC_FILTER_HPP_

#include "geo.hpp"

namespace geo {
namespace pc {
namespace filter {

	class PointClassFilter : public PointFilter {
	public:
		std::vector<int> classes;
		PointClassFilter(int cls) {
			classes.push_back(cls);
		}
		PointClassFilter(const std::vector<int>& classes) :
			classes(classes) {
		}
		bool keep(const geo::pc::Point& pt) const {
			int cls = pt.classId();
			for(size_t i = 0; i < classes.size(); ++i) {
				if(classes[i] == cls)
					return true;
			}
			return false;
		}
		void print() const {
			std::stringstream ss;
			for(int c : classes)
				ss << c << ", ";
			g_debug("Class filter: " << ss.str());
		}
	};

	class PointScanAngleFilter : public PointFilter  {
	public:
		double minAngle;
		double maxAngle;
		PointScanAngleFilter(double minAngle = -geo::maxvalue<double>(), double maxAngle = geo::maxvalue<double>()) :
			minAngle(minAngle), maxAngle(maxAngle) {
		}
		bool keep(const geo::pc::Point& pt) const {
			double a = pt.scanAngle();
			if(a < minAngle || a > maxAngle)
				return false;
			return true;
		}
		void print() const {
			g_debug("Scan angle filter: " << minAngle << ", " << maxAngle);
		}
	};

	class PointEdgeFilter : public PointFilter  {
	public:
		bool keepEdge;
		PointEdgeFilter(bool keepEdge = false) :
			keepEdge(keepEdge) {
		}
		bool keep(const geo::pc::Point& pt) const {
			return keepEdge || !pt.isEdge();
		}
		void print() const {
			g_debug("Keep edge filter: " << keepEdge);
		}
	};

	class PointZRangeFilter : public PointFilter  {
	public:
		double minZ;
		double maxZ;
		PointZRangeFilter(double minZ = -geo::maxvalue<double>(), double maxZ = geo::maxvalue<double>()) :
			minZ(minZ), maxZ(maxZ) {
		}
		bool keep(const geo::pc::Point& pt) const {
			double z = pt.z();
			return z >= minZ && z <= maxZ;
		}
		void print() const {
			g_debug("Z range filter: " << minZ << ", " << maxZ);
		}
	};

	class PointIntensityFilter : public PointFilter  {
	public:
		double minIntensity;
		double maxIntensity;
		PointIntensityFilter(double minIntensity = -geo::maxvalue<double>(), double maxIntensity = geo::maxvalue<double>()) :
			minIntensity(minIntensity), maxIntensity(maxIntensity) {
		}
		bool keep(const geo::pc::Point& pt) const {
			double i = pt.intensity();
			return i >= minIntensity && i <= maxIntensity;
		}
		void print() const {
			g_debug("Intensity filter: " << minIntensity << ", " << maxIntensity);
		}
	};

	class PointFirstOnlyFilter : public PointFilter  {
	public:
		bool keepFirst;
		PointFirstOnlyFilter(bool keepFirst = false) :
			keepFirst(keepFirst) {
		}
		bool keep(const geo::pc::Point& pt) const {
			return !keepFirst || pt.isFirst();
		}
		void print() const {
			g_debug("First only filter: " << keepFirst);
		}
	};

	class PointLastOnlyFilter : public PointFilter  {
	public:
		bool keepLast;
		PointLastOnlyFilter(bool keepLast = false) :
			keepLast(keepLast) {
		}
		bool keep(const geo::pc::Point& pt) const {
			return !keepLast || pt.isLast();
		}
		void print() const {
			g_debug("Last only filter: " << keepLast);
		}
	};

} // filter
} // pc
} // geo



#endif /* LIBGEO_INCLUDE_PC_FILTER_HPP_ */
