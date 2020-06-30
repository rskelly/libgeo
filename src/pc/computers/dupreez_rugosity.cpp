/*
 * dupreez_rugosity.cpp
 *
 *  Created on: Apr 3, 2018
 *      Author: rob
 */

#include <vector>

#include <CGAL/Plane_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_2_algorithms.h>

#include "pointcloud.hpp"
#include "pc_computer.hpp"

using namespace geo::pc::compute;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K> Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3 Point_3;
typedef K::Plane_3 Plane_3;
typedef Delaunay::Finite_faces_iterator Finite_faces_iterator;
typedef Delaunay::Face Face;

/**
 * Compute planar area.
 */
double computePArea(double x1, double y1, double z1, double x2, double y2,
		double z2, double x3, double y3, double z3) {
	double side0 = std::sqrt(std::pow(x1 - x2, 2.0) + std::pow(y1 - y2, 2.0) + std::pow(z1 - z2, 2.0));
	double side1 = std::sqrt(std::pow(x2 - x3, 2.0) + std::pow(y2 - y3, 2.0) + std::pow(z2 - z3, 2.0));
	double side2 = std::sqrt(std::pow(x3 - x1, 2.0) + std::pow(y3 - y1, 2.0) + std::pow(z3 - z1, 2.0));
	double s = (side0 + side1 + side2) / 2.0;
	return std::sqrt(s * (s - side0) * (s - side1) * (s - side2));
}

/**
 * Compute the area of a face.
 */
double computeFArea(const Face& face) {
	Point_3 p1 = face.vertex(0)->point();
	Point_3 p2 = face.vertex(1)->point();
	Point_3 p3 = face.vertex(2)->point();
	return computePArea(p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z(), p3.x(),
			p3.y(), p3.z());
}

double toPlane(const Point_3& p, const Plane_3& plane, const Point_3&) {
	return (p.x() * plane.a() + p.y() * plane.b() + plane.d()) / -plane.c();
}

double polyArea(const std::list<Point_3>& hull, const Plane_3& plane,
		const Point_3& centroid) {
	double area = 0.0;
	auto it0 = hull.begin();
	auto it1 = hull.begin();
	it1++;
	do {
		double z0 = toPlane(*it0, plane, centroid);
		double z1 = toPlane(*it1, plane, centroid);
		area += computePArea(it0->x(), it0->y(), z0, it1->x(), it1->y(), z1,
				centroid.x(), centroid.y(), centroid.z());
		it0++;
		it1++;
		if (it1 == hull.end())
			it1 = hull.begin();
	} while (it0 != hull.end());
	return area;
}

int RugosityComputer::compute(double, double, const std::vector<geo::pc::Point>&, const std::vector<geo::pc::Point>& filtered, double radius, std::vector<double>& out) {

	if(filtered.size() >= 3) {

		double area = radius * radius * M_PI;
		double density = filtered.size() / area;

		std::list<Point_3> verts;
		for (const geo::pc::Point& pt : filtered)
			verts.emplace_back(pt.x(), pt.y(), pt.value());

		// Convex hull and POBF
		std::list<Point_3> hull;
		Plane_3 plane;
		Point_3 centroid;
		CGAL::convex_hull_2(verts.begin(), verts.end(), std::back_inserter(hull), Gt());
		CGAL::linear_least_squares_fitting_3(hull.begin(), hull.end(), plane, centroid, CGAL::Dimension_tag<0>());

		// POBF surface area.
		double parea = polyArea(hull, plane, centroid);

		// If the poly area is zero, quit.
		if(parea <= 0) {
			out.push_back(std::nan(""));
		} else {

			// Delaunay 3D surface area.
			double tarea = 0.0;
			Delaunay dt(verts.begin(), verts.end());
			for (Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it)
				tarea += computeFArea(*it);

			// TODO: This is an attempt at modelling the relationship between the ACR and
			// density. The fractal dimension is involved. Should be redone and documented.
			double densityFactor = 1.0 / (2.49127261 + 9.01659384 * std::sqrt(density * 32.65748276));

			double acr = (tarea / parea) * densityFactor;

			out.push_back(acr);
		}
	} else {
		out.push_back(std::nan(""));
	}
	return 1;
}

int RugosityComputer::bandCount() const {
	return 1;
}

std::vector<std::pair<std::string, std::string>> RugosityComputer::bandMeta() const {
	return {{"dupreez_rugosity", "ACR Rugosity (DuPreez 2014)"}};
}

