/*
 * perturb.cpp
 *
 *  Created on: Feb 16, 2018
 *      Author: rob
 */

#include <string>
#include <fstream>
#include <random>

#include <iostream>
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>

#include <liblas/liblas.hpp>

int main(int argc, char** argv) {

	std::string infile1(argv[1]);
	std::string infile2(argv[2]);

	liblas::ReaderFactory rf;

	pcl::PointCloud<pcl::PointXYZ>::Ptr pc1(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr pc2(new pcl::PointCloud<pcl::PointXYZ>);

	{
		std::ifstream instr(infile1, std::ios::binary|std::ios::in);
		liblas::Reader rdr = rf.CreateWithStream(instr);
		liblas::Header hdr(rdr.GetHeader());
		pc1->width = hdr.GetPointRecordsCount();
		pc1->height = 1;
		pc1->is_dense = true;
		// pc1.sensor_orientation_
		// pc1.sensor_origin_
		size_t i = 0;
		while(rdr.ReadNextPoint()) {
			const liblas::Point& pt = rdr.GetPoint();
			pcl::PointXYZ& p0 = pc1->points[i++];
			p0.x = pt.GetX();
			p0.y = pt.GetY();
			p0.z = pt.GetZ();
		}
	}

	{
		std::ifstream instr(infile2, std::ios::binary|std::ios::in);
		liblas::Reader rdr = rf.CreateWithStream(instr);
		liblas::Header hdr(rdr.GetHeader());
		pc2->width = hdr.GetPointRecordsCount();
		pc2->height = 1;
		pc2->is_dense = true;
		// pc1.sensor_orientation_
		// pc1.sensor_origin_
		size_t i = 0;
		while(rdr.ReadNextPoint()) {
			const liblas::Point& pt = rdr.GetPoint();
			pcl::PointXYZ& p0 = pc2->points[i++];
			p0.x = pt.GetX();
			p0.y = pt.GetY();
			p0.z = pt.GetZ();
		}
	}

	pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
	icp.setInputSource(pc1);
	icp.setInputTarget(pc2);

	pcl::PointCloud<pcl::PointXYZ> final;
	icp.align(final);

	std::cout << "Converged: " << icp.hasConverged() << "\n";
	std::cout << "Fitness: " << icp.getFitnessScore() << "\n";
	std::cout << "Trans: " << icp.getFinalTransformation() << "\n";

}
