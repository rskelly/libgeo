/*
 * sbet2csv.cpp
 *
 *  Created on: Feb 14, 2018
 *      Author: rob
 */


#include <grid.hpp>
#include <proj_api.h>
#include <fstream>

#include "pc_trajectory.hpp"

using namespace geo::raster;
using namespace geo::pc;

double toDeg(double c) {
	return c * 180.0 / M_PI;
}

void usage() {
	std::cerr << "Usage: <sbet file> <output srid>\n";
}


void printRec(SBETRecord& rec) {
	std::cout << std::setprecision(12)
		<< rec.time << ","
		<< rec.latitude << "," << rec.longitude << "," << rec.altitude << ","
		<< rec.xVelocity << "," << rec.yVelocity << "," << rec.zVelocity << ","
		<< rec.roll << "," << rec.pitch << "," << rec.heading << ","
		<< rec.wanderAngle << ","
		<< rec.xAcceleration << "," << rec.yAcceleration << "," << rec.zAcceleration << ","
		<< rec.xAngularRate << "," << rec.yAngularRate << "," << rec.zAngularRate << "\n";
}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	std::string sbetfile(argv[1]);
	int srid = atoi(argv[2]);

	projPJ sp = pj_init_plus("+init=epsg:4326");

	std::stringstream os;
	os << "+init=epsg:" << srid;
	projPJ op = pj_init_plus(os.str().c_str());

	SBETRecord rec;
	std::ifstream str(sbetfile, std::ios::binary|std::ios::in);
	do {
		str.read((char*) &rec, sizeof(SBETRecord));
		pj_transform(sp, op, 1, 1, &rec.longitude, &rec.latitude, NULL);
		printRec(rec);
	} while(true);

}


