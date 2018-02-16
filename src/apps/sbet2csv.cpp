/*
 * sbet2csv.cpp
 *
 *  Created on: Feb 14, 2018
 *      Author: rob
 */


#include <proj_api.h>
#include <fstream>

#include "raster.hpp"

using namespace geo::raster;

double toDeg(double c) {
	return c * 180.0 / M_PI;
}

void usage() {
	std::cerr << "Usage: <sbet file> <output srid>\n";
}

typedef struct {
	double time; //	seconds	double
	double latitude; //	radians	double
	double longitude; //	radians	double
	double altitude; //	meters	double
	double xVelocity; //	meters/second	double
	double yVelocity;	//meters/second	double
	double zVelocity; //	meters/second	double
	double roll; //	radians	double
	double pitch; //	radians	double
	double heading; //	radians	double
	double wanderAngle;//	radians	double
	double xAcceleration;//	meters/second2	double
	double yAcceleration;//	meters/second2	double
	double zAcceleration; //	meters/second2	double
	double xAngularRate; //	radians/second	double
	double yAngularRate; //	radians/second	double
	double zAngularRate; //	radians/second	double
} Record;

void printRec(Record& rec) {
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

	Record rec;
	std::ifstream str(sbetfile, std::ios::binary|std::ios::in);
	do {
		str.read((char*) &rec, sizeof(Record));
		pj_transform(sp, op, 1, 1, &rec.longitude, &rec.latitude, NULL);
		printRec(rec);
	} while(true);

}


