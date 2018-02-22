/*
 * perturb.cpp
 *
 *  Created on: Feb 16, 2018
 *      Author: rob
 */

#include <string>
#include <fstream>
#include <random>

#include <liblas/liblas.hpp>

void usage() {
	std::cerr << "Usage: pcperturb <infile> <outfile> <xshift> <yshift> <zshift> <sigma> \n"
			<< " Sigma is the std. devation of the distribution by which\n"
			<< " the shifts will be multiplied.\n";
}

int main(int argc, char** argv) {

	if(argc < 7) {
		usage();
		return 1;
	}

	std::string infile(argv[1]);
	std::string outfile(argv[2]);
	double xshift = atof(argv[3]);
	double yshift = atof(argv[4]);
	double zshift = atof(argv[5]);
	double sig = atof(argv[6]);

	std::ifstream instr(infile, std::ios::binary|std::ios::in);
	std::ofstream ostr(outfile, std::ios::binary|std::ios::out);
	liblas::ReaderFactory rf;
	liblas::Reader rdr = rf.CreateWithStream(instr);
	liblas::Header hdr(rdr.GetHeader());

	std::default_random_engine gen;
	std::normal_distribution<double> nd(0, sig);

	double minx = std::numeric_limits<double>::max();
	double maxx = std::numeric_limits<double>::lowest();
	double minz = std::numeric_limits<double>::max();
	double maxz = std::numeric_limits<double>::lowest();
	double miny = std::numeric_limits<double>::max();
	double maxy = std::numeric_limits<double>::lowest();

	liblas::Writer wtr(ostr, hdr);
	while(rdr.ReadNextPoint()) {
		liblas::Point pt(rdr.GetPoint());
		pt.SetX(pt.GetX() + xshift + nd(gen));
		pt.SetY(pt.GetY() + yshift + nd(gen));
		pt.SetZ(pt.GetZ() + zshift + nd(gen));
		minx = std::min(minx, pt.GetX());
		maxx = std::max(maxx, pt.GetX());
		miny = std::min(miny, pt.GetY());
		maxy = std::max(maxy, pt.GetY());
		minz = std::min(minz, pt.GetZ());
		maxz = std::max(maxz, pt.GetZ());
		wtr.WritePoint(pt);
	}
	hdr.SetMin(minx, miny, minz);
	hdr.SetMax(maxx, maxy, maxz);
	wtr.SetHeader(hdr);
}
