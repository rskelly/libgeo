/*
 * externalmergesort.cpp
 *
 *  Created on: Apr 8, 2018
 *      Author: rob
 */

#include "externalmergesort.hpp"

int main(int argc, char** argv) {

	std::string inFile = argv[1];
	std::string tmpDir = argv[2];
	std::string outFile = argv[3];

	std::vector<std::string> inFiles;
	inFiles.push_back(inFile);

	ExternalMergeSort em(outFile, inFiles, tmpDir, 1024 * 1024 , 8);

	em.sort(false);

	geo::pc::Point pt;
	double x, y;
	while(em.next(pt)) {
		//std::cerr << "pt: " << pt.x() << ", " << pt.y() << "\n";
		x = pt.x();
		y = pt.y();
	}

}
