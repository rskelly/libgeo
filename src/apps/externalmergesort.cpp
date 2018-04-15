/*
 * externalmergesort.cpp
 *
 *  Created on: Apr 8, 2018
 *      Author: rob
 */

#include "externalmergesort.hpp"

int main(int argc, char** argv) {

	/*
	std::string tmpDir = argv[1];
	std::string outFile = argv[2];
	std::vector<std::string> inFiles;
	for(int i = 3; i < argc; ++i)
		inFiles.push_back(argv[i]);

	ExternalMergeSort em(outFile, inFiles, tmpDir, 1024 * 1024 , 8);

	em.sort(false);

	geo::pc::Point pt;
	double bounds[4] = {99999999.0, 99999999.0, -99999999.0, -99999999.0};
	int i = 0;
	while(em.next(pt)) {
		//std::cerr << "pt: " << pt.x() << ", " << pt.y() << "\n";
		double x = pt.x();
		double y = pt.y();
		if(x < bounds[0]) bounds[0] = x;
		if(y < bounds[1]) bounds[1] = y;
		if(x > bounds[2]) bounds[2] = x;
		if(y > bounds[3]) bounds[3] = y;
		if(++i % 1000 == 0)
			std::cerr << "bounds: " << bounds[0] << ", " << bounds[1] << ", " << bounds[2] << ", " << bounds[3] << "\n";
	}
	std::cerr << "total: " << i;
	*/
}
