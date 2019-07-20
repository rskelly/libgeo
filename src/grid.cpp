/*
 * grid.cpp
 *
 *  Created on: Jun 16, 2019
 *      Author: rob
 */

#include <string>
#include <algorithm>

#include "grid.hpp"

using namespace geo::grid;

Interleave geo::grid::interleaveFromString(const std::string& str) {
	std::string lower;
	std::transform(str.begin(), str.end(), lower.begin(), ::tolower);
	if(lower == "band") {
		return Interleave::BIL;
	} else {
		return Interleave::BIP;
	}
}

