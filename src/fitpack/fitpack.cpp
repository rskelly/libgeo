/*
 * fitpack_dummy.cpp
 *
 * Implementation for if the fortran compiler isn't used in util.
 *
 *  Created on: May 6, 2020
 *      Author: rob
 */

#include "fitpack.hpp"

void surfit_(int*, int*, const double*, const double*, const double*, const double*,
		double*, double*, double*, double*, int*, int*,
		double*, int*, int*, int*, double*,
		int*, double*, int*, double*, double*, double*,
		double*, int*, double*, int*, int*, int*,
		int*) {}

void surev_(int*, double*, int*, double*, int*,
		double*, const double*, int*, const double*, int*, double*, int*,
		double*, int*, int*, int*, int*) {}

void curfit_(int*, int*, const double*, const double*, const double*,
		double*, double*, int*, double*, int*, int*, double*, double*,
		double*, double*, int*, int*, int*) {}

void curev_(int*, double*, int*, double*, int*, int*,
		const double*, int*, double*, int*, int*) {}
