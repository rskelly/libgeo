/*
 * fitpack_dummy.cpp
 *
 * Implementation for if the fortran compiler isn't used in util.
 *
 *  Created on: May 6, 2020
 *      Author: rob
 */

/*
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
*/

void surfit_(int* iopt, int* m, const double* x, const double* y, const double* z, const double* w,
		double* xb, double* xe, double* yb, double* ye, int* kx, int* ky,
		double* s, int* nxest, int* nyest, int* nmax, double* eps,
		int* nx, double* tx, int* ny, double* ty, double* c, double* fp,
		double* wrk1, int* lwrk1, double* wrk2, int* lwrk2, int* iwrk, int* kwrk,
		int* ier) {}

void surev_(int* idim, double* tu, int* nu, double* tv, int* nv,
		double* c, const double* u, int* mu, const double* v, int* mv, double* f, int* mf,
		double* wrk, int* lwrk, int* iwrk, int* kwrk, int* ier) {}

void curfit_(int* iopt, int* m, const double* x, const double* y, const double* w,
		double* xb, double* xe, int* k, double* s, int* nest, int* n, double* t, double* c,
		double* fp, double* wrk, int* lwrk, int* iwrk, int* ier) {}

void curev_(int* idim, double* t, int* n, double* c, int* nc, int* k,
		const double* u, int* m, double* x, int* mx, int* ier) {}
