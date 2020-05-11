/**
 * This file provides an interface to selected fitpack routines in the src/fitpack folder.
 * There is a dummy c file there, and fortran files which are selected depending on compiler
 * arguments.
 */

#include "geo.hpp"

#ifdef __cplusplus
extern "C" {
#endif

	// Using const for arrays passed directly in from c++ call.

	G_DLL_EXPORT void surfit_(int* iopt, int* m, const double* x, const double* y, const double* z, const double* w,
		double* xb, double* xe, double* yb, double* ye, int* kx, int* ky,
		double* s, int* nxest, int* nyest, int* nmax, double* eps,
		int* nx, double* tx, int* ny, double* ty, double* c, double* fp,
		double* wrk1, int* lwrk1, double* wrk2, int* lwrk2, int* iwrk, int* kwrk,
		int* ier);

	G_DLL_EXPORT void surev_(int* idim, double* tu, int* nu, double* tv, int* nv,
		double* c, const double* u, int* mu, const double* v, int* mv, double* f, int* mf,
		double* wrk, int* lwrk, int* iwrk, int* kwrk, int* ier);

	G_DLL_EXPORT void curfit_(int* iopt, int* m, const double* x, const double* y, const double* w,
		double* xb, double* xe, int* k, double* s, int* nest, int* n, double* t, double* c,
		double* fp, double* wrk, int* lwrk, int* iwrk, int* ier);

	G_DLL_EXPORT void curev_(int* idim, double* t, int* n, double* c, int* nc, int* k,
		const double* u, int* m, double* x, int* mx, int* ier);

#ifdef __cplusplus
}
#endif
