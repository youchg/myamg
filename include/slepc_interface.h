#ifndef __SLEPC_INTERFACE_H__
#define __SLEPC_INTERFACE_H__

#ifdef WITH_SLEPC

#include "matrix.h"
void Eigen_solver_slepc(dmatcsr *A, dmatcsr *M, int nev, double *eval, double **evec);
#endif

#endif
