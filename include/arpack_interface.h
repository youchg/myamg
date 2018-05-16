#ifndef __ARPACK_INTERFACE_H__
#define __ARPACK_INTERFACE_H__

#include "matrix.h"
#include "multigrid.h"
#include "amg_param.h"
//void DSEigenSolver(dmatcsr *A,  dmatcsr *M, int nev, int indictor, double tau, 
//                   double *eval, double **evec, bool change);
void Eigen_solver_arpack_dn(dmatcsr *A, dmatcsr *M, int nev, double *eval, double **evec);

void Eigen_solver_arpack_dn_amg(multigrid *amg, int level, int nev, double *eval, double **evec, amg_param param);
#endif
