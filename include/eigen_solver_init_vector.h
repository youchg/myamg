#ifndef __EIGEN_SOLVER_INIT_VECTOR__
#define __EIGEN_SOLVER_INIT_VECTOR__

#include "matrix.h"
#include "multigrid.h"
#include "amg_param.h"

void Get_eigen_solver_init_vector(multigrid *amg, int nev, double *eval, double **evec, amg_param param);

#endif
