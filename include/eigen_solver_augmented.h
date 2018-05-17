#ifndef __AMG_EIGEN_SOLVER_AUGMENTED_H__
#define __AMG_EIGEN_SOLVER_AUGMENTED_H__

#include "multigrid.h"
#include "amg_param.h"

void Eigen_solver_amg_augmented(multigrid *amg, 
                                 int nev, double *eval, double **evec, 
                                 amg_param param);
#endif
