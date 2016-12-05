#ifndef __AMG_EIGEN_SOLVER_H__
#define __AMG_EIGEN_SOLVER_H__

#include "multigrid.h"
#include "amg_param.h"
void Eigen_solver_amg(multigrid *amg, 
                      int nev, double *eval, 
                      double **evec, int init_approx_level, int is_init_approx_level_correction, 
                      amg_param param);

void Eigen_solver_amg_nested(multigrid *amg, 
                             int nev, double *eval, double **evec,
                             amg_param param);
#endif
