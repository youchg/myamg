#ifndef __PAR_AMG_EIGEN_SOLVER_H__
#define __PAR_AMG_EIGEN_SOLVER_H__

#ifdef WITH_MPI

#include "par_multigrid.h"
#include "amg_param.h"
void Eigen_solver_par_amg(par_multigrid *pamg, 
                          int nev, double *eval, 
                          par_dvec **evec, int init_approx_level, int is_init_approx_level_correction, 
                          amg_param param);

void Eigen_solver_par_amg_nested(par_multigrid *pamg, 
                                 int nev, double *eval, par_dvec **evec,
                                 amg_param param);
#endif
#endif
