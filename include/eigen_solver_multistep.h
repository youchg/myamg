#ifndef __AMG_EIGEN_SOLVER_MULTISTEP_H__
#define __AMG_EIGEN_SOLVER_MULTISTEP_H__

#include "multigrid.h"
#include "amg_param.h"
void Eigen_solver_amg_multistep(multigrid *amg, 
                                int nev, double *eval, double **evec, 
				int init_approx_level, int is_init_approx_level_correction, 
				int multistep,
				amg_param param);
#endif
