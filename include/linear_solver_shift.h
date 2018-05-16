#ifndef __LINEAR_SOLVER_SHIFT_H__
#define __LINEAR_SOLVER_SHIFT_H__

#include "preprocess.h"
#include "matrix.h"
#include "multigrid.h"
#include "amg_param.h"
#include "linear_solver.h"


double Linear_solver_shift_amgcycle(multigrid *amg, double shift, int current_level,
                                     double *b, double *x, 
			             amg_param param);

void Linear_solver_shift_amg(multigrid *amg, double shift, int current_level,
			     double *b, double *x, 
			     amg_param param, 
			     double *resi_norm, int *ncycle);

void Linear_solver_shift_pcg_amg(multigrid *amg, double shift, int current_level, 
			         double *b,  double *x, 
			         amg_param param, 
			         double *resi, double *resi_norm, int *nit,
			         int *ncycle_amg_total);
#endif
