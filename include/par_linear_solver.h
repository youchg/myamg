#ifndef __PAR_LINEAR_SOLVER_H__
#define __PAR_LINEAR_SOLVER_H__

#include "preprocess.h"
#include "amg_param.h"
#include "matrix.h"
#include "multigrid.h"
#include "par_matrix_vector.h"
#include "par_multigrid.h"

#if WITH_MPI
double Linear_solver_par_amgcycle(par_multigrid *pamg, int current_level,
                                  par_dvec *b,         par_dvec *x, 
                                  amg_param param);

#if 0
void Get_residual(dmatcsr *A, double *b, double *x, 
	          double *r, double *rnorm);


void Linear_solver_cg(dmatcsr *A, double *b,  double *x, 
	              double tol, int max_iter, 
	              double *resi, double *resi_norm, int *niter,
	              int print_level);



void Linear_solver_amg(multigrid *amg, int current_level,
	               double *b, double *x, 
	               amg_param param, 
	               double *resi_norm, int *ncycle);


void Linear_solver_pcg_amg(multigrid *amg, int current_level, 
                           double *b,  double *x, 
                           amg_param param, 
                           double *resi, double *resi_norm, int *nit,
                           int *ncycle_amg_total);

#endif

#endif
#endif
