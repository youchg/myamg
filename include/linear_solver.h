#ifndef __LINEAR_SOLVER_H__
#define __LINEAR_SOLVER_H__

#include "preprocess.h"
#include "matrix.h"
#include "multigrid.h"
#include "amg_param.h"

void Get_residual(dmatcsr *A, double *b, double *x, 
	          double *r, double *rnorm);

/* A = D - L - U */
void Split_dmatcsr(dmatcsr *A, dmatcsr *D, dmatcsr *L, dmatcsr *U);
void Get_diag_dmatcsr(dmatcsr *A, double *d);
void Zoom_dmatcsr(dmatcsr *A, double *d);

void Linear_solver_cg(dmatcsr *A, double *b,  double *x, 
	              double tol, int max_iter, 
	              double *resi, double *resi_norm, int *niter,
	              int print_level);

/*
void Linear_solver_gs(dmatcsr *A, dmatcsr *D, dmatcsr *L, dmatcsr *U, 
                      double *b,  double *x, 
	              double tol, int max_iter, 
	              double *resi, double *resi_norm, int *niter,
	              int print_level);
		      */
void Linear_solver_gs(dmatcsr *A, double *d,
                      double *b,  double *x, 
	              double tol, int max_iter, 
	              double *resi, double *resi_norm, int *niter,
	              int print_level);

double Linear_solver_amgcycle(multigrid *amg, int current_level,
                               double *b, double *x, 
                               amg_param param);

void Linear_solver_amg(multigrid *amg, int current_level,
	               double *b, double *x, 
	               amg_param param, 
	               double *resi_norm, int *ncycle);


void Linear_solver_pcg_amg(multigrid *amg, int current_level, 
                           double *b,  double *x, 
                           amg_param param, 
                           double *resi, double *resi_norm, int *nit,
                           int *ncycle_amg_total);

#ifdef WITH_UMFPACK
void Linear_solver_direct(dmatcsr *A, double *b, double *x);
#endif

#endif
