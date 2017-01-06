#ifndef __PAR_LINEAR_ALGEBRA_H__
#define __PAR_LINEAR_ALGEBRA_H__

#include "par_matrix_vector.h"

#if WITH_MPI

void Multi_par_dmatcsr_dvec(par_dmatcsr *A, par_dvec *x, par_dvec *y);

double Get_par_dvec_2norm(par_dvec *x);

void Multi_par_dmatcsr_dmatcsr(par_dmatcsr *A, par_dmatcsr *B, par_dmatcsr *C);

void Transpose_par_dmatcsr(par_dmatcsr *A, par_dmatcsr *AT);

#endif

#endif
