#ifndef __PAR_LINEAR_ALGEBRA_H__
#define __PAR_LINEAR_ALGEBRA_H__

#include "par_matrix_vector.h"

#ifdef WITH_MPI

void Multi_par_dmatcsr_dvec(par_dmatcsr *A, par_dvec *x, par_dvec *y);

double Get_par_dvec_2norm(par_dvec *x);

void Multi_par_dmatcsr_dmatcsr(par_dmatcsr *A, par_dmatcsr *B, par_dmatcsr *C);

void Transpose_par_dmatcsr(par_dmatcsr *A, par_dmatcsr *AT);

void Sumself_par_dvec_axpby(par_dvec *x, double a, par_dvec *y, double b);

void Copy_par_dvec(par_dvec *x, par_dvec *x_copy);

double Multi_par_dvec_dvec(par_dvec *x, par_dvec *y);

void Scale_par_dvec(par_dvec *x, double a);
void Normalize_par_dvec(par_dvec *x);
#endif

#endif
