#ifndef __LINEAR_ALGEBRA_H__
#define __LINEAR_ALGEBRA_H__

#include "matrix.h"

void Transpose_dmatcsr       (dmatcsr *A, dmatcsr *AT);
void Transpose_imatcsr       (imatcsr *A, imatcsr *AT);
void Transpose_imatcsr_struct(imatcsr *A, imatcsr *AT);

void Remove_zero_dmatcsr(dmatcsr *A);

dmatcsr *Expand_dmatcsr_struct(dmatcsr *A, int expan);
void     Expand_dmatcsr       (dmatcsr *A, int expan, double **vec, double **mat);

void Multi_dmatcsr_dmatcsr             (dmatcsr *A, dmatcsr *B, dmatcsr *C);
void Multi_dmatcsr_dmatcsr_unsorted    (dmatcsr *A, dmatcsr *B, dmatcsr *C);
void Multi_dmatcsr_dmatcsr_quick_sorted(dmatcsr *A, dmatcsr *B, dmatcsr *C);

void Multi_dmatcsr_dmatcsr_dmatcsr                 (dmatcsr *A, dmatcsr *B, dmatcsr *C, dmatcsr *D);
void Multi_dmatcsr_dmatcsr_dmatcsr_unsorted        (dmatcsr *A, dmatcsr *B, dmatcsr *C, dmatcsr *D);
void Multi_dmatcsr_dmatcsr_dmatcsr_quick_sorted    (dmatcsr *A, dmatcsr *B, dmatcsr *C, dmatcsr *D);
void Multi_dmatcsr_dmatcsr_dmatcsr_insertion_sorted(dmatcsr *A, dmatcsr *B, dmatcsr *C, dmatcsr *D);

void Multi_dmatcsr_dvec(dmatcsr *A, double *x, double *y);

double Multi_dvec_dvec   (double *x, double *y, int length);
void   Sum_dvec_axpby    (double *x, double a, double *y, double b, double *z, int length);
void   Sumself_dvec_axpby(double *x, double a, double *y, double b, int length);
void   Scale_dvec        (double *x, double a, int length);
void   Normalize_dvec    (double *x, int length);
double Get_dvec_2norm    (double *x, int length);

int Max_abs_dvec(double *x, int length, double *max, int *pos);
int Max_dvec    (double *x, int length, double *max, int *pos);
int Max_ivec    (int    *x, int length, int    *max, int *pos);

int Equal_dmatcsr_struct(dmatcsr *A, dmatcsr *M);
dmatcsr *Sum_dmatcsr_mApnB(double m, dmatcsr *A, double n, dmatcsr *B);
#endif
