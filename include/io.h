#ifndef __IO_H__
#define __IO_H__

#include "preprocess.h"
#include "amg_param.h"
#include "matrix.h"
#include "multigrid.h"

dmatcsr *Read_dmatcsr(const char *filename);
imatcsr *Read_imatcsr(const char *filename);

void Print_dmatcsr(dmatcsr *A);
void Print_imatcsr(imatcsr *A);

void Print_ivec(int *vec, int length);

void Write_dvec(double *x, int length, const char *filename);
void Read_dvec (double *x, int length, const char *filename);

void Write_dmatcsr_csr(dmatcsr *A, const char *filename);
void Write_imatcsr_csr(imatcsr *A, const char *filename);
void Write_dmatcsr_ps (dmatcsr *A, const char *filename);
void Write_imatcsr_ps (imatcsr *A, const char *filename);

#if WITH_BMP
#include "bmp.h"
void Write_dmatcsr_bmp(dmatcsr *A, const char *filename, const int width, const int height, RGB_QUAD(*ColorMap)(double));
void Write_imatcsr_bmp(imatcsr *A, const char *filename, const int width, const int height, RGB_QUAD(*ColorMap)(double));
#endif

void Input_amg_param(const char *filename, amg_param *param, dmatcsr **A, dmatcsr **M, double **b);
void Print_amg_param(amg_param param);
void Init_amg_param_argv(int argc, char* argv[], amg_param *param, dmatcsr **A, dmatcsr **M, double **b);
void Init_nev_argv(int argc, char* argv[], int *nev, int *nb, int *ne);

void Print_amg(multigrid *amg);

#endif
