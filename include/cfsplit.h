#ifndef __CFSPLIT_H__
#define __CFSPLIT_H__

#include "matrix.h"
#include "amg_param.h"

void Reset_dof(dmatcsr *A);

int Generate_strong_coupling_set(dmatcsr *A, imatcsr *S, amg_param param);

int Pre_split     (dmatcsr *A, imatcsr *S, int *dof);
int Pre_split_fasp(dmatcsr *A, imatcsr *S, int *dof);
int Post_split(imatcsr *S, int *dof);

int CLJP_split(dmatcsr *A, imatcsr *S, int *dof);

int Generate_sparsity_P_dir(imatcsr *S, int *dof, dmatcsr *P);
int Generate_P_dir(dmatcsr *A, imatcsr *S, int *dof, dmatcsr *P);

int Generate_sparsity_P_std(imatcsr *S, int *dof, dmatcsr *P);
int Generate_P_std(dmatcsr *A, imatcsr *S, int *dof, dmatcsr *P);

void Truncate_P(dmatcsr *P, amg_param param);

#endif
