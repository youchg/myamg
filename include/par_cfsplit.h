#ifndef __PAR_CFSPLIT_H__
#define __PAR_CFSPLIT_H__

#include "matrix.h"
#include "par_matrix_vector.h"
#include "amg_param.h"
#include "cfsplit.h"

int Generate_par_strong_coupling_set(par_dmatcsr *A, par_imatcsr *S, amg_param param);

//int Split_par_CLJP(par_dmatcsr *A, par_imatcsr *S);
int Split_par_CLJP(par_dmatcsr *A, par_imatcsr *S, par_ivec *dof);

int Get_par_interpolation_direct(par_dmatcsr *A, par_imatcsr *S, par_ivec *dof, par_dmatcsr *P);
int Generate_par_sparsity_P_dir(par_dmatcsr *A, par_imatcsr *S, 
	                        int *cfmarker, int *cfmarker_offd, 
				par_dmatcsr *P);

int Generate_par_P_dir         (par_dmatcsr *A,      par_imatcsr *S, 
	                        int *cfmarker,       int *cfmarker_offd, 
		                int *cfmarker_index, int *cfmarker_offd_index, 
		                par_dmatcsr *P);

void Truncate_par_P(par_dmatcsr *P, amg_param param);

int Get_nonnegtive_ivec_all_value_ascend(int *x, int length, int max_val, int *nvalue, int *value);
imatcsr *Get_S_ext(par_dmatcsr *A, par_imatcsr *S);

int *Insert_ascend_ivec_to_ivec(int *a, int length_a, int *b, int length_b, int *position);
int Get_nonnegtive_ivec_all_value_ascend(int *x, int length, int max_val, int *nvalue, int *value);

#endif
