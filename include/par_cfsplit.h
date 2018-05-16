#ifndef __PAR_CFSPLIT_H__
#define __PAR_CFSPLIT_H__

#include "par_matrix_vector.h"

#ifdef WITH_MPI

#include "matrix.h"
#include "amg_param.h"
#include "cfsplit.h"

int Generate_par_strong_coupling_set(par_dmatcsr *A, par_imatcsr *S, amg_param param);

//int Split_par_CLJP(par_dmatcsr *A, par_imatcsr *S);
int Split_par_CLJP(par_dmatcsr *A, par_imatcsr *S, par_ivec *dof);

//int Get_par_interpolation_direct(par_dmatcsr *A, par_imatcsr *S, par_ivec *dof, par_dmatcsr *P);
int Get_par_interpolation_direct(par_dmatcsr *A, par_imatcsr *S, par_ivec *dof, par_dmatcsr *P, int *ncpt_proc);
int Generate_par_sparsity_P_dir(par_dmatcsr *A, par_imatcsr *S, 
	                        int *cfmarker, int *cfmarker_offd, 
				par_dmatcsr *P);

int Generate_par_P_dir         (par_dmatcsr *A,      par_imatcsr *S, 
	                        int *cfmarker,       int *cfmarker_offd, 
		                int *cfmarker_index, int *cfmarker_offd_index, 
		                par_dmatcsr *P);

void Truncate_par_P(par_dmatcsr *P, amg_param param);

imatcsr *Get_S_ext(par_dmatcsr *A, par_imatcsr *S);

void Reset_par_dof(void);

#endif
#endif
