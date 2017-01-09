#ifndef __PAR_MULTIGRID_H__
#define __PAR_MULTIGRID_H__

#include "par_matrix_vector.h"

#if WITH_MPI

typedef struct PAR_MULTIGRID_
{
    int max_level;
    int actual_level;
    
    MPI_Comm  *comm;
    MPI_Group *group;

    par_dmatcsr **A;
    par_dmatcsr **P;
    par_dmatcsr **R;
    par_dmatcsr **M;
} par_multigrid;

// type = 1   boundary problem or eigenvalue problem
// type = 2   generalized eigenvalue problem
par_multigrid *Init_par_amg(int max_level, int type);
//par_multigrid *Build_par_amg(par_dmatcsr *A, par_dmatcsr *M, int max_level);
par_multigrid *Build_par_amg(par_dmatcsr *A, par_dmatcsr *M, int max_level, MPI_Comm comm, MPI_Group group);

void Free_par_multigrid(par_multigrid *par_amg);

//void ProlongCoarse2Fine (multigrid *amg, int coarse_level, int fine_level,   double *coarse_vec, double *fine_vec);
//void RestrictFine2Coarse(multigrid *amg, int fine_level,   int coarse_level, double *fine_vec,   double *coarse_vec);

#endif
#endif
