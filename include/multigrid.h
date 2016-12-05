#ifndef __MULTIGRID_H__
#define __MULTIGRID_H__

#include "matrix.h"

typedef struct MULTIGRID_
{
    int max_level;
    int actual_level;
    
    dmatcsr **A;
    dmatcsr **P;
    dmatcsr **R;
    dmatcsr **M;
} multigrid;

// type = 1   boundary problem or eigenvalue problem
// type = 2   generalized eigenvalue problem
multigrid *Init_amg(int max_level, int type);
multigrid *Build_amg(dmatcsr *A, dmatcsr *M, int max_level);

void Free_multigrid(multigrid *amg);

void ProlongCoarse2Fine (multigrid *amg, int coarse_level, int fine_level,   double *coarse_vec, double *fine_vec);
void RestrictFine2Coarse(multigrid *amg, int fine_level,   int coarse_level, double *fine_vec,   double *coarse_vec);

#endif
