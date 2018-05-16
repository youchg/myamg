#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "linear_algebra.h"
#include "multigrid.h"

// type = 1   boundary problem or eigenvalue problem
// type = 2   generalized eigenvalue problem
multigrid *Init_amg(int max_level, int type)
{
    int i;
    multigrid *amg = (multigrid*)malloc(sizeof(multigrid));
    
    amg->max_level    = max_level;
    amg->actual_level = 0;
    amg->A = (dmatcsr**)malloc(max_level*sizeof(dmatcsr*));
    amg->P = (dmatcsr**)malloc(max_level*sizeof(dmatcsr*));
    amg->R = (dmatcsr**)malloc(max_level*sizeof(dmatcsr*));
    
    for(i=0; i<max_level; i++)
    {
        amg->A[i] = (dmatcsr*)malloc(sizeof(dmatcsr));
        amg->P[i] = (dmatcsr*)malloc(sizeof(dmatcsr));
        amg->R[i] = (dmatcsr*)malloc(sizeof(dmatcsr));
    }
    
    switch( type )
    {
        case 1: amg->M = NULL; break;
        case 2: amg->M = (dmatcsr**)malloc(max_level*sizeof(dmatcsr*));
                for(i=0; i<max_level; i++)
                    amg->M[i] = (dmatcsr*)malloc(sizeof(dmatcsr));
                break;
        default: printf("Cannot initialize multi-grid: unknown type!\n");
    }

    return amg;
}

multigrid *Build_amg(dmatcsr *A, dmatcsr *M, int max_level)
{
    int type = (NULL==M) ? 1 : 2;
    multigrid *amg = Init_amg(max_level, type);
  
    amg->actual_level = 1;
    free(amg->A[0]);
    amg->A[0] = Copy_dmatcsr(A);
    if(NULL != M)
    {
        free(amg->M[0]);
        amg->M[0] = Copy_dmatcsr(M);
    }
    return amg;
}

void Free_multigrid(multigrid *amg)
{
    int max_level    = amg->max_level;
    int actual_level = amg->actual_level;
    int i;
    
    if(NULL != amg->A)
    {
        for(i=0; i<actual_level; i++)
            Free_dmatcsr(amg->A[i]);
        for(i=actual_level; i<max_level; i++)
        {
            free(amg->A[i]);
            amg->A[i] = NULL;
        }
    }
    free(amg->A);
    amg->A = NULL;
    
    if(NULL != amg->P)
    {
        for(i=0; i<actual_level-1; i++)
            Free_dmatcsr(amg->P[i]);
        for(i=actual_level-1; i<max_level; i++)
        {
            free(amg->P[i]);
            amg->P[i] = NULL;
        }
    }
    free(amg->P);
    amg->P = NULL;
    
    if(NULL != amg->R)
    {
        for(i=0; i<actual_level-1; i++)
            Free_dmatcsr(amg->R[i]);
        for(i=actual_level-1; i<max_level; i++)
        {
            free(amg->R[i]);
            amg->R[i] = NULL;
        }
    }
    free(amg->R);
    amg->R = NULL;
    
    if(NULL != amg->M)
    {
        for(i=0; i<actual_level; i++)
            Free_dmatcsr(amg->M[i]);
        for(i=actual_level; i<max_level; i++)
        {
            free(amg->M[i]);
            amg->M[i] = NULL;
        }
    }
    free(amg->M);
    amg->M = NULL;

    amg->max_level    = 0;
    amg->actual_level = 0;
    
    free(amg);
    amg = NULL;
}


void RestrictFine2Coarse(multigrid *amg, int fine_level, int coarse_level, double *fine_vec, double *coarse_vec)
{
    if(1 < coarse_level-fine_level)
    {
        int i;
        double *work1 = (double*)calloc(amg->A[fine_level]->nr, sizeof(double));
        double *work2 = (double*)calloc(amg->A[fine_level]->nr, sizeof(double));
        double *tmp;
        Multi_dmatcsr_dvec(amg->R[fine_level], fine_vec, work1);
        for(i=fine_level+1; i<coarse_level-1; i++)
        {
            Multi_dmatcsr_dvec(amg->R[i], work1, work2);
            tmp   = work1;
            work1 = work2;
            work2 = tmp;
        }
        Multi_dmatcsr_dvec(amg->R[coarse_level-1], work1, coarse_vec);
        free(work1);
        free(work2);
    }
    else
    {
        Multi_dmatcsr_dvec(amg->R[fine_level], fine_vec, coarse_vec);
    }
}

void ProlongCoarse2Fine(multigrid *amg, int coarse_level, int fine_level, double *coarse_vec, double *fine_vec)
{
    if( 1 < coarse_level-fine_level)
    {
        int i;
        double *work1 = (double*)calloc(amg->A[fine_level]->nr, sizeof(double));
        double *work2 = (double*)calloc(amg->A[fine_level]->nr, sizeof(double));
        double *tmp;
        Multi_dmatcsr_dvec(amg->P[coarse_level-1], coarse_vec, work1);
        for(i=coarse_level-2; i>fine_level; i--)
        {
            Multi_dmatcsr_dvec(amg->P[i], work1, work2);
            tmp   = work1;
            work1 = work2;
            work2 = tmp;
        }
        Multi_dmatcsr_dvec(amg->P[fine_level], work1, fine_vec);
        free(work1);
        free(work2);
    }
    else
    {
        Multi_dmatcsr_dvec(amg->P[fine_level], coarse_vec, fine_vec);
    }
}
