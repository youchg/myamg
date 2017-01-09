#include "par_multigrid.h"

#if WITH_MPI

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include "linear_algebra.h"

// type = 1   boundary problem or eigenvalue problem
// type = 2   generalized eigenvalue problem
par_multigrid *Init_par_amg(int max_level, int type)
{
    int i;
    par_multigrid *par_amg = (par_multigrid*)malloc(sizeof(par_multigrid));
    
    par_amg->max_level    = max_level;
    par_amg->actual_level = 0;
    par_amg->comm  = (MPI_Comm*)    malloc(max_level*sizeof(MPI_Comm));
    par_amg->group = (MPI_Group*)   malloc(max_level*sizeof(MPI_Group));
    par_amg->A     = (par_dmatcsr**)malloc(max_level*sizeof(par_dmatcsr*));
    par_amg->P     = (par_dmatcsr**)malloc(max_level*sizeof(par_dmatcsr*));
    par_amg->R     = (par_dmatcsr**)malloc(max_level*sizeof(par_dmatcsr*));
    
    for(i=0; i<max_level; i++)
    {
	par_amg->comm[i]  = MPI_COMM_NULL;
	par_amg->group[i] = MPI_GROUP_NULL;
        par_amg->A[i] = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
        par_amg->P[i] = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
        par_amg->R[i] = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
    }
    
    switch( type )
    {
        case 1: par_amg->M = NULL; break;
        case 2: par_amg->M = (par_dmatcsr**)malloc(max_level*sizeof(par_dmatcsr*));
                for(i=0; i<max_level; i++)
                    par_amg->M[i] = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
                break;
        default: printf("Cannot initialize par-multi-grid: unknown type!\n");
    }

    return par_amg;
}

par_multigrid *Build_par_amg(par_dmatcsr *A, par_dmatcsr *M, int max_level, MPI_Comm comm, MPI_Group group)
{
    int type = (NULL==M) ? 1 : 2;
    par_multigrid *par_amg = Init_par_amg(max_level, type);
  
    par_amg->actual_level = 1;

    MPI_Comm_dup(comm, &par_amg->comm[0]);
    MPI_Group_excl(group, 0, NULL, &par_amg->group[0]);

    //int group_coarse_size;
    //int group_coarse_rank;
    //MPI_Group_size(par_amg->group[0], &group_coarse_size);
    //MPI_Group_rank(par_amg->group[0], &group_coarse_rank);
    //printf("group_coarse_size = %d, group_coarse_rank %d.\n", group_coarse_size, group_coarse_rank);

    free(par_amg->A[0]);
    par_amg->A[0] = Copy_par_dmatcsr(A);

    if(NULL != M)
    {
        free(par_amg->M[0]);
        par_amg->M[0] = Copy_par_dmatcsr(M);
    }
    return par_amg;
}

void Free_par_multigrid(par_multigrid *par_amg)
{
    int max_level    = par_amg->max_level;
    int actual_level = par_amg->actual_level;
    int i;
    
    if(NULL != par_amg->A)
    {
        for(i=0; i<actual_level; i++)
            Free_par_dmatcsr(par_amg->A[i]);
        for(i=actual_level; i<max_level; i++)
        {
            free(par_amg->A[i]);
            par_amg->A[i] = NULL;
        }
    }
    free(par_amg->A);
    par_amg->A = NULL;
    
    if(NULL != par_amg->P)
    {
        for(i=0; i<actual_level-1; i++)
            Free_par_dmatcsr(par_amg->P[i]);
        for(i=actual_level-1; i<max_level; i++)
        {
            free(par_amg->P[i]);
            par_amg->P[i] = NULL;
        }
    }
    free(par_amg->P);
    par_amg->P = NULL;
    
    if(NULL != par_amg->R)
    {
        for(i=0; i<actual_level-1; i++)
            Free_par_dmatcsr(par_amg->R[i]);
        for(i=actual_level-1; i<max_level; i++)
        {
            free(par_amg->R[i]);
            par_amg->R[i] = NULL;
        }
    }
    free(par_amg->R);
    par_amg->R = NULL;
    
    if(NULL != par_amg->M)
    {
        for(i=0; i<actual_level; i++)
            Free_par_dmatcsr(par_amg->M[i]);
        for(i=actual_level; i<max_level; i++)
        {
            free(par_amg->M[i]);
            par_amg->M[i] = NULL;
        }
    }
    free(par_amg->M);
    par_amg->M = NULL;

    par_amg->max_level    = 0;
    par_amg->actual_level = 0;

    if(NULL != par_amg->group)
    {
        for(i=0; i<actual_level; i++)
	{
	    if(MPI_GROUP_NULL!=par_amg->group[i])
		MPI_Group_free(&par_amg->group[i]);
	}
    }
    free(par_amg->group);
    par_amg->group = NULL;

    if(NULL != par_amg->comm)
    {
        for(i=0; i<actual_level; i++)
	{
	    if((MPI_COMM_WORLD!=par_amg->comm[i]) && (MPI_COMM_NULL!=par_amg->comm[i]))
		MPI_Comm_free(&par_amg->comm[i]);
	}
    }
    free(par_amg->comm);
    par_amg->comm = NULL;
    
    free(par_amg);
    par_amg = NULL;
}


#if 0
void Restrict_par_f2c(par_multigrid *pamg, int fine_level, int coarse_level, par_dvec *fine_vec, par_dvec *coarse_vec)
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
#endif

#endif
