#ifdef WITH_MPI

#include "par_multigrid.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include "par_matrix_vector.h"
#include "par_linear_algebra.h"

// type = 1   boundary problem or eigenvalue problem
// type = 2   generalized eigenvalue problem
par_multigrid *Init_par_amg(int max_level, int type)
{
    int i;
    par_multigrid *par_amg = (par_multigrid*)malloc(sizeof(par_multigrid));
    
    par_amg->max_level           = max_level;
    par_amg->actual_level_global = 0;
    par_amg->actual_level        = 0;
    par_amg->comm  = (MPI_Comm*)    malloc(max_level*sizeof(MPI_Comm));
    par_amg->group = (MPI_Group*)   malloc(max_level*sizeof(MPI_Group));
    par_amg->A     = (par_dmatcsr**)malloc(max_level*sizeof(par_dmatcsr*));
    par_amg->P     = (par_dmatcsr**)malloc(max_level*sizeof(par_dmatcsr*));
    par_amg->R     = (par_dmatcsr**)malloc(max_level*sizeof(par_dmatcsr*));
    
    for(i=0; i<max_level; i++)
    {
	par_amg->comm[i]  = MPI_COMM_NULL;
	par_amg->group[i] = MPI_GROUP_NULL;
        par_amg->A[i] = Init_empty_par_dmatcsr();
        par_amg->P[i] = Init_empty_par_dmatcsr();
        par_amg->R[i] = Init_empty_par_dmatcsr();
        //par_amg->A[i] = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
        //par_amg->P[i] = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
        //par_amg->R[i] = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
    }
    
    switch( type )
    {
        case 1: par_amg->M = NULL; break;
        case 2: par_amg->M = (par_dmatcsr**)malloc(max_level*sizeof(par_dmatcsr*));
                for(i=0; i<max_level; i++)
                    par_amg->M[i] = Init_empty_par_dmatcsr();
                break;
        default: printf("Cannot initialize par-multi-grid: unknown type!\n");
    }

    return par_amg;
}

par_multigrid *Build_par_amg(par_dmatcsr *A, par_dmatcsr *M, int max_level, MPI_Comm comm, MPI_Group group)
{
    int type = (NULL==M) ? 1 : 2;
    par_multigrid *par_amg = Init_par_amg(max_level, type);
  
    par_amg->actual_level        = 1;
    par_amg->actual_level_global = 1;

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

void Restrict_par_f2c(par_multigrid *pamg, int fine_level, int coarse_level, par_dvec *fine_vec, par_dvec *coarse_vec)
{
    //printf("coarse level = %d, fine_level = %d\n", coarse_level, fine_level);
    if(1 < coarse_level-fine_level)
    {
        int i;
        par_dvec *work1 = Init_par_dvec_mv(pamg->A[fine_level]);
        par_dvec *work2 = Init_par_dvec_mv(pamg->A[fine_level]);
        par_dvec *tmp;
	if(MPI_COMM_NULL != pamg->comm[fine_level])
	    Multi_par_dmatcsr_dvec(pamg->R[fine_level], fine_vec, work1);
        for(i=fine_level+1; i<coarse_level-1; i++)
        {
	    if(MPI_COMM_NULL != pamg->comm[i])
		Multi_par_dmatcsr_dvec(pamg->R[i], work1, work2);
            tmp   = work1;
            work1 = work2;
            work2 = tmp;
        }
	if(MPI_COMM_NULL != pamg->comm[coarse_level-1])
	    Multi_par_dmatcsr_dvec(pamg->R[coarse_level-1], work1, coarse_vec);
        Free_par_dvec(work1);
        Free_par_dvec(work2);
    }
    else
    {
	//printf("*************************************\n");
	if(MPI_COMM_NULL != pamg->comm[fine_level])
	    Multi_par_dmatcsr_dvec(pamg->R[fine_level], fine_vec, coarse_vec);
	//printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
    }
}

void Prolong_par_c2f(par_multigrid *pamg, int coarse_level, int fine_level, par_dvec *coarse_vec, par_dvec *fine_vec)
{
    if(1 < coarse_level-fine_level)
    {
        int i;
        par_dvec *work1 = Init_par_dvec_mv(pamg->A[fine_level]);
        par_dvec *work2 = Init_par_dvec_mv(pamg->A[fine_level]);
        par_dvec *tmp;
	if(MPI_COMM_NULL != pamg->comm[coarse_level-1])
	    Multi_par_dmatcsr_dvec(pamg->P[coarse_level-1], coarse_vec, work1);
        for(i=coarse_level-2; i>fine_level; i--)
        {
	    if(MPI_COMM_NULL != pamg->comm[i])
		Multi_par_dmatcsr_dvec(pamg->P[i], work1, work2);
            tmp   = work1;
            work1 = work2;
            work2 = tmp;
        }
	if(MPI_COMM_NULL != pamg->comm[fine_level])
	    Multi_par_dmatcsr_dvec(pamg->P[fine_level], work1, fine_vec);
        Free_par_dvec(work1);
        Free_par_dvec(work2);
    }
    else
    {
	if(MPI_COMM_NULL != pamg->comm[fine_level])
	    Multi_par_dmatcsr_dvec(pamg->P[fine_level], coarse_vec, fine_vec);
    }
}

void Print_par_amg(par_multigrid *par_amg)
{
    int  myrank;
    MPI_Comm comm = par_amg->comm[0];
    MPI_Comm_rank(comm, &myrank);

    int i;

    //其实 nr_global 可以通过某层任意进程直接获得
    int actual_level_global = par_amg->actual_level_global;
    int *nr_nn_self   = (int*)calloc(actual_level_global*2, sizeof(int));
    int *nr_nn_global = (int*)calloc(actual_level_global*2, sizeof(int));
    for(i=0; i<actual_level_global; i++)
    {
	if(MPI_COMM_NULL != par_amg->comm[i])
	{
	    nr_nn_self[2*i]   = par_amg->A[i]->diag->nr;
	    nr_nn_self[2*i+1] = par_amg->A[i]->diag->nn + par_amg->A[i]->offd->nn;
	}
    }
    MPI_Allreduce(nr_nn_self, nr_nn_global, actual_level_global*2, MPI_INT, MPI_SUM, par_amg->comm[0]);

    if(0 == myrank)
    {
	double gc = 0.0;
	double oc = 0.0;
	printf("========================= multigrid =========================\n");
	printf("            level    nrow      nnz        sparse \n"); 
	for(i=0; i<actual_level_global; i++)
	{
	    int nr_g = nr_nn_global[2*i];
	    int nn_g = nr_nn_global[2*i+1];
	    double sp = (double)nn_g/(double)nr_g;
	    printf("             %2d   %7d   %7d    %9.6f\n", i, nr_g, nn_g, sp);
	    gc += (double)nr_g;
	    oc += (double)nn_g;
	}
	printf("-------------------------------------------------------------\n");
	printf("grid complexity = %f, operator complexity = %f\n", 
	       gc/(double)nr_nn_global[0], oc/(double)nr_nn_global[1]);
	printf("=============================================================\n");
    }

    free(nr_nn_global); nr_nn_global = NULL;
    free(nr_nn_self); nr_nn_self = NULL;
}
#endif
