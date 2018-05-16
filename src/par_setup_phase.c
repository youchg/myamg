#ifdef WITH_MPI

#include "par_setup_phase.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "preprocess.h"
#include "tool.h"
#include "amg_param.h"
#include "par_matrix_vector.h"
#include "par_linear_algebra.h"
#include "par_cfsplit.h"
#include "par_multigrid.h"

int is_print = 0;

void Setup_par_phase(par_multigrid *par_amg, amg_param param)
{
    int max_level    = par_amg->max_level;
    int i, j;
    if(NULL != par_amg->M)
    {
        for(i=0; i<max_level-1; i++)
        {
	    if(MPI_COMM_NULL != par_amg->comm[i])
	    {
		int  nproc_global;
		MPI_Comm_size(par_amg->comm[i], &nproc_global);

		int *ncpt_proc = (int*)calloc(nproc_global, sizeof(int));
		if(!Coarsen_par(par_amg->A[i],   par_amg->M[i],   par_amg->P[i],   par_amg->R[i], 
			    par_amg->A[i+1], par_amg->M[i+1], ncpt_proc,       param))
		{
		    free(ncpt_proc); ncpt_proc = NULL;
		    break;
		}

		int nproc_cpt = 0;
		for(j=0; j<nproc_global; j++) if(ncpt_proc[j] > 0) nproc_cpt++;
		int *proc_coarse = (int*)malloc(nproc_cpt* sizeof(int));
		nproc_cpt = 0;
		for(j=0; j<nproc_global; j++) if(ncpt_proc[j] > 0) proc_coarse[nproc_cpt++] = j;

		MPI_Group_incl(par_amg->group[i], nproc_cpt, proc_coarse, &par_amg->group[i+1]);
		MPI_Comm_create(par_amg->comm[i], par_amg->group[i+1], &par_amg->comm[i+1]);
		if(MPI_COMM_NULL != par_amg->comm[i+1])
		{
		    int *map_proc_f2c = (int*)malloc(nproc_global * sizeof(int));
		    for(j=0; j<nproc_global; j++) map_proc_f2c[j] = -1;
		    int index = 0;
		    for(j=0; j<nproc_global; j++) if(ncpt_proc[j] > 0) map_proc_f2c[j] = index++;
		    assert(index == nproc_cpt);

		    par_comm_info *AH_comm_info = par_amg->A[i+1]->comm_info;
		    for(j=0; j<AH_comm_info->nproc_neighbor; j++)
		    {
			assert(map_proc_f2c[AH_comm_info->proc_neighbor[j]] >= 0);
			AH_comm_info->proc_neighbor[j] = map_proc_f2c[AH_comm_info->proc_neighbor[j]];
		    }
		    index = 0;
		    for(j=0; j<nproc_global; j++)
		    {
			if(ncpt_proc[j] > 0)
			{
			    par_amg->A[i+1]->row_start[index]   = par_amg->A[i+1]->row_start[j];
			    par_amg->A[i+1]->row_start[index+1] = par_amg->A[i+1]->row_start[j+1];
			    par_amg->A[i+1]->col_start[index]   = par_amg->A[i+1]->col_start[j];
			    par_amg->A[i+1]->col_start[index+1] = par_amg->A[i+1]->col_start[j+1];
			    index++;
			}
		    }
		    par_amg->A[i+1]->comm = par_amg->comm[i+1];

		    par_comm_info *MH_comm_info = par_amg->M[i+1]->comm_info;
		    for(j=0; j<MH_comm_info->nproc_neighbor; j++)
		    {
			assert(map_proc_f2c[MH_comm_info->proc_neighbor[j]] >= 0);
			MH_comm_info->proc_neighbor[j] = map_proc_f2c[MH_comm_info->proc_neighbor[j]];
		    }
		    index = 0;
		    for(j=0; j<nproc_global; j++)
		    {
			if(ncpt_proc[j] > 0)
			{
			    par_amg->M[i+1]->row_start[index]   = par_amg->M[i+1]->row_start[j];
			    par_amg->M[i+1]->row_start[index+1] = par_amg->M[i+1]->row_start[j+1];
			    par_amg->M[i+1]->col_start[index]   = par_amg->M[i+1]->col_start[j];
			    par_amg->M[i+1]->col_start[index+1] = par_amg->M[i+1]->col_start[j+1];
			    index++;
			}
		    }
		    par_amg->M[i+1]->comm = par_amg->comm[i+1];

		    free(map_proc_f2c); map_proc_f2c = NULL;
		}

		free(proc_coarse); proc_coarse = NULL;
		free(ncpt_proc); ncpt_proc = NULL;

		par_amg->actual_level++;
	    }
	    else
	    {
		break;
	    }
        }
    }
    else
    {
        for(i=0; i<max_level-1; i++)
        {
	    //if(i == 6) is_print = 1;
	    if(MPI_COMM_NULL != par_amg->comm[i])
	    {
		int  nproc_global;
		MPI_Comm_size(par_amg->comm[i], &nproc_global);

		int *ncpt_proc = (int*)calloc(nproc_global, sizeof(int));
		if(!Coarsen_par(par_amg->A[i],   NULL, par_amg->P[i],   par_amg->R[i], 
			        par_amg->A[i+1], NULL, ncpt_proc,       param))
		{
		    free(ncpt_proc); ncpt_proc = NULL;
		    //printf("multigrid actual level = %d\n", actual_level);
		    break;
		}

		int nproc_cpt = 0;
		for(j=0; j<nproc_global; j++) if(ncpt_proc[j] > 0) nproc_cpt++;
		int *proc_coarse = (int*)malloc(nproc_cpt* sizeof(int));
		nproc_cpt = 0;
		for(j=0; j<nproc_global; j++) if(ncpt_proc[j] > 0) proc_coarse[nproc_cpt++] = j;

		assert(MPI_GROUP_NULL != par_amg->group[i]);

		MPI_Group_incl(par_amg->group[i], nproc_cpt, proc_coarse, &par_amg->group[i+1]);
		MPI_Comm_create(par_amg->comm[i], par_amg->group[i+1], &par_amg->comm[i+1]);
		if(MPI_COMM_NULL != par_amg->comm[i+1])
		{
		    int *map_proc_f2c = (int*)malloc(nproc_global * sizeof(int));
		    for(j=0; j<nproc_global; j++) map_proc_f2c[j] = -1;
		    int index = 0;
		    for(j=0; j<nproc_global; j++) if(ncpt_proc[j] > 0) map_proc_f2c[j] = index++;
		    assert(index == nproc_cpt);

		    par_comm_info *AH_comm_info = par_amg->A[i+1]->comm_info;
		    for(j=0; j<AH_comm_info->nproc_neighbor; j++)
		    {
			assert(map_proc_f2c[AH_comm_info->proc_neighbor[j]] >= 0);
			AH_comm_info->proc_neighbor[j] = map_proc_f2c[AH_comm_info->proc_neighbor[j]];
		    }
		    index = 0;
		    for(j=0; j<nproc_global; j++)
		    {
			if(ncpt_proc[j] > 0)
			{
			    par_amg->A[i+1]->row_start[index]   = par_amg->A[i+1]->row_start[j];
			    par_amg->A[i+1]->row_start[index+1] = par_amg->A[i+1]->row_start[j+1];
			    par_amg->A[i+1]->col_start[index]   = par_amg->A[i+1]->col_start[j];
			    par_amg->A[i+1]->col_start[index+1] = par_amg->A[i+1]->col_start[j+1];
			    index++;
			}
		    }
		    par_amg->A[i+1]->comm = par_amg->comm[i+1];

		    free(map_proc_f2c); map_proc_f2c = NULL;
		}

		free(proc_coarse); proc_coarse = NULL;
		free(ncpt_proc); ncpt_proc = NULL;

		par_amg->actual_level++;
	    }
	    else
	    {
		break;
	    }
        }
    }

    MPI_Allreduce(&par_amg->actual_level, &par_amg->actual_level_global, 1, MPI_INT, MPI_MAX, par_amg->comm[0]);
}

int Coarsen_par(par_dmatcsr *A,  par_dmatcsr *M, 
	        par_dmatcsr *P,  par_dmatcsr *R, 
		par_dmatcsr *AH, par_dmatcsr *MH, 
		int *ncpt_proc,  amg_param param)
{
    par_ivec *dof = Init_par_ivec_length_comm(A->diag->nr, A->comm);
    par_imatcsr *S = (par_imatcsr*)malloc(sizeof(par_imatcsr));
    Generate_par_strong_coupling_set(A, S, param);
    int ncpt_total = Split_par_CLJP(A, S, dof);

    if(ncpt_total < param.max_coarsest_dof)
    {
	Free_par_imatcsr(S);
	Free_par_ivec(dof);
	return FAIL;
    }

    //if(is_print) printf("get interpolation...\n");
    Get_par_interpolation_direct(A, S, dof, P, ncpt_proc);
    //if(is_print) printf("get prolongation...\n");
    Transpose_par_dmatcsr(P, R);
    //if(is_print) Write_par_dmatcsr_csr(A, "../output/A6.par.dat", 0);
    //if(is_print) Write_par_dmatcsr_csr(P, "../output/P6.par.dat", 0);
    //if(is_print) Write_par_dmatcsr_csr(R, "../output/R6.par.dat", 0);

    //if(is_print) printf("get coarse matrix...\n");
    //if(is_print) printf("doing first multi...\n");
    par_dmatcsr *AP = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
    Multi_par_dmatcsr_dmatcsr(A, P, AP);
    //if(is_print) Write_par_dmatcsr_csr(AP, "../output/AP.par.dat", 0);
    Remove_par_dmatcsr_extra_proc_neighbor(AP);
    //if(is_print) Print_par_dmatcsr(AP, 3);
    //if(is_print) printf("doing second multi...\n");
    MPI_Barrier(A->comm);
    MPI_Barrier(A->comm);
    MPI_Barrier(A->comm);
    Multi_par_dmatcsr_dmatcsr(R, AP, AH);
    Remove_par_dmatcsr_extra_proc_neighbor(AH);
    //if(is_print) Write_par_dmatcsr_csr(AH, "../output/AH.par.dat", 0);
    //if(is_print) printf("Done with coarse matrix...\n");
    
    //Truncate_P(P, param);
    
    if(NULL != M)
    {
	par_dmatcsr *MP = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
	Multi_par_dmatcsr_dmatcsr(M, P, MP);
	Remove_par_dmatcsr_extra_proc_neighbor(MP);
	Multi_par_dmatcsr_dmatcsr(R, MP, MH);
	Remove_par_dmatcsr_extra_proc_neighbor(MH);
	Free_par_dmatcsr(MP);
    }

    Free_par_dmatcsr(AP);
    Free_par_imatcsr(S);
    Free_par_ivec(dof);

    //printf("coarsening dond...\n");
    return SUCCESS;
}


#endif
