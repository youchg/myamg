#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "preprocess.h"
#include "matrix.h"
#include "io.h"
#include "linear_algebra.h"
#include "tool.h"
#include "par_matrix_vector.h"
#include "par_cfsplit.h"
#include "par_linear_algebra.h"

#include "mpi.h"

int print_rank = 0;

int my_CLJP_split(imatcsr *S);

int print_num1 = 0;
int print_num2 = 0;

int main(int argc, char *argv[])
{
    amg_param param;
    Init_amg_param(&param);

    MPI_Init(&argc, &argv);

    int  myrank;
    int nproc_global;
    //int nproc_global, myname_len;
    //char myname[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_size(MPI_COMM_WORLD, &nproc_global);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    //MPI_Get_processor_name(myname, &myname_len);

    //char file[256] = "../../dat/fem3d/hydrogen-stiff-4913.dat";
    //char file[256] = "../../dat/fem2d_poisson_square/gmg_A_refine5.m";
    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine5.m";
    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine4.m";
    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine3.m";
    //char file[256] = "../../dat/fdm2d9pt/A_fdm9pt_6400x6400.dat";
    char file[256] = "../dat/fdm2d9pt/A_fdm9pt_49x49.dat";
    par_dmatcsr *A = Read_par_dmatcsr(file, MPI_COMM_WORLD);
    //Print_par_dmatcsr(A, 3);

    par_imatcsr *S = (par_imatcsr*)malloc(sizeof(par_imatcsr));
    Generate_par_strong_coupling_set(A, S, param);
    par_ivec *dof = Init_par_ivec_length_comm(A->diag->nr, A->comm);
    Split_par_CLJP(A, S, dof);

    int *ncpt_proc = (int*)calloc(nproc_global, sizeof(int));
    par_dmatcsr *P = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
    Get_par_interpolation_direct(A, S, dof, P, ncpt_proc);
    Write_par_dmatcsr_csr(P, "../output/P.par.dat", 0);
    //Print_par_dmatcsr(P, 3);
    int i;
    int nproc_cpt = 0;
    for(i=0; i<nproc_global; i++) if(ncpt_proc[i] > 0) nproc_cpt++;
    int *proc_coarse = (int*)malloc(nproc_cpt* sizeof(int));
    nproc_cpt = 0;
    for(i=0; i<nproc_global; i++) if(ncpt_proc[i] > 0) proc_coarse[nproc_cpt++] = i;

    par_dmatcsr *R = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
    Transpose_par_dmatcsr(P, R);
    Print_par_dmatcsr(R, 3);
    Write_par_dmatcsr_csr(R, "../output/R.par.dat", 0);

    if(myrank == print_rank) printf("Begin to perfom Multi_par A*P...\n");
    par_dmatcsr *AP = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
    MPI_Barrier(MPI_COMM_WORLD);
    double tb_multi_par = MPI_Wtime();
    Multi_par_dmatcsr_dmatcsr(A, P, AP);
    MPI_Barrier(MPI_COMM_WORLD);
    double te_multi_par = MPI_Wtime();
    if(myrank == print_rank) printf("Multi_par A*P time: %f\n", te_multi_par-tb_multi_par);
    if(myrank == print_rank) Write_dmatcsr_csr(AP->offd, "../output/APoffd.dat");
    //Print_par_dmatcsr(AP, 3);
    MPI_Barrier(MPI_COMM_WORLD);
    Remove_par_dmatcsr_extra_proc_neighbor(AP);
    //Print_par_dmatcsr(AP, 3);
    //Write_par_dmatcsr_csr(AP, "../output/AP.par.dat", 0);
    //Write_par_dmatcsr_csr(P, "../output/P.par.dat", 0);
    
    par_dmatcsr *AH = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
    MPI_Barrier(MPI_COMM_WORLD);
    tb_multi_par = MPI_Wtime();
    Multi_par_dmatcsr_dmatcsr(R, AP, AH);
    MPI_Barrier(MPI_COMM_WORLD);
    te_multi_par = MPI_Wtime();
    Write_par_dmatcsr_csr(AH, "../output/AH.par.dat", 0);
    if(myrank == print_rank) printf("Multi_par R*A*P time: %f\n", te_multi_par-tb_multi_par);
    MPI_Barrier(MPI_COMM_WORLD);
    Remove_par_dmatcsr_extra_proc_neighbor(AH);
    //Print_par_dmatcsr(AH, 3);
    
    MPI_Group group_global, group_coarse;
    MPI_Comm_group(MPI_COMM_WORLD, &group_global);
    MPI_Group_incl(group_global, nproc_cpt, proc_coarse, &group_coarse);

    int group_coarse_size;
    int group_coarse_rank;
    MPI_Group_size(group_coarse, &group_coarse_size);
    MPI_Group_rank(group_coarse, &group_coarse_rank);
#if 0
    int is_proc_coarse = ncpt_proc[myrank];
    if(is_proc_coarse > 0)
	printf("group_coarse_size = %d, group_coarse_rank %d.\n", group_coarse_size, group_coarse_rank);
#endif

    MPI_Comm comm_coarse;
    MPI_Comm_create(MPI_COMM_WORLD, group_coarse, &comm_coarse);
    int nproc_coarse, myrank_coarse;
    if(MPI_COMM_NULL != comm_coarse)
    {
	MPI_Comm_size(comm_coarse, &nproc_coarse);
	MPI_Comm_rank(comm_coarse, &myrank_coarse);
	//printf("myrank_coarse %d nproc_coarse %d\n", myrank_coarse, nproc_coarse);

	int *map_proc_f2c = (int*)malloc(nproc_global * sizeof(int));
	for(i=0; i<nproc_global; i++) map_proc_f2c[i] = -1;
	int index = 0;
	for(i=0; i<nproc_global; i++) if(ncpt_proc[i] > 0) map_proc_f2c[i] = index++;
	assert(index == nproc_cpt);

	//if(myrank_coarse == print_rank) Print_ivec(map_proc_f2c, nproc_global);

	par_comm_info *AH_comm_info = AH->comm_info;
	for(i=0; i<AH_comm_info->nproc_neighbor; i++)
	{
	    assert(map_proc_f2c[AH_comm_info->proc_neighbor[i]] >= 0);
	    AH_comm_info->proc_neighbor[i] = map_proc_f2c[AH_comm_info->proc_neighbor[i]];
	}

	index = 0;
	for(i=0; i<nproc_global; i++)
	{
	    if(ncpt_proc[i] > 0)
	    {
		AH->row_start[index]   = AH->row_start[i];
		AH->row_start[index+1] = AH->row_start[i+1];
		AH->col_start[index]   = AH->col_start[i];
		AH->col_start[index+1] = AH->col_start[i+1];
		index++;
	    }
	}
	AH->comm = comm_coarse;
	//Print_par_dmatcsr(AH, 3);
	
	par_imatcsr *SH = (par_imatcsr*)malloc(sizeof(par_imatcsr));
	Generate_par_strong_coupling_set(AH, SH, param);
	par_ivec *dofH = Init_par_ivec_length_comm(AH->diag->nr, AH->comm);
	Split_par_CLJP(AH, SH, dofH);

	int *ncpt_proc_coarse = (int*)calloc(nproc_coarse, sizeof(int));
	par_dmatcsr *PH = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
	Get_par_interpolation_direct(AH, SH, dofH, PH, ncpt_proc_coarse);

	free(ncpt_proc_coarse); ncpt_proc_coarse = NULL;
	Free_par_dmatcsr(PH);
	Free_par_ivec(dofH);
	Free_par_imatcsr(SH);
	free(map_proc_f2c); map_proc_f2c = NULL;
    }

    if(MPI_COMM_NULL != comm_coarse) MPI_Comm_free(&comm_coarse);
    MPI_Group_free(&group_coarse);
    MPI_Group_free(&group_global);


    Free_par_dmatcsr(AH);

    Free_par_dmatcsr(R);
    Free_par_dmatcsr(AP);

    free(proc_coarse); proc_coarse = NULL;
    free(ncpt_proc); ncpt_proc = NULL;
    Free_par_dmatcsr(P);
    Free_par_ivec(dof);
    Free_par_imatcsr(S);

    print_num1 = A->col_start[myrank];
    print_num2 = A->col_start[myrank+1];

    Free_par_dmatcsr(A);

#if 0
    if(myrank == print_rank)
    {
	fflush(stdout);
	//printf("print from %d to %d\n", print_num1, print_num2-1);

	dmatcsr *AA = Read_dmatcsr(file);
	imatcsr *SS = (imatcsr*)malloc(sizeof(imatcsr));
	int *dof = (int*)calloc(AA->nr, sizeof(int));
	Generate_strong_coupling_set(AA, SS, param);
	CLJP_split(AA, SS, dof);
	dmatcsr *PP = (dmatcsr*)malloc(sizeof(dmatcsr));
	Generate_sparsity_P_dir(SS, dof, PP);
	Generate_P_dir(AA, SS, dof, PP);

	dmatcsr *RR = (dmatcsr*)malloc(sizeof(dmatcsr));
	Transpose_dmatcsr(PP, RR);

	//Write_dmatcsr_csr(AA, "../output/A.dat");
	Write_dmatcsr_csr(PP, "../output/P.dat");
	Write_dmatcsr_csr(RR, "../output/R.dat");
	dmatcsr *AAPP = (dmatcsr*)malloc(sizeof(dmatcsr));

	double tb_multi = Get_time();
	Multi_dmatcsr_dmatcsr(AA, PP, AAPP);
	double te_multi = Get_time();
	printf("Multi A*P time: %f\n", te_multi-tb_multi);

	dmatcsr *RRAAPP = (dmatcsr*)malloc(sizeof(dmatcsr));
	tb_multi = Get_time();
	Multi_dmatcsr_dmatcsr(RR, AAPP, RRAAPP);
	te_multi = Get_time();
	printf("Multi R*A*P time: %f\n", te_multi-tb_multi);
	Write_dmatcsr_csr(RRAAPP, "../output/RRAAPP.dat");

	Write_dmatcsr_csr(AAPP, "../output/AAPP.dat");
	//Write_imatcsr_csr(SS, "../output/S.dat");
	//my_CLJP_split(SS);
	free(dof);
	Free_dmatcsr(RRAAPP);
	Free_dmatcsr(RR);
	Free_dmatcsr(PP);
	Free_imatcsr(SS);
	Free_dmatcsr(AA);
	Free_dmatcsr(AAPP);
    }
#endif

    MPI_Finalize();

    return 0;
}

#if 0
int my_CLJP_split(imatcsr *S)
{
    //printf("CLJP spliting...\n");

    int i;
    
    imatcsr *ST = (imatcsr*)malloc(sizeof(imatcsr));
    Transpose_imatcsr_struct(S, ST);
    
    int  ST_nr = ST->nr;
    int *ST_ia = ST->ia;
    
    srand(1);
    int *lambda_ST = (int*)malloc(ST_nr*sizeof(int));
    for(i=0; i<ST_nr; i++) lambda_ST[i] = ST_ia[i+1]-ST_ia[i];

    //printf("print from %d to %d\n", print_num1, print_num2-1);
    //Print_ivec(lambda_ST+print_num1, print_num2-print_num1);


    free(lambda_ST);
    Free_imatcsr(ST);

    return 1;
}
#endif
