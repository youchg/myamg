#ifdef WITH_MPI

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "par_matrix_vector.h"
#include "par_linear_algebra.h"
#include "par_multigrid.h"
#include "par_eigen_solver.h"
#include "arpack_interface.h"
#include "par_linear_solver.h"
#include "linear_algebra.h"
#include "amg_param.h"
#include "tool.h"
#include "io.h"

extern int print_rank;

/* using 'dvec':
 * 1. to form    rhs 
 * 2. as initial approximation 
 * 3. as final   approximation
 * */
static double Correction_par_solve_linear(par_multigrid *pamg, int current_level, int n, double *dval, par_dvec **dvec, amg_param param);

static double Correction_par_expand_matrix_RAhV_VTAhV(par_multigrid *pamg, int current_level, int coarsest_level, 
				                  int n, par_dvec **V, int proc_root_coarsest, 
						  par_dmatcsr *Af, par_dmatcsr *Ac, 
						  dmatcsr *AH, dmatcsr *AL);

static double Correction_par_get_new_evec(par_multigrid *pamg, int current_level, int coarsest_level,
	                                  int n, par_dvec **V,  
					  par_dvec **evec_coarsest, double **evec_seq_coarsest_span,
					  par_dvec **evec);

/* 当 init_approx_level >= 0 时，eval, evec 必须给定初值，并且不能为 0;
   当 init_approx_level == 0 时，is_init_approx_level_correction 必须大于 0
*/
//Problem: when already have a finer result such as 0-level approximation, but
//want to correct it from a coarser level, there maybe some problem
//----20160301
void Eigen_solver_par_amg(par_multigrid *pamg, 
                          int nev, double *eval, 
                          par_dvec **evec, int init_approx_level, int is_init_approx_level_correction, 
                          amg_param param)
{
    int myrank;
    MPI_Comm_rank(pamg->comm[0], &myrank);

    double t1, t2, t3, t4;
    double eigen_time    = 0;
    double linear_time   = 0;
    double expand_time   = 0;
    double new_evec_time = 0;
    t1 = MPI_Wtime();

    int status = 0;

    int i, j, k, m;
    int finest_level   = 0;
    int nlevel         = pamg->actual_level_global;
    int niter_outer    = param.amgeigen_nouter_iter;
    int coarsest_level = param.amgeigen_coarsest_level;
    if(coarsest_level <= 0) coarsest_level += nlevel-1;

    /*
     * 将 coarsest level 汇总到 group[coarsest_level] 进程组的 root 进程
     */
    dmatcsr *AH = NULL;
    dmatcsr *MH = NULL;
    dmatcsr *Alarge = NULL;
    dmatcsr *Mlarge = NULL;

    // 最粗层：以 0 号进程为根进程；
    // 需要找到最粗层0号进程对应最细层的进程号，然后作为最细层的根进程
    int proc_root            = -1;
    int is_proc_root         = -1;
    int proc_root_current    = -1;
    int is_proc_root_current = -1;
    int proc_root_coarsest   = 0;
    if(MPI_COMM_NULL != pamg->comm[coarsest_level])
    {
	AH = (dmatcsr*)malloc(sizeof(dmatcsr));
	MH = (dmatcsr*)malloc(sizeof(dmatcsr));
	Get_dmatcsr_from_par_dmatcsr(pamg->A[coarsest_level], AH, proc_root_coarsest);
	Get_dmatcsr_from_par_dmatcsr(pamg->M[coarsest_level], MH, proc_root_coarsest);
    }

    double **evec_seq_coarsest_span = (double**)malloc(nev * sizeof(double*));
    for(j=0; j<nev; j++) evec_seq_coarsest_span[j] = (double*)calloc(nev, sizeof(double));

    par_dvec **evec_coarsest = (par_dvec**)malloc(nev * sizeof(par_dvec*));
    for(j=0; j<nev; j++) evec_coarsest[j] = NULL;

    double   **evec_seq_coarsest         = NULL;
    double   **evec_seq_coarsest_expand  = NULL;
    int       *evec_coarsest_comm_num    = NULL;
    int       *evec_coarsest_comm_displs = NULL;
    int        myrank_coarsest           = -1;
    int        nproc_coarsest            = 0;

    MPI_Comm comm_coarsest = pamg->comm[coarsest_level];
    if(MPI_COMM_NULL != comm_coarsest)
    {
	MPI_Comm_rank(comm_coarsest, &myrank_coarsest);
	MPI_Comm_size(comm_coarsest, &nproc_coarsest);
	if(myrank_coarsest == proc_root_coarsest)
	{
	    is_proc_root = myrank;
	}

	for(j=0; j<nev; j++) evec_coarsest[j] = Init_par_dvec_mv(pamg->A[coarsest_level]);

	evec_seq_coarsest = (double**)malloc(nev * sizeof(double*));
	for(j=0; j<nev; j++) evec_seq_coarsest[j] = NULL;

	evec_seq_coarsest_expand = (double**)malloc(nev * sizeof(double*));
	for(i=0; i<nev; i++) evec_seq_coarsest_expand[i] = NULL;

	if(myrank_coarsest == proc_root_coarsest)
	{
	    Alarge = Expand_dmatcsr_struct(AH, nev);
	    Mlarge = Expand_dmatcsr_struct(MH, nev);
	    for(i=0; i<nev; i++) evec_seq_coarsest_expand[i] = (double*)calloc(Alarge->nr, sizeof(double));

	    evec_coarsest_comm_num = (int*)calloc(nproc_coarsest, sizeof(int));
	    for(j=0; j<nproc_coarsest; j++) 
		evec_coarsest_comm_num[j] = pamg->A[coarsest_level]->row_start[j+1] - pamg->A[coarsest_level]->row_start[j];

	    evec_coarsest_comm_displs = (int*)calloc(nproc_coarsest, sizeof(int));
	    for(j=0; j<nproc_coarsest-1; j++) 
		evec_coarsest_comm_displs[j+1] = evec_coarsest_comm_displs[j] + evec_coarsest_comm_num[j];
	}
    }

    if(init_approx_level < 0)
    {
	if(MPI_COMM_NULL != comm_coarsest)
	{
	    if(myrank_coarsest == proc_root_coarsest)
	    {
		Write_dmatcsr_csr(AH, "../output/AH.dat");
		Write_dmatcsr_csr(MH, "../output/MH.dat");

		for(i=0; i<nev; i++) evec_seq_coarsest[i] = (double*)calloc(AH->nr, sizeof(double));

		for(i=0; i<nev; i++) eval[i] = 0.0;
		t3 = MPI_Wtime();
		Eigen_solver_arpack_dn(AH, MH, nev, eval, evec_seq_coarsest);
		for(j=0; j<nev; j++) Normalize_dvec(evec_seq_coarsest[j], AH->nr);
		t4 = MPI_Wtime();
		eigen_time += t4-t3;
		if(1)
		{
		    printf("Solving eigenvalue problem by Eigen_solver_arpack_dn() on the coarest level...\n");
		    printf("========= coarest eigenvalue on level %d =========\n", coarsest_level);
		    for(i=0; i<nev; i++) printf("%d: %15.12f\n", i, eval[i]);
		    printf("==================================================\n");
		}
	    }

	    //将 evec_seq_coarsest 散发到最粗层网格的各个进程
	    for(j=0; j<nev; j++)
	    {
		MPI_Scatterv(evec_seq_coarsest[j], evec_coarsest_comm_num, evec_coarsest_comm_displs, MPI_DOUBLE, 
			     evec[j]->value,       pamg->A[coarsest_level]->diag->nr,                 MPI_DOUBLE, 
			     proc_root_coarsest,   comm_coarsest);
	    }
	}
    }

    MPI_Allreduce(&is_proc_root, &proc_root, 1, MPI_INT, MPI_MAX, pamg->comm[0]);
    MPI_Bcast(eval, nev, MPI_DOUBLE, proc_root, pamg->comm[0]);
    //printf("rank %d: proc_root = %d\n", myrank, proc_root);
    //=====================================================================
    int amg_correction_begin_level;
    if(init_approx_level < 0)
	amg_correction_begin_level = coarsest_level -1;
    else if(is_init_approx_level_correction > 0)
	amg_correction_begin_level = init_approx_level;
    else
	amg_correction_begin_level = init_approx_level - 1;
    /*排除：在第0层上给定初值，但是不需要矫正。以后可以改成这种情况直接返回。*/
    assert(init_approx_level!=0 || is_init_approx_level_correction > 0);
    
    //================= amg method for eigenvalue problem =================
    par_dvec **dvec_amg = (par_dvec**)malloc(nev * sizeof(par_dvec*));
    for(i=0; i<nev; i++) dvec_amg[i] = Init_par_dvec_mv(pamg->A[finest_level]);

    //if(status) printf("Begin amg cycle...\n");
    for(i=amg_correction_begin_level; i>=finest_level; i--)
    {
	//if(status) printf("level = %d, nlevel = %d, coarsest_level = %d\n", i, nlevel, coarsest_level);
	for(m=0; m<niter_outer; m++)
	{
	    int myrank_current = -1;
	    if(MPI_COMM_NULL != pamg->comm[i])
	    {
		MPI_Barrier(pamg->comm[i]);
		MPI_Comm_rank(pamg->comm[i], &myrank_current);
		int nrank;
		MPI_Comm_size(pamg->comm[i], &nrank);
		//if(status) printf("correction: m = %d\n", m);
		//if(status) printf("Solving linear system on current level %d ...\n", i);
		for(j=0; j<nev; j++)
		{
		    if(m == 0)
		    {
			if(init_approx_level < 0)
			{
			    Multi_par_dmatcsr_dvec(pamg->P[i], evec[j], dvec_amg[j]);
			}
			else
			{
			    if(i==amg_correction_begin_level && is_init_approx_level_correction>0)
				Copy_par_dvec(evec[j], dvec_amg[j]);
			    else
				Multi_par_dmatcsr_dvec(pamg->P[i], evec[j], dvec_amg[j]);
			}
		    }
		    else
		    {
			Copy_par_dvec(evec[j], dvec_amg[j]);
		    }
		}

		//printf("rank %2d in %2d: solve linear start...\n", myrank, nrank);
		linear_time += Correction_par_solve_linear(pamg, i, nev, eval, dvec_amg, param);
		//printf("rank %2d in %2d: solve linear end...\n", myrank, nrank);
		//printf("rank %2d in %2d: expand matrix A...\n", myrank, nrank);
		expand_time += Correction_par_expand_matrix_RAhV_VTAhV(pamg, i, coarsest_level, nev, dvec_amg, proc_root_coarsest, pamg->A[i], pamg->A[coarsest_level], AH, Alarge);
		//printf("rank %2d in %2d: expand matrix M...\n", myrank, nrank);
		expand_time += Correction_par_expand_matrix_RAhV_VTAhV(pamg, i, coarsest_level, nev, dvec_amg, proc_root_coarsest, pamg->M[i], pamg->M[coarsest_level], MH, Mlarge);
		//printf("rank %d: expand matrix done...\n", myrank);

		if(myrank == proc_root)
		{
		    //if(status) printf("Solving corrected eigenvalue problem...\n");
		    printf("Solving corrected eigenvalue problem...\n");
		    t3 = MPI_Wtime(); 
		    memset(eval, 0, nev*sizeof(double));
		    for(j=0; j<nev; j++) memset(evec_seq_coarsest_expand[j], 0, Alarge->nr*sizeof(double));
		    if(i==0)
		    {
		    Write_dmatcsr_csr(Alarge, "../output/Alarge.dat");
		    Write_dmatcsr_csr(Mlarge, "../output/Mlarge.dat");
		    //exit(-1);
		    }

		    Eigen_solver_arpack_dn(Alarge, Mlarge, nev, eval, evec_seq_coarsest_expand);
		    //Print_dvec(eval, nev);
		    //Print_dvec(evec_seq_coarsest_expand[0], Alarge->nr);
		    Insertion_ascend_sort_dvec_dvecvec(eval, evec_seq_coarsest_expand, 0, nev-1);
		    t4 = MPI_Wtime();
		    eigen_time += t4-t3;

		    for(j=0; j<nev; j++)
		    {
			for(k=0; k<nev; k++)
			    evec_seq_coarsest_span[j][k] = evec_seq_coarsest_expand[j][k+AH->nr];
		    }
		    
		    if(status)
		    //if(1)
		    {
			printf("========= corrected eigenvalue on level %d ===========\n", i);
			for(j=0; j<nev; j++) printf("%15.12f\n", eval[j]);
			//printf("\n");
			printf("======================================================\n");
		    }

		    is_proc_root_current = myrank_current;
		}

		if(MPI_COMM_NULL != pamg->comm[coarsest_level])
		{
		    for(j=0; j<nev; j++)
		    {
			MPI_Scatterv(evec_seq_coarsest_expand[j], evec_coarsest_comm_num, evec_coarsest_comm_displs, MPI_DOUBLE, 
				     evec_coarsest[j]->value,     pamg->A[coarsest_level]->diag->nr,                 MPI_DOUBLE, 
				     proc_root_coarsest,          pamg->comm[coarsest_level]);
		    }
		}


	    MPI_Allreduce(&is_proc_root_current, &proc_root_current, 1, MPI_INT, MPI_MAX, pamg->comm[i]);
	    is_proc_root_current = -1;

	    for(j=0; j<nev; j++)
		MPI_Bcast(evec_seq_coarsest_span[j], nev, MPI_DOUBLE, proc_root_current, pamg->comm[i]); 
	    new_evec_time += Correction_par_get_new_evec(pamg, i, coarsest_level, nev, dvec_amg, evec_coarsest, evec_seq_coarsest_span, evec);
	    }
	    MPI_Bcast(eval, nev, MPI_DOUBLE, proc_root, pamg->comm[0]);
	}
    }


    if(status) printf("amg cycle end.\n");
    //=====================================================================
    
    for(i=0; i<nev; i++) {Free_par_dvec(dvec_amg[i]); dvec_amg[i] = NULL; }
    free(dvec_amg); dvec_amg = NULL;

    if(myrank_coarsest == proc_root_coarsest)
    {
	for(i=0; i<nev; i++) {free(evec_seq_coarsest_expand[i]); evec_seq_coarsest_expand[i] = NULL;}
	free(evec_seq_coarsest_expand); evec_seq_coarsest_expand = NULL;

	free(evec_coarsest_comm_num); evec_coarsest_comm_num = NULL;
	free(evec_coarsest_comm_displs); evec_coarsest_comm_displs = NULL;

	Free_dmatcsr(Alarge);
	Free_dmatcsr(Mlarge);

	Free_dmatcsr(AH);
	Free_dmatcsr(MH);
    }

    if(MPI_COMM_NULL != pamg->comm[coarsest_level])
    {
	for(i=0; i<nev; i++) {free(evec_seq_coarsest[i]); evec_seq_coarsest[i] = NULL;}
	free(evec_seq_coarsest); evec_seq_coarsest = NULL;

	for(j=0; j<nev; j++) Free_par_dvec(evec_coarsest[j]);

	if(myrank_coarsest != proc_root_coarsest)
	{
	    free(evec_seq_coarsest_expand); evec_seq_coarsest_expand = NULL;
	    free(AH); AH = NULL;
	    free(MH); MH = NULL;
	}
    }

    free(evec_coarsest); evec_coarsest = NULL;

    for(i=0; i<nev; i++) {free(evec_seq_coarsest_span[i]); evec_seq_coarsest_span[i] = NULL; }
    free(evec_seq_coarsest_span); evec_seq_coarsest_span = NULL;

    t2 = MPI_Wtime();
    if(param.amgeigen_print_level>0 && myrank==proc_root)
    {
	printf("direct eigen     time = %f\n", eigen_time);
        printf("amg linear solve time = %f\n", linear_time);
        printf("expand matrix    time = %f\n", expand_time);
        printf("get new evec     time = %f\n", new_evec_time);
        printf("correction total time = %f\n", t2-t1);
    }
}


void Eigen_solver_par_amg_nested(par_multigrid *pamg, 
                                 int nev, double *eval, par_dvec **evec,
                                 amg_param param)
{
    Eigen_solver_par_amg(pamg, nev, eval, evec, -1, 1, param);
}

static double Correction_par_solve_linear(par_multigrid *pamg, int current_level, int n, double *dval, par_dvec **dvec, amg_param param)
{
    double tb = MPI_Wtime();
    if(MPI_COMM_NULL != pamg->comm[current_level])
    {
	int myrank;
	MPI_Comm_rank(pamg->comm[current_level], &myrank);
	double resi_norm = -1.0;
	int    ncycle    = -1;
	if(0 == myrank)
	{
	    printf(".............. Par Eigen -- par linear solver amg .............\n");
	    printf("level   j   amgcycle            rn               time \n");
	}

	par_dmatcsr *A = pamg->A[current_level];
	par_dmatcsr *M = pamg->M[current_level];

	par_dvec *rhs = Init_par_dvec_mv(A);

	int j;
	for(j=0; j<n; j++)
	{
	    //memset(rhs, 0, M->nr*sizeof(double));
	    double tb = MPI_Wtime();
	    Multi_par_dmatcsr_dvec(M, dvec[j], rhs);
#if 0
	    int print_j = 1;
#if 0
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==0 && j==print_j) Print_dvec(rhs->value, M->diag->nr);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==1 && j==print_j) Print_dvec(rhs->value, M->diag->nr);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==2 && j==print_j) Print_dvec(rhs->value, M->diag->nr);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==3 && j==print_j) Print_dvec(rhs->value, M->diag->nr);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(j==print_j) printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
#endif
#if 1
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==0 && j==print_j) Print_dvec(dvec[j]->value, M->diag->nr);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==1 && j==print_j) Print_dvec(dvec[j]->value, M->diag->nr);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==2 && j==print_j) Print_dvec(dvec[j]->value, M->diag->nr);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==3 && j==print_j) Print_dvec(dvec[j]->value, M->diag->nr);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==0 && j == print_j) printf("rank %d: dval[j] = %f\n", myrank, dval[j]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==1 && j == print_j) printf("rank %d: dval[j] = %f\n", myrank, dval[j]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==2 && j == print_j) printf("rank %d: dval[j] = %f\n", myrank, dval[j]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==3 && j == print_j) printf("rank %d: dval[j] = %f\n", myrank, dval[j]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
#endif
#endif
	    if(MABS(dval[j]) > EPS) Scale_par_dvec(dvec[j], 1.0/dval[j]);
	    Linear_solver_par_amg(pamg, current_level, rhs, dvec[j], param, &resi_norm, &ncycle);
	    double te = MPI_Wtime();
	    if(myrank == 0) printf(" %2d    %2d     %3d       %18.15f     %f\n", current_level, j, ncycle, resi_norm, te-tb);
	}
	if(0 == myrank)
	    printf(".......................................................\n");

	Free_par_dvec(rhs);
    }

    return MPI_Wtime() - tb;
}

//expand_time += Correction_par_expand_matrix_RAhV_VTAhV(pamg, i, coarsest_level, nev, dvec_amg, proc_root_coarsest, pamg->A[i], pamg->A[coarsest_level], AH, Alarge);
static double Correction_par_expand_matrix_RAhV_VTAhV(par_multigrid *pamg, int current_level, int coarsest_level, 
				                  int n, par_dvec **V, int proc_root_coarsest, 
						  par_dmatcsr *Af, par_dmatcsr *Ac, 
						  dmatcsr *AH, dmatcsr *AL)
{
    int myrank_0, nrank_0;
    MPI_Comm_rank(pamg->comm[current_level], &myrank_0);
    MPI_Comm_size(pamg->comm[current_level], &nrank_0);

    double tb = MPI_Wtime();

    int j, k;

    int length = Af->diag->nr;

    //printf("rank %2d in %2d: begin to expand matrix...\n", myrank_0, nrank_0);
    par_dvec **AhV = (par_dvec**)malloc(n*sizeof(par_dvec*));
    for(j=0; j<n; j++) AhV[j] = Init_par_dvec_mv(Af);
    for(j=0; j<n; j++) Multi_par_dmatcsr_dvec(Af, V[j], AhV[j]); 
    //printf("rank %2d in %2d: multi par AhV...\n", myrank_0, nrank_0);
    //MPI_Barrier(pamg->comm[current_level]);

    par_dvec **RAhV = (par_dvec**)malloc(n*sizeof(par_dvec*));
    for(j=0; j<n; j++) RAhV[j] = NULL;
    //printf("rank %2d in %2d: malloc...\n", myrank_0, nrank_0);

    //printf("rank %2d in %2d: multi VTAhV...\n", myrank_0, nrank_0);
    double **VTAhV = (double**)malloc(n*sizeof(double*));
    for(j=0; j<n; j++) VTAhV[j] = (double*)calloc(n, sizeof(double));
    for(j=0; j<n; j++) for(k=0; k<n; k++) VTAhV[j][k] = Multi_par_dvec_dvec_length(V[j], AhV[k], length, pamg->comm[current_level]);

    double **RAhV_gatherv        = NULL;
    int     *RAhV_gatherv_num    = NULL;
    int     *RAhV_gatherv_displs = NULL;
    int      myrank_coarsest     = -1;
    int      nproc_coarsest      = -1;
    if(MPI_COMM_NULL != pamg->comm[current_level])
    {
	for(j=0; j<n; j++) RAhV[j] = Init_par_dvec_mv(Af);
    //printf("rank %d: expanding matrix 1111111111111 ...\n", myrank_0);
	for(j=0; j<n; j++) Restrict_par_f2c(pamg, current_level, coarsest_level, AhV[j], RAhV[j]);
    //printf("rank %d: expanding matrix 2222222222222 ...\n", myrank_0);
	//printf("rank %d: A->nr = %d\n", myrank_coarsest, pamg->A[coarsest_level]->diag->nr);
    }

    if(MPI_COMM_NULL != pamg->comm[coarsest_level])
    {
	MPI_Comm_rank(pamg->comm[coarsest_level], &myrank_coarsest);
	MPI_Comm_size(pamg->comm[coarsest_level], &nproc_coarsest);

        RAhV_gatherv        = (double**)malloc(n*sizeof(double*));
	for(j=0; j<n; j++) RAhV_gatherv[j] = NULL;

	if(myrank_coarsest == proc_root_coarsest)
	{
	    for(j=0; j<n; j++) RAhV_gatherv[j] = (double*)calloc(AH->nr, sizeof(double));

	    RAhV_gatherv_num = (int*)calloc(nproc_coarsest, sizeof(int));
	    for(j=0; j<nproc_coarsest; j++) 
		RAhV_gatherv_num[j] = Ac->row_start[j+1] - Ac->row_start[j];

	    RAhV_gatherv_displs = (int*)calloc(nproc_coarsest, sizeof(int));
	    for(j=0; j<nproc_coarsest-1; j++) 
		RAhV_gatherv_displs[j+1] = RAhV_gatherv_displs[j] + RAhV_gatherv_num[j];

	    //Print_ivec(RAhV_gatherv_num, nproc_coarsest);
	    //Print_ivec(RAhV_gatherv_displs, nproc_coarsest);
	}

	for(j=0; j<n; j++)
	{
	    MPI_Gatherv(RAhV[j]->value,     Ac->diag->nr,     MPI_DOUBLE, 
		        RAhV_gatherv[j],    RAhV_gatherv_num, RAhV_gatherv_displs, MPI_DOUBLE, 
			proc_root_coarsest, pamg->comm[coarsest_level]);
	}

	if(myrank_coarsest == proc_root_coarsest)
	{
	    Expand_dmatcsr(AL, n, RAhV_gatherv, VTAhV);

	    free(RAhV_gatherv_displs); RAhV_gatherv_displs = NULL;
	    free(RAhV_gatherv_num);    RAhV_gatherv_num    = NULL;
	}

	for(j=0; j<n; j++)
	{
	    free(RAhV_gatherv[j]);
	    RAhV_gatherv[j] = NULL;
	}
	free(RAhV_gatherv); 
	RAhV_gatherv = NULL;
    }

    if(MPI_COMM_NULL != pamg->comm[current_level])
	for(j=0; j<n; j++) Free_par_dvec(RAhV[j]); 
    free(RAhV); RAhV = NULL;
    for(j=0; j<n; j++) {free(VTAhV[j]); VTAhV[j] = NULL;} free(VTAhV); VTAhV = NULL;
    for(j=0; j<n; j++) Free_par_dvec( AhV[j]); 
    free( AhV);  AhV = NULL;

    //printf("rank %2d in %2d: done with expand matrix...\n", myrank_0, nrank_0);
    return MPI_Wtime() - tb;
}

static double Correction_par_get_new_evec(par_multigrid *pamg, int current_level, int coarsest_level,
	                                  int n, par_dvec **V,  
					  par_dvec **evec_coarsest, double **evec_seq_coarsest_span,
					  par_dvec **evec)
{
    double tb = MPI_Wtime();

    int j, k;
    if(MPI_COMM_NULL != pamg->comm[current_level])
    {
	int nr  = pamg->A[current_level]->diag->nr;
	int myrank;
	MPI_Comm_rank(pamg->comm[current_level], &myrank);
	int nrank;
	MPI_Comm_size(pamg->comm[current_level], &nrank);
	//int nrc = pamg->A[coarsest_level]->diag->nr;

	par_dvec **PV = (par_dvec**)malloc(n * sizeof(par_dvec*));
	for(j=0; j<n; j++) PV[j] = Init_par_dvec_mv(pamg->A[current_level]);

	for(j=0; j<n; j++)
	{
	    for(k=0; k<nr; k++) evec[j]->value[k] = 0.0;
#if 0
	    for(k=0; k<nrank; k++)
	    {
		MPI_Barrier(pamg->comm[current_level]);
		MPI_Barrier(pamg->comm[current_level]);
		MPI_Barrier(pamg->comm[current_level]);
		MPI_Barrier(pamg->comm[current_level]);
		if(myrank == k) if(evec_coarsest[j]!=NULL) Print_dvec(evec_coarsest[j]->value, nrc);
		MPI_Barrier(pamg->comm[current_level]);
		MPI_Barrier(pamg->comm[current_level]);
		MPI_Barrier(pamg->comm[current_level]);
		MPI_Barrier(pamg->comm[current_level]);
	    }
#endif
	    Prolong_par_c2f(pamg, coarsest_level, current_level, evec_coarsest[j], PV[j]);
	    for(k=0; k<n; k++) Sumself_par_dvec_axpby_length(V[k], evec_seq_coarsest_span[j][k], evec[j], 1.0, nr);
	    Sumself_par_dvec_axpby_length(PV[j], 1.0, evec[j], 1.0, nr);
	    //Print_dvec(evec[j]->value, nr);
#if 0
	    int print_j = 1;
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==0 && j==print_j) Print_dvec(evec[j]->value, nr);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==1 && j==print_j) Print_dvec(evec[j]->value, nr);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==2 && j==print_j) Print_dvec(evec[j]->value, nr);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    if(myrank==3 && j==print_j) Print_dvec(evec[j]->value, nr);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
	    MPI_Barrier(pamg->comm[current_level]);
#endif
	    Normalize_par_dvec_length(evec[j], nr, pamg->comm[current_level]);
	}
	for(j=0; j<n; j++) Free_par_dvec(PV[j]); 
	free(PV); PV = NULL;
    }

    return MPI_Wtime() - tb;
}


#endif
