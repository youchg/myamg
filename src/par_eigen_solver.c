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

#define par_amg_eigen_linear_solver 1 //1:amgsolver 2:pcg_amg

extern int print_rank;

/* using 'dvec':
 * 1. to form    rhs 
 * 2. as initial approximation 
 * 3. as final   approximation
 * */
static double Correction_par_solve_linear(par_multigrid *pamg, int current_level, int n, double *dval, par_dvec **dvec, amg_param param);

static double Correction_par_expand_matrix_RAhV_VTAhV(par_multigrid *pamg, int current_level, int coarsest_level, 
				                  int n, par_dvec **V, int proc_root_coarsest, 
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

    int proc_root_coarsest = 0;
    if(MPI_COMM_NULL != pamg->comm[coarsest_level])
    {
	AH = (dmatcsr*)malloc(sizeof(dmatcsr));
	MH = (dmatcsr*)malloc(sizeof(dmatcsr));
	Get_dmatcsr_from_par_dmatcsr(pamg->A[coarsest_level], AH, proc_root_coarsest);
	Get_dmatcsr_from_par_dmatcsr(pamg->M[coarsest_level], MH, proc_root_coarsest);
    }

    double **evec_seq_coarsest_span = NULL;
    for(j=0; j<nev; j++) evec_seq_coarsest_span[j] = (double*)calloc(nev, sizeof(double));

    par_dvec **evec_coarsest             = NULL;
    double   **evec_seq_coarsest         = NULL;
    double   **evec_seq_coarsest_expand  = NULL;
    int       *evec_coarsest_comm_num    = NULL;
    int       *evec_coarsest_comm_displs = NULL;
    int        myrank_coarsest           = -1;
    int        nproc_coarsest            = 0;
    if(init_approx_level < 0)
    {
	if(MPI_COMM_NULL != pamg->comm[coarsest_level])
	{
	    MPI_Comm comm_coarsest = pamg->comm[coarsest_level];
	    MPI_Comm_rank(comm_coarsest, &myrank_coarsest);
	    MPI_Comm_size(comm_coarsest, &nproc_coarsest);

	    for(j=0; j<nev; j++) evec_coarsest[j] = Init_par_dvec_mv(pamg->A[coarsest_level]);

	    if(myrank_coarsest == proc_root_coarsest)
	    {
		evec_seq_coarsest = (double**)malloc(nev * sizeof(double*));
		for(i=0; i<nev; i++) evec_seq_coarsest[i] = (double*)calloc(AH->nr, sizeof(double));

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

		Alarge = Expand_dmatcsr_struct(AH, nev);
		Mlarge = Expand_dmatcsr_struct(MH, nev);
		evec_seq_coarsest_expand = (double**)malloc(nev * sizeof(double*));
		for(i=0; i<nev; i++) evec_seq_coarsest_expand[i] = (double*)calloc(Alarge->nr, sizeof(double));

		evec_coarsest_comm_num = (int*)calloc(nproc_coarsest, sizeof(int));
		for(j=0; j<nproc_coarsest; j++) 
		    evec_coarsest_comm_num[j] = pamg->A[coarsest_level]->row_start[j+1] - pamg->A[coarsest_level]->row_start[j];

		evec_coarsest_comm_displs = (int*)calloc(nproc_coarsest, sizeof(int));
		for(j=0; j<nproc_coarsest-1; j++) 
		    evec_coarsest_comm_displs[j+1] = evec_coarsest_comm_displs[j] + evec_coarsest_comm_num[j];
	    }

	    //将 evec_seq_coarsest 散发到最粗层网格的各个进程
	    for(j=0; j<nev; j++)
	    {
		MPI_Scatterv(evec_seq_coarsest[j], evec_coarsest_comm_num, evec_coarsest_comm_displs, MPI_DOUBLE, 
			     evec,                 pamg->A[coarsest_level]->diag->nr,                    MPI_DOUBLE, 
			     myrank_coarsest,      comm_coarsest);
	    }

	}
    }

    int proc_root = -1;
    if(myrank_coarsest == proc_root_coarsest) proc_root = myrank;
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
	    if(MPI_COMM_NULL != pamg->comm[i])
	    {
		MPI_Barrier(pamg->comm[i]);
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

		linear_time += Correction_par_solve_linear(pamg, i, nev, eval, dvec_amg, param);
		expand_time += Correction_par_expand_matrix_RAhV_VTAhV(pamg, i, coarsest_level, nev, dvec_amg, proc_root_coarsest, AH, Alarge);
		expand_time += Correction_par_expand_matrix_RAhV_VTAhV(pamg, i, coarsest_level, nev, dvec_amg, proc_root_coarsest, MH, Mlarge);
		//Write_dmatcsr_csr(Mlarge, "../output/Mlarge.dat");

		if(myrank == proc_root)
		{
		    if(status) printf("Solving corrected eigenvalue problem...\n");
		    t3 = MPI_Wtime(); 
		    memset(eval, 0, nev*sizeof(double));
		    for(j=0; j<nev; j++) memset(evec_seq_coarsest_expand[j], 0, Alarge->nr*sizeof(double));
		    Eigen_solver_arpack_dn(Alarge, Mlarge, nev, eval, evec_seq_coarsest_expand);
		    Insertion_ascend_sort_dvec_dvecvec(eval, evec_seq_coarsest_expand, 0, nev-1);
		    t4 = MPI_Wtime();
		    eigen_time += t4-t3;

		    for(j=0; j<nev; j++)
		    {
			for(k=0; k<nev; k++)
			    evec_seq_coarsest_span[j][k] = evec_seq_coarsest_expand[j][k+AH->nr];
		    }
		    
		    if(status)
		    {
			printf("========= corrected eigenvalue on level %d ===========\n", i);
			for(j=0; j<nev; j++) printf("%15.12f\n", eval[j]);
			printf("\n");
			printf("======================================================\n");
		    }
		}

		if(MPI_COMM_NULL != pamg->comm[coarsest_level])
		{
		    for(j=0; j<nev; j++)
		    {
			MPI_Scatterv(evec_seq_coarsest_expand[j], evec_coarsest_comm_num, evec_coarsest_comm_displs, MPI_DOUBLE, 
				     evec_coarsest,               pamg->A[coarsest_level]->diag->nr,                 MPI_DOUBLE, 
				     proc_root_coarsest,          pamg->comm[coarsest_level]);
		    }
		}

	    }

	    for(j=0; j<nev; j++)
		MPI_Bcast(evec_seq_coarsest_span[j], nev, MPI_DOUBLE, proc_root, pamg->comm[i]); 
	    new_evec_time += Correction_par_get_new_evec(pamg, i, coarsest_level, nev, dvec_amg, evec_coarsest, evec_seq_coarsest_span, evec);
	}
    }

    MPI_Bcast(eval, nev, MPI_DOUBLE, proc_root, pamg->comm[0]);

    if(status) printf("amg cycle end.\n");
    //=====================================================================
    
    for(i=0; i<nev; i++) {Free_par_dvec(dvec_amg[i]); dvec_amg[i] = NULL; }
    free(dvec_amg); dvec_amg = NULL;

    if(myrank_coarsest == proc_root_coarsest)
    {
	for(i=0; i<nev; i++) {free(evec_seq_coarsest[i]); evec_seq_coarsest[i] = NULL;}
	free(evec_seq_coarsest); evec_seq_coarsest = NULL;

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
	for(j=0; j<nev; j++) Free_par_dvec(evec_coarsest[j]);
	free(evec_coarsest); evec_coarsest = NULL;
    }

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
#if par_amg_eigen_linear_solver == 1
    //double resi_norm_amg      = -1;
    //int    ncycle             =  0;
    //printf(".............. Par Eigen -- par linear solver amg .............\n");
    //printf("level   j   amgcycle            rn               time \n");
#elif par_amg_eigen_linear_solver == 2
    double resi_norm_pcg      = -1;
    int    ncycle             =  0;
    int    niter              =  0;
    printf(".................. Par Eigen -- par linear solver pcg ..................\n");
    printf("level   j   pcgiter   amgcycle            rn              time \n");
#endif

    par_dmatcsr *A = pamg->A[current_level];
    par_dmatcsr *M = pamg->M[current_level];

    par_dvec *rhs = Init_par_dvec_mv(A);

    //double t1, t2;
    int j;
    for(j=0; j<n; j++)
    {
	//t1 = Get_time();
	//memset(rhs, 0, M->nr*sizeof(double));
	Multi_par_dmatcsr_dvec(M, dvec[j], rhs);
#if 0
	Scale_par_dvec(rhs, dval[j]);
#else
	if(MABS(dval[j]) > EPS) Scale_par_dvec(dvec[j], 1.0/dval[j]);
#endif

#if   par_amg_eigen_linear_solver == 1
	Linear_solver_par_amg(pamg, current_level, rhs, dvec[j], param, NULL, NULL);
	//Linear_solver_par_amg(pamg, current_level, rhs, dvec[j], param, &resi_norm_amg, &ncycle);
	//t2 = Get_time();
	//printf(" %2d    %2d     %3d       %18.15f     %f\n", current_level, j, ncycle, resi_norm_amg, t2-t1);
#elif par_amg_eigen_linear_solver == 2
	Linear_solver_par_pcg_amg(pamg, current_level, rhs, dvec[j], param, NULL, &resi_norm_pcg, &niter, &ncycle);
	//t2 = Get_time();
	//printf(" %2d    %2d     %3d      %3d       %18.15f     %f\n", current_level, j, niter, ncycle, resi_norm_pcg, t2-t1);
#endif
    }

#if   par_amg_eigen_linear_solver == 1
    //printf(".......................................................\n");
#elif par_amg_eigen_linear_solver == 2
    //printf("................................................................\n");
#endif

    Free_par_dvec(rhs);

    return MPI_Wtime() - tb;
}

static double Correction_par_expand_matrix_RAhV_VTAhV(par_multigrid *pamg, int current_level, int coarsest_level, 
				                  int n, par_dvec **V, int proc_root_coarsest, 
						  dmatcsr *AH, dmatcsr *AL)
{
    double tb = MPI_Wtime();

    int j, k;

    par_dvec **AhV = (par_dvec**)malloc(n*sizeof(par_dvec*));
    for(j=0; j<n; j++) AhV[j] = Init_par_dvec_mv(pamg->A[current_level]);
    for(j=0; j<n; j++) Multi_par_dmatcsr_dvec(pamg->A[current_level], V[j], AhV[j]); 

    par_dvec **RAhV = (par_dvec**)malloc(n*sizeof(par_dvec*));
    for(j=0; j<n; j++) RAhV[j] = Init_par_dvec_mv(pamg->A[coarsest_level]);
    for(j=0; j<n; j++) Restrict_par_f2c(pamg, current_level, coarsest_level, AhV[j], RAhV[j]);

    double **VTAhV = (double**)malloc(n*sizeof(double*));
    for(j=0; j<n; j++) VTAhV[j] = (double*)calloc(n, sizeof(double));
    for(j=0; j<n; j++) for(k=0; k<n; k++) VTAhV[j][k] = Multi_par_dvec_dvec(V[j], AhV[k]);

    double **RAhV_gatherv        = NULL;
    int     *RAhV_gatherv_num    = NULL;
    int     *RAhV_gatherv_displs = NULL;
    int      myrank_coarsest     = -1;
    int      nproc_coarsest      = -1;
    if(MPI_COMM_NULL != pamg->comm[coarsest_level])
    {
	MPI_Comm_rank(pamg->comm[coarsest_level], &myrank_coarsest);
	MPI_Comm_size(pamg->comm[coarsest_level], &nproc_coarsest);

	if(myrank_coarsest == proc_root_coarsest)
	{
	    RAhV_gatherv = (double**)malloc(n*sizeof(double*));
	    for(j=0; j<n; j++) RAhV_gatherv[j] = (double*)calloc(AH->nr, sizeof(double));

	    RAhV_gatherv_num = (int*)calloc(nproc_coarsest, sizeof(int));
	    for(j=0; j<nproc_coarsest; j++) 
		RAhV_gatherv_num[j] = pamg->A[coarsest_level]->row_start[j+1] - pamg->A[coarsest_level]->row_start[j];

	    RAhV_gatherv_displs = (int*)calloc(nproc_coarsest, sizeof(int));
	    for(j=0; j<nproc_coarsest-1; j++) 
		RAhV_gatherv_displs[j+1] = RAhV_gatherv_displs[j] + RAhV_gatherv_num[j];
	}

	for(j=0; j<n; j++)
	{
	    MPI_Gatherv(RAhV[j],            pamg->A[coarsest_level]->diag->nr,     MPI_DOUBLE, 
		        RAhV_gatherv[j],    RAhV_gatherv_num, RAhV_gatherv_displs, MPI_DOUBLE, 
			proc_root_coarsest, pamg->comm[coarsest_level]);
	}

	if(myrank_coarsest == proc_root_coarsest)
	{
	    Expand_dmatcsr(AL, n, RAhV_gatherv, VTAhV);

	    for(j=0; j<nproc_coarsest; j++)
	    {
		free(RAhV_gatherv[j]);        RAhV_gatherv[j]        = NULL;
	    }
	    free(RAhV_gatherv_displs); RAhV_gatherv_displs = NULL;
	    free(RAhV_gatherv_num);    RAhV_gatherv_num    = NULL;
	    free(RAhV_gatherv);        RAhV_gatherv        = NULL;
	}
    }

    for(j=0; j<n; j++) {free(VTAhV[j]); VTAhV[j] = NULL;} free(VTAhV); VTAhV = NULL;

    for(j=0; j<n; j++) Free_par_dvec(RAhV[j]); free(RAhV); RAhV = NULL;
    for(j=0; j<n; j++) Free_par_dvec( AhV[j]); free( AhV);  AhV = NULL;

    return MPI_Wtime() - tb;
}

static double Correction_par_get_new_evec(par_multigrid *pamg, int current_level, int coarsest_level,
	                                  int n, par_dvec **V,  
					  par_dvec **evec_coarsest, double **evec_seq_coarsest_span,
					  par_dvec **evec)
{
    double tb = MPI_Wtime();

    int j, k;
    //int nr  = pamg->A[current_level]->diag->nr;
    int nrc = pamg->A[coarsest_level]->diag->nr;

    par_dvec **PV = (par_dvec**)malloc(n * sizeof(par_dvec*));
    for(j=0; j<n; j++) PV[j] = Init_par_dvec_mv(pamg->A[current_level]);

    for(j=0; j<n; j++)
    {
	for(k=0; k<evec[j]->length; k++) evec[j]->value[k] = 0.0;
	Prolong_par_c2f(pamg, coarsest_level, current_level, evec_coarsest[j], PV[j]);
	for(k=0; k<n; k++) Sumself_par_dvec_axpby(V[k], evec_seq_coarsest_span[j][nrc+k], evec[j], 1.0);
	Sumself_par_dvec_axpby(PV[j], 1.0, evec[j], 1.0);
        Normalize_par_dvec(evec[j]);
    }
    for(j=0; j<n; j++) Free_par_dvec(PV[j]); 
    free(PV); PV = NULL;

    return MPI_Wtime() - tb;
}


#endif
