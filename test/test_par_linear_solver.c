#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "preprocess.h"
#include "io.h"
#include "tool.h"

#include "linear_solver.h"
#include "linear_algebra.h"
#include "matrix.h"
#include "multigrid.h"
#include "setup_phase.h"

#include "par_matrix_vector.h"
#include "par_linear_algebra.h"
#include "par_linear_solver.h"
#include "par_multigrid.h"
#include "par_setup_phase.h"

#define cg        0
#define amgcycle  1
#define amgpcg    0

int print_rank = 0;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int  myrank;
    int nproc_global;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc_global);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    amg_param param;
    Init_amg_param(&param);
    param.linear_solver_base_type = CG;
    param.coarsening_type  = CLJP;
    param.cg_max_iter      = 20;
    param.amgsolver_tol    = 1e-08;
    param.max_coarsest_dof = 1;
    param.amgsolver_print_level = 1;
    //param.amgsolver_max_cycle = 1;

    //char file[256] = "../../dat/fem2d_poisson_square/gmg_A_refine8.m";
    //char file[256] = "../../dat/fem3d/hydrogen-stiff-4913.dat";
    //char file[256] = "../../dat/fem2d_poisson_square/gmg_A_refine5.m";
    char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine5.m";
    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine4.m";
    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine3.m";
    //char file[256] = "../../dat/fdm2d9pt/A_fdm9pt_122500x122500.dat";
    //char file[256] = "../../dat/fdm2d9pt/A_fdm9pt_6400x6400.dat";
    //char file[256] = "../dat/fdm2d9pt/A_fdm9pt_49x49.dat";

    int i;

    double tb_ini = MPI_Wtime();
    par_dmatcsr *A = Read_par_dmatcsr(file, MPI_COMM_WORLD);

    par_dvec *sol   = Init_par_dvec_mv(A);
    par_dvec *b     = Init_par_dvec_mv(A);
    par_dvec *x_ini = Init_par_dvec_mv(A);

    srand((unsigned)time(0));
    for(i=0; i<sol->length; i++) sol->value[i] = 1.0;//((double)rand())/((double)RAND_MAX);
    Multi_par_dmatcsr_dvec(A, sol, b);
    double te_ini = MPI_Wtime();
    if(myrank == print_rank) printf("Initial and read time = %f\n", te_ini-tb_ini);
    
#if cg
    /*------------------ cg ------------------*/
    if(myrank == print_rank)
	printf("\ncalling cg method...\n");
    par_dvec *x_cg = Init_par_dvec_mv(A);
    Copy_par_dvec(x_ini, x_cg);
    double tol_cg       = param.cg_tol;
    int    max_iter_cg  = param.cg_max_iter;
    double resi_norm_cg = 9999.9;
    int    niter_cg     = -1;
    double tb_cg = MPI_Wtime();
    Linear_solver_par_cg(A, b, x_cg, tol_cg, max_iter_cg, NULL, &resi_norm_cg, &niter_cg, 1);
    double te_cg = MPI_Wtime();
    Get_par_residual(A, b, x_cg, NULL, &resi_norm_cg);
    if(myrank == print_rank)
	printf("CG method: resi = %18.15f, niter = %3d, time = %f\n", resi_norm_cg, niter_cg, te_cg-tb_cg);
    Free_par_dvec(x_cg);
    /*----------------------------------------*/
#endif

#if amgcycle
    /*------------------ MYAMG ------------------*/
    if(myrank == print_rank)
	printf("\ncalling AMG method...\n");
    par_dvec *x_amg = Init_par_dvec_mv(A);
    Copy_par_dvec(x_ini, x_amg);
    double resi_norm_amg = 9999.9;
    int    ncycle_amg    = -1;
    
    amg_param param_myamg = param;

    MPI_Group mpi_group_world;
    MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_world);
    par_multigrid *pamg = Build_par_amg(A, NULL, param_myamg.max_level, MPI_COMM_WORLD, mpi_group_world);
    //par_multigrid *pamg = Build_par_amg(A, NULL, 9, MPI_COMM_WORLD, mpi_group_world);
    double tb_setup = MPI_Wtime();
    Setup_par_phase(pamg, param_myamg);
    double te_setup = MPI_Wtime();
    if(myrank == print_rank)
	printf("setup phase time: %f\n\n", te_setup-tb_setup);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //Print_par_dmatcsr(pamg->A[3], 0);
    //Write_par_dmatcsr_csr(pamg->A[0], "../output/A0_4.par.dat", 0);
    //Write_par_dmatcsr_csr(pamg->A[1], "../output/A1_4.par.dat", 0);
    //Write_par_dmatcsr_csr(pamg->A[2], "../output/A2_4.par.dat", 0);
    //Write_par_dmatcsr_csr(pamg->A[3], "../output/A3_4.par.dat", 0);
    //Write_par_dmatcsr_csr(pamg->P[0], "../output/P0_4.par.dat", 0);
    //Write_par_dmatcsr_csr(pamg->P[1], "../output/P1_4.par.dat", 0);
    //Write_par_dmatcsr_csr(pamg->P[2], "../output/P2_4.par.dat", 0);
    //Write_par_dmatcsr_csr(pamg->R[0], "../output/R0_4.par.dat", 0);
    //Write_par_dmatcsr_csr(pamg->R[1], "../output/R1_4.par.dat", 0);
    //Write_par_dmatcsr_csr(pamg->R[2], "../output/R2_4.par.dat", 0);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    Print_par_amg(pamg);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //if(MPI_COMM_NULL != pamg->comm[pamg->actual_level-1]) Write_par_dmatcsr_csr(pamg->A[pamg->actual_level-1], "../output/AH.par.dat", 0);

    //if(MPI_COMM_NULL != pamg->comm[pamg->actual_level-1]) printf("rank %d: coarsest nr = %d\n", myrank, pamg->A[pamg->actual_level-1]->nr_global);

    double tb_amg = MPI_Wtime();
    Linear_solver_par_amg(pamg, 0, b, x_amg, param_myamg, &resi_norm_amg, &ncycle_amg);
    double te_amg = MPI_Wtime();
    if(myrank == print_rank)
	printf("amg method: \nresi = %18.15f, ncycle = %3d, \ntime = %f, setup time = %f\n", resi_norm_amg, ncycle_amg, te_amg-tb_amg, te_setup-tb_setup);
    Free_par_multigrid(pamg);
    Free_par_dvec(x_amg);
    /*-----------------------------------------*/
#endif

#if amgpcg
    /*------------------ PCG+AMG ------------------*/
    printf("\ncalling pcg amg method...\n");
    double *x_pcg_amg = (double*)calloc(size, sizeof(double));
    memcpy(x_pcg_amg, x_ini, size*sizeof(double));
    
    double resi_norm_pcg_amg = 9999.9;
    int    niter_pcg         = -1;
    int    ncycle            = -1;
    
    amg_param param_pcg = param;
    param_pcg.amgsolver_max_cycle = 5;
    
    double tb_pcg_amg_setup = Get_time();
    multigrid *pcg_amg = Build_amg(A, NULL, param_pcg.max_level);
    Setup_phase(pcg_amg, param_pcg);
    double te_pcg_amg_setup = Get_time();
    Print_amg(pcg_amg);
    
    double tb_pcg_amg = Get_time();
    Linear_solver_pcg_amg(pcg_amg, 0, b, x_pcg_amg, param_pcg, NULL, &resi_norm_pcg_amg, &niter_pcg, &ncycle);
    double te_pcg_amg = Get_time();
    printf("*************** pcg amg method ***************\n");
    printf("resi       = %g\n", resi_norm_pcg_amg);
    printf("pcg iter   = %d\n", niter_pcg);
    printf("amg cycle  = %d\n", ncycle);
    printf("time       = %g\n", te_pcg_amg-tb_pcg_amg);
    printf("setup time = %g\n", te_pcg_amg_setup-tb_pcg_amg_setup);
    printf("**********************************************\n");
    Free_multigrid(pcg_amg);
    free(x_pcg_amg);
    /*-----------------------------------------*/
#endif

#if 1
    if(myrank == 0)
    {
	dmatcsr *A_seq = Read_dmatcsr(file);

	int size = A_seq->nr;
	double *sol_seq   = (double*)calloc(size, sizeof(double));
	double *b_seq     = (double*)calloc(size, sizeof(double));
	double *x_ini_seq = (double*)calloc(size, sizeof(double));
	srand((unsigned)time(0));
	for(i=0; i<size; i++) sol_seq[i] = 1.0;//((double)rand())/((double)RAND_MAX);
	Multi_dmatcsr_dvec(A_seq, sol_seq, b_seq);


#if cg
	/*------------------ cg ------------------*/
	printf("\ncalling cg method...\n");
	double *x_cg_seq = (double*)calloc(size, sizeof(double));
	memcpy(x_cg_seq, x_ini_seq, size*sizeof(double));
	double tol_cg_seq       = param.cg_tol;
	int    max_iter_cg_seq  = param.cg_max_iter;
	double resi_norm_cg_seq = 9999.9;
	int    niter_cg_seq     = -1;
	double tb_cg = Get_time();
	Linear_solver_cg(A_seq, b_seq, x_cg_seq, tol_cg_seq, max_iter_cg_seq, NULL, &resi_norm_cg_seq, &niter_cg_seq, 1);
	double te_cg = Get_time();
	Get_residual(A_seq, b_seq, x_cg_seq, NULL, &resi_norm_cg_seq);
	printf("CG method: resi = %18.15f, niter = %3d, time = %f\n", resi_norm_cg_seq, niter_cg_seq, te_cg-tb_cg);
	free(x_cg_seq);
	/*----------------------------------------*/
#endif
#if amgcycle
	/*------------------ MYAMG ------------------*/
	printf("\ncalling AMG method...\n");
	double *x_amg_seq = (double*)calloc(size, sizeof(double));
	memcpy(x_amg_seq, x_ini_seq, size*sizeof(double));
	double resi_norm_amg_seq = 9999.9;
	int    ncycle_amg_seq    = -1;
	
	amg_param param_myamg_seq = param;
	param_myamg_seq.amgsolver_tol = 1e-08;

	multigrid *amg_seq = Build_amg(A_seq, NULL, param_myamg_seq.max_level);
	double tb_setup = Get_time();
	Setup_phase(amg_seq, param_myamg_seq);
	double te_setup = Get_time();
	Print_amg(amg_seq);
	printf("setup phase time: %f\n", te_setup-tb_setup);

	//Write_dmatcsr_csr(amg_seq->A[3], "../output/A3.dat");
	//Write_dmatcsr_csr(amg_seq->A[2], "../output/A2.dat");

	double tb_amg = Get_time();
	Linear_solver_amg(amg_seq, 0, b_seq, x_amg_seq, param_myamg_seq, &resi_norm_amg_seq, &ncycle_amg_seq);
	double te_amg = Get_time();
	printf("amg method: \nresi = %18.15f, ncycle = %3d, \ntime = %f, setup time = %f\n", resi_norm_amg_seq, ncycle_amg_seq, te_amg-tb_amg, te_setup-tb_setup);

	Free_multigrid(amg_seq);
	free(x_amg_seq);
	/*-----------------------------------------*/
#endif

	free(x_ini_seq);
	free(b_seq);
	free(sol_seq);

	Free_dmatcsr(A_seq);
    }
#endif

    Free_par_dvec(b);
    Free_par_dvec(x_ini);
    Free_par_dvec(sol);
    Free_par_dmatcsr(A);

    MPI_Group_free(&mpi_group_world);
    MPI_Finalize();


    return 0;
}

