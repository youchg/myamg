#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "io.h"
#include "matrix.h"
#include "multigrid.h"
#include "setup_phase.h"
#include "tool.h"
#include "par_matrix_vector.h"
#include "par_multigrid.h"
#include "par_setup_phase.h"
#include "par_linear_algebra.h"

#include "mpi.h"

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
    param.max_coarsest_dof = 1;

    //char file[256] = "../../dat/fem2d_poisson_square/gmg_A_refine9.m";
    //char file[256] = "../../dat/fem3d/hydrogen-stiff-4913.dat";
    //char file[256] = "../../dat/fem2d_poisson_square/gmg_A_refine5.m";
    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine5.m";
    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine4.m";
    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine3.m";
    //char file[256] = "../../dat/fdm2d9pt/A_fdm9pt_122500x122500.dat";
    //char file[256] = "../../dat/fdm2d9pt/A_fdm9pt_6400x6400.dat";
    char file[256] = "../dat/fdm2d9pt/A_fdm9pt_49x49.dat";
    par_dmatcsr *A = Read_par_dmatcsr(file, MPI_COMM_WORLD);

    MPI_Group mpi_group_world;
    MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_world);

    par_multigrid *pamg = Build_par_amg(A, NULL, 20, MPI_COMM_WORLD, mpi_group_world);
    Setup_par_phase(pamg, param);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank == print_rank) printf("\n");
    //if(myrank == print_rank) printf("test function Get_par_dvec_2norm()\n");

    int i;

    par_dvec *x0   = Init_par_dvec_mv(pamg->A[0]);
    par_dvec *x3   = Init_par_dvec_mv(pamg->A[3]);
    par_dvec *A0x0 = Init_par_dvec_mv(pamg->A[0]);
    par_dvec *P0x3 = Init_par_dvec_mv(pamg->A[0]);
    for(i=0; i<x0->length; i++) x0->value[i] = (double)(pamg->A[0]->row_start[myrank] + i + 1)/pamg->A[0]->nc_global/pamg->A[0]->nc_global;

    Multi_par_dmatcsr_dvec(pamg->A[0], x0, A0x0);


    Restrict_par_f2c(pamg, 0, 3, A0x0, x3);
    Prolong_par_c2f(pamg, 3, 0, x3, P0x3);

    double norm_x0 = Get_par_dvec_2norm(x0);
    if(myrank == print_rank) printf("norm(x0) = %15.12f, length = %d, length_global = %d\n", norm_x0, x0->length, x0->length_global);

    double norm_A0x0 = Get_par_dvec_2norm(A0x0);
    if(myrank == print_rank) printf("norm(A0x0) = %15.12f\n", norm_A0x0);

    if(MPI_COMM_NULL != pamg->comm[3])
    {
	double norm_x3 = Get_par_dvec_2norm(x3);
	//printf("norm(x3) = %15.12f, length = %d, length_global = %d\n", norm_x3, x3->length, x3->length_global);

	int myrank3;
	MPI_Comm_rank(pamg->comm[3], &myrank3);
	if(myrank3 == print_rank) printf("norm(x3) = %15.12f, length = %d, length_global = %d\n", norm_x3, x3->length, x3->length_global);
    }

    double norm_P0x3 = Get_par_dvec_2norm(P0x3);
    if(myrank == print_rank) printf("norm(P0x3) = %18.15f\n", norm_P0x3);

    Write_par_dmatcsr_csr(pamg->A[0], "../output/A0.par.dat", 0);
    Write_par_dmatcsr_csr(pamg->P[0], "../output/P0.par.dat", 0);
    Write_par_dmatcsr_csr(pamg->R[0], "../output/R0.par.dat", 0);
    if(MPI_COMM_NULL != pamg->comm[1])
    {
	Write_par_dmatcsr_csr(pamg->R[1], "../output/R1.par.dat", 0);
	Write_par_dmatcsr_csr(pamg->P[1], "../output/P1.par.dat", 0);
    }
    if(MPI_COMM_NULL != pamg->comm[2])
    {
	Write_par_dmatcsr_csr(pamg->R[2], "../output/R2.par.dat", 0);
	Write_par_dmatcsr_csr(pamg->P[2], "../output/P2.par.dat", 0);
    }

    Free_par_dvec(P0x3);
    Free_par_dvec(x0);
    Free_par_dvec(x3);
    Free_par_dvec(A0x0);

#if 0

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank == print_rank) printf("\n");
    //if(myrank == print_rank) printf("test function Get_par_dvec_2norm()\n");

    par_dvec *x0 = Init_par_dvec_mv(pamg->A[0]);

    int i;
    for(i=0; i<x0->length; i++) x0->value[i] = (double)(A->row_start[myrank] + i + 1)/A->nc_global/A->nc_global;
    double norm_x0 = Get_par_dvec_2norm(x0);
    if(myrank == print_rank) printf("norm(x0) = %15.12f, length = %d, length_global = %d\n", norm_x0, x0->length, x0->length_global);

    par_dvec *A0x0 = Init_par_dvec_mv(pamg->A[0]);
    Multi_par_dmatcsr_dvec(pamg->A[0], x0, A0x0);
    double norm_A0x0 = Get_par_dvec_2norm(A0x0);
    if(myrank == print_rank) printf("norm(A0x0) = %15.12f\n", norm_A0x0);

    //MPI_Barrier(MPI_COMM_WORLD);
    //Print_par_dmatcsr(pamg->R[0], 3);
    par_dvec *R0A0x0 = Init_par_dvec_mv(pamg->A[0]);
    Multi_par_dmatcsr_dvec(pamg->R[0], A0x0, R0A0x0);
    double norm_R0A0x0 = Get_par_dvec_2norm(R0A0x0);
    if(myrank == print_rank) printf("norm(R0A0x0) = %15.12f\n", norm_R0A0x0);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //Print_par_dmatcsr(pamg->P[0], 3);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    par_dvec *P0R0A0x0 = Init_par_dvec_mv(pamg->A[0]);
    Multi_par_dmatcsr_dvec(pamg->P[0], R0A0x0, P0R0A0x0);
    double norm_P0R0A0x0 = Get_par_dvec_2norm(P0R0A0x0);
    if(myrank == print_rank) printf("norm(P0R0A0x0) = %15.12f\n", norm_P0R0A0x0);

    Write_par_dmatcsr_csr(pamg->A[0], "../output/A0.par.dat", 0);
    Write_par_dmatcsr_csr(pamg->P[0], "../output/P0.par.dat", 0);
    Write_par_dmatcsr_csr(pamg->R[0], "../output/R0.par.dat", 0);

    Free_par_dvec(P0R0A0x0);
    Free_par_dvec(R0A0x0);
    Free_par_dvec(A0x0);
    Free_par_dvec(x0);
#endif

    //if(MPI_COMM_NULL != pamg->comm[1]) Print_par_dmatcsr(pamg->A[1], 3);

    MPI_Group_free(&mpi_group_world);
    Free_par_multigrid(pamg);
    Free_par_dmatcsr(A);

    MPI_Finalize();
    return 0;
}
