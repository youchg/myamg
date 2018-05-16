#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "preprocess.h"
#include "io.h"
#include "tool.h"

#include "par_linear_algebra.h"
#include "par_linear_solver.h"
#include "par_matrix_vector.h"

/* cg convergence factor 
   \sigma = \frac{\sqrt(k)-1}{\sqrt(k)+1}
   k = \frac{\lamda_{max}}{\lambda_{min}}
   */

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
    param.cg_max_iter      = 20;

    //char file[256] = "../../dat/fem2d_poisson_square/gmg_A_refine9.m";
    //char file[256] = "../../dat/fem3d/hydrogen-stiff-4913.dat";
    //char file[256] = "../../dat/fem2d_poisson_square/gmg_A_refine5.m";
    char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine5.m";
    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine4.m";
    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine3.m";
    //char file[256] = "../../dat/fdm2d9pt/A_fdm9pt_122500x122500.dat";
    //char file[256] = "../../dat/fdm2d9pt/A_fdm9pt_6400x6400.dat";
    //char file[256] = "../dat/fdm2d9pt/A_fdm9pt_49x49.dat";
    par_dmatcsr *A = Read_par_dmatcsr(file, MPI_COMM_WORLD);

    par_dvec *sol = Init_par_dvec_mv(A);
    par_dvec *b   = Init_par_dvec_mv(A);
    
    int i;

    srand((unsigned)time(0));
    for(i=0; i<sol->length; i++) sol->value[i] = 1.0;//((double)rand())/((double)RAND_MAX);
    Multi_par_dmatcsr_dvec(A, sol, b);
    
    if(myrank == print_rank) printf("\ncalling cg method...\n");
    par_dvec *x_cg = Init_par_dvec_mv(A);

    double tol_cg       = param.cg_tol;
    int    max_iter_cg  = param.cg_max_iter;
    double resi_norm_cg = 9999.9;
    int    niter_cg     = -1;

    double tb_cg = MPI_Wtime();
    Linear_solver_par_cg(A, b, x_cg, tol_cg, max_iter_cg, NULL, &resi_norm_cg, &niter_cg, 1);
    double te_cg = MPI_Wtime();
    Get_par_residual(A, b, x_cg, NULL, &resi_norm_cg);
    
    par_dvec *error = Init_par_dvec_mv(A);
    for(i=0; i<error->length; i++) error->value[i] = sol->value[i] - x_cg->value[i];
    double enorm = Get_par_dvec_2norm(error);
    
    if(myrank == print_rank)
	printf("CG method: resi = %18.15f, error = %18.15f, niter = %3d, time = %f\n", resi_norm_cg, enorm, niter_cg, te_cg-tb_cg);
    
    Free_par_dvec(error);
    Free_par_dvec(x_cg);
    Free_par_dvec(b);
    Free_par_dvec(sol);
    Free_par_dmatcsr(A);

    MPI_Finalize();
    return 0;
}
