#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "preprocess.h"
#include "io.h"
#include "linear_algebra.h"
#include "linear_solver.h"
#include "matrix.h"
#include "tool.h"


/* cg convergence factor 
   \sigma = \frac{\sqrt(k)-1}{\sqrt(k)+1}
   k = \frac{\lamda_{max}}{\lambda_{min}}
   */
int main(int argc, char* argv[])
{
    dmatcsr *A;

    amg_param param;
    Init_amg_param(&param);
    param.cg_max_iter      = 20;

    int i;

    double tb_ini = Get_time();
    
    //Init_amg_param_argv(argc, argv, &param, &A, NULL, NULL);
    char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine5.m";
    A = Read_dmatcsr(file);
    Print_dmatcsr(A);
    //Print_amg_param(param);
    
    int     size = A->nr;
    double *sol  = (double*)calloc(size, sizeof(double));
    double *b    = (double*)calloc(size, sizeof(double));
    
    srand((unsigned)time(0));
    for(i=0; i<size; i++) sol[i] = 1.0;//((double)rand())/((double)RAND_MAX);
    Multi_dmatcsr_dvec(A, sol, b);
    
    double te_ini = Get_time();
    printf("Initial and read time = %f\n", te_ini-tb_ini);

    printf("\ncalling cg method...\n");
    double *x_cg = (double*)calloc(size, sizeof(double));
    double tol_cg       = param.cg_tol;
    int    max_iter_cg  = param.cg_max_iter;
    double resi_norm_cg = 9999.9;
    int    niter_cg     = -1;
    double tb_cg = Get_time();
    Linear_solver_cg(A, b, x_cg, tol_cg, max_iter_cg, NULL, &resi_norm_cg, &niter_cg, 1);
    double te_cg = Get_time();
    Get_residual(A, b, x_cg, NULL, &resi_norm_cg);
    
    double *error = (double*)calloc(size, sizeof(double));
    for(i=0; i<size; i++) error[i] = sol[i] - x_cg[i];
    double enorm = Get_dvec_2norm(error, size);
    
    printf("CG method: resi = %18.15f, error = %18.15f, niter = %3d, time = %f\n", resi_norm_cg, enorm, niter_cg, te_cg-tb_cg);
    
    free(error);
    free(x_cg);
    free(b);
    free(sol);
    Free_dmatcsr(A);
    return 0;
}
