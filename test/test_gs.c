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

int main(int argc, char* argv[])
{
    dmatcsr *A;
    amg_param param;
    int i;

    double tb_ini = Get_time();
    
    Init_argv(argc, argv, &param, &A, NULL, NULL);
    Remove_zero_dmatcsr(A);
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

    //dmatcsr *D = Create_dmatcsr(A->nr, A->nc, A->nr);
    //dmatcsr *L = Create_dmatcsr(A->nr, A->nc, A->nn);
    //dmatcsr *U = Create_dmatcsr(A->nr, A->nc, A->nn);
    
    //Split_dmatcsr(A, D, L, U);
    
    //Print_dmatcsr(D);
    //Print_dmatcsr(L);
    //Print_dmatcsr(U);
    //Write_dmatcsr_csr(A, "../output/A.dmatcsr");
    //Write_dmatcsr_csr(D, "../output/D.dmatcsr");
    //Write_dmatcsr_csr(L, "../output/L.dmatcsr");
    //Write_dmatcsr_csr(U, "../output/U.dmatcsr");
    double *d = (double*)calloc(A->nr, sizeof(double));
    Get_diag_dmatcsr(A, d);
    
    double *x_gs = (double*)calloc(size, sizeof(double));
    double rnorm_gs = 9999.9;
    int    niter_gs = -1;
    
    double tb_gs = Get_time();
    //Linear_solver_gs(A, D, L, U, b, x_gs, 1e-14, 20, NULL, &rnorm_gs, &niter_gs, 1);
    Linear_solver_gs(A, d, b, x_gs, 1e-14, 20, NULL, &rnorm_gs, &niter_gs, 1);
    double te_gs = Get_time();
    
    //Get_residual(A, b, x_gs, NULL, &rnorm_gs);
    double *error = (double*)calloc(size, sizeof(double));
    for(i=0; i<size; i++) error[i] = sol[i] - x_gs[i];
    double enorm = Get_dvec_2norm(error, size);
    
    printf("GS method: resi = %18.15f, error = %18.15f, niter = %3d, time = %f\n", rnorm_gs, enorm, niter_gs, te_gs-tb_gs);
    
    free(error);
    free(x_gs);
    //Free_dmatcsr(U);
    //Free_dmatcsr(L);
    //Free_dmatcsr(D);
    free(d);
    free(b);
    free(sol);
    Free_dmatcsr(A);
    return 0;
}
