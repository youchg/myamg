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
#include "multigrid.h"
#include "setup_phase.h"
#include "fasp_interface.h"
#include "tool.h"

#define cg        0
#define umfpack   0
#define fasp      0
#define amgcycle  1
#define amgpcg    0

/* cg convergence factor 
   \sigma = \frac{\sqrt(k)-1}{\sqrt(k)+1}
   k = \frac{\lamda_{max}}{\lambda_{min}}
   */
int main(int argc, char* argv[])
{
    dmatcsr *A;
    amg_param param;
    int i;

    double tb_ini = Get_time();
    Init_amg_param_argv(argc, argv, &param, &A, NULL, NULL);
    //Remove_zero_dmatcsr(A);
    Print_dmatcsr(A);
    int size = A->nr;
    Print_amg_param(param);
    double *sol   = (double*)calloc(size, sizeof(double));
    double *b     = (double*)calloc(size, sizeof(double));
    double *x_ini = (double*)calloc(size, sizeof(double));
    srand((unsigned)time(0));
    for(i=0; i<size; i++) sol[i] = 1.0;//((double)rand())/((double)RAND_MAX);
    Multi_dmatcsr_dvec(A, sol, b);
    double te_ini = Get_time();
    printf("Initial and read time = %f\n", te_ini-tb_ini);
    
    //Write_dmatcsr_fasp(A, "../output/csrmat_FE.dat");
    //Write_dvec(b, size, "../output/rhs_FE.dat");

#if cg
    /*------------------ cg ------------------*/
    printf("\ncalling cg method...\n");
    double *x_cg = (double*)calloc(size, sizeof(double));
    memcpy(x_cg, x_ini, size*sizeof(double));
    double tol_cg       = param.cg_tol;
    int    max_iter_cg  = param.cg_max_iter;
    double resi_norm_cg = 9999.9;
    int    niter_cg     = -1;
    double tb_cg = Get_time();
    Linear_solver_cg(A, b, x_cg, tol_cg, max_iter_cg, NULL, &resi_norm_cg, &niter_cg, 0);
    double te_cg = Get_time();
    Get_residual(A, b, x_cg, NULL, &resi_norm_cg);
    printf("CG method: resi = %18.15f, niter = %3d, time = %f\n", resi_norm_cg, niter_cg, te_cg-tb_cg);
    free(x_cg);
    /*----------------------------------------*/
#endif

#if umfpack
    /*------------------ UMFPACK ------------------*/
    printf("\ncalling UMFPACK direct method...\n");
    double *x_umfpack = (double*)calloc(size, sizeof(double));
    memcpy(x_umfpack, x_ini, size*sizeof(double));
    double rn_umfpack;
    double tb_umfpack = Get_time();
    Linear_solver_direct(A, b, x_umfpack);
    double te_umfpack = Get_time();
    Get_residual(A, b, x_umfpack, NULL, &rn_umfpack);
    printf("umfpack method: rn = %18.15f, time = %f\n", rn_umfpack, te_umfpack-tb_umfpack);
    free(x_umfpack);
    /*---------------------------------------------*/
#endif

#if fasp
    /*------------------ FASP AMG ------------------*/
    printf("\ncalling FASP AMG method...\n");
    double *x_faspamg = (double*)calloc(size, sizeof(double));
    memcpy(x_faspamg, x_ini, size*sizeof(double));
    double resi_norm_faspamg = 9999.9;
    int    ncycle_fasp       = -1;
    
    amg_param param_fasp = param;
    param_fasp.amgsolver_tol = 1e-08;

    double tb_faspsetup = Get_time();
    multigrid *faspamg = Build_amg(A, NULL, param_fasp.max_level);
    Setup_phase_fasp(faspamg);
    double te_faspsetup = Get_time();
    Print_amg(faspamg);

    double tb_faspamg = Get_time();
    Linear_solver_amg(faspamg, 0, b, x_faspamg, param_fasp, &resi_norm_faspamg, &ncycle_fasp);
    double te_faspamg = Get_time();
    printf("faspamg method: \nresi = %18.15f, ncycle = %3d, \ntime = %f, setup time = %f\n", resi_norm_faspamg, ncycle_fasp, te_faspamg-tb_faspamg, te_faspsetup-tb_faspsetup);
    Free_multigrid(faspamg);
    free(x_faspamg);
    /*-----------------------------------------*/
#endif

#if amgcycle
    /*------------------ MYAMG ------------------*/
    printf("\ncalling AMG method...\n");
    double *x_amg = (double*)calloc(size, sizeof(double));
    memcpy(x_amg, x_ini, size*sizeof(double));
    double resi_norm_amg = 9999.9;
    int    ncycle_amg    = -1;
    
    amg_param param_myamg = param;
    param_myamg.amgsolver_tol = 1e-08;

    multigrid *amg = Build_amg(A, NULL, param_myamg.max_level);
    double tb_setup = Get_time();
    Setup_phase(amg, param_myamg);
    double te_setup = Get_time();
    Print_amg(amg);
    printf("setup phase time: %f\n", te_setup-tb_setup);
    //Write_dmatcsr_csr(amg->A[0], "../output/A0.dmatcsr");
    //Write_dmatcsr_csr(amg->P[0], "../output/P0.dmatcsr");
    //Write_dmatcsr_csr(amg->R[0], "../output/R0.dmatcsr");
    //Write_dmatcsr_csr(amg->A[1], "../output/A1.dmatcsr");

    double tb_amg = Get_time();
    Linear_solver_amg(amg, 0, b, x_amg, param_myamg, &resi_norm_amg, &ncycle_amg);
    double te_amg = Get_time();
    printf("amg method: \nresi = %18.15f, ncycle = %3d, \ntime = %f, setup time = %f\n", resi_norm_amg, ncycle_amg, te_amg-tb_amg, te_setup-tb_setup);
    Free_multigrid(amg);
    free(x_amg);
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

    free(b);
    free(x_ini);
    free(sol);
    Free_dmatcsr(A);
    return 0;
}
