#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "preprocess.h"
#include "io.h"
#include "linear_algebra.h"
#include "linear_solver.h"
#include "eigen_solver.h"
#include "matrix.h"
#include "multigrid.h"
#include "setup_phase.h"
#include "fasp_interface.h"
#include "arpack_interface.h"
#include "tool.h"

#define nmax_correction   100
#define tol_correction    1e-09

int main(int argc, char* argv[])
{
    if(argc < 2) {printf("too few arguments!\n"); exit(0);}

    int nev = 1;
    //int error_nev_b = 0;
    //int error_nev_e = 0;
    //Init_nev_argv(argc, argv, &nev, &error_nev_b, &error_nev_e);
    //printf("nev = %d, nb = %d, ne = %d\n", nev, error_nev_b, error_nev_e);


    double tb_init = Get_time();
    dmatcsr *A, *M;
    amg_param param;
    Init_amg_param_argv(argc, argv, &param, &A, &M, NULL);
    Print_dmatcsr(A);
    Print_dmatcsr(M);
    //Print_amg_param(param);
    double te_init = Get_time();
    printf("\ninit time: %f\n", te_init-tb_init);
    

    multigrid *amg = Build_amg(A, M, param.max_level);
    double tb_setup = Get_time();
    Setup_phase(amg, param);
    double te_setup = Get_time();
    Print_amg(amg);
    dmatcsr *AH = amg->A[amg->actual_level-1];
    dmatcsr *MH = amg->M[amg->actual_level-1];
    Write_dmatcsr_csr(AH, "../output/AH.dmatcsr");
    Write_dmatcsr_csr(MH, "../output/MH.dmatcsr");
    printf("setup phase time: %f\n", te_setup-tb_setup);

#if 1
    int i, j;
    double  *corre_time  = (double*)calloc(nmax_correction, sizeof(double));

    /* solve eigenvalue problem */
    double tb_correction_amg = Get_time();
    double  *eval_amg = (double*) calloc(nev, sizeof(double));
    double **evec_amg = (double**)malloc(nev* sizeof(double*));
    for(i=0; i<nev; i++) evec_amg[i] = (double*)calloc(A->nc, sizeof(double));
    double tb_amg, te_amg;
    printf("=============== 0 ===============\n");
    tb_amg = Get_time();
    amg_param param_eigen = param;
    param_eigen.amgeigen_nouter_iter = 1;
    param_eigen.amgsolver_max_cycle  = 1;
    param_eigen.pcg_amg_max_iter     = 1;
    Eigen_solver_amg_nested(amg, nev, eval_amg, evec_amg, param_eigen);
    te_amg = Get_time();
    printf("* 0 * approximate eigenvalue: \n");/* show the result */
    corre_time[0] = te_amg - tb_amg;
    for(j=0; j<nev; j++) printf("%2d: %12.9f\n", j, eval_amg[j]);
    printf("correction %2d time : %20.15f\n", 0, te_amg - tb_amg);
    printf("***************************************************\n");
    printf("begin to correct eigenpair on the finest level...\n");
    
    int ncorrection = 0;
    for(i=1; i<nmax_correction; i++)
    {
	printf("=============== %d ===============\n", i);
	tb_amg = Get_time();
        param_eigen.amgsolver_max_cycle  = 1;
        param_eigen.pcg_amg_max_iter     = 1;
	Eigen_solver_amg(amg, nev, eval_amg, evec_amg, 0, 1, param_eigen);
	te_amg = Get_time();
        corre_time[i] = te_amg - tb_amg;
	for(j=0; j<nev; j++) printf("%2d: %12.9f\n", j, eval_amg[j]);
	printf("correction %2d time : %20.15f\n", i, corre_time[i]);
	ncorrection++;
    }
    double te_correction_amg = Get_time();

    printf("***************************************************\n");
    printf("******** whole correction time: %f *********\n", te_correction_amg - tb_correction_amg);
    printf("***************************************************\n");

    free(corre_time);
    for(i=0; i<nev; i++) free(evec_amg[i]);
    free(evec_amg);
    free(eval_amg);
#endif

    Free_multigrid(amg);


    Free_dmatcsr(M);
    Free_dmatcsr(A);
    return 0;
}

