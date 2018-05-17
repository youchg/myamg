#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>
#include "preprocess.h"
#include "io.h"
#include "linear_algebra.h"
#include "linear_solver.h"
#include "eigen_solver.h"
#include "eigen_solver_augmented.h"
#include "eigen_solver_augmented_one.h"
#include "matrix.h"
#include "multigrid.h"
#include "setup_phase.h"
#include "fasp_interface.h"
#include "arpack_interface.h"
#include "tool.h"

#define eigenpair_given   1
#define direct_method_all 0
#define amg_method        1
#define direct_method_amg 0
#define amg_eigen         1

#define precondition      0

#define direct_nev        32

#define tol_correction    1e-09

#define nmax_correction   20

int print_rank = 0;
int main(int argc, char* argv[])
{
    int nev = 0;
    int error_nev_b = 0;
    int error_nev_e = 0;
    int augmented_index = 12;
    Init_nev_argv(argc, argv, &nev, &error_nev_b, &error_nev_e);
    printf("nev = %d, nb = %d, ne = %d\n", nev, error_nev_b, error_nev_e);
    error_nev_b = augmented_index;
    error_nev_e = augmented_index;
    assert(nev >= augmented_index+1);

#if !(eigenpair_given || direct_method_all || direct_method_amg)
#undef  direct_method_amg
#define direct_method_amg 1
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("!! Using direct_method_amg by force. !!\n");
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
#endif
    if(argc < 2)
    {
	printf("too few arguments!\n");
	exit(0);
    }
    double tb_init = Get_time();
    dmatcsr *A;
    dmatcsr *M;
    amg_param param;
    Init_amg_param_argv(argc, argv, &param, &A, &M, NULL);
    Print_dmatcsr(A);
    Print_dmatcsr(M);
    Remove_zero_dmatcsr(A);
    Remove_zero_dmatcsr(M);
    Print_dmatcsr(A);
    Print_dmatcsr(M);

#if precondition
    double *diag = (double*)calloc(A->nr, sizeof(double));
    Get_diag_dmatcsr(A, diag);
    //Write_dvec(diag, A->nr, "../output/diag.dat");
    Zoom_dmatcsr(A, diag);
    Zoom_dmatcsr(M, diag);
    //Write_dmatcsr_csr(A, "../output/zoomA.dmatcsr");
    //Write_dmatcsr_csr(M, "../output/zoomM.dmatcsr");
#endif

    Print_amg_param(param);
    double te_init = Get_time();
    printf("\ninit time: %f\n", te_init-tb_init);
    
    int i;
#if eigenpair_given
    double eval_given[direct_nev] = {
#if 0
	 19.739351152446577,
	 49.348398772994052,
	 49.348426661655587,
	 78.957543954641395,
	 98.696926072478092,
	 98.696926072787591,
	128.306055963593963,
	128.306290898228610,
	167.784999753374791,
	167.785014924843836,
	177.654996436879685,
	197.394417265854258,
	197.394417273111628,
	246.743050204326266,
	246.743972243036808,
	256.612819086151660,
	256.612819086431614,
	286.222396771513957,
	286.222479299660563,
	315.832405282174534,
	335.571880346205262,
	335.571880413890710,
	365.180563509384399,
	365.180577155566198,
	394.790444748240247,
	394.790444752357189,
	404.659991649028711,
	404.662532642826989,
	444.140187371669356,
	444.140438372226185,
	493.488516543021433,
	493.488516544614299,
#else
         19.742181652312652, 
         49.360802349344290,
         49.367944185095041,
         79.004391701675502,
         98.754512911503113,
         98.754533209294024,
        128.394168556935540,
        128.454367342076296,
        167.940430957701437,
        167.944318594906576,
        177.893344633664412,
        197.653679181721060,
        197.654155507478691,
        247.074076746612491,
        247.310545615300498,
        256.969593559649127,
        256.969611929320706,
        286.724147309196951,
        286.745315764746238,
        316.584651924283833,
        336.361880577667876,
        336.366327347913796,
        365.888058939655650,
        365.891559701172412,
        395.720359990433792,
        395.720630471374477,
        405.555049749392197,
        406.207452347479659,
        445.426637813072546,
        445.491139460087368,
        494.768409981298817,
        494.768519621335997,

#endif

    };
#endif

#if !eigenpair_given && direct_method_all
    double   tb_direct_all          = Get_time();
    double  *eval_direct_all = (double*) calloc(direct_nev, sizeof(double));
    double **evec_direct_all = (double**)malloc(direct_nev* sizeof(double*));
    for(i=0; i<direct_nev; i++) evec_direct_all[i] = (double*)calloc(A->nc, sizeof(double));
    printf("\ncalling direct method all...\n");
    Eigen_solver_arpack_dn(A, M, direct_nev, eval_direct_all, evec_direct_all);
    printf("================= direct all result ===================\n");
    for(i=0; i<direct_nev; i++) printf("%2d: %20.15f\n", i, eval_direct_all[i]);
    printf("===================================================\n");
    double te_direct_all = Get_time();
    printf("direct eigen time: %f\n", te_direct_all - tb_direct_all);
#endif

#if amg_method
    multigrid *amg = Build_amg(A, M, param.max_level);
    double tb_setup = Get_time();
    Setup_phase(amg, param);
    double te_setup = Get_time();
    Print_amg(amg);
    printf("setup phase time: %f\n", te_setup-tb_setup);

#if !eigenpair_given && direct_method_amg
    //printf("Direct eigen solver begin: memory use (MB): %f\n", Get_memory());

    double   tb_direct_amg          = Get_time();
    double  *eval_direct_amg = (double*) calloc(direct_nev, sizeof(double));
    double **evec_direct_amg = (double**)malloc(direct_nev* sizeof(double*));

    for(i=0; i<direct_nev; i++) evec_direct_amg[i] = (double*)calloc(A->nc, sizeof(double));
    printf("calling direct method amg...\n");
    Eigen_solver_arpack_dn_amg(amg, 0, direct_nev, eval_direct_amg, evec_direct_amg, param);
    printf("================= direct amg result ===================\n");
    for(i=0; i<direct_nev; i++) printf("%2d: %20.15f\n", i, eval_direct_amg[i]);
    printf("===================================================\n");
    double te_direct_amg = Get_time();
    printf("direct eigen amg time: %f\n", te_direct_amg - tb_direct_amg);

    //printf("Direct eigen solver end: memory use (MB): %f\n", Get_memory());
#endif

#if amg_eigen
    //printf("AMG eigen solver begin: memory use (MB): %f\n", Get_memory());

    double  *total_error = (double*)calloc(nmax_correction, sizeof(double));
    double  *corre_time  = (double*)calloc(nmax_correction, sizeof(double));

    /* solve eigenvalue problem */
    double tb_correction_amg = Get_time();
    double  *eval_amg = (double*) calloc(nev, sizeof(double));
    double **evec_amg        = (double**)malloc(nev* sizeof(double*));
    for(i=0; i<nev; i++) evec_amg[i] = (double*)calloc(A->nc, sizeof(double));
    srand(1);
    int j = 0;
    // random initial vector
    for(i = 0; i < nev; ++i) {
      for(j = 0; j < A->nc; ++j) {
        evec_amg[i][j] = (double)(rand()) / (double)((unsigned long)RAND_MAX + 1);
      }
    }
    double tb_amg, te_amg;
    printf("=============== 0 ===============\n");
    tb_amg = Get_time();
    amg_param param_eigen = param;
    param_eigen.amgeigen_nouter_iter = 1;
    param_eigen.amgsolver_max_cycle  = 1;
    param_eigen.pcg_amg_max_iter     = 1;
    Eigen_solver_amg_augmented(amg, nev, eval_amg, evec_amg, param_eigen);
    te_amg = Get_time();
    printf("* 0 * approximate eigenvalue: \n");/* show the result */
    corre_time[0] = te_amg - tb_amg;
#if eigenpair_given 
    //for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_given[i]) / eval_given[i];
    for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_given[i]);
    //for(i=0; i<nev; i++) printf("%2d: %20.15f %20.15f\n", i, eval_amg[i], fabs(eval_amg[i]-eval_given[i]));
    for(i=error_nev_b; i<=error_nev_e; i++) printf("%2d: %20.15f %20.15f\n", i, eval_amg[i], eval_amg[i]-eval_given[i]);
#elif direct_method_all
    //for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_direct_all[i]) / eval_direct_all[i];
    for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_direct_all[i]);
    //for(i=0; i<nev; i++) printf("%2d: %20.15f %20.15f\n", i, eval_amg[i], fabs(eval_amg[i]-eval_direct_all[i]));
    for(i=error_nev_b; i<=error_nev_e; i++) printf("%2d: %20.15f %20.15f\n", i, eval_amg[i], eval_amg[i]-eval_direct_all[i]);
#elif direct_method_amg
    //for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_direct_amg[i]) / eval_direct_amg[i];
    for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_direct_amg[i]);
    for(i=0; i<nev; i++) printf("%2d: %20.15f %20.15f\n", i, eval_amg[i], fabs(eval_amg[i]-eval_direct_amg[i]));
#endif
    printf("correction %2d time : %20.15f\n", 0, te_amg - tb_amg);
    printf("correction %2d error: %20.15f\n", 0, total_error[0]);
    printf("***************************************************\n");
    printf("***************************************************\n");
    printf("begin to correct eigenpair on the finest level...\n");
    
    int ncorrection = 0;
    for(i=1; i<nmax_correction; i++)
    {
	printf("=============== %d ===============\n", i);
	tb_amg = Get_time();
        param_eigen.amgsolver_max_cycle  = 1;
        param_eigen.pcg_amg_max_iter     = 1;
        assert(nev >= augmented_index);
	Eigen_solver_amg_augmented_one(amg, augmented_index, &eval_amg[augmented_index], evec_amg[augmented_index], param_eigen);
	int j;
	te_amg = Get_time();
        corre_time[i] = te_amg - tb_amg;
#if eigenpair_given
	//for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_given[j])      / eval_given[j];
	for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_given[j]);
	//for(j=0; j<nev; j++) printf("%2d: %20.15f %20.15f\n", j, eval_amg[j], fabs(eval_amg[j]-eval_given[j]));
        for(j=error_nev_b; j<=error_nev_e; j++) printf("%2d: %20.15f %20.15f\n", j, eval_amg[j], eval_amg[j]-eval_given[j]);
#elif direct_method_all
	//for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_direct_all[j]) / eval_direct_all[j];
	for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_direct_all[j]);
	//for(j=0; j<nev; j++) printf("%2d: %20.15f %20.15f\n", j, eval_amg[j], fabs(eval_amg[j]-eval_direct_all[j]));
        for(j=error_nev_b; j<=error_nev_e; j++) printf("%2d: %20.15f %20.15f\n", j, eval_amg[j], eval_amg[j]-eval_direct_all[j]);
#elif direct_method_amg
	//for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_direct_amg[j]) / eval_direct_amg[j];
	for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_direct_amg[j]);
	for(j=0; j<nev; j++) printf("%2d: %20.15f %20.15f\n", j, eval_amg[j], fabs(eval_amg[j]-eval_direct_amg[j]));
#endif
	printf("correction %2d time : %20.15f\n", i, corre_time[i]);
	printf("correction %2d error: %20.15f\n", i, total_error[i]);
	ncorrection++;
	if(total_error[i] < tol_correction) break;
    }
    double te_correction_amg = Get_time();
    printf("==================================\n");

    printf("=============== correction information ===============\n");
    printf("correction             error            ratio        time\n");
    printf("    %2d       %20.15f     %s     %f\n", 0, total_error[0], "--------", corre_time[0]);
    for(i=1; i<=ncorrection; i++) 
	printf("    %2d       %20.15f     %f     %f\n", i, total_error[i], total_error[i]/total_error[i-1], corre_time[i]);
    printf("======================================================\n");

    printf("***************************************************\n");
    printf("******** whole correction time: %f *********\n", te_correction_amg - tb_correction_amg);
    printf("***************************************************\n");

    //printf("AMG eigen solver end: memory use (MB): %f\n", Get_memory());

    free(corre_time);
    free(total_error);
    for(i=0; i<nev; i++) free(evec_amg[i]);
    free(evec_amg);
    free(eval_amg);
#endif

#if !eigenpair_given && direct_method_amg
    for(i=0; i<direct_nev; i++) free(evec_direct_amg[i]);
    free(evec_direct_amg);
    free(eval_direct_amg);
#endif

    Free_multigrid(amg);
#endif //amg method

#if !eigenpair_given && direct_method_all
    for(i=0; i<direct_nev; i++) free(evec_direct_all[i]);
    free(evec_direct_all);
    free(eval_direct_all);
#endif

#if precondition
    free(diag);
#endif
    Free_dmatcsr(M);
    Free_dmatcsr(A);
    return 0;
}

