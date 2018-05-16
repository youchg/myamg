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

#define eigenpair_given   1
#define direct_method_all 0
#define direct_method_amg 0
#define amg_method        1

#define precondition      0

#define direct_nev        30

#define nmax_correction   20
//#define nev               13
//#define error_nev_b       0
//#define error_nev_e       nev-1

#define tol_correction    1e-09

int print_rank = 0;
int main(int argc, char* argv[])
{
    int nev = 0;
    int error_nev_b = 0;
    int error_nev_e = 0;
    Init_nev_argv(argc, argv, &nev, &error_nev_b, &error_nev_e);

#if !(eigenpair_given || direct_method_all || direct_method_amg)
#undef  direct_method_amg
#define direct_method_amg 1
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("!! Using direct_method_amg by force. !!\n");
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
#endif
    if(argc < 3)
    {
	printf("Too few arguments, exit...\n");
	exit(-1);
    }
    else if(argc == 3)
    {
	printf("Too few arguments, set nev = %d, nb = %d, ne = %d\n", direct_nev, 0, direct_nev-1);
	nev = direct_nev;
	error_nev_b  = 0;
	error_nev_e  = direct_nev - 1;
    }
    else if(argc == 5)
    {
	printf("Too few arguments, nev = %d, set nb = %d, ne = %d\n", nev, 0, nev-1);
	error_nev_b  = 0;
	error_nev_e  = nev - 1;
    }
    else if(argc != 9)
    {
	printf("Unresolved command line parameters, exit...\n");
	exit(-1);
    }
    printf("nev = %d, nb = %d, ne = %d\n", nev, error_nev_b, error_nev_e);
    assert(nev > 0);
    assert(error_nev_b >= 0);
    assert(error_nev_e >= error_nev_b);

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
	            5.567636362502870, 
	            6.450275478656486,
	            7.133449569066343,
	            7.805653861754617,
		    8.437714879048290,
		    9.073708666780679,
		    9.694558679214822,
		   10.323221000756817,
		   10.946226316489714,
		   11.578538113884502,
		   12.209959940792181,
		   12.851373263735441,
		   13.494677839558323,
		   14.148317005625978,
		   14.805630470017009,
		   15.473457016348881,
		   16.146178226559016,
		   16.829504824188092,
		   17.518602607650379,
		   18.218350388322200,
		   18.924522729212772,
		   19.641362350947009,
		   20.365128157582124,
		   21.099562606661241,
		   21.297344444147694,
		   21.841317726052672,
		   22.593733594323009,
		   23.353786603403247,
		   23.409302262446829,
		   24.124487245224749,
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

    double  *total_error = (double*)calloc(nmax_correction, sizeof(double));
    double  *corre_time  = (double*)calloc(nmax_correction, sizeof(double));

#if !eigenpair_given && direct_method_amg
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
#endif

    /* solve eigenvalue problem */
    double tb_correction_amg = Get_time();
    double  *eval_amg = (double*) calloc(nev, sizeof(double));
    double **evec_amg        = (double**)malloc(nev* sizeof(double*));
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
#if eigenpair_given 
    //for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_given[i]) / eval_given[i];
    for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_given[i]);
    for(i=0; i<nev; i++) printf("%2d: %20.15f %20.15f\n", i, eval_amg[i], fabs(eval_amg[i]-eval_given[i]));
#elif direct_method_all
    //for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_direct_all[i]) / eval_direct_all[i];
    for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_direct_all[i]);
    for(i=0; i<nev; i++) printf("%2d: %20.15f %20.15f\n", i, eval_amg[i], fabs(eval_amg[i]-eval_direct_all[i]));
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
	Eigen_solver_amg(amg, nev, eval_amg, evec_amg, 0, 1, param_eigen);
	int j;
	te_amg = Get_time();
        corre_time[i] = te_amg - tb_amg;
#if eigenpair_given
	//for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_given[j])      / eval_given[j];
	for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_given[j]);
	for(j=0; j<nev; j++) printf("%2d: %20.15f %20.15f\n", j, eval_amg[j], fabs(eval_amg[j]-eval_given[j]));
#elif direct_method_all
	//for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_direct_all[j]) / eval_direct_all[j];
	for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_direct_all[j]);
	for(j=0; j<nev; j++) printf("%2d: %20.15f %20.15f\n", j, eval_amg[j], fabs(eval_amg[j]-eval_direct_all[j]));
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

    free(corre_time);
    free(total_error);
    for(i=0; i<nev; i++) free(evec_amg[i]);
    free(evec_amg);
    free(eval_amg);

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

