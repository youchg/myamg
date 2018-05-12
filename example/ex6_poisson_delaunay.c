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

#define eigenpair_given   0
#define direct_method_all 0
#define direct_method_amg 1
#define amg_method        1

#define precondition      0

#define direct_nev        13//32

#define nmax_correction   20

#define tol_correction    1e-9

int print_rank = 0;
int main(int argc, char* argv[])
{
    int nev = 0;
    int error_nev_b = 0;
    int error_nev_e = 0;
    Init_nev_argv(argc, argv, &nev, &error_nev_b, &error_nev_e);
    printf("nev = %d, nb = %d, ne = %d\n", nev, error_nev_b, error_nev_e);

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
    //Remove_zero_dmatcsr(A);
    //Remove_zero_dmatcsr(M);
    //Print_dmatcsr(A);
    //Print_dmatcsr(M);

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
	 19.739203676169186,
 	 49.348096273196418,
	 49.348104775244764,
	 78.957103216475232,
	 98.696500327977517,
	 98.696514274521434,
	128.305687568375049,
	128.305698861836106,
	167.784777146442906,
	167.784794337525170,
	177.654594254445982,
	197.394228908278592,
	197.394231580738392,
	246.743540132765531,
	246.743580128812198,
	256.613445582319628,
	256.613479238986315,
	286.223219893151509,
	286.223246191436033,
	315.833123356659939,
	335.573111897455533,
	335.573122649725292,
	365.183159433302706,
	365.183186179283837,
	394.793336178531888,
	394.793345641497467,
	404.663432747543879,
	404.663468426061002,
	444.143882437289960,
	444.143906862065705,
	493.494685952597990,
	493.494755401723012,
#endif


#if 1
19.739321565904365,
49.348322495822067,
49.348330785288965,
78.957363349172596,
98.696719175791429,
98.696750921254704,
128.305792390932851,
128.305819662987005,
167.784613898406633,
167.784636869244764,
177.654333411710297,
197.393753105679366,
197.393809934040291,
246.742400450661620,
246.742438206683801,
256.612109705594946,
256.612152226199271,
286.221353338815788,
286.221365497927081,
315.830656416798377,
335.570090183798243,
335.570172602859145,
365.179390794456253,
365.179456915746186,
394.788733284832688,
394.788788676746776,
404.658519996561949,
404.658597342113922,
444.137698569071631,
444.137729397037845,
493.486633007975229,
493.486698661234016,
#endif

#if 0
	 19.739819119754731,
	 49.351830760406258,
	 49.351847952831143,
	 78.966604821783420,
	 98.711275970328018,
	 98.711352131047505,
        128.330640103796640,
	128.330736738368500,
        167.827310862382632,
        167.827457003635061,
	177.702360492539327,
	197.453088659956876,
	197.453293336271969,
	246.835491574125143,
	246.835710334150320,
	256.712668300863527,
	256.713108689211310,
	286.346700458569842,
	286.347086929429793,
	315.983929661519198,
	335.742843293247290,
	335.743324422616695,
	365.383749217224135,
	365.384764423248384,
	395.027839004284715,
	395.028514106624755,
	404.910429824543996,
	404.911407550863544,
	444.441011249880546,
	444.442143172008514,
	493.860971959731387,
	493.862475347726843,
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
    for(i=0; i<nev; i++) printf("%2d: %20.15f\n", i, eval_amg[i]);
    printf("correction %2d time : %20.15f\n", 0, te_amg - tb_amg);
    corre_time[0] = te_amg - tb_amg;
#if eigenpair_given 
    for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_given[i]);
#elif direct_method_all
    for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_direct_all[i]);
#elif direct_method_amg
    for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_direct_amg[i]);
#endif
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
        param_eigen.pcg_amg_max_iter     = 2;
	Eigen_solver_amg(amg, nev, eval_amg, evec_amg, 0, 1, param_eigen);
	int j;
	te_amg = Get_time();
	for(j=0; j<nev; j++) printf("%2d: %20.15f\n", j, eval_amg[j]);
	printf("correction %2d time : %20.15f\n", i, te_amg - tb_amg);
        corre_time[i] = te_amg - tb_amg;
#if eigenpair_given
	for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_given[j]);
#elif direct_method_all
	for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_direct_all[j]);
#elif direct_method_amg
	for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_direct_amg[j]);
#endif
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
