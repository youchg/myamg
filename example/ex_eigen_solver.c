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
#include "eigen_solver_shift.h"
#include "eigen_solver_multistep.h"
#include "eigen_solver_multistep2.h"
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

#define direct_nev        14

#define nmax_correction   10

#define tol_correction    1e-10

int print_rank = -1; //useless
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
#if 1
	 19.742181652312645,
	 49.360802349344354,
	 49.367944185094970,
	 79.004391701675800,
	 98.754512911503582,
         98.754533209293612,
	128.394168556935057,
	128.454367342076466,
	167.940430957700670,
	167.944318594906093,
	177.893344633664498,
	197.653679181721998,
	197.654155507479231,
	247.074076746614253,
#endif

#if 0
	    10.857195080016409,
	    17.335104047922837,
	    26.671610677454293,
	    28.505230862366687,
	    41.804525628594597,
	    43.820240665791161,
	    52.515200028290948,
	    54.115874724164939,
	    64.824306320110182,
	    71.379471726937822,
	    79.380690379853291,
	    81.126795449851784,
	    90.502801250074597,
	    94.221835539270273,
#endif

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

    Write_dmatcsr_csr(amg->A[amg->actual_level-1], "../output/AH.dat");
    Write_dmatcsr_csr(amg->M[amg->actual_level-1], "../output/MH.dat");

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
	int j;
	tb_amg = Get_time();
        param_eigen.amgsolver_max_cycle  = 10;
        param_eigen.pcg_amg_max_iter     = 10;
        param_eigen.amgeigen_nouter_iter = 4;
	//Eigen_solver_amg(amg, nev, eval_amg, evec_amg, 0, 1, param_eigen);
	//int is_shift = 1;
	double *shift = (double*)calloc(nev, sizeof(double));
	shift[0] = - 9.54;
	shift[1] = -10.06;
	shift[2] = -10.06;
	//Eigen_solver_shift_amg(amg, nev, eval_amg, evec_amg, 0, 1, is_shift, shift, param_eigen);
	//Eigen_solver_amg_multistep(amg, nev, eval_amg, evec_amg, 0, 1, 2, param_eigen);
	Eigen_solver_amg_multistep2(amg, nev, eval_amg, evec_amg, param_eigen);
	te_amg = Get_time();
	free(shift);
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
