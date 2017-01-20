#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "preprocess.h"
#include "io.h"
#include "fasp_interface.h"
#include "arpack_interface.h"
#include "tool.h"
#include "par_linear_algebra.h"
#include "par_linear_solver.h"
#include "par_eigen_solver.h"
#include "par_matrix_vector.h"
#include "par_multigrid.h"
#include "par_setup_phase.h"

#define eigenpair_given   0
#define direct_method_all 0
#define direct_method_amg 1
#define amg_method        1

#define precondition      0

#define direct_nev        14

#define nmax_correction   2

#define tol_correction    1e-10

int print_rank = 0;
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int  myrank;
    int nproc_global;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc_global);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    MPI_Group mpi_group_world;
    MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_world);

    int nev = 0;
    int error_nev_b = 0;
    int error_nev_e = 0;
    Init_nev_argv(argc, argv, &nev, &error_nev_b, &error_nev_e);
    if(myrank == print_rank) printf("nev = %d, nb = %d, ne = %d\n", nev, error_nev_b, error_nev_e);

#if !(eigenpair_given || direct_method_all || direct_method_amg)
#undef  direct_method_amg
#define direct_method_amg 1
    if(myrank == print_rank)
    {
	printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	printf("!! Using direct_method_amg by force. !!\n");
	printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    }
#endif
    if(argc < 2)
    {
	if(myrank == print_rank) printf("too few arguments!\n");
	exit(0);
    }

    double tb_init = MPI_Wtime();
    par_dmatcsr *A;
    par_dmatcsr *M;
    amg_param param;
    Init_par_amg_param_argv(argc, argv, &param, &A, &M, MPI_COMM_WORLD);
    //Print_par_dmatcsr(A, 0);
    //Print_par_dmatcsr(M, 0);

    if(myrank == print_rank)
	Print_amg_param(param);
    double te_init = MPI_Wtime();
    if(myrank == print_rank)
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
			 };
#endif

#if amg_method
    par_multigrid *pamg = Build_par_amg(A, M, param.max_level, MPI_COMM_WORLD, mpi_group_world);
    double tb_setup = MPI_Wtime();
    Setup_par_phase(pamg, param);
    double te_setup = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    Print_par_amg(pamg);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank == print_rank)
	printf("setup phase time: %f\n", te_setup-tb_setup);

    double  *total_error = (double*)calloc(nmax_correction, sizeof(double));
    double  *corre_time  = (double*)calloc(nmax_correction, sizeof(double));

    /* solve eigenvalue problem */
    double tb_correction_amg = MPI_Wtime();
    double    *eval_amg = (double*)   calloc(nev, sizeof(double));
    par_dvec **evec_amg = (par_dvec**)malloc(nev* sizeof(par_dvec*));
    for(i=0; i<nev; i++) evec_amg[i] = Init_par_dvec_mv(A);
    double tb_amg, te_amg;
    if(myrank == print_rank)
	printf("=============== 0 ===============\n");
    tb_amg = MPI_Wtime();
    amg_param param_eigen = param;
    param_eigen.amgeigen_nouter_iter = 1;
    param_eigen.amgsolver_max_cycle  = 1;
    param_eigen.pcg_amg_max_iter     = 1;
    Eigen_solver_par_amg_nested(pamg, nev, eval_amg, evec_amg, param_eigen);
    te_amg = MPI_Wtime();
    if(myrank == print_rank)
    {
	printf("* 0 * approximate eigenvalue: \n");/* show the result */
	for(i=0; i<nev; i++) printf("%2d: %20.15f\n", i, eval_amg[i]);
	printf("correction %2d time : %20.15f\n", 0, te_amg - tb_amg);
    }
    corre_time[0] = te_amg - tb_amg;
#if eigenpair_given 
    for(i=error_nev_b; i<=error_nev_e; i++) total_error[0] += fabs(eval_amg[i] - eval_given[i]);
#endif
    if(myrank == print_rank)
    {
	printf("correction %2d error: %20.15f\n", 0, total_error[0]);
	printf("***************************************************\n");
	printf("***************************************************\n");
	printf("begin to correct eigenpair on the finest level...\n");
    }
    
    int ncorrection = 0;
    for(i=1; i<nmax_correction; i++)
    {
	if(myrank == print_rank)
	    printf("=============== %d ===============\n", i);
	int j;
	tb_amg = MPI_Wtime();
        param_eigen.amgsolver_max_cycle  = 10;
        param_eigen.pcg_amg_max_iter     = 10;
        param_eigen.amgeigen_nouter_iter = 4;
	Eigen_solver_par_amg(pamg, nev, eval_amg, evec_amg, 0, 1, param_eigen);
	te_amg = MPI_Wtime();
	if(myrank == print_rank)
	{
	    for(j=0; j<nev; j++) printf("%2d: %20.15f\n", j, eval_amg[j]);
	    printf("correction %2d time : %20.15f\n", i, te_amg - tb_amg);
	}
        corre_time[i] = te_amg - tb_amg;
#if eigenpair_given
	for(j=error_nev_b; j<=error_nev_e; j++) total_error[i] += fabs(eval_amg[j] - eval_given[j]);
#endif
	if(myrank == print_rank)
	    printf("correction %2d error: %20.15f\n", i, total_error[i]);
	ncorrection++;
	if(total_error[i] < tol_correction) break;
    }
    double te_correction_amg = MPI_Wtime();
    if(myrank == print_rank)
	printf("==================================\n");

    if(myrank == print_rank)
    {
	printf("=============== correction information ===============\n");
	printf("correction             error            ratio        time\n");
	printf("    %2d       %20.15f     %s     %f\n", 0, total_error[0], "--------", corre_time[0]);
	for(i=1; i<=ncorrection; i++) 
	    printf("    %2d       %20.15f     %f     %f\n", i, total_error[i], total_error[i]/total_error[i-1], corre_time[i]);
	printf("======================================================\n");

	printf("***************************************************\n");
	printf("******** whole correction time: %f *********\n", te_correction_amg - tb_correction_amg);
	printf("***************************************************\n");
    }

    free(corre_time);
    free(total_error);
    for(i=0; i<nev; i++) Free_par_dvec(evec_amg[i]);
    free(evec_amg); evec_amg = NULL;
    free(eval_amg); eval_amg = NULL;

    Free_par_multigrid(pamg);
#endif //amg method

    Free_par_dmatcsr(M);
    Free_par_dmatcsr(A);

    MPI_Group_free(&mpi_group_world);
    MPI_Finalize();

    return 0;
}
