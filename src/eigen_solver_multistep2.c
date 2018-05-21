#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "matrix.h"
#include "linear_algebra.h"
#include "multigrid.h"
#include "eigen_solver_multistep.h"
#include "arpack_interface.h"
#include "linear_solver.h"
#include "amg_param.h"
#include "tool.h"
#include "io.h"

#define amg_eigen_linear_solver 1 //1:amgsolver 2:pcg_amg

/* using 'dvec':
 * 1. to form rhsdvec' 
 * 2. as initial approximation 
 * 3. as final   approximation
 * */
static double Correction_solve_linear(multigrid *amg, int current_level, int n, double *dval, double **dvec, amg_param param);


static double Correction_expand_matrix_RAhV_VTAhV_multistep(dmatcsr *AL, dmatcsr *Ah, dmatcsr *AH, 
	                                                    multigrid *amg, int current_level, int coarsest_level, 
				                            int n, double **V, int length);

static double Correction_get_new_evec_multistep(multigrid *amg, int current_level, int coarsest_level,
	                                        int n, double **V, int length, double **evec_expand, double **evec);


/* 当 init_approx_level >= 0 时，eval, evec 必须给定初值，并且不能为 0;
   当 init_approx_level == 0 时，is_init_approx_level_correction 必须大于 0
*/
//Problem: when already have a finer result such as 0-level approximation, but
//want to correct it from a coarser level, there maybe some problem
//----20160301
void Eigen_solver_amg_multistep2(multigrid *amg, int nev, double *eval, double **evec, 
                                 amg_param param)
{
    double t1, t2, t3, t4;
    double eigen_time    = 0;
    double linear_time   = 0;
    double expand_time   = 0;
    double new_evec_time = 0;
    t1 = Get_time();

    int status = 1;

    int i, j, m;
    int nlevel         = amg->actual_level;
    int niter_outer    = param.amgeigen_nouter_iter;
    int finest_level   = 0;
    int coarsest_level = param.amgeigen_coarsest_level;
    if(coarsest_level <= 0) coarsest_level += nlevel-1;

    dmatcsr *matA = amg->A[finest_level];
    dmatcsr *matM = amg->M[finest_level];
    //dmatcsr *matP = amg->P[finest_level];
    dmatcsr *AH   = amg->A[coarsest_level];
    dmatcsr *MH   = amg->M[coarsest_level]; 
    //=====================================================================

    int multistep = 3;
    double **dvec_expand = (double**)malloc(multistep*nev * sizeof(double*));
    for(i=0; i<multistep*nev; i++) dvec_expand[i] = (double*)NULL;

    double **dvec_amg = (double**)malloc(nev * sizeof(double*));
    for(i=0; i<nev; i++) dvec_amg[i] = (double*)calloc(matA->nr, sizeof(double));


    //================= amg method for eigenvalue problem =================
    if(status) printf("Begin amg cycle...\n");


    niter_outer = 3;
    int nev_new = nev;
    double *eval_new = (double *)calloc(nev_new, sizeof(double));
    for(m=0; m<niter_outer; m++)
    {
	if(status) printf("correction: m = %d\n", m);
       
	for(j=0; j<nev; j++) memcpy(dvec_amg[j], evec[j], matA->nr*sizeof(double));
	linear_time += Correction_solve_linear(amg, finest_level, nev, eval, dvec_amg, param);

	for(j=0; j<nev; j++) dvec_expand[j+2*nev] = dvec_expand[j+nev];
	for(j=0; j<nev; j++) dvec_expand[j+nev]   = evec[j];
	for(j=0; j<nev; j++) dvec_expand[j]       = dvec_amg[j];

	int length = (m == 0) ? multistep-1: multistep;
	//int length = multistep;
	printf("length = %d\n", length);
	dmatcsr *Alarge = Expand_dmatcsr_struct(AH, nev*length);
	dmatcsr *Mlarge = Expand_dmatcsr_struct(MH, nev*length);
	expand_time += Correction_expand_matrix_RAhV_VTAhV_multistep(Alarge, matA, AH, amg, finest_level, coarsest_level, nev, dvec_expand, length); 
	expand_time += Correction_expand_matrix_RAhV_VTAhV_multistep(Mlarge, matM, MH, amg, finest_level, coarsest_level, nev, dvec_expand, length); 
	Write_dmatcsr_csr(Alarge, "../output/Alarge.dat");
	Write_dmatcsr_csr(Mlarge, "../output/Mlarge.dat");
	dmatcsr *Alarge_new = Copy_dmatcsr(Alarge);
	dmatcsr *Mlarge_new = Copy_dmatcsr(Mlarge);
	Write_dmatcsr_csr(Alarge_new, "../output/Alarge_new.dat");
	Write_dmatcsr_csr(Mlarge_new, "../output/Mlarge_new.dat");
	dmatcsr *Alarge_new_new = Copy_dmatcsr(Alarge);
	dmatcsr *Mlarge_new_new = Copy_dmatcsr(Mlarge);

	double **evec_expand     = (double**)malloc(nev     * sizeof(double*));
	double **evec_expand_new = (double**)malloc(nev_new * sizeof(double*));
	for(j=0; j<nev;     j++) evec_expand[j]     = (double*)calloc(Alarge->nr, sizeof(double));
	for(j=0; j<nev_new; j++) evec_expand_new[j] = (double*)calloc(Alarge->nr, sizeof(double));


	if(status) printf("Solving corrected eigenvalue problem...\n");
	t3 = Get_time(); 
	memset(eval,     0, nev    *sizeof(double));
	memset(eval_new, 0, nev_new*sizeof(double));
	//for(j=0; j<nev; j++) memset(evec_expand[j], 0, Alarge->nr*sizeof(double));
	Eigen_solver_arpack_dn(Alarge,     Mlarge,     nev,     eval,     evec_expand);
	Eigen_solver_arpack_dn(Alarge_new, Mlarge_new, nev_new, eval_new, evec_expand_new);
	//Insertion_ascend_sort_dvec_dvecvec(eval,     evec_expand,     0, nev-1);
	//Insertion_ascend_sort_dvec_dvecvec(eval_new, evec_expand_new, 0, nev_new-1);
	int output_num = Alarge->nr;
	printf("========= corrected eigenvalue ===========\n");
	for(j=0; j<output_num; j++)
	{
	    for(i=0; i<nev; i++)
		printf("%8.5f %8.5f", evec_expand[i][j], evec_expand_new[i][j]);
	    printf("\n");
	}
	/*
	printf("========= NEW corrected eigenvalue ===========\n");
	for(j=0; j<output_num; j++)
	{
	    for(i=0; i<nev_new; i++)
		printf("%8.5f ", evec_expand_new[i][j]);
	    printf("\n");
	}
	*/
	t4 = Get_time();
	eigen_time += t4-t3;
	
	if(status)
	{
	    printf("========= corrected eigenvalue ===========\n");
	    for(j=0; j<nev; j++) printf("%15.12f\n", eval[j]);
	    printf("======================================================\n");

	    printf("========= NEW corrected eigenvalue ===========\n");
	    for(j=0; j<nev_new; j++) printf("%15.12f\n", eval_new[j]);
	    printf("======================================================\n");
	}

	//=============================================================================================
	//=============================================================================================
	//=============================================================================================
	//=============================================================================================
	//=============================================================================================
	//=============================================================================================
	//=============================================================================================
	//=============================================================================================
	//=============================================================================================
	//=============================================================================================
	int nev_new2 = nev;
	double **evec_expand_new2 = (double**)malloc(nev_new2 * sizeof(double*));
	for(j=0; j<nev_new2; j++) evec_expand_new2[j] = (double*)calloc(Alarge->nr, sizeof(double));
	double *eval_new2 = (double *)calloc(nev_new2, sizeof(double));
	Eigen_solver_arpack_dn(Alarge_new_new, Mlarge_new_new, nev_new2, eval_new2, evec_expand_new2);
	Insertion_ascend_sort_dvec_dvecvec(eval_new2, evec_expand_new2, 0, nev_new2-1);
	
	printf("========= NEW NEW corrected eigenvalue ===========\n");
	for(j=0; j<nev_new2; j++) printf("%15.12f\n", eval_new2[j]);
	printf("======================================================\n");
	free(eval_new2);
	for(j=0; j<nev_new2; j++) free(evec_expand_new2[j]);
	free(evec_expand_new2);
	Free_dmatcsr(Alarge_new_new);
	Free_dmatcsr(Mlarge_new_new);
	//=============================================================================================
	//=============================================================================================
	//=============================================================================================
	//=============================================================================================

	new_evec_time += Correction_get_new_evec_multistep(amg, finest_level, coarsest_level, nev, dvec_expand, length, evec_expand, evec);

	for(j=0; j<nev;     j++) {free(evec_expand[j]);     evec_expand[j]     = NULL;}
	for(j=0; j<nev_new; j++) {free(evec_expand_new[j]); evec_expand_new[j] = NULL;};
	free(evec_expand);     evec_expand     = NULL;
	free(evec_expand_new); evec_expand_new = NULL;

	Free_dmatcsr(Alarge); Alarge = NULL;
	Free_dmatcsr(Mlarge); Mlarge = NULL;

	Free_dmatcsr(Alarge_new); Alarge_new = NULL;
	Free_dmatcsr(Mlarge_new); Mlarge_new = NULL;
	
    }
	free(eval_new); eval_new = NULL;
    if(status) printf("amg cycle end.\n");
    //=====================================================================
    for(i=0; i<nev; i++) {free(dvec_amg[i]); dvec_amg[i] = NULL;}
    free(dvec_amg); dvec_amg = NULL;

    free(dvec_expand); dvec_expand = NULL;

    t2 = Get_time();
    if(param.amgeigen_print_level > 0)
    {
	printf("direct eigen     time = %f\n", eigen_time);
        printf("amg linear solve time = %f\n", linear_time);
        printf("expand matrix    time = %f\n", expand_time);
        printf("get new evec     time = %f\n", new_evec_time);
        printf("correction total time = %f\n", t2-t1);
    }
}

static double Correction_solve_linear(multigrid *amg, int current_level, int n, double *dval, double **dvec, amg_param param)
{
    double tb = Get_time();

#if   amg_eigen_linear_solver == 1
    double resi_norm_amg      = -1;
    int    ncycle             =  0;
    printf(".............. Eigen -- linear solver amg .............\n");
    printf("level   j   amgcycle            rn               time \n");
#elif amg_eigen_linear_solver == 2
    double resi_norm_pcg      = -1;
    int    ncycle             =  0;
    int    niter              =  0;
    printf(".................. Eigen -- linear solver pcg ..................\n");
    printf("level   j   pcgiter   amgcycle            rn              time \n");
#endif

    dmatcsr *A = amg->A[current_level];
    dmatcsr *M = amg->M[current_level];

    double *rhs = (double*)calloc(A->nr, sizeof(double));

    double t1, t2;
    int j;
    for(j=0; j<n; j++)
    {
	t1 = Get_time();
	memset(rhs, 0, M->nr*sizeof(double));
	Multi_dmatcsr_dvec(M, dvec[j], rhs);
#if 0
	Scale_dvec(rhs, dval[j], M->nr);
#else
	if(MABS(dval[j]) > MYAMGEPS) Scale_dvec(dvec[j], 1.0/dval[j], M->nr);
#endif

#if   amg_eigen_linear_solver == 1
	Linear_solver_amg(amg, current_level, rhs, dvec[j], param, &resi_norm_amg, &ncycle);
	t2 = Get_time();
	printf(" %2d    %2d     %3d       %18.15f     %f\n", current_level, j, ncycle, resi_norm_amg, t2-t1);
#elif amg_eigen_linear_solver == 2
	Linear_solver_pcg_amg(amg, current_level, rhs, dvec[j], param, NULL, &resi_norm_pcg, &niter, &ncycle);
	t2 = Get_time();
	printf(" %2d    %2d     %3d      %3d       %18.15f     %f\n", current_level, j, niter, ncycle, resi_norm_pcg, t2-t1);
#endif
    }

#if   amg_eigen_linear_solver == 1
    printf(".......................................................\n");
#elif amg_eigen_linear_solver == 2
    printf("................................................................\n");
#endif

    free(rhs); rhs = NULL;

    return Get_time() - tb;
}

static double Correction_get_new_evec_multistep(multigrid *amg, int current_level, int coarsest_level,
	                              int n, double **V, int length, double **evec_expand, double **evec)
{
    double tb = Get_time();
    int j, k;
    int nr  = amg->A[current_level]->nr;
    int nrc = amg->A[coarsest_level]->nr;

    double           **PV    = (double**)malloc(n*  sizeof(double*));
    for(j=0; j<n; j++) PV[j] = (double*) calloc(nr, sizeof(double));

    for(j=0; j<n; j++)
    {
	memset(evec[j], 0, nr*sizeof(double));
	ProlongCoarse2Fine(amg, coarsest_level, current_level, evec_expand[j], PV[j]);
	for(k=0; k<n*length; k++) Sumself_dvec_axpby(V[k], evec_expand[j][nrc+k], evec[j], 1.0, nr);
	Sumself_dvec_axpby(PV[j], 1.0, evec[j], 1.0, nr);
        Normalize_dvec(evec[j], nr);
    }
    for(j=0; j<n; j++) { free(PV[j]); PV[j] = NULL; } free(PV); PV = NULL;

    return Get_time() - tb;
}

static double Correction_expand_matrix_RAhV_VTAhV_multistep(dmatcsr *AL, dmatcsr *Ah, dmatcsr *AH, 
	                                                    multigrid *amg, int current_level, int coarsest_level, 
				                            int n, double **V, int length)
{
    double tb = Get_time();

    int j, k;

    int size = n * length;
    double **  AhV = (double**)malloc(size*sizeof(double*));
    double ** RAhV = (double**)malloc(size*sizeof(double*));
    double **VTAhV = (double**)malloc(size*sizeof(double*));
    for(j=0; j<size; j++)   AhV[j] = (double*)calloc(Ah->nr, sizeof(double));
    for(j=0; j<size; j++)  RAhV[j] = (double*)calloc(AH->nr, sizeof(double));
    for(j=0; j<size; j++) VTAhV[j] = (double*)calloc(size,   sizeof(double));

    for(j=0; j<size; j++) Multi_dmatcsr_dvec(Ah, V[j], AhV[j]); 
    for(j=0; j<size; j++) RestrictFine2Coarse(amg, current_level, coarsest_level, AhV[j], RAhV[j]);
    for(j=0; j<size; j++) for(k=0; k<size; k++) VTAhV[j][k] = Multi_dvec_dvec(V[j], AhV[k], Ah->nr);

    Expand_dmatcsr(AL, size, RAhV, VTAhV);

    for(j=0; j<size; j++) { free(VTAhV[j]); VTAhV[j] = NULL; } free(VTAhV); VTAhV = NULL;
    for(j=0; j<size; j++) { free( RAhV[j]);  RAhV[j] = NULL; } free( RAhV);  RAhV = NULL;
    for(j=0; j<size; j++) { free(  AhV[j]);   AhV[j] = NULL; } free(  AhV);   AhV = NULL;

    return Get_time() - tb;
}
