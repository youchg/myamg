#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "matrix.h"
#include "linear_algebra.h"
#include "multigrid.h"
#include "eigen_solver.h"
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


static double Correction_expand_matrix_RAhV_VTAhV(dmatcsr *AL, dmatcsr *Ah, dmatcsr *AH, 
	                                          multigrid *amg, int current_level, int coarsest_level, 
				                  int n, double **V);

static double Correction_get_new_evec(multigrid *amg, int current_level, int coarsest_level,
	                              int n, double **V, double **evec_expand, double **evec);


void Eigen_solver_amg_augmented(multigrid *amg, 
                                 int nev, double *eval, double **evec, 
                                 amg_param param)
{
  //printf("Begin Eigen_solver_amg(): memory use (MB): %f\n", Get_memory());
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

  dmatcsr *A      = amg->A[finest_level];
  dmatcsr *M      = amg->M[finest_level];
  dmatcsr *AH     = amg->A[coarsest_level];
  dmatcsr *MH     = amg->M[coarsest_level]; 

  //=====================================================================
  dmatcsr *Alarge = Expand_dmatcsr_struct(AH, nev);
  dmatcsr *Mlarge = Expand_dmatcsr_struct(MH, nev);

  double **dvec_amg    = (double**)malloc(nev * sizeof(double*));
  for(i=0; i<nev; i++) dvec_amg[i]    = (double*)calloc(     A->nr, sizeof(double));
  double **evec_expand = (double**)malloc(nev * sizeof(double*));
  for(i=0; i<nev; i++) evec_expand[i] = (double*)calloc(Alarge->nr, sizeof(double));

  //================= amg method for eigenvalue problem =================
  if(status) printf("Doing amg cycle...\n");
  for(m=0; m<niter_outer; m++)
  {
    if(status) printf("\ncorrection: m = %d\n", m);

    if(status) printf("Constructing the argumented eigenvalue problem...\n");
    for(j=0; j<nev; j++) memcpy(dvec_amg[j], evec[j], A->nr*sizeof(double));
    expand_time += Correction_expand_matrix_RAhV_VTAhV(Alarge, A, AH, amg, 0, coarsest_level, nev, dvec_amg); 
    expand_time += Correction_expand_matrix_RAhV_VTAhV(Mlarge, M, MH, amg, 0, coarsest_level, nev, dvec_amg); 
    if(status) printf("Constructing the argumented eigenvalue problem... done.\n");

    if(status) printf("Solving the argumented eigenvalue problem...\n");
    t3 = Get_time(); 
    memset(eval, 0, nev*sizeof(double));
    for(j=0; j<nev; j++) memset(evec_expand[j], 0, Alarge->nr*sizeof(double));
    Eigen_solver_arpack_dn(Alarge, Mlarge, nev, eval, evec_expand);
    Insertion_ascend_sort_dvec_dvecvec(eval, evec_expand, 0, nev-1);
    t4 = Get_time();
    eigen_time += t4-t3;
    if(status) printf("Solving the argumented eigenvalue problem... done.\n");

    if(status) printf("Updating approximate eigenvectors...\n");
    new_evec_time += Correction_get_new_evec(amg, 0, coarsest_level, nev, dvec_amg, evec_expand, evec);
    if(status) printf("Updating approximate eigenvectors... done.\n");

    if(status) printf("Solving linear system...\n");
    linear_time += Correction_solve_linear(amg, 0, nev, eval, evec, param);
    if(status) printf("Solving linear system... done.\n");

    if(1)
    {
      printf("========= corrected eigenvalue on level %d ===========\n", finest_level);
      for(j=0; j<nev; j++) printf(" %02d  %20.15f\n", j, eval[j]);
      printf("\n");
      printf("======================================================\n");
    }
  }
  if(status) printf("Doing amg cycle... done.\n");
  //=====================================================================

  for(i=0; i<nev; i++) { free(evec_expand[i]); evec_expand[i] = NULL; }
  for(i=0; i<nev; i++) { free(dvec_amg[i]);    dvec_amg[i]    = NULL; }
  free(evec_expand);    evec_expand = NULL;
  free(dvec_amg);       dvec_amg    = NULL;
  Free_dmatcsr(Alarge); Alarge      = NULL;
  Free_dmatcsr(Mlarge); Mlarge      = NULL;

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
	//int print_j = 1;
	//if(j == print_j) Print_dvec(rhs, M->nr);
#if 0
	Scale_dvec(rhs, dval[j], M->nr);
#else
	//if(j == print_j) Print_dvec(dvec[j], M->nr);
	//if(j == print_j) Print_dvec(&dval[j], 1);
	if(MABS(dval[j]) > EPS) Scale_dvec(dvec[j], 1.0/dval[j], M->nr);
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

static double Correction_get_new_evec(multigrid *amg, int current_level, int coarsest_level,
	                              int n, double **V, double **evec_expand, double **evec)
{
    double tb = Get_time();
    int j, k;
    int nr  = amg->A[current_level]->nr;
    int nrc = amg->A[coarsest_level]->nr;

    double           **PV    = (double**)malloc(n*   sizeof(double*));
    for(j=0; j<n; j++) PV[j] = (double*) calloc(nr,  sizeof(double));

    for(j=0; j<n; j++)
    {
	memset(evec[j], 0, nr*sizeof(double));
	ProlongCoarse2Fine(amg, coarsest_level, current_level, evec_expand[j], PV[j]);
	for(k=0; k<n; k++) Sumself_dvec_axpby(V[k], evec_expand[j][nrc+k], evec[j], 1.0, nr);
	Sumself_dvec_axpby(PV[j], 1.0, evec[j], 1.0, nr);
	//if(j == 1) Print_dvec(evec[j], nr);
        Normalize_dvec(evec[j], nr);
    }
    for(j=0; j<n; j++) { free(PV[j]); PV[j] = NULL; } free(PV); PV = NULL;

    return Get_time() - tb;
}

static double Correction_expand_matrix_RAhV_VTAhV(dmatcsr *AL, dmatcsr *Ah, dmatcsr *AH, 
	                                          multigrid *amg, int current_level, int coarsest_level, 
				                  int n, double **V)
{
    double tb = Get_time();

    int j, k;

    double **  AhV = (double**)malloc(n*sizeof(double*));
    double ** RAhV = (double**)malloc(n*sizeof(double*));
    double **VTAhV = (double**)malloc(n*sizeof(double*));
    for(j=0; j<n; j++)   AhV[j] = (double*)calloc(Ah->nr, sizeof(double));
    for(j=0; j<n; j++)  RAhV[j] = (double*)calloc(AH->nr, sizeof(double));
    for(j=0; j<n; j++) VTAhV[j] = (double*)calloc(n,      sizeof(double));

    for(j=0; j<n; j++) Multi_dmatcsr_dvec(Ah, V[j], AhV[j]); 
    for(j=0; j<n; j++) RestrictFine2Coarse(amg, current_level, coarsest_level, AhV[j], RAhV[j]);
    for(j=0; j<n; j++) for(k=0; k<n; k++) VTAhV[j][k] = Multi_dvec_dvec(V[j], AhV[k], Ah->nr);

    Expand_dmatcsr(AL, n, RAhV, VTAhV);

    for(j=0; j<n; j++) { free(VTAhV[j]); VTAhV[j] = NULL; } free(VTAhV); VTAhV = NULL;
    for(j=0; j<n; j++) { free( RAhV[j]);  RAhV[j] = NULL; } free( RAhV);  RAhV = NULL;
    for(j=0; j<n; j++) { free(  AhV[j]);   AhV[j] = NULL; } free(  AhV);   AhV = NULL;

    return Get_time() - tb;
}
