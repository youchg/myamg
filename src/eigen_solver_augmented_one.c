#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "matrix.h"
#include "linear_algebra.h"
#include "multigrid.h"
#include "eigen_solver.h"
#include "arpack_interface.h"
#include "slepc_interface.h"
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
static double Correction_solve_linear(multigrid *amg, int current_level, double dval, double *dvec, amg_param param);


static double Correction_expand_matrix_RAhV_VTAhV(dmatcsr *AL, dmatcsr *Ah, dmatcsr *AH, 
	                                          multigrid *amg, int current_level, int coarsest_level, 
				                  int n, double **V, double ***RAhV0, double ***VTAhV0);

static double Correction_get_new_evec(multigrid *amg, int current_level, int coarsest_level,
	                              int index, double *V, double **evec_expand, double *evec);

static double Get_vector_A_norm(double *vec, dmatcsr *A);

/*
 * @brief Augmented subspace method for one eigenvector
 *
 * @param amg        algebraic multigrid
 * @param index  eigenvalue index
 * @param eval       eigenvalue
 * @param evec       eigenvector
 * @param param      solving parameters
 */
void Eigen_solver_amg_augmented_one(multigrid *amg, 
                                    int index, double *eval, double *evec, 
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
  dmatcsr *Alarge = Expand_dmatcsr_struct(AH, 1);
  dmatcsr *Mlarge = Expand_dmatcsr_struct(MH, 1);

  double *dvec_amg    = (double*)calloc(A->nr, sizeof(double));

  int      augmented_nev  = index + 1;
  double  *augmented_eval = (double* )calloc(augmented_nev,  sizeof(double));
  double **augmented_evec = (double**)malloc(augmented_nev * sizeof(double*));
  for(i=0; i<augmented_nev; i++) augmented_evec[i] = (double*)calloc(Alarge->nr, sizeof(double));

  //================= amg method for eigenvalue problem =================
  if(status) printf("Doing amg cycle...\n");
  for(m=0; m<niter_outer; m++)
  {
    if(status) printf("\ncorrection: m = %d eig = %f\n", m, *eval);

    if(status) printf("Constructing the argumented eigenvalue problem...\n");
    memcpy(dvec_amg, evec, A->nr*sizeof(double));
    // expand matrix on finest level (0) only 1 dimension
    double **b    = NULL;
    double **beta = NULL;
    expand_time += Correction_expand_matrix_RAhV_VTAhV(Alarge, A, AH, amg, 0, coarsest_level, 1, &dvec_amg, NULL, NULL); 
    expand_time += Correction_expand_matrix_RAhV_VTAhV(Mlarge, M, MH, amg, 0, coarsest_level, 1, &dvec_amg, &b, &beta); 
    if(status) printf("Constructing the argumented eigenvalue problem... done.\n");

    if(status) printf("Solving the argumented eigenvalue problem...\n");
    t3 = Get_time(); 
    memset(augmented_eval, 0, augmented_nev*sizeof(double));
    for(j=0; j<augmented_nev; j++) memset(augmented_evec[j], 0, Alarge->nr*sizeof(double));
    //Eigen_solver_arpack_dn(Alarge, Mlarge, augmented_nev, augmented_eval, augmented_evec);
    Eigen_solver_slepc(Alarge, Mlarge, augmented_nev, augmented_eval, augmented_evec);
    Insertion_ascend_sort_dvec_dvecvec(augmented_eval, augmented_evec, 0, augmented_nev-1);
    t4 = Get_time();
    eigen_time += t4-t3;
    if(status) printf("Solving the argumented eigenvalue problem... done.\n");

    if(status) printf("Updating approximate eigenvectors...\n");
    // pick the target eigenvector
    int    max_proj_index = -1;
    double max_proj       = 0.0;
    int    k              = -1;
    double *proj = (double*)calloc(augmented_nev, sizeof(double));
    for(j = 0; j < augmented_nev; ++j) {
      for(k = 0; k < AH->nr; ++k) {
        proj[j] += augmented_evec[j][k] * b[0][k];
      }
      printf("@@ proj = %20.15f + %20.15f * %20.15f = %20.15f\n", proj[j], augmented_evec[j][k], beta[0][0], proj[j] + augmented_evec[j][k] * beta[0][0]);
      proj[j] += augmented_evec[j][k] * beta[0][0];
      if(fabs(proj[j]) > max_proj) {
        max_proj       = fabs(proj[j]);
        max_proj_index = j;
      }
    }

    free(b[0]); b[0] = NULL;
    free(b);    b    = NULL;
    free(beta[0]); beta[0] = NULL;
    free(beta);    beta    = NULL;
    new_evec_time += Correction_get_new_evec(amg, 0, coarsest_level, max_proj_index, dvec_amg, augmented_evec, evec);
    *eval = augmented_eval[max_proj_index];
    if(status) printf("Updating approximate eigenvectors... done.\n");

    if(status) printf("Solving linear system...\n");
    linear_time += Correction_solve_linear(amg, 0, *eval, evec, param);
    double scale = Get_vector_A_norm(evec, M);
    assert(fabs(scale) > 1.0e-15);
    Scale_dvec(evec, 1.0/scale, M->nr);
    if(status) printf("Solving linear system... done.\n");

    if(1)
    {
      printf("========= corrected eigenvalue on level %d ===========\n", finest_level);
      for(j = 0; j < augmented_nev; ++j) {
        printf("  %2d  %20.15f  %20.15f", j, augmented_eval[j], proj[j]);
        if(j == max_proj_index) {
          printf(" (max proj) ");
        }
        printf("\n");
      }
      printf("\n");
      printf("======================================================\n");
    }

    free(proj); proj = NULL;
  }
  if(status) printf("Doing amg cycle... done.\n");
  //=====================================================================

  for(i=0; i<augmented_nev; i++) { free(augmented_evec[i]); augmented_evec[i] = NULL; }
  free(augmented_evec);    augmented_evec = NULL;
  free(augmented_eval);    augmented_eval = NULL;
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

static double Correction_solve_linear(multigrid *amg, int current_level, double dval, double *dvec, amg_param param)
{
    double tb = Get_time();

#if   amg_eigen_linear_solver == 1
    double resi_norm_amg      = -1;
    int    ncycle             =  0;
    printf(".............. Eigen -- linear solver amg .............\n");
    printf("level     amgcycle            rn               time \n");
#elif amg_eigen_linear_solver == 2
    double resi_norm_pcg      = -1;
    int    ncycle             =  0;
    int    niter              =  0;
    printf(".................. Eigen -- linear solver pcg ..................\n");
    printf("level     pcgiter   amgcycle            rn              time \n");
#endif

    dmatcsr *A = amg->A[current_level];
    dmatcsr *M = amg->M[current_level];

    double *rhs = (double*)calloc(A->nr, sizeof(double));

    double t1, t2;
	t1 = Get_time();
	memset(rhs, 0, M->nr*sizeof(double));
	Multi_dmatcsr_dvec(M, dvec, rhs);
	if(MABS(dval) > MYAMGEPS) Scale_dvec(dvec, 1.0/dval, M->nr);

#if   amg_eigen_linear_solver == 1
	Linear_solver_amg(amg, current_level, rhs, dvec, param, &resi_norm_amg, &ncycle);
	t2 = Get_time();
	printf(" %2d        %3d       %18.15f     %f\n", current_level, ncycle, resi_norm_amg, t2-t1);
#elif amg_eigen_linear_solver == 2
	Linear_solver_pcg_amg(amg, current_level, rhs, dvec, param, NULL, &resi_norm_pcg, &niter, &ncycle);
	t2 = Get_time();
	printf(" %2d        %3d      %3d       %18.15f     %f\n", current_level, niter, ncycle, resi_norm_pcg, t2-t1);
#endif

#if   amg_eigen_linear_solver == 1
    printf(".......................................................\n");
#elif amg_eigen_linear_solver == 2
    printf("................................................................\n");
#endif

    free(rhs); rhs = NULL;

    return Get_time() - tb;
}

// evec = P evec_expand(0:nrc-1) + evec_expand(nrc) * V
static double Correction_get_new_evec(multigrid *amg, int current_level, int coarsest_level,
	                              int index, double *V, double **evec_expand, double *evec)
{
    double tb = Get_time();
    int nr  = amg->A[current_level]->nr;
    int nrc = amg->A[coarsest_level]->nr;

    double *PV = (double*)calloc(nr,  sizeof(double));

    memset(evec, 0, nr*sizeof(double));
    ProlongCoarse2Fine(amg, coarsest_level, current_level, evec_expand[index], PV);
    Sumself_dvec_axpby(V, evec_expand[index][nrc], evec, 1.0, nr);
    Sumself_dvec_axpby(PV, 1.0, evec, 1.0, nr);
    //Normalize_dvec(evec, nr);
    double scale = Get_vector_A_norm(evec, amg->M[current_level]);
    assert(fabs(scale) > 1.0e-15);
    Scale_dvec(evec, 1.0/scale, nr);

    free(PV); PV = NULL;

    return Get_time() - tb;
}

static double Correction_expand_matrix_RAhV_VTAhV(dmatcsr *AL, dmatcsr *Ah, dmatcsr *AH, 
	                                          multigrid *amg, int current_level, int coarsest_level, 
				                  int n, double **V, double ***RAhV0, double ***VTAhV0)
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

    if(NULL == VTAhV0) {
      for(j=0; j<n; j++) { free(VTAhV[j]); VTAhV[j] = NULL; } free(VTAhV); VTAhV = NULL;
    } else {
      *VTAhV0 = VTAhV;
    }
    if(NULL == RAhV0) {
      for(j=0; j<n; j++) { free( RAhV[j]);  RAhV[j] = NULL; } free( RAhV);  RAhV = NULL;
    } else {
      *RAhV0 = RAhV;
    }
    for(j=0; j<n; j++) { free(  AhV[j]);   AhV[j] = NULL; } free(  AhV);   AhV = NULL;

    return Get_time() - tb;
}

static double Get_vector_A_norm(double *vec, dmatcsr *A)
{
  double *tmp = (double*)calloc(A->nr, sizeof(double));
  Multi_dmatcsr_dvec(A, vec, tmp);
  double norm = Multi_dvec_dvec(vec, tmp, A->nr);
  free(tmp); tmp = NULL;
  return sqrt(norm);
}
