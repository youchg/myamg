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
#define amg_eigen_expand_type   1 //1:RAhV+VTAhV  2:RAhV+VTRTAHRV  3:AHRV+VTAhV   4:AHRV+VTRTAHRV

/* using 'dvec':
 * 1. to form rhsdvec' 
 * 2. as initial approximation 
 * 3. as final   approximation
 * */
static double Correction_solve_linear(multigrid *amg, int current_level, int n, double *dval, double **dvec, amg_param param);


#if   amg_eigen_expand_type == 1
static double Correction_expand_matrix_RAhV_VTAhV(dmatcsr *AL, dmatcsr *Ah, dmatcsr *AH, 
	                                          multigrid *amg, int current_level, int coarsest_level, 
				                  int n, double **V);
#elif amg_eigen_expand_type == 2
static double Correction_expand_matrix_RAhV_VTRTAHRV(dmatcsr *AL, dmatcsr *Ah, dmatcsr *AH, 
	                                             multigrid *amg, int current_level, int coarsest_level, 
				                     int n, double **V);
#elif amg_eigen_expand_type == 3
static double Correction_expand_matrix_AHRV_VTAhV(dmatcsr *AL, dmatcsr *Ah, dmatcsr *AH, 
	                                          multigrid *amg, int current_level, int coarsest_level, 
				                  int n, double **V);
#elif amg_eigen_expand_type == 4
static double Correction_expand_matrix_AHRV_VTRTAHRV(dmatcsr *AL, dmatcsr *Ah, dmatcsr *AH, 
	                                             multigrid *amg, int current_level, int coarsest_level, 
				                     int n, double **V);
#endif

static double Correction_get_new_evec(multigrid *amg, int current_level, int coarsest_level,
	                              int n, double **V, double **evec_expand, double **evec);


/* 当 init_approx_level >= 0 时，eval, evec 必须给定初值，并且不能为 0;
   当 init_approx_level == 0 时，is_init_approx_level_correction 必须大于 0
*/
//Problem: when already have a finer result such as 0-level approximation, but
//want to correct it from a coarser level, there maybe some problem
//----20160301
void Eigen_solver_amg(multigrid *amg, 
                      int nev, double *eval, 
                      double **evec, int init_approx_level, int is_init_approx_level_correction, 
                      amg_param param)
{
    double t1, t2, t3, t4;
    double eigen_time    = 0;
    double linear_time   = 0;
    double expand_time   = 0;
    double new_evec_time = 0;
    t1 = Get_time();

    int status = 0;

    int i, j, m;
    int nlevel         = amg->actual_level;
    int niter_outer    = param.amgeigen_nouter_iter;
    int finest_level   = 0;
    int coarsest_level = param.amgeigen_coarsest_level;
    if(coarsest_level <= 0) coarsest_level += nlevel-1;

    dmatcsr *matA;
    dmatcsr *matM;
    dmatcsr *matP;

    dmatcsr *A      = amg->A[finest_level];
    dmatcsr *AH     = amg->A[coarsest_level];
    dmatcsr *MH     = amg->M[coarsest_level]; 
    //=====================================================================
    int amg_correction_begin_level;
    if(init_approx_level < 0)
    {
        t3 = Get_time();

        Eigen_solver_arpack_dn(AH, MH, nev, eval, evec);
        //for(j=0; j<nev; j++) Normalize_dvec(evec[j], AH->nr);

        t4 = Get_time();
        eigen_time += t4-t3;
        if(1)
        {
	    printf("Solving eigenvalue problem by Eigen_solver_arpack_dn() on the coarest level...\n");
            printf("========= coarest eigenvalue on level %d =========\n", coarsest_level);
            for(i=0; i<nev; i++) printf("%d: %15.12f\n", i, eval[i]);
	    printf("==================================================\n");

#if 0
            for(i=0; i<nev; i++) 
                printf("norm evec[%d] = %f\n", i, Get_dvec_2norm(evec[i], AH->nr));
            for(j=0; j<AH->nr; j++) 
                printf("evec[0][%d] = %f\n", j, evec[0][j]);

            double *Atmp = (double*)calloc(AH->nr, sizeof(double));
            double *Mtmp = (double*)calloc(AH->nr, sizeof(double));
            for(i=0; i<nev; i++) {
              Multi_dmatcsr_dvec(MH, evec[i], Mtmp);
              for(j=0; j<nev; j++) {
                printf("%f ", Multi_dvec_dvec(evec[j], Mtmp, AH->nr));
              }
              printf("\n");
            }
            printf("\n");
            for(i=0; i<nev; i++) {
              Multi_dmatcsr_dvec(AH, evec[i], Atmp);
              for(j=0; j<nev; j++) {
                printf("%f ", Multi_dvec_dvec(evec[j], Atmp, AH->nr));
              }
              printf("\n");
            }
            for(i=0; i<nev; i++) {
              Multi_dmatcsr_dvec(AH, evec[i], Atmp);
              Multi_dmatcsr_dvec(MH, evec[i], Mtmp);
                printf("%f ", Multi_dvec_dvec(evec[i], Atmp, AH->nr) / Multi_dvec_dvec(evec[i], Mtmp, AH->nr));
              printf("\n");
            }
            free(Atmp);
            free(Mtmp);
#endif
        }
    }

    //amg_correction_begin_level = (init_approx_level < 0) ? coarsest_level-1 : init_approx_level;
    if(init_approx_level < 0)
	amg_correction_begin_level = coarsest_level -1;
    else if(is_init_approx_level_correction > 0)
	amg_correction_begin_level = init_approx_level;
    else
	amg_correction_begin_level = init_approx_level - 1;
    /*排除：在第0层上给定初值，但是不需要矫正。以后可以改成这种情况直接返回。*/
    assert(init_approx_level!=0 || is_init_approx_level_correction > 0);

    dmatcsr *Alarge = Expand_dmatcsr_struct(AH, nev);
    dmatcsr *Mlarge = Expand_dmatcsr_struct(MH, nev);
    
    double **dvec_amg    = (double**)malloc(nev * sizeof(double*));
    double **evec_expand = (double**)malloc(nev * sizeof(double*));
    for(i=0; i<nev; i++) dvec_amg[i]    = (double*)calloc(     A->nr, sizeof(double));
    for(i=0; i<nev; i++) evec_expand[i] = (double*)calloc(Alarge->nr, sizeof(double));

    //================= amg method for eigenvalue problem =================
    if(status) printf("Begin amg cycle...\n");

    for(i=amg_correction_begin_level; i>=finest_level; i--)
    {
        if(status) printf("level = %d, nlevel = %d, coarsest_level = %d\n", i, nlevel, coarsest_level);
        
	matA = amg->A[i];
	matM = amg->M[i];
	matP = amg->P[i];

        for( m=0; m<niter_outer; m++ )
        {
            if(status) printf("correction: m = %d\n", m);
	    if(status) printf("Solving linear system on current level %d ...\n", i);

	    for( j=0; j<nev; j++ )
            {
		if(m == 0)
		{
		    if(init_approx_level < 0)
		    {
			Multi_dmatcsr_dvec(matP, evec[j], dvec_amg[j]);
                //printf("norm dvec_amg[%d] = %f\n", j, Get_dvec_2norm(dvec_amg[j], matP->nr));
		    }
		    else
		    {
			if(i==amg_correction_begin_level && is_init_approx_level_correction>0)
			    memcpy(dvec_amg[j], evec[j], matA->nr*sizeof(double));
			else
			    Multi_dmatcsr_dvec(matP, evec[j], dvec_amg[j]);
		    }
		}
		else
		{
		    memcpy(dvec_amg[j], evec[j], matA->nr*sizeof(double));
		}

	    }

	    linear_time += Correction_solve_linear(amg, i, nev, eval, dvec_amg, param);
            for(j=0; j<nev; j++)
                printf("after linear norm dvec_amg[%d] = %f\n", j, Get_dvec_2norm(dvec_amg[j], matA->nr));
#if   amg_eigen_expand_type == 1
	    expand_time += Correction_expand_matrix_RAhV_VTAhV(Alarge, matA, AH, amg, i, coarsest_level, nev, dvec_amg); 
	    expand_time += Correction_expand_matrix_RAhV_VTAhV(Mlarge, matM, MH, amg, i, coarsest_level, nev, dvec_amg); 
#elif amg_eigen_expand_type == 2
	    expand_time += Correction_expand_matrix_RAhV_VTRTAHRV(Alarge, matA, AH, amg, i, coarsest_level, nev, dvec_amg); 
	    expand_time += Correction_expand_matrix_RAhV_VTRTAHRV(Mlarge, matM, MH, amg, i, coarsest_level, nev, dvec_amg); 
#elif amg_eigen_expand_type == 3
	    expand_time += Correction_expand_matrix_AHRV_VTAhV(Alarge, matA, AH, amg, i, coarsest_level, nev, dvec_amg); 
	    expand_time += Correction_expand_matrix_AHRV_VTAhV(Mlarge, matM, MH, amg, i, coarsest_level, nev, dvec_amg); 
#elif amg_eigen_expand_type == 4
	    expand_time += Correction_expand_matrix_AHRV_VTRTAHRV(Alarge, matA, AH, amg, i, coarsest_level, nev, dvec_amg); 
	    expand_time += Correction_expand_matrix_AHRV_VTRTAHRV(Mlarge, matM, MH, amg, i, coarsest_level, nev, dvec_amg); 
#endif
	    //if(i==0)
	    {
	    Write_dmatcsr_csr(Alarge, "../output/Alarge_seq.dat");
	    Write_dmatcsr_csr(Mlarge, "../output/Mlarge_seq.dat");
	    //exit(-1);
	    }
   
	    if(status) printf("Solving corrected eigenvalue problem...\n");
	    t3 = Get_time(); 
            memset(eval, 0, nev*sizeof(double));
            for(j=0; j<nev; j++) memset(evec_expand[j], 0, Alarge->nr*sizeof(double));
	    Eigen_solver_arpack_dn(Alarge, Mlarge, nev, eval, evec_expand);
	    Insertion_ascend_sort_dvec_dvecvec(eval, evec_expand, 0, nev-1);
	    t4 = Get_time();
	    eigen_time += t4-t3;
	    
	    //if(status)
	    if(1)
	    {
	        printf("========= corrected eigenvalue on level %d ===========\n", i);
                for(j=0; j<nev; j++) printf("%15.12f\n", eval[j]);
                printf("\n");
                printf("======================================================\n");
	    }

	    new_evec_time += Correction_get_new_evec(amg, i, coarsest_level, nev, dvec_amg, evec_expand, evec);
            exit(-1);
        }
    }

    if(status) printf("amg cycle end.\n");
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

void Eigen_solver_amg_nested(multigrid *amg, 
                             int nev, double *eval, double **evec,
                             amg_param param)
{
      Eigen_solver_amg(amg, nev, eval, evec, -1, 1, param);
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

#if   amg_eigen_expand_type == 1
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
    printf("\n\n");
    for(j=0; j<n; j++) for(k=0; k<n; k++) printf("VTAhV[%d][%d] = %f\n", j, k, VTAhV[j][k]);
    printf("\n\n");

    printf("\n\n RAhV = \n");
    for(k=0; k<AH->nr; k++) {
      for(j=0; j<n; j++) {
        printf("%g ", RAhV[j][k]);
      }
      printf("\n");
    }
    printf("\n\n");

    Expand_dmatcsr(AL, n, RAhV, VTAhV);

    for(j=0; j<n; j++) { free(VTAhV[j]); VTAhV[j] = NULL; } free(VTAhV); VTAhV = NULL;
    for(j=0; j<n; j++) { free( RAhV[j]);  RAhV[j] = NULL; } free( RAhV);  RAhV = NULL;
    for(j=0; j<n; j++) { free(  AhV[j]);   AhV[j] = NULL; } free(  AhV);   AhV = NULL;

    return Get_time() - tb;
}
#elif amg_eigen_expand_type == 2
static double Correction_expand_matrix_RAhV_VTRTAHRV(dmatcsr *AL, dmatcsr *Ah, dmatcsr *AH, 
	                                             multigrid *amg, int current_level, int coarsest_level, 
				                     int n, double **V)
{
    double tb = Get_time();

    int j, k;

    double **      RV = (double**)malloc(n*sizeof(double*));
    double **     AhV = (double**)malloc(n*sizeof(double*));
    double **    RAhV = (double**)malloc(n*sizeof(double*));
    double **    AHRV = (double**)malloc(n*sizeof(double*));
    double **VTRTAHRV = (double**)malloc(n*sizeof(double*));
    for(j=0; j<n; j++)       RV[j] = (double*)calloc(AH->nr, sizeof(double));
    for(j=0; j<n; j++)      AhV[j] = (double*)calloc(Ah->nr, sizeof(double));
    for(j=0; j<n; j++)     RAhV[j] = (double*)calloc(AH->nr, sizeof(double));
    for(j=0; j<n; j++)     AHRV[j] = (double*)calloc(AH->nr, sizeof(double));
    for(j=0; j<n; j++) VTRTAHRV[j] = (double*)calloc(n,      sizeof(double));

    for(j=0; j<n; j++) Multi_dmatcsr_dvec(Ah, V[j], AhV[j]); 
    for(j=0; j<n; j++) RestrictFine2Coarse(amg, current_level, coarsest_level, AhV[j], RAhV[j]); 
    for(j=0; j<n; j++) RestrictFine2Coarse(amg, current_level, coarsest_level,   V[j],   RV[j]); 
    for(j=0; j<n; j++) Multi_dmatcsr_dvec(AH, RV[j], AHRV[j]); 
    for(j=0; j<n; j++) for(k=0; k<n; k++) VTRTAHRV[j][k] = Multi_dvec_dvec(RV[j], AHRV[k], AH->nr);

    Expand_dmatcsr(AL, n, RAhV, VTRTAHRV);

    for(j=0; j<n; j++) { free(VTRTAHRV[j]); VTRTAHRV[j] = NULL; } free(VTRTAHRV); VTRTAHRV = NULL;
    for(j=0; j<n; j++) { free(    AHRV[j]);     AHRV[j] = NULL; } free(    AHRV);     AHRV = NULL;
    for(j=0; j<n; j++) { free(    RAhV[j]);     RAhV[j] = NULL; } free(    RAhV);     RAhV = NULL;
    for(j=0; j<n; j++) { free(     AhV[j]);      AhV[j] = NULL; } free(     AhV);      AhV = NULL;
    for(j=0; j<n; j++) { free(      RV[j]);       RV[j] = NULL; } free(      RV);       RV = NULL;

    return Get_time() - tb;
}
#elif amg_eigen_expand_type == 3
static double Correction_expand_matrix_AHRV_VTAhV(dmatcsr *AL, dmatcsr *Ah, dmatcsr *AH, 
	                                          multigrid *amg, int current_level, int coarsest_level, 
				                  int n, double **V)
{
    double tb = Get_time();

    int j, k;

    double **  AhV = (double**)malloc(n*sizeof(double*));
    double **   RV = (double**)malloc(n*sizeof(double*));
    double ** AHRV = (double**)malloc(n*sizeof(double*));
    double **VTAhV = (double**)malloc(n*sizeof(double*));
    for(j=0; j<n; j++)   AhV[j] = (double*)calloc(Ah->nr, sizeof(double));
    for(j=0; j<n; j++)    RV[j] = (double*)calloc(AH->nr, sizeof(double));
    for(j=0; j<n; j++)  AHRV[j] = (double*)calloc(AH->nr, sizeof(double));
    for(j=0; j<n; j++) VTAhV[j] = (double*)calloc(n,      sizeof(double));

    for(j=0; j<n; j++) Multi_dmatcsr_dvec(Ah, V[j], AhV[j]); 
    for(j=0; j<n; j++) RestrictFine2Coarse(amg, current_level, coarsest_level, V[j], RV[j]); 
    for(j=0; j<n; j++) Multi_dmatcsr_dvec(AH, RV[j], AHRV[j]); 
    for(j=0; j<n; j++) for(k=0; k<n; k++) VTAhV[j][k] = Multi_dvec_dvec(V[j], AhV[k], Ah->nr);
    //for(j=0; j<n; j++) for(k=0; k<n; k++) VTAV[j][k] = Multi_dvec_dvec(RV[j], AHRV[k], AH->nr);

    Expand_dmatcsr(AL, n, AHRV, VTAhV);

    for(j=0; j<n; j++) { free(VTAhV[j]); VTAhV[j] = NULL; } free(VTAhV); VTAhV = NULL;
    for(j=0; j<n; j++) { free( AHRV[j]);  AHRV[j] = NULL; } free( AHRV);  AHRV = NULL;
    for(j=0; j<n; j++) { free(   RV[j]);    RV[j] = NULL; } free(   RV);    RV = NULL;
    for(j=0; j<n; j++) { free(  AhV[j]);   AhV[j] = NULL; } free(  AhV);   AhV = NULL;

    return Get_time() - tb;
}
#elif amg_eigen_expand_type == 4
static double Correction_expand_matrix_AHRV_VTRTAHRV(dmatcsr *AL, dmatcsr *Ah, dmatcsr *AH, 
	                                             multigrid *amg, int current_level, int coarsest_level, 
				                     int n, double **V)
{
    double tb = Get_time();

    int j, k;

    double **     AhV = (double**)malloc(n*sizeof(double*));
    double **      RV = (double**)malloc(n*sizeof(double*));
    double **    AHRV = (double**)malloc(n*sizeof(double*));
    double **VTRTAHRV = (double**)malloc(n*sizeof(double*));
    for(j=0; j<n; j++)      AhV[j] = (double*)calloc(Ah->nr, sizeof(double));
    for(j=0; j<n; j++)       RV[j] = (double*)calloc(AH->nr, sizeof(double));
    for(j=0; j<n; j++)     AHRV[j] = (double*)calloc(AH->nr, sizeof(double));
    for(j=0; j<n; j++) VTRTAHRV[j] = (double*)calloc(n,      sizeof(double));

    for(j=0; j<n; j++) Multi_dmatcsr_dvec(Ah, V[j], AhV[j]); 
    for(j=0; j<n; j++) RestrictFine2Coarse(amg, current_level, coarsest_level, V[j], RV[j]); 
    for(j=0; j<n; j++) Multi_dmatcsr_dvec(AH, RV[j], AHRV[j]); 
    for(j=0; j<n; j++) for(k=0; k<n; k++) VTRTAHRV[j][k] = Multi_dvec_dvec(RV[j], AHRV[k], AH->nr);

    Expand_dmatcsr(AL, n, AHRV, VTRTAHRV);

    for(j=0; j<n; j++) { free(VTRTAHRV[j]); VTRTAHRV[j] = NULL; } free(VTRTAHRV); VTRTAHRV = NULL;
    for(j=0; j<n; j++) { free(    AHRV[j]);     AHRV[j] = NULL; } free(    AHRV);     AHRV = NULL;
    for(j=0; j<n; j++) { free(      RV[j]);       RV[j] = NULL; } free(      RV);       RV = NULL;
    for(j=0; j<n; j++) { free(     AhV[j]);      AhV[j] = NULL; } free(     AhV);      AhV = NULL;

    return Get_time() - tb;
}
#endif
