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

#define direct_nev        32
#define tol_correction    1e-09
#define nmax_correction   2

#define filename_prefix     "/home/ycg/Software/youchg/dat/fem2d_poisson_square/"
#define gmg_finest_level    7
#define gmg_coarsest_level  4

int print_rank = 0;
int main(int argc, char* argv[])
{
  int nev = 5;
  int augmented_index = 3;
  assert(nev >= augmented_index+1);

  int i = 0;

  double eval_given[direct_nev] = {
#if gmg_finest_level == 7
         19.742181652312652, // 00
         49.360802349344290, // 01
         49.367944185095041, // 02
         79.004391701675502, // 03
         98.754512911503113, // 04
         98.754533209294024, // 05
        128.394168556935540, // 06
        128.454367342076296, // 07
        167.940430957701437, // 08
        167.944318594906576, // 09
        177.893344633664412, // 10
        197.653679181721060, // 11
        197.654155507478691, // 12
        247.074076746612491, // 13
        247.310545615300498, // 14
        256.969593559649127, // 15
        256.969611929320706, // 16
        286.724147309196951, // 17
        286.745315764746238, // 18
        316.584651924283833, // 19
        336.361880577667876, // 20
        336.366327347913796, // 21
        365.888058939655650, // 22
        365.891559701172412, // 23
        395.720359990433792, // 24
        395.720630471374477, // 25
        405.555049749392197, // 26
        406.207452347479659, // 27
        445.426637813072546, // 28
        445.491139460087368, // 29
        494.768409981298817, // 30
        494.768519621335997, // 31
#elif gmg_finest_level == 11
         19.739351152446577, // 00
         49.348398772994052, // 01
         49.348426661655587, // 02
         78.957543954641395, // 03
         98.696926072478092, // 04
         98.696926072787591, // 05
        128.306055963593963, // 06
        128.306290898228610, // 07
        167.784999753374791, // 08
        167.785014924843836, // 09
        177.654996436879685, // 10
        197.394417265854258, // 11
        197.394417273111628, // 12
        246.743050204326266, // 13
        246.743972243036808, // 14
        256.612819086151660, // 15
        256.612819086431614, // 16
        286.222396771513957, // 17
        286.222479299660563, // 18
        315.832405282174534, // 19
        335.571880346205262, // 20
        335.571880413890710, // 21
        365.180563509384399, // 22
        365.180577155566198, // 23
        394.790444748240247, // 24
        394.790444752357189, // 25
        404.659991649028711, // 26
        404.662532642826989, // 27
        444.140187371669356, // 28
        444.140438372226185, // 29
        493.488516543021433, // 30
        493.488516544614299, // 31
#endif
  };

#if 1
  if(argc < 2) {
    printf("Too few arguments!\n");
    exit(0);
  }
#endif

  /* ==========================================================================
   * ====================  Initial Parameters  ================================
   * ========================================================================== */
  double tb_init = Get_time();

  amg_param param;
  Init_amg_param_argv(argc, argv, &param, NULL, NULL, NULL);
  Print_amg_param(param);

#if 1
  int actual_level = gmg_finest_level - gmg_coarsest_level + 1;
  multigrid *amg = Init_amg(actual_level, 2); // 2: generalized eigenvalue problem
  {
    for(i = 0; i < actual_level; i++) {
      char str[10];
      char Afile[256];
      strcpy(Afile, filename_prefix"gmg_A_refine");
      sprintf(str, "%d", gmg_finest_level - i);
      strcat(Afile, str);
      strcat(Afile, ".dat");
      printf("Read A[%2d] from %s\n", i, Afile);
      free(amg->A[i]);
      amg->A[i] = Read_dmatcsr(Afile);
      Print_dmatcsr(amg->A[i]);
      //Remove_zero_dmatcsr(A);
      //Print_dmatcsr(A);
    }
    printf("\n");
    for(i = 0; i < actual_level; i++) {
      char str[10];
      char Mfile[256];
      strcpy(Mfile, filename_prefix"gmg_M_refine");
      sprintf(str, "%d", gmg_finest_level - i);
      strcat(Mfile, str);
      strcat(Mfile, ".dat");
      printf("Read M[%2d] from %s\n", i, Mfile);
      free(amg->M[i]);
      amg->M[i] = Read_dmatcsr(Mfile);
      Print_dmatcsr(amg->M[i]);
    }
    printf("\n");
    for(i = 0; i < actual_level - 1; i++) {
      char str[10];
      char Pfile[256];
      strcpy(Pfile, filename_prefix"gmg_P_refine");
      sprintf(str, "%d", gmg_finest_level - i);
      strcat(Pfile, str);
      strcat(Pfile, ".dat");
      printf("Read P[%2d] from %s\n", i, Pfile);
      free(amg->P[i]);
      amg->P[i] = Read_dmatcsr(Pfile);
      Print_dmatcsr(amg->P[i]);
    }
    printf("\n");
    for(i = 0; i < actual_level - 1; i++) {
      char str[10];
      char Rfile[256];
      strcpy(Rfile, filename_prefix"gmg_R_refine");
      sprintf(str, "%d", gmg_finest_level - i);
      strcat(Rfile, str);
      strcat(Rfile, ".dat");
      printf("Read R[%2d] from %s\n", i, Rfile);
      free(amg->R[i]);
      amg->R[i] = Read_dmatcsr(Rfile);
      Print_dmatcsr(amg->R[i]);
    }
  }
  amg->actual_level = actual_level;
  Print_amg(amg);
#else
  char str[10];
  char Afile[256];
  char Mfile[256];
  strcpy(Afile, filename_prefix"gmg_A_refine");
  strcpy(Mfile, filename_prefix"gmg_M_refine");
  sprintf(str, "%d", gmg_finest_level - i);
  strcat(Afile, str);
  strcat(Mfile, str);
  strcat(Afile, ".dat");
  strcat(Mfile, ".dat");
  printf("Read A[%2d] from %s\n", i, Afile);
  printf("Read M[%2d] from %s\n", i, Mfile);
  dmatcsr *AA = Read_dmatcsr(Afile);
  dmatcsr *MM = Read_dmatcsr(Mfile);
  Print_dmatcsr(AA);
  Print_dmatcsr(MM);
  Remove_zero_dmatcsr(AA);
  Remove_zero_dmatcsr(MM);
  Print_dmatcsr(AA);
  Print_dmatcsr(MM);
  multigrid *amg = Build_amg(AA, MM, param.max_level);
  Free_dmatcsr(AA);
  Free_dmatcsr(MM);
  double tb_setup = Get_time();
  Setup_phase(amg, param);
  double te_setup = Get_time();
  Print_amg(amg);
  printf("setup phase time: %f\n", te_setup-tb_setup);
#endif

  double te_init = Get_time();
  printf("\ninit time: %f\n", te_init-tb_init);

#if 1
  double   tb_direct_all   = Get_time();
  double  *eval_direct_all = (double*) calloc(direct_nev, sizeof(double));
  double **evec_direct_all = (double**)malloc(direct_nev* sizeof(double*));
  int      direct_level    = amg->actual_level - 1;
  for(i=0; i<direct_nev; i++) evec_direct_all[i] = (double*)calloc(amg->A[direct_level]->nc, sizeof(double));
  printf("\ncalling direct method all...\n");
  Eigen_solver_arpack_dn(amg->A[direct_level], amg->M[direct_level], direct_nev, eval_direct_all, evec_direct_all);
  printf("================= direct all result ===================\n");
  for(i=0; i<direct_nev; i++) printf("%2d: %20.15f\n", i, eval_direct_all[i]);
  printf("===================================================\n");
  for(i=0; i<direct_nev; i++) { free(evec_direct_all[i]); } 
  free(evec_direct_all); evec_direct_all = NULL;
  free(eval_direct_all); eval_direct_all = NULL;
  double te_direct_all = Get_time();
  printf("direct eigen time: %f\n", te_direct_all - tb_direct_all);
#endif

  dmatcsr *A = amg->A[0];

  double  *total_error = (double*)calloc(nmax_correction, sizeof(double));
  double  *corre_time  = (double*)calloc(nmax_correction, sizeof(double));

  /* solve eigenvalue problem */
  double tb_correction_amg = Get_time();
  double  *eval_amg = (double*) calloc(nev, sizeof(double));
  double **evec_amg = (double**)malloc(nev* sizeof(double*));
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
  param_eigen.amgsolver_max_cycle  = 100;
  param_eigen.pcg_amg_max_iter     = 1;
  //Eigen_solver_amg_augmented(amg, nev, eval_amg, evec_amg, param_eigen);
  Eigen_solver_amg_nested(amg, nev, eval_amg, evec_amg, param_eigen);
  te_amg = Get_time();
  printf("* 0 * approximate eigenvalue: \n");/* show the result */
  corre_time[0] = te_amg - tb_amg;

  total_error[0] = eval_amg[augmented_index] - eval_given[augmented_index];
  printf("%2d: %20.15f %20.15f\n", augmented_index, eval_amg[augmented_index], total_error[0]);
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
    //Eigen_solver_amg_augmented(amg, nev, eval_amg, evec_amg, param_eigen);
    Eigen_solver_amg_augmented_one(amg, augmented_index, &eval_amg[augmented_index], evec_amg[augmented_index], param_eigen);
    te_amg = Get_time();
    corre_time[i] = te_amg - tb_amg;
    total_error[i] = eval_amg[augmented_index] - eval_given[augmented_index];
    printf("%2d: %20.15f %20.15f\n", augmented_index, eval_amg[augmented_index], total_error[i]);
    printf("correction %2d time : %20.15f\n", i, corre_time[i]);
    printf("correction %2d error: %20.15f\n", i, total_error[i]);
    ncorrection++;
    if(fabs(total_error[i]) < tol_correction) break;
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

  Free_multigrid(amg);
  return 0;
}

