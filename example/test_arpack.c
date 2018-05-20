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
#include "matrix.h"
#include "arpack_interface.h"
#include "tool.h"

int main()
{
  dmatcsr *A = Read_dmatcsr("../output/Alarge_seq.dat");
  dmatcsr *M = Read_dmatcsr("../output/Mlarge_seq.dat");
  //dmatcsr *A = Read_dmatcsr("../../dat/fem2d_poisson_square/gmg_A_refine7.dat");
  //dmatcsr *M = Read_dmatcsr("../../dat/fem2d_poisson_square/gmg_M_refine7.dat");
  int nev = 13;
  double * eval = (double*) calloc(nev, sizeof(double));
  double **evec = (double**)calloc(nev, sizeof(double*));
  int i = -1;
  for(i = 0; i < nev; ++i) {
    evec[i] = (double*)calloc(A->nr, sizeof(double));
  }
  Eigen_solver_arpack_dn(A, M, nev, eval, evec);
  for(i = 0; i < nev; ++i) {
    printf("eval[%02d] = %20.15f\n", i, eval[i]);
  }
  for(i = 0; i < nev; ++i) {
    free(evec[i]); evec[i] = NULL;
  }
  free(evec); evec = NULL;
  free(eval); eval = NULL;
  Free_dmatcsr(M);
  Free_dmatcsr(A);
  return 0;
}
