#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "io.h"
#include "slepc_interface.h"

int print_rank = 0;
int main()
{
  char Afile[256] = "../../dat/tmp/Alarge_seq.dat";
  char Mfile[256] = "../../dat/tmp/Mlarge_seq.dat";
  dmatcsr *A = Read_dmatcsr(Afile);
  dmatcsr *M = Read_dmatcsr(Mfile);
  Print_dmatcsr(A);

  int i;
  int nev = 3;
  double  *eval = (double *)calloc(nev,  sizeof(double));
  double **evec = (double**)malloc(nev * sizeof(double*));
  for(i=0; i<nev; i++) evec[i] = (double*)calloc(A->nr, sizeof(double));

  Eigen_solver_slepc(A, M, nev, eval, evec);

  printf("\neigenvalues: \n");
  for(i=0; i<nev; i++) printf("%02d  %20.15f\n", i, eval[i]);

  for(i=0; i<nev; i++) {
    free(evec[i]); evec[i] = NULL;
  }
  free(evec); evec = NULL;

  free(eval); eval = NULL;

  Free_dmatcsr(M);
  Free_dmatcsr(A);

  printf("\n");
  return 0;
}
