#ifdef WITH_SLEPC

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "preprocess.h"
#include "matrix.h"
#include "tool.h"
#include <slepceps.h>

static void Transform_dmatcsr_petscmat(dmatcsr *A, Mat *AA);

void Eigen_solver_slepc(dmatcsr *A, dmatcsr *M, int nev, double *eval, double **evec)
{
  //SlepcInitialize(NULL, NULL, NULL, NULL);
  Mat AA, MM;
  Transform_dmatcsr_petscmat(A, &AA);
  Transform_dmatcsr_petscmat(M, &MM);

  EPS eps;
  ST  st;
  EPSCreate(PETSC_COMM_WORLD, &eps);
  EPSSetOperators(eps, AA, MM);
  EPSSetProblemType(eps, EPS_GNHEP);
  EPSSetDimensions(eps, nev, PETSC_DEFAULT, PETSC_DEFAULT);
  EPSGetST(eps, &st);
  STSetType(st, "sinvert");
  EPSSetTarget(eps, 0.0);
  EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE);
  EPSSolve(eps);

  int i = -1;
  Vec vv;
  MatCreateVecs(AA, NULL, &vv);
  int *idx = (int*)calloc(A->nr, sizeof(int));
  for(i = 0; i < A->nr; ++i) {
    idx[i] = i;
  }

  int nconv;
  EPSGetConverged(eps, &nconv);
  assert(nconv >= nev);
  for (i = 0; i < nev; i++) {
    PetscScalar kr, ki;
    EPSGetEigenpair(eps, i, &kr, &ki, vv, NULL);
    assert(fabs(ki) < MYAMGEPS);
    VecGetValues(vv, A->nr, idx, evec[i]);
    eval[i] = kr;

#if 0
    PetscReal   errorr, errora;
    EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &errorr);
    EPSComputeError(eps, i, EPS_ERROR_ABSOLUTE, &errora);

    if (ki != 0.0) {
      PetscPrintf(PETSC_COMM_WORLD, " %02d %.16e%+.16ei %.6e %.6e\n",                i, (double)kr, (double)ki, (double)errorr, (double)errora);
    } else {
      PetscPrintf(PETSC_COMM_WORLD, " %02d      % .16e        %.6e          %.6e\n", i, (double)kr, (double)errorr, (double)errora);
    }
#endif
  }
  

  free(idx); idx = NULL;
  VecDestroy(&vv);
  MatDestroy(&AA);
  MatDestroy(&MM);
  EPSDestroy(&eps);
  //SlepcFinalize();
}

static void Transform_dmatcsr_petscmat(dmatcsr *A, Mat *AA)
{
  int i = -1;
  MatCreateSeqAIJ(PETSC_COMM_WORLD, A->nr, A->nc, PETSC_DEFAULT, NULL, AA);
  MatSetUp(*AA);
  for(i = 0; i < A->nr; ++i) {
    MatSetValues(*AA, 1, &i, A->ia[i+1] - A->ia[i], A->ja + A->ia[i], A->va + A->ia[i], INSERT_VALUES);
  }
  MatAssemblyBegin(*AA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*AA,MAT_FINAL_ASSEMBLY);
}

#endif
