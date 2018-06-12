#include "eigen_solver_init_vector.h"

#include <stdlib.h>
#include <string.h>
#include "arpack_interface.h"
#include "slepc_interface.h"

void Get_eigen_solver_init_vector(multigrid *amg, int nev, double *eval, double **evec, amg_param param)
{
  int nlevel         = amg->actual_level;
  int finest_level   = 0;
  int coarsest_level = param.amgeigen_coarsest_level;
  if(coarsest_level <= 0) coarsest_level += nlevel-1;

  dmatcsr *AH     = amg->A[coarsest_level];
  dmatcsr *MH     = amg->M[coarsest_level]; 

  int i;
  double **evec_coarsest = (double**)malloc(nev * sizeof(double*));
  for(i=0; i<nev; i++) evec_coarsest[i] = (double*)calloc(AH->nr, sizeof(double));

  Eigen_solver_arpack_dn(AH, MH, nev, eval, evec_coarsest);

  for(i=0; i<nev; i++) {
    memset(evec[i], 0, amg->A[finest_level]->nr*sizeof(double));
    ProlongCoarse2Fine(amg, coarsest_level, finest_level, evec_coarsest[i], evec[i]);
  }

  for(i=0; i<nev; i++) { free(evec_coarsest[i]); evec_coarsest[i] = NULL; }
}

