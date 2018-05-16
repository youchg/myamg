#ifndef __SETUP_PHASE_H__
#define __SETUP_PHASE_H__
#include "multigrid.h"
#include "amg_param.h"

void Setup_phase(multigrid *amg, amg_param param);
int Coarsen(dmatcsr *A, dmatcsr *M, dmatcsr *P, dmatcsr *R, dmatcsr *AH, dmatcsr *MH, amg_param param);

#endif
