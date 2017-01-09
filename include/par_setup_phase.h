#ifndef __PAR_SETUP_PHASE_H__
#define __PAR_SETUP_PHASE_H__

#include "par_multigrid.h"
#include "amg_param.h"

#if WITH_MPI
void Setup_par_phase(par_multigrid *par_amg, amg_param param);
int Coarsen_par(par_dmatcsr *A,  par_dmatcsr *M,  par_dmatcsr *P,  par_dmatcsr *R, 
		par_dmatcsr *AH, par_dmatcsr *MH, int *ncpt_proc,  amg_param param);
#endif

#endif
