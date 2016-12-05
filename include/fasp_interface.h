#ifndef __FASP_INTERFACE_H__
#define __FASP_INTERFACE_H__

#include "matrix.h"
#include "multigrid.h"
#include "fasp.h"
#include "fasp_functs.h"

void Setup_phase_fasp(multigrid *amg);
void Init_fasp_amg_param(AMG_param *amgparam);
void Free_fasp_data(AMG_data *mgl);
void Copy_dmatcsr_to_fasp  (dmatcsr *A, dCSRmat *B);
void Copy_dmatcsr_from_fasp(dCSRmat *B, dmatcsr *A);

void Write_dmatcsr_fasp(dmatcsr *A, const char *filename);

#endif
