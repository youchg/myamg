#include <stdio.h>
#include <stdlib.h>
#include "fasp_interface.h"
#include "multigrid.h"
#include "matrix.h"
#include "linear_algebra.h"
#include "fasp.h"
#include "fasp_functs.h"

void Setup_phase_fasp(multigrid *amg)
{
    AMG_param amgparam;
    Init_fasp_amg_param(&amgparam);
    amgparam.max_levels = amg->max_level;

    AMG_data *mgl = fasp_amg_data_create(amg->max_level);
    Copy_dmatcsr_to_fasp(amg->A[0], &mgl[0].A);
    mgl[0].b = fasp_dvec_create(mgl[0].A.col);
    mgl[0].x = fasp_dvec_create(mgl[0].A.col);

    const SHORT prtlvl   = amgparam.print_level;
    const SHORT amg_type = amgparam.AMG_type;
    switch (amg_type)
    {
        case SA_AMG: // Smoothed Aggregation AMG setup
	    if(prtlvl > PRINT_NONE) printf("\nCalling SA AMG ...\n");
	    fasp_amg_setup_sa(mgl, &amgparam); break;

	case UA_AMG: // Unsmoothed Aggregation AMG setup
	    if(prtlvl > PRINT_NONE) printf("\nCalling UA AMG ...\n");
	    fasp_amg_setup_ua(mgl, &amgparam); break;

	default: // Classical AMG setup
	    if(prtlvl > PRINT_NONE) printf("\nCalling classical AMG ...\n");
	    fasp_amg_setup_rs(mgl, &amgparam);
    }

    int actual_level  = mgl[0].num_levels;
    amg->actual_level = mgl[0].num_levels;

    int i;
    for(i=0; i<actual_level; i++)
    {
        if(i < actual_level-1)
        {
            Copy_dmatcsr_from_fasp(&mgl[i].P, amg->P[i]);
            Copy_dmatcsr_from_fasp(&mgl[i].R, amg->R[i]);
        }
        if(i > 0)
        {
            Copy_dmatcsr_from_fasp(&mgl[i].A, amg->A[i]);
            if(amg->M != NULL) Multi_dmatcsr_dmatcsr_dmatcsr(amg->R[i-1], amg->M[i-1], amg->P[i-1], amg->M[i]);
        }
    }
    
    Free_fasp_data(mgl); 
}

void Init_fasp_amg_param(AMG_param *amgparam)
{
    fasp_param_amg_init(amgparam);
    
    amgparam->AMG_type             = CLASSIC_AMG;// CLASSIC_AMG 
                                                 // SA_AMG 
                                                 // UA_AMG   
    amgparam->max_levels           = 20;
    amgparam->coarse_dof           = 100;//max number of coarse dof
    amgparam->cycle_type           = V_CYCLE;

    amgparam->coarsening_type      = 1;// 1 Modified RS
                                       // 3 Compatible Relaxation
                                       // 4 Aggressive 
    amgparam->interpolation_type   = 1;// 1 Direct | 2 Standard | 3 Energy-min
    amgparam->strong_threshold     = 0.25;
    amgparam->truncation_threshold = 0.2;
    amgparam->max_row_sum          = 0.9;
    
    // Aggregation AMG specific
    amgparam->aggregation_type     = PAIRWISE;
    amgparam->pair_number          = 2;
    amgparam->strong_coupled       = 0.08;
    amgparam->max_aggregation      = 20;
    amgparam->tentative_smooth     = 0.67;
    amgparam->smooth_filter        = OFF;
    amgparam->quality_bound        = 8.0;

    amgparam->print_level          = 0;

    amgparam->maxit                = 100;
    amgparam->tol                  = 1e-8;
    
    // ILU smoother parameters
    amgparam->ILU_type             = ILUk;
    amgparam->ILU_levels           = 0;
    amgparam->ILU_lfil             = 0;
    amgparam->ILU_droptol          = 0.001;
    amgparam->ILU_relax            = 0;
    
    // Schwarz smoother parameters
    amgparam->Schwarz_levels       = 0;
    amgparam->Schwarz_mmsize       = 200;
    amgparam->Schwarz_maxlvl       = 2;
    amgparam->Schwarz_type         = 1;
}

void Free_fasp_data(AMG_data *mgl)
{    
    const INT max_levels = mgl[0].max_levels;
    INT i;

    for (i=0; i<max_levels; ++i) {
        //if (&mgl[i].A) { fasp_dcsr_free(&mgl[i].A); }
        //if (&mgl[i].P) { fasp_dcsr_free(&mgl[i].P); }
        //if (&mgl[i].R) { fasp_dcsr_free(&mgl[i].R); }
        if (&mgl[i].b) { fasp_dvec_free(&mgl[i].b); }
        if (&mgl[i].x) { fasp_dvec_free(&mgl[i].x); }
        if (&mgl[i].w) { fasp_dvec_free(&mgl[i].w); }
        if (&mgl[i].cfmark) { fasp_ivec_free(&mgl[i].cfmark); }
        if (&mgl[i].LU) { fasp_ilu_data_free(&mgl[i].LU); }
        
        if (&mgl[i].Schwarz) {fasp_Schwarz_data_free (&mgl[i].Schwarz);}
        //if (&mgl[i].schwarz) {fasp_schwarz_data_free (&mgl[i].schwarz);}
    }
    
    for (i=0; i<mgl->near_kernel_dim; ++i) {
        fasp_mem_free(mgl->near_kernel_basis[i]); 
        mgl->near_kernel_basis[i]=NULL;
    }
    
    fasp_mem_free(mgl->near_kernel_basis); mgl->near_kernel_basis = NULL;
    fasp_mem_free(mgl); mgl = NULL;
}

void Copy_dmatcsr_to_fasp(dmatcsr *A, dCSRmat *B)
{
    B->row = A->nr;
    B->col = A->nc;
    B->nnz = A->nn;
    
    B->IA  = A->ia;
    B->JA  = A->ja;
    B->val = A->va;
}
void Copy_dmatcsr_from_fasp(dCSRmat *B, dmatcsr *A)
{
    A->nr = B->row;
    A->nc = B->col;
    A->nn = B->nnz;
    
    A->ia = B->IA;
    A->ja = B->JA;
    A->va = B->val;
}

void Write_dmatcsr_fasp(dmatcsr *A, const char *filename)
{
    FILE *file = fopen(filename,"w");
    if(!file)
    {
        printf("\nError: Cannot open %s!\n", filename);
	exit(0);
    }
    
    int i;
    
    fprintf(file, "%d\n", A->nr);
    fprintf(file, "\n");
    
    int nr = A->nr+1; // nr plus 1
    int nn = A->nn;
    
    for(i=0; i<nr; i++) fprintf(file, "%d\n",      A->ia[i]+1);
    fprintf(file, "\n");
    for(i=0; i<nn; i++) fprintf(file, "%d\n",      A->ja[i]+1);
    fprintf(file, "\n");
    for(i=0; i<nn; i++) fprintf(file, "%15.12f\n", A->va[i]);
    
    fclose(file);
}
