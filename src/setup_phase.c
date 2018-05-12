#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "preprocess.h"
#include "matrix.h"
#include "linear_algebra.h"
#include "tool.h"
#include "io.h"
#include "cfsplit.h"
#include "multigrid.h"
#include "amg_param.h"
#include "setup_phase.h"

#define DETAIL_TIME  0

void Setup_phase(multigrid *amg, amg_param param)
{
#if DETAIL_TIME > 0
    double t1, t2;
    double tc[amg->max_level];
#endif
    int max_level = amg->max_level;
    int i;
    if(NULL != amg->M)
    {
        for(i=0; i<max_level-1; i++)
        {
#if DETAIL_TIME > 0
            t1 = Get_time();
#endif
            if(!Coarsen(amg->A[i], amg->M[i], amg->P[i], amg->R[i], amg->A[i+1], amg->M[i+1], param))
		break;
#if DETAIL_TIME > 0
            t2 = Get_time(); tc[i] = t2 - t1;
#endif
            amg->actual_level++;
        }
    }
    else
    {
        for(i=0; i<max_level-1; i++)
        {
#if DETAIL_TIME > 0
            t1 = Get_time();
#endif
            if(!Coarsen(amg->A[i], NULL, amg->P[i], amg->R[i], amg->A[i+1], NULL, param))
		    break;
#if DETAIL_TIME > 0
            t2 = Get_time(); tc[i] = t2 - t1;
#endif
            amg->actual_level++;
        }
    }
#if DETAIL_TIME > 0
    for(i=0; i<amg->actual_level-1; i++)
        printf("coarsen --%d-- time: %f\n", i, tc[i]);
#endif
    //assert(amg->actual_level > 1);
}

int Coarsen(dmatcsr *A, dmatcsr *M, dmatcsr *P, dmatcsr *R, dmatcsr *AH, dmatcsr *MH, amg_param param)
{
#if DETAIL_TIME > 3
    double t1, t2, tb, te;
    double time_gen_S;
    double time_pre, time_post;
    double time_CLJP;
    double time_gen_spa_P, time_gen_P;
    double time_gen_R;
    double time_gen_AH;
    double time_total;
    tb = Get_time();
#endif

    int ncpt = 0;
    int     *dof = (int *)   calloc(A->nr, sizeof(int));
    imatcsr *S   = (imatcsr*)malloc(sizeof(imatcsr));
    Reset_dof(A);

#if DETAIL_TIME > 3    
    t1 = Get_time();
#endif

    Generate_strong_coupling_set(A, S, param);

#if DETAIL_TIME > 3
    t2 = Get_time(); time_gen_S = t2 - t1; 
    printf("gen  S       : %f\n", time_gen_S);
    t1 = Get_time();
#endif

    if(STD_RS == param.coarsening_type)
    {
	ncpt = Pre_split(A, S, dof);
	//ncpt = Pre_split_fasp(A, S, dof);
#if DETAIL_TIME > 3
	t2 = Get_time(); time_pre = t2 - t1; 
	printf("pre  spliting: %f\n", time_pre);
	t1 = Get_time();
#endif

	ncpt = Post_split(S, dof);
#if DETAIL_TIME > 3
	t2 = Get_time(); time_post = t2 - t1; 
	printf("post spliting: %f\n", time_post);
#endif
    }
    else if(CLJP == param.coarsening_type)
    {
	//printf("coarsening type: CLJP\n");
	ncpt = CLJP_split(A, S, dof);
#if DETAIL_TIME > 3
	t2 = Get_time(); time_CLJP = t2 - t1; 
	printf("CLJP spliting: %f\n", time_CLJP);
#endif
    }

#if DETAIL_TIME > 3
	t1 = Get_time();
#endif

    if(ncpt < param.max_coarsest_dof)
    {
	Free_imatcsr(S);
	free(dof);
	return FAIL;
    }
    
    if(param.interpolation_type == DIR)
    {
        Generate_sparsity_P_dir(S, dof, P);
    }
    else if(param.interpolation_type == STD)
    {
        Generate_sparsity_P_std(S, dof, P);
    }

#if DETAIL_TIME > 3    
    t2 = Get_time(); time_gen_spa_P = t2 - t1;
    printf("gen  spa P   : %f\n", time_gen_spa_P);
    t1 = Get_time();
#endif

    if(param.interpolation_type == DIR)
    {
        Generate_P_dir(A, S, dof, P);
    }
    else if(param.interpolation_type == STD)
    {
        Generate_P_std(A, S, dof, P);
    }
    Truncate_P(P, param);
    
#if DETAIL_TIME > 3
    t2 = Get_time(); time_gen_P = t2 - t1; 
    printf("gen  P       : %f\n", time_gen_P);
    t1 = Get_time();
#endif

    Transpose_dmatcsr(P, R);
    
#if DETAIL_TIME > 3
    t2 = Get_time(); time_gen_R = t2 - t1; 
    printf("gen  R       : %f\n", time_gen_R);
    t1 = Get_time();
#endif
    
    Multi_dmatcsr_dmatcsr_dmatcsr(R, A, P, AH);

#if DETAIL_TIME > 3
    t2 = Get_time(); time_gen_AH = t2 - t1;
    printf("gen  AH      : %f\n", time_gen_AH);
#endif    

#if 0
    dmatcsr *AH2 = malloc(sizeof(dmatcsr));
    dmatcsr *RA  = malloc(sizeof(dmatcsr));
    double tabcb2 = Get_time();
    Multi_dmatcsr_dmatcsr(R, A, RA);
    Multi_dmatcsr_dmatcsr(RA, P, AH2);
    double tabce2 = Get_time();
    printf("gen AH2     :%f\n", tabce2-tabcb2);
#endif

    if(NULL != M) Multi_dmatcsr_dmatcsr_dmatcsr(R, M, P, MH);

    Free_imatcsr(S);
    free(dof);
#if DETAIL_TIME > 3
    te = Get_time(); time_total = te - tb;
    printf("total        : %f\n", time_total);
    printf("****************************************\n\n");
#endif
    return SUCCESS;
}
