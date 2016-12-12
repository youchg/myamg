#include "amg_param.h"
#include "preprocess.h"
void Init_amg_param(amg_param *param)
{
    /*---------- amg setup phase ----------*/
    param->strong_connection_threshold          = 0.25;
    param->strong_diagonally_dominant_threshold = 0.1;
    param->truncation_threshold                 = 0.2;
    param->positive_connected                   = NO;
    param->coarsening_type                      = STD_RS;
    param->interpolation_type                   = DIR;
    param->max_level                            = 20;
    param->max_coarsest_dof                     = 100;
    param->setup_phase_print_level              = 0;

    /*---------- linear solver base type ----------*/
    param->linear_solver_base_type = GS;
    
    /*---------- linear solver cg ----------*/
    param->cg_max_iter    = 1000;
    param->cg_tol         = 1e-16;
    
    /*---------- linear solver gs ----------*/
    param->gs_max_iter    = 1000;
    param->gs_tol         = 1e-16;

    /*---------- linear solver amg cycle ----------*/
    param->amgcycle_tol                 = 1e-16;
    param->amgcycle_coarsest_tol        = 1e-18;
    param->amgcycle_mu                  = 1;
    param->amgcycle_pre_post_smooth     = 4;
    param->amgcycle_coarsest_level      = 0;
    param->amgcycle_max_coarsest_smooth = 1000;
    param->amgcycle_print_level         = 0;

    /*---------- amgsolver ----------*/
    param->amgsolver_tol                     = 1e-16;
    param->amgsolver_max_cycle               = 100;
    param->amgsolver_max_convergence_factor  = 0.99;
    param->amgsolver_nmax_convergence_factor = 1;
    param->amgsolver_print_level             = 0;

    /*---------- amgpcg ----------*/
    param->pcg_amg_tol         = 1e-16;
    param->pcg_amg_max_iter    = 100;
    param->pcg_amg_print_level = 0;

    /*---------- eigen solver amg ----------*/
    param->amgeigen_coarsest_level = 0;
    param->amgeigen_nouter_iter    = 1;
    param->amgeigen_print_level    = 0;
}
