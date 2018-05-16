#ifndef __AMG_PARAM_H__
#define __AMG_PARAM_H__

typedef struct AMG_PARAMATER_
{
    /*---------- amg setup phase ----------*/
    double strong_connection_threshold;
    double strong_diagonally_dominant_threshold;
    double truncation_threshold;
    int    positive_connected;
    int    coarsening_type;
    int    interpolation_type;
    int    max_level;
    int    max_coarsest_dof;
    int    setup_phase_print_level;

    /*---------- linear solver base type ----------*/
    int    linear_solver_base_type;
    
    /*---------- linear solver cg ----------*/
    int    cg_max_iter;
    double cg_tol;
    
    /*---------- linear solver gs ----------*/
    int    gs_max_iter;
    double gs_tol;

    /*---------- linear solver amg cycle ----------*/
    double amgcycle_tol;
    double amgcycle_coarsest_tol;
    int    amgcycle_mu;
    int    amgcycle_pre_post_smooth;
    int    amgcycle_coarsest_level;
    int    amgcycle_max_coarsest_smooth;
    int    amgcycle_print_level;

    /*---------- linear solver amg ----------*/
    double amgsolver_tol;
    int    amgsolver_max_cycle;
    double amgsolver_max_convergence_factor;
    int    amgsolver_nmax_convergence_factor;
    int    amgsolver_print_level;

    /*---------- linear solver pcg amg ----------*/
    double pcg_amg_tol;
    int    pcg_amg_max_iter;
    int    pcg_amg_print_level;

    /*---------- eigen solver amg ----------*/
    int amgeigen_coarsest_level;
    int amgeigen_nouter_iter;
    int amgeigen_print_level;

}amg_param;

void Init_amg_param(amg_param *param);

#endif
