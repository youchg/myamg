#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "par_linear_solver.h"

#if WITH_MPI

#include "preprocess.h"
#include "tool.h"
#include "amg_param.h"
#include "matrix.h"
#include "par_matrix_vector.h"
#include "par_linear_algebra.h"
#include "par_multigrid.h"

#if 1
double Linear_solver_par_amgcycle(par_multigrid *pamg, int current_level,
                                  par_dvec *b,         par_dvec *x, 
                                  amg_param param)
{
    double tol            = param.amgcycle_tol;
    int    mu             = param.amgcycle_mu;
    int    nsmooth        = param.amgcycle_pre_post_smooth;
    int    coarsest_level = param.amgcycle_coarsest_level;
    if(coarsest_level <= 0) coarsest_level += pamg->actual_level-1;
    
    int print_level = param.amgcycle_print_level;
    double rn = 99999.9;
    if(current_level == coarsest_level)
    {
	//if(MPI_COMM_NULL != pamg->comm[coarsest_level])
        //printf("multigrid solver accurate solve on level %d\n", fine_level);
        double coarsest_tol     = param.amgcycle_coarsest_tol;
        int max_coarsest_smooth = param.amgcycle_max_coarsest_smooth;
        par_dmatcsr *A = pamg->A[coarsest_level];
        if(1 == A->nr_global)
	{
	    //printf("b[0] = %15.12f\n", b->value[0]);
	    //printf("A[0] = %15.12f\n", A->diag->va[0]);
	    //printf("\n\n");
	    x->value[0] = b->value[0] / A->diag->va[0];
	    rn = 0.0;
	}
	else if(param.linear_solver_base_type == CG)
            Linear_solver_par_cg(A, b, x, coarsest_tol, max_coarsest_smooth, NULL, &rn, NULL, 0);
	else
	{
	    printf("ERROR: Unknown linear solver base type!\n");
	}
        return rn;
    }
    else
    {   
        //printf("multigrid solver pre-smoothing  on level %d\n", fine_level);
        int fine_level   = current_level;
        int coarse_level = current_level + 1;
        
        par_dmatcsr *A    = pamg->A[fine_level];
        par_dvec    *resi = Init_par_dvec_mv(A);
        double rn_pre;
        //pre-smooth
        if(param.linear_solver_base_type == CG)
        {
            Linear_solver_par_cg(A, b, x, tol, nsmooth, resi, &rn_pre, NULL, 0);
        }
	else
	{
	    printf("ERROR: Unknown linear solver base type!\n");
	}
        
        par_dvec *coarse_x    = Init_par_dvec_mv(pamg->A[coarse_level]);
        par_dvec *coarse_resi = Init_par_dvec_mv(pamg->A[coarse_level]);
        Restrict_par_f2c(pamg, current_level, coarse_level, resi, coarse_resi);

        //double te_pre = Get_time();
        //if(print_level > 7) printf("pre smoothing time = %f, rn = %18.15f\n", te_pre-tb_pre, rn_pre);
        
        //double tb_cycle = Get_time();
        int i;
        for( i=0; i<mu; i++ )
	    if(MPI_COMM_NULL != pamg->comm[coarse_level])
		Linear_solver_par_amgcycle(pamg, coarse_level, coarse_resi, coarse_x, param);
        //double te_cycle = Get_time();
        //if(print_level > 7) printf("recursive cycling time = %f, ", te_cycle-tb_cycle);
        
        //double tb_post = Get_time();
        //printf("multigrid solver post-smoothing on level %d\n", fine_level);
        par_dvec *fine_x = Init_par_dvec_mv(A);
        Prolong_par_c2f(pamg, coarse_level, current_level, coarse_x, fine_x);
        Sumself_par_dvec_axpby(fine_x, 1.0, x, 1.0);
        if(print_level > 7)
        {
            Get_par_residual(A, b, x, NULL, &rn);
            printf("rn = %18.15f\n", rn);
        }
        //post-smooth
        if(param.linear_solver_base_type == CG)
        {
            Linear_solver_par_cg(A, b, x, tol, nsmooth, resi, &rn, NULL, 0);//post-smooth
        }
        
        //double te_post = Get_time();
        //if(print_level > 7) printf("post smoothing time = %f, rn = %18.15f\n", te_post-tb_post, rn);
        //if(print_level > 7) printf("........................................\n");
        
        Free_par_dvec(resi);
        Free_par_dvec(coarse_resi);
        Free_par_dvec(coarse_x);
        Free_par_dvec(fine_x);
        
        return rn;
    }
}
#endif


void Get_par_residual(par_dmatcsr *A, par_dvec *b, par_dvec *x, par_dvec *r, double *rnorm)
{
    par_dvec *resi = (NULL==r) ? Init_par_dvec_mv(A) : r;
    Multi_par_dmatcsr_dvec(A, x, resi);
    Sumself_par_dvec_axpby(b, 1.0, resi, -1.0);
    if(NULL != rnorm) *rnorm = Get_par_dvec_2norm(resi);
    if(NULL == r)     Free_par_dvec(resi);
}


/*
 * A(x_k + \alpha p_k) != A x_k + \alpha A p_k
 * and this will cause the output "resi_norm" 
 * is not exactly the same with 
 * the result obtained by calling Get_residual 
 * outside this function.
*/
void Linear_solver_par_cg(par_dmatcsr *A, par_dvec *b, par_dvec *x, double tol, int max_iter, 
		          par_dvec *resi, double *resi_norm, int *niter, int print_level)
{
    int  myrank, nproc_total;
    MPI_Comm comm = A->comm;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nproc_total);

    double rnorm;
    par_dvec *r = (NULL==resi) ? Init_par_dvec_mv(A) : resi;
    Get_par_residual(A, b, x, r, &rnorm);
    int nit = 0;
    if(rnorm < tol)
    {
        if(NULL == resi)       Free_par_dvec(r);
        if(NULL != resi_norm) *resi_norm = tol;
        if(NULL != niter)     *niter     = 0;
        return;
    }
    par_dvec *p  = Init_par_dvec_mv(A);
    par_dvec *Ap = Init_par_dvec_mv(A);
    double alpha, beta;
    double tmp = rnorm;
    while((rnorm > tol) && (nit < max_iter))
    {
	nit++;
	if(nit == 1)
	{
	    Copy_par_dvec(r, p);
	}
	else
	{
	    beta = rnorm*rnorm / (tmp*tmp);
	    Sumself_par_dvec_axpby(r, 1.0, p, beta);
	}
        Multi_par_dmatcsr_dvec(A, p, Ap);
        alpha = rnorm*rnorm / Multi_par_dvec_dvec(p, Ap);
        Sumself_par_dvec_axpby(p,   alpha, x, 1.0);
        Sumself_par_dvec_axpby(Ap, -alpha, r, 1.0);
	tmp = rnorm;
        rnorm = Get_par_dvec_2norm(r);
        if((print_level>0) && (myrank==0)) printf("cg iter = %3d, rnorm = %18.15f, ratio = %f\n", nit, rnorm, rnorm/tmp);
    }
    Free_par_dvec(p);
    Free_par_dvec(Ap);
    if(NULL == resi)       Free_par_dvec(r);
    if(NULL != resi_norm) *resi_norm = rnorm;
    if(NULL != niter)     *niter     = nit;
    if((print_level>0) && (myrank==0))
    {
        printf("---------------------------------------------\n");
        printf("cg iter = %3d, rnorm = %18.15f, tol = %18.15f\n", nit, rnorm, tol);
    }
}


	       
void Linear_solver_par_amg(par_multigrid *pamg, int current_level,
			   par_dvec *b, par_dvec *x, 
			   amg_param param, 
			   double *resi_norm, int *ncycle)
{
    int  myrank, nproc_total;
    MPI_Comm comm = pamg->comm[0];
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nproc_total);

    double tol                     = param.amgsolver_tol;
    double max_convergence_factor  = param.amgsolver_max_convergence_factor;
    int    max_cycle               = param.amgsolver_max_cycle;
    int    nmax_convergence_factor = param.amgsolver_nmax_convergence_factor;

    int print_level = param.amgsolver_print_level;
    
    double tb, te;
    double rn1 = 9999999.0;//residual norm
    double rn0 = 9999999.0;
    int    nslow = 0;
    int    iter  = 0;
    if((print_level>0) && (myrank==0)) printf("~~~~~~~~~~~~~~~~~~~~ linear solver amg ~~~~~~~~~~~~~~~~~~~~\n");
    //if((print_level>0) && (myrank==0)) printf("niter            rnorm            ratio\n");
    if((print_level>0) && (myrank==0)) printf("niter            rnorm            ratio      time\n");
    //Get_par_residual(pamg->A[current_level], b, x, NULL, &rn1);
    //if((print_level>0) && (myrank==0)) printf("%3d    %22.15f   %s   %s\n", 0, rn1, "  ----  ", "  ----  ");
    while(iter<max_cycle && rn1 > tol)
    {
	iter++;
	tb = MPI_Wtime();
	rn1 = Linear_solver_par_amgcycle(pamg, current_level, b, x, param);
	te = MPI_Wtime();
	//if((print_level>0) && (myrank==0)) printf("%3d    %22.15f   %f\n", iter, rn1, rn1/rn0);
	if((print_level>0) && (myrank==0)) printf("%3d    %22.15f   %f   %f\n", iter, rn1, rn1/rn0, te-tb);
	if(rn1 / rn0 > max_convergence_factor)
	{
	    nslow++;
	    if(nslow >= nmax_convergence_factor)
	    {
		if((print_level>0) && (myrank==0)) 
		    printf("warning: the convergence of amg linear solver may be slowed down.\n");
		break;
	    }
	}
	else
	{
	    nslow = 0;
	}
	rn0 = rn1;
    }
    if((print_level>0) && (myrank==0)) printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    if(NULL != resi_norm) *resi_norm = rn1;
    if(NULL != ncycle)    *ncycle    = iter;
}

#if 0
void Linear_solver_pcg_amg(multigrid *amg, int current_level, 
                           double *b,  double *x, 
                           amg_param param, 
                           double *resi, double *resi_norm, int *nit,
                           int *ncycle_amg_total)
{
    amg_param param_local = param;
    double tol            = param_local.pcg_amg_tol;
    double max_iter       = param_local.pcg_amg_max_iter;
    int    print_level    = param_local.pcg_amg_print_level;
    if(print_level <= 1)    param_local.amgsolver_print_level = 0;
    
    dmatcsr *A = amg->A[current_level];
    int size   = A->nc;
    int niter  = 0;
    double rnorm;
    //double bnorm;
    double t1, t2;
    
    double *r  = (NULL==resi)?(double*)calloc(size, sizeof(double)):resi;
    t1    = Get_time();
    Get_residual(A, b, x, r, NULL);
    rnorm  = Get_dvec_2norm(r, size);
    //bnorm  = Get_dvec_2norm(b, size);
    t2     = Get_time();
    if(print_level > 0) printf("************************* linear solver pcg amg *************************\n"); 
    if(print_level > 0) printf("niter    amgcycle               rnorm             time          ratio\n");
    if(print_level > 0) printf("-------------------------------------------------------------------------\n");
    if(print_level > 0) printf(" %2d         %3d      %22.15f      %f        ----\n", 0, 0, rnorm, t2-t1);
    if(print_level > 0) printf("-------------------------------------------------------------------------\n");
    //if(rnorm <= tol*bnorm)
    if(rnorm <= tol)
    {
        if(NULL == resi)       free(r);
        if(NULL != resi_norm) *resi_norm = rnorm;
	if(NULL != ncycle_amg_total) *ncycle_amg_total = 0;
        if(NULL != nit)       *nit       = 0;
        return;
    }
    
    double *p    = (double*)calloc(size, sizeof(double));
    double *Ap   = (double*)calloc(size, sizeof(double));
    double *z    = (double*)calloc(size, sizeof(double));
    double alpha = 0.0;
    double beta  = 0.0;
    double tmp1  = 1.0;
    double tmp2  = 1.0;
    
    int ncycle_amgsolver = 0;
    int ncycle_amgsolver_total = 0;
    
    double rnorm2 = rnorm;
    //while((rnorm > tol*bnorm) && (niter < max_iter))
    while((rnorm > tol) && (niter < max_iter))
    {
	t1 = Get_time();
	memset(z, 0, size*sizeof(double));
	Linear_solver_amg(amg, current_level, r, z, param_local, NULL, &ncycle_amgsolver);
	ncycle_amgsolver_total += ncycle_amgsolver;
        niter++;
        if(niter == 1)
        {
            memcpy(p, z, size*sizeof(double));
            tmp1 = Multi_dvec_dvec(r, z, size);
        }
        else
        {
            tmp2 = tmp1;
            tmp1 = Multi_dvec_dvec(r, z, size);
            beta = tmp1/tmp2;
            Sumself_dvec_axpby(z, 1.0, p, beta, size);
        }
        Multi_dmatcsr_dvec(A, p, Ap);
        tmp2 = Multi_dvec_dvec(p, Ap, size);
        alpha = tmp1/tmp2;
        Sumself_dvec_axpby(p,   alpha, x, 1.0, size);
        Sumself_dvec_axpby(Ap, -alpha, r, 1.0, size);
        rnorm = Get_dvec_2norm(r, size);
	t2 = Get_time();
        //if(print_level > 0) printf("pcg: niter = %2d, amgcycle = %3d, rnorm = %18.15f, r/b = %g\n", niter, ncycle_amgsolver, rnorm, rnorm/bnorm);
        if(print_level > 0) printf(" %2d         %3d      %22.15f      %f     %f\n", 
		                    niter, ncycle_amgsolver, rnorm, t2-t1, rnorm/rnorm2);
	rnorm2 = rnorm;
    }
    if(print_level > 0) printf("-------------------------------------------------------------------------\n");
    if(print_level > 0) printf("*************************************************************************\n");
    
    if(NULL != resi_norm)        *resi_norm        = rnorm;
    if(NULL != nit)              *nit              = niter;
    if(NULL != ncycle_amg_total) *ncycle_amg_total = ncycle_amgsolver_total;
    if(NULL == resi)              free(r);
    free(p);
    free(Ap);
    free(z);
}

#endif
#endif
