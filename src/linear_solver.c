#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "preprocess.h"
#include "matrix.h"
#include "linear_algebra.h"
#include "linear_solver.h"
#include "multigrid.h"
#include "tool.h"
#include "amg_param.h"
#include "umfpack.h"

void Get_residual(dmatcsr *A, double *b, double *x, 
	          double *r, double *rnorm)
{
    double *resi = (NULL==r) ? (double*)calloc(A->nr, sizeof(double)) : r;
    //memset(resi, 0, A->nr*sizeof(double));
    Multi_dmatcsr_dvec(A, x, resi);
    Sumself_dvec_axpby(b, 1.0, resi, -1.0, A->nr);
    if(NULL != rnorm) *rnorm = Get_dvec_2norm(resi, A->nr);
    if(NULL == r)     free(resi);
}

/* A = D - L - U */
void Split_dmatcsr(dmatcsr *A, dmatcsr *D, dmatcsr *L, dmatcsr *U)
{
    //double tb = Get_time();
    assert(A->nr == A->nc);
    
    int       nr = A->nr;
    int    *A_ia = A->ia;
    int    *A_ja = A->ja;
    double *A_va = A->va;
    
    int    *D_ia = D->ia;
    int    *D_ja = D->ja;
    double *D_va = D->va;
    
    int    *L_ia = L->ia;
    int    *L_ja = L->ja;
    double *L_va = L->va;
    
    int    *U_ia = U->ia;
    int    *U_ja = U->ja;
    double *U_va = U->va;
    
    int Dnn = 0;
    int Lnn = 0;
    int Unn = 0;
    
    int i, j, k, rb, re, ip1;
    for(i=0; i<nr; i++)
    {
        ip1 = i + 1;
        rb = A_ia[i];
        re = A_ia[ip1];
        
        D_ia[ip1] = D_ia[i];
        L_ia[ip1] = L_ia[i];
        U_ia[ip1] = U_ia[i];
        
        for(j=rb; j<re; j++)
        {
            k = A_ja[j];
            if(k < i)
            {
                L_ia[ip1]++;
                L_ja[Lnn] = k;
                L_va[Lnn] = -A_va[j];
                Lnn++;
            }
            else if(k > i)
            {
                U_ia[ip1]++;
                U_ja[Unn] = k;
                U_va[Unn] = -A_va[j];
                Unn++;
            }
            else
            {
                D_ia[ip1]++;
                D_ja[Dnn] = k;
                D_va[Dnn] = A_va[j];
                Dnn++;
            }
        }
    }
    
    assert(Dnn = D->nn);
    assert(Dnn + Lnn + Unn == A->nn);
    
    int    *ja_tmp;
    double *va_tmp;
    
    ja_tmp = realloc(L_ja, Lnn*sizeof(int));
    va_tmp = realloc(L_va, Lnn*sizeof(double));
    if((NULL != ja_tmp) && (NULL != va_tmp))
    {
        L->nn = Lnn;
        L->ja = ja_tmp;
        L->va = va_tmp;
    }
    else
    {
        printf("Error when calling \"Split_dmatcsr\": realloc L\n");
    }
    
    ja_tmp = realloc(U_ja, Unn*sizeof(int));
    va_tmp = realloc(U_va, Unn*sizeof(double));
    if((NULL != ja_tmp) && (NULL != va_tmp))
    {
        U->nn = Unn;
        U->ja = ja_tmp;
        U->va = va_tmp;
    }
    else
    {
        printf("Error when calling \"Split_dmatcsr\": realloc U\n");
    }
    
    //double te = Get_time();
    //printf("split matrix %7d time: %f\n", A->nr, te-tb);
}

void Get_diag_dmatcsr(dmatcsr *A, double *d)
{
    assert(A->nr == A->nc);
    
    int       nr = A->nr;
    int    *A_ia = A->ia;
    int    *A_ja = A->ja;
    double *A_va = A->va;
    
    int i, j, k, rb, re;
    for(i=0; i<nr; i++)
    {
        rb = A_ia[i];
        re = A_ia[i+1];
        for(j=rb; j<re; j++)
        {
            k = A_ja[j];
            if(k == i) d[i] = A_va[j];
        }
    }
}

void Zoom_dmatcsr(dmatcsr *A, double *d)
{
    int       nr = A->nr;
    int    *A_ia = A->ia;
    int    *A_ja = A->ja;
    double *A_va = A->va;
    
    int i, j, k, rb, re;
    for(i=0; i<nr; i++)
    {
        rb = A_ia[i];
        re = A_ia[i+1];
        for(j=rb; j<re; j++)
        {
            k = A_ja[j];
	    A_va[j] /= sqrt(d[i] * d[k]);
        }
    }
}

/*
 * A(x_k + \alpha p_k) != A x_k + \alpha A p_k
 * and this will cause the output "resi_norm" 
 * is not exactly the same with 
 * the result obtained by calling Get_residual 
 * outside this function.
*/
void Linear_solver_cg(dmatcsr *A, double *b,  double *x, 
	              double tol, int max_iter, 
		      double *resi, double *resi_norm, int *niter, 
		      int print_level)
{
    int size = A->nc;
    double rnorm;
    double *r = (NULL==resi) ? (double*)malloc(size*sizeof(double)) : resi;
    Get_residual(A, b, x, r, &rnorm);
    int nit = 0;
    if(rnorm < tol)
    {
        if(NULL == resi)       free(r);
        if(NULL != resi_norm) *resi_norm = tol;
        if(NULL != niter)     *niter     = 0;
        return;
    }
    double *p  = (double*)malloc(size*sizeof(double));  
    double *Ap = (double*)malloc(size*sizeof(double));
    double alpha, beta;
    double tmp = rnorm;
    while((rnorm > tol) && (nit < max_iter))
    {
	nit++;
	if(nit == 1)
	{
            memcpy(p, r, size*sizeof(double));
	}
	else
	{
	    beta = rnorm*rnorm / (tmp*tmp);
	    Sumself_dvec_axpby(r, 1.0, p, beta, size);
	}
        Multi_dmatcsr_dvec(A, p, Ap);
        alpha = rnorm*rnorm / Multi_dvec_dvec(p, Ap, size);
        Sumself_dvec_axpby(p,   alpha, x, 1.0, size);
        Sumself_dvec_axpby(Ap, -alpha, r, 1.0, size);
	tmp = rnorm;
        rnorm = Get_dvec_2norm(r, size);
        if(print_level > 0) printf("cg iter = %3d, rnorm = %18.15f, ratio = %f\n", nit, rnorm, rnorm/tmp);
    }
    free(p);
    free(Ap);
    if(NULL == resi)       free(r);
    if(NULL != resi_norm) *resi_norm = rnorm;
    if(NULL != niter)     *niter     = nit;
    if(print_level > 0)
    {
        printf("---------------------------------------------\n");
        printf("cg iter = %3d, rnorm = %18.15f, tol = %18.15f\n", nit, rnorm, tol);
    }
}

void Linear_solver_gs(dmatcsr *A, double *d, 
                      double *b,  double *x, 
	              double tol, int max_iter, 
	              double *resi, double *resi_norm, int *niter,
	              int print_level)
{
    int size   = A->nc;
    int *ia    = A->ia;
    int *ja    = A->ja;
    double *va = A->va;

    double rnorm, rnorm_prev;
    double *r = (NULL==resi) ? (double*)malloc(size*sizeof(double)) : resi;
    Get_residual(A, b, x, r, &rnorm);
    
    int i, j, k, rb, re;
    double val = 0.0;
    
    rnorm_prev = rnorm;
    int nit = 0;
    //while((rnorm > tol) && (nit < max_iter))
    while(nit < max_iter)
    {
	nit++;
	
	for(i=0; i<size; i++)
	{
	    rb = ia[i];
	    re = ia[i+1];
	    
	    val = 0.0;
	    for(j=rb; j<re; j++)
	    {
	        k = ja[j];
		val += va[j] * x[k];
	    }
	    
	    //assert(fabs(d[i]) > EPS);
	    x[i] += (b[i] - val) / d[i];
	}
	
        if(print_level > 0) 
        {
            Get_residual(A, b, x, r, &rnorm);
            printf("gs iter = %3d, rnorm = %18.15f, ratio = %f\n", nit, rnorm, rnorm/rnorm_prev);
            rnorm_prev = rnorm;
        }
    }
    
    if(NULL != resi || NULL != resi_norm) Get_residual(A, b, x, r, &rnorm);
    if(NULL != resi_norm) *resi_norm = rnorm;
    if(NULL != niter)     *niter     = nit;
    if(NULL == resi)       free(r);
    if(print_level > 0)
    {
        printf("---------------------------------------------\n");
        printf("gs iter = %3d, rnorm = %18.15f, tol = %18.15f\n", nit, rnorm, tol);
    }
}


double Linear_solver_amgcycle(multigrid *amg, int current_level,
                              double *b, double *x, 
                              amg_param param)
{
    double tol            = param.amgcycle_tol;
    int    mu             = param.amgcycle_mu;
    int    nsmooth        = param.amgcycle_pre_post_smooth;
    int    coarsest_level = param.amgcycle_coarsest_level;
    if(coarsest_level <= 0) coarsest_level += amg->actual_level-1;
    
    int print_level = param.amgcycle_print_level;
    double rn = 99999.9;
    if(current_level == coarsest_level)
    {
        //printf("multigrid solver accurate solve on level %d\n", fine_level);
#if WITH_UMFPACK
	Linear_solver_direct(amg->A[coarsest_level], b, x);
	Get_residual(amg->A[coarsest_level], b, x, NULL, &rn);
#else
        double coarsest_tol     = param.amgcycle_coarsest_tol;
        int max_coarsest_smooth = param.amgcycle_max_coarsest_smooth;
        dmatcsr *A = amg->A[coarsest_level];
        if(param.linear_solver_base_type == CG)
        {
            Linear_solver_cg(A, b, x, coarsest_tol, max_coarsest_smooth, NULL, &rn, NULL, 0);
        }
        else if(param.linear_solver_base_type == GS)
        {
            //Remove_zero_dmatcsr(A);
            //dmatcsr *D = Create_dmatcsr(A->nr, A->nc, A->nr);
            //dmatcsr *L = Create_dmatcsr(A->nr, A->nc, A->nn);
            //dmatcsr *U = Create_dmatcsr(A->nr, A->nc, A->nn);
            //Split_dmatcsr(A, D, L, U);
	    double *d = (double*)calloc(A->nr, sizeof(double));
	    Get_diag_dmatcsr(A, d);
            //Linear_solver_gs(A, D, L, U, b, x, coarsest_tol, max_coarsest_smooth, NULL, &rn, NULL, 0);
            Linear_solver_gs(A, d, b, x, coarsest_tol, max_coarsest_smooth, NULL, &rn, NULL, 0);
            //Free_dmatcsr(U);
            //Free_dmatcsr(L);
            //Free_dmatcsr(D);
	    free(d);
        }
#endif
        return rn;
    }
    else
    {   
        if(print_level > 7) printf("........................................\n");
        double tb_pre = Get_time();
        //printf("multigrid solver pre-smoothing  on level %d\n", fine_level);
        int fine_level   = current_level;
        int coarse_level = current_level + 1;
        
        dmatcsr *A   = amg->A[fine_level];
        double *resi = calloc(A->nr, sizeof(double));
        double rn_pre;
	double *d = (double*)calloc(A->nr, sizeof(double));
        //dmatcsr *D, *L, *U;
        //D = L = U = NULL;
        //pre-smooth
        if(param.linear_solver_base_type == CG)
        {
            Linear_solver_cg(A, b, x, tol, nsmooth, resi, &rn_pre, NULL, 0);
        }
        else if(param.linear_solver_base_type == GS)
        {
            //Remove_zero_dmatcsr(A);
            //D = Create_dmatcsr(A->nr, A->nc, A->nr);
            //L = Create_dmatcsr(A->nr, A->nc, A->nn);
            //U = Create_dmatcsr(A->nr, A->nc, A->nn);
            //Split_dmatcsr(A, D, L, U);
	    Get_diag_dmatcsr(A, d);
            //Linear_solver_gs(A, D, L, U, b, x, tol, nsmooth, resi, &rn_pre, NULL, 0);
            Linear_solver_gs(A, d, b, x, tol, nsmooth, resi, &rn_pre, NULL, 0);
        }
        
        double *coarse_x    = calloc(amg->A[coarse_level]->nr, sizeof(double));
        double *coarse_resi = calloc(amg->A[coarse_level]->nr, sizeof(double));
        RestrictFine2Coarse(amg, current_level, coarse_level, resi, coarse_resi);

        double te_pre = Get_time();
        if(print_level > 7) printf("pre smoothing time = %f, rn = %18.15f\n", te_pre-tb_pre, rn_pre);
        
        double tb_cycle = Get_time();
        int i;
        for( i=0; i<mu; i++ )
            Linear_solver_amgcycle(amg, coarse_level, coarse_resi, coarse_x, param);
        double te_cycle = Get_time();
        if(print_level > 7) printf("recursive cycling time = %f, ", te_cycle-tb_cycle);
        
        double tb_post = Get_time();
        //printf("multigrid solver post-smoothing on level %d\n", fine_level);
        double *fine_x = calloc(A->nc, sizeof(double));
        ProlongCoarse2Fine(amg, coarse_level, current_level, coarse_x, fine_x);
        Sumself_dvec_axpby(fine_x, 1.0, x, 1.0, A->nc);
        if(print_level > 7)
        {
            Get_residual(A, b, x, NULL, &rn);
            printf("rn = %18.15f\n", rn);
        }
        //post-smooth
        if(param.linear_solver_base_type == CG)
        {
            Linear_solver_cg(A, b, x, tol, nsmooth, resi, &rn, NULL, 0);//post-smooth
        }
        else if(param.linear_solver_base_type == GS)
        {
            Linear_solver_gs(A, d, b, x, tol, nsmooth, resi, &rn, NULL, 0);
            //Linear_solver_gs(A, D, L, U, b, x, tol, nsmooth, resi, &rn, NULL, 0);
            //Free_dmatcsr(U);
            //Free_dmatcsr(L);
            //Free_dmatcsr(D);
        }
        
        double te_post = Get_time();
        if(print_level > 7) printf("post smoothing time = %f, rn = %18.15f\n", te_post-tb_post, rn);
        if(print_level > 7) printf("........................................\n");
        
	free(d);
        free(resi);
        free(coarse_resi);
        free(coarse_x);
        free(fine_x);
        
        return rn;
    }
}
	       
void Linear_solver_amg(multigrid *amg, int current_level,
	               double *b, double *x, 
	               amg_param param, 
	               double *resi_norm, int *ncycle)
{
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
    if(print_level > 0) printf("~~~~~~~~~~~~~~~~~~~~ linear solver amg ~~~~~~~~~~~~~~~~~~~~\n");
    if(print_level > 0) printf("niter            rnorm            ratio      time\n");
    //Get_residual(amg->A[current_level], b, x, NULL, &rn1);
    //if(print_level > 0) printf("%3d    %22.15f   %s   %s\n", 0, rn1, "  ----  ", "  ----  ");
    while(iter<max_cycle && rn1 > tol)
    {
	iter++;
	tb = Get_time();
	rn1 = Linear_solver_amgcycle(amg, current_level, b, x, param);
	te = Get_time();
	if(print_level > 0) printf("%3d    %22.15f   %f   %f\n", iter, rn1, rn1/rn0, te-tb);
	if(rn1 / rn0 > max_convergence_factor)
	{
	    nslow++;
	    if(nslow >= nmax_convergence_factor)
	    {
		if(print_level > 1) printf("warning: the convergence of amg linear solver may be slowed down.\n");
		break;
	    }
	}
	else
	{
	    nslow = 0;
	}
	rn0 = rn1;
    }
    if(print_level > 0) printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    if(NULL != resi_norm) *resi_norm = rn1;
    if(NULL != ncycle)    *ncycle    = iter;
}

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

/* UMFPACK_A   ---- A   x = b
 * UMFPACK_At  ---- A^H x = b * UMFPACK_Aat ---- A^T x = b */
void Linear_solver_direct(dmatcsr *A, double *b, double *x)
{
    void *Symbolic, *Numeric;
    int status;
    status = umfpack_di_symbolic(A->nr, A->nc, A->ia, A->ja, A->va, &Symbolic, NULL, NULL);
    if(status) { printf("Error in umfpack_di_symbolic: %d!\n", status); exit(-1); }
    status = umfpack_di_numeric(A->ia, A->ja, A->va, Symbolic, &Numeric, NULL, NULL);
    if(status) { printf("Error in umfpack_di_numeric: %d!\n", status); exit(-1); }
    umfpack_di_free_symbolic(&Symbolic);
    status = umfpack_di_solve(UMFPACK_Aat, A->ia, A->ja, A->va, x, b, Numeric, NULL, NULL);
    if(status) { printf("Error in umfpack_di_solve: %d!\n", status); exit(-1); }
    umfpack_di_free_numeric(&Numeric);
}
