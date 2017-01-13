#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "preprocess.h"
#include "matrix.h"
#include "linear_algebra.h"
#include "linear_solver.h"
#include "multigrid.h"
#include "amg_param.h"
#include "tool.h"

//the subrutine for the non-symmetric matrix of the double precision	  		  
extern void dnaupd_(int *ido, const char *bmat, int *n, const char *which, 
                    int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, 
                    int *iparam, int *ipntr, 
                    double *workd, double *workl, int *lworkl, int *info);
	      
extern void dneupd_(int *rvec, char *hownmy, int *select, double *dr, double *di, 
	            double *z, int *ldz, double *sigmar, double *sigmai, double *workev, 
	            const char *bmat, int *n, const char *which, int *nev, double *tol, 
	            double *resid, int *ncv, double *v, int *ldv, int *iparam, 
	            int *ipntr, double *workd, double *workl, int *lworkl, int *ierr);

//dndrv4
//对复特征值，nev 可能会增加 1
#define CG      1

#ifdef WITH_UMFPACK
#define UMFPACK 2
#define ARPACK_LINEAR_SOLVER 2
#else
#define ARPACK_LINEAR_SOLVER 1
#endif

void Eigen_solver_arpack_dn(dmatcsr *A, dmatcsr *M, int nev, double *eval, double **evec)
{
    int    nev_local = nev;
    /* parameters for calling dnaupd_ */
    int     ido      = 0;
    char    bmat     = (NULL==M)? 'I': 'G';
    int     n        = A->nr;//matrix dimension
    char   *which    = "LM";//LM or SM
    double  tol      = 0.0;
    double *resid    = calloc(n, sizeof(double)); 
    int     ncv      = 4*nev_local;//Arnoldi vectors number
    if(ncv > n) ncv  = n;
    int     ldv      = n;
    double *v        = calloc(ldv*ncv, sizeof(double));
    int    *iparam   = calloc(11, sizeof(int)); 
    iparam[0]        = 1;//selecting the implicit shifts
    iparam[2]        = n;//maximum number of Arnoldi iterations
    iparam[6]        = (NULL==M)? 1: 3;//eigenproblem type
    int    *ipntr    = calloc(14, sizeof(int));//output, pointer
    double *workd    = calloc(3*n, sizeof(double)); 
    int     lworkl   = 3*ncv*ncv + 6*ncv;
    double *workl    = calloc(lworkl, sizeof(double)); 
    int     info     = 0;//random initial vector
    
    /* tmp variables */
    double *rhs = calloc(n, sizeof(double));
    int i, j;
   
    /* parameters for calling dneupd_ */
    int     rvec   = 1;//compute the Ritz or Schur vectors
    char    howmny = 'A';//compute nev Ritz vectors
    int    *select = calloc(ncv, sizeof(int));//internal workspace
    double *dr     = calloc(ncv, sizeof(double));//real      part of the Ritz value
    double *di     = calloc(ncv, sizeof(double));//imaginary part of the Ritz value
    double  sigmar = 0.0;//real      part of the shift
    double  sigmai = 0.0;//imaginary part of the shift
    double *workev = (double*)calloc(3*ncv, sizeof(double)); 
    int     ierr;//info
    
#if ARPACK_LINEAR_SOLVER == CG
    /*粗略测试，假设cg iteration足够多，tol = 10^-n
      特征值的有效位数为 n-3 */
    double cg_tol = 1e-15;
    int    cg_ite = 1000;
    double rnorm  = 9999.9;
    int    ncgit  = -1;

    //printf("linear solver type: CG\n");
#else
    //printf("linear solver type: direct\n");
#endif

    //=================================================================
    // start reverse communivation process
    //=================================================================
    do
    {
        dnaupd_(&ido, &bmat, &n, which, &nev_local, &tol, resid, 
	        &ncv, v, &ldv, iparam, ipntr, workd, workl,
	        &lworkl, &info);
	/*
        printf("info = %d\n", info);
        printf("Arnoldi size = %d\n",iparam[4]);
        printf("iparam[6] = %d\n", iparam[6]);
	*/

        if(iparam[6] == 3)
        {
            //printf("ido = %d\n", ido);
            switch(ido)
            {	
	        case -1:
	            /* Perform y <-- OP(x) = inv[A-sigma*M]*M*x to 
	            force starting vector into the range of OP.
	            User should supply 
	            a matrix vector multiplication routine [input: workd(ipntr(1))] and
	            a linear system solver [output: workd(ipntr(2))]. */
	            Multi_dmatcsr_dvec(M, workd+ipntr[0]-1, rhs);
#if ARPACK_LINEAR_SOLVER == CG
	            Linear_solver_cg(A, rhs, workd+ipntr[1]-1, cg_tol, cg_ite, NULL, &rnorm, &ncgit, 0);
		    //printf("DNEigenSolver (ido = -1): residual norm = %18.16f, niter = %d\n", rnorm, ncgit);
#else
		    Linear_solver_direct(A, rhs, workd+ipntr[1]-1);
		    //printf("DNEigenSolver (ido =  -1)\n");
#endif
	            break;
	        case 1:
	            /* Perform y <-- OP(x) = inv[A-sigma*M]*M*x where
	            M*x has been saved in workd(ipntr(3)).
	            User only need 
	            a linear system solver [input: workd(ipntr(3))] [output: workd(ipntr(2))]. */
#if ARPACK_LINEAR_SOLVER == CG
	            Linear_solver_cg(A, workd+ipntr[2]-1, workd+ipntr[1]-1, cg_tol, cg_ite, NULL, &rnorm, &ncgit, 0);
		    //printf("DNEigenSolver (ido =  1): residual norm = %18.16f, niter = %d\n", rnorm, ncgit);
#else
		    Linear_solver_direct(A, workd+ipntr[2]-1, workd+ipntr[1]-1);
		    //printf("DNEigenSolver (ido =  1)\n");
#endif
                    break;
                case 2:
                    /* Perform y <-- M*x.
	            User only need 
	            matrix vector multiplication routine
	            [input: workd(ipntr(1))] [output: workd(ipntr(2))]. */
	            Multi_dmatcsr_dvec(M, workd+ipntr[0]-1, workd+ipntr[1]-1);
		    //printf("DNEigenSolver (ido =  2)\n");
	            break;
	        default:
	            break;
            }
        }
        else if(iparam[6] == 1)
        {
            switch(ido)
            {	
	        case -1:
	        case  1:
	            /* Perform y <-- A*x.
	            User only need 
	            matrix vector multiplication routine
	            [input: workd(ipntr(1))] [output: workd(ipntr(2))]. */
	            Multi_dmatcsr_dvec(A, workd+ipntr[0]-1, workd+ipntr[1]-1);
	            break;
	        default:
	            break;
	    }
	}
    }while((ido==1) || (ido==-1) || (ido==2));
    
    if(info != 0)
    {
        printf("Error with dnaupd_, info = %d\n",info);
        printf("Converged Ritz value number = %d\n", iparam[4]);
        printf("Check documentation in dnaupd_\n\n");
        assert(1 == 0);
        exit(-1);
    } 

    //=================================================================
    // compute eigenvector
    //=================================================================
    dneupd_(&rvec, &howmny, select, dr, di, v, &ldv, &sigmar, &sigmai, workev,
	    &bmat, &n, which, &nev_local, &tol, resid, &ncv, v, &ldv, 
	    iparam, ipntr, workd, workl, &lworkl, &ierr);
    if (ierr != 0) 
    {
        printf("Error with zneupd, info = %d\n",ierr);
        printf("Check the documentation of dneupd_.\n\n");
        assert(1 == 0);
        exit(-1);
    }

    if(nev_local != nev)
    {
	assert(nev_local == nev + 1);
	printf("Warning: there may exist complex eigenvalue(s).\n");
    }
    
    /* for output */
    //for (i=0; i<nev; i++) printf("eval %d: %f\n", i, dr[i]);
    for (i=0; i<nev; i++) eval[i] = dr[i];

    if(NULL != evec)
    {
        for (i=0; i<nev; i++)
        {
	    for (j=0; j<n; j++)
	        evec[i][j] = v[i*n+j];
        }
    }   

    free(workev);
    free(di);
    free(dr);
    free(select);
    free(rhs);
    free(workl);
    free(workd);
    free(ipntr);
    free(iparam);
    free(v);
    free(resid);
    
    Insertion_ascend_sort_dvec_dvecvec(eval, evec, 0, nev-1);

#if 0
    printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
    double *vvv = (double*)calloc(nev+1, sizeof(double));
    for(i=0; i<nev; i++)
    {
	for(j=0; j<n-nev; j++)
	    vvv[i] += fabs(evec[i][j]);
    }
    for(j=0; j<n-nev; j++) vvv[nev] += fabs(v[nev*n+j]);
    printf("print nev+1 leading coefficients...\n");
    printf("%2d: ", -1);
    for(i=0; i<nev+1; i++)
	printf("%11.6f ", vvv[i]);
    printf("\n");
    for(j=0; j<nev; j++)
    {
	printf("%2d: ", j);
	for(i=0; i<nev; i++)
	    printf("%11.6f ", evec[i][n-nev+j]);
	printf("%11.6f ", v[nev*n+n-nev+j]);
	printf("\n");
    }
    printf("..........................................\n\n");
    double *vv2 = (double*)calloc(nev+1, sizeof(double));
    double *leading = (double*)calloc(nev+1, sizeof(double));
    for(i=0; i<nev; i++)
    {
	for(j=0; j<nev; j++)
	{
	    if(j != i)
		vv2[i] += fabs(evec[i][n-nev+j]);
	    else
		leading[i] = evec[i][n-nev+j];
	}
    }
    for(j=0; j<nev; j++) vv2[nev] += fabs(v[nev*n+n-nev+j]);
    leading[nev] = v[nev*(n+1)];
    for(i=0; i<nev+1; i++)
	printf("%2d: %20.15f %20.15f %20.15f\n", i, leading[i], vv2[i], vvv[i]);
    free(vvv);
    free(vv2);
    free(leading);
    printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n");
#endif
}

void Eigen_solver_arpack_dn_amg(multigrid *amg, int level, int nev, double *eval, double **evec, amg_param param)
{
    dmatcsr *A = amg->A[level];
    dmatcsr *M = amg->M[level];

    amg_param param_local = param;
    param_local.amgsolver_tol       = 1e-08;
    param_local.pcg_amg_tol         = 1e-16;
    param_local.pcg_amg_max_iter    = 1000;
    //param_local.amgsolver_max_cycle = 10;
    
    int nev_local = nev;
    /* parameters for calling dnaupd_ */
    int     ido      = 0;
    char    bmat     = (NULL==M)? 'I': 'G';
    int     n        = A->nr;//matrix dimension
    char   *which    = "LM";//LM or SM
    double  tol      = 0.0;
    double *resid    = calloc(n, sizeof(double)); 
    int     ncv      = 4*nev_local;//Arnoldi vectors number
    if(ncv > n) ncv  = n;
    int     ldv      = n;
    double *v        = calloc(ldv*ncv, sizeof(double));
    int    *iparam   = calloc(11, sizeof(int)); 
    iparam[0]        = 1;//selecting the implicit shifts
    iparam[2]        = n;//maximum number of Arnoldi iterations
    iparam[6]        = (NULL==M)? 1: 3;//eigenproblem type
    int    *ipntr    = calloc(14, sizeof(int));//output, pointer
    double *workd    = calloc(3*n, sizeof(double)); 
    int     lworkl   = 3*ncv*ncv + 6*ncv;
    double *workl    = calloc(lworkl, sizeof(double)); 
    int     info     = 0;//random initial vector
    
    /* tmp variables */
    double *rhs = calloc(n, sizeof(double));
    int i, j;
   
    /* parameters for calling dneupd_ */
    int     rvec   = 1;//compute the Ritz or Schur vectors
    char    howmny = 'A';//compute nev Ritz vectors
    int    *select = calloc(ncv, sizeof(int));//internal workspace
    double *dr     = calloc(ncv, sizeof(double));//real      part of the Ritz value
    double *di     = calloc(ncv, sizeof(double));//imaginary part of the Ritz value
    double  sigmar = 0.0;//real      part of the shift
    double  sigmai = 0.0;//imaginary part of the shift
    double *workev = (double*)calloc(3*ncv, sizeof(double)); 
    int     ierr;//info
    
    /* for amg linear solver */
    double resi_norm_amg      = 9999.9;
    int    niter              = -1;
    int    ncycle             = -1;
    
    //=================================================================
    // start reverse communivation process
    //=================================================================
    do
    {
        dnaupd_(&ido, &bmat, &n, which, &nev_local, &tol, resid, 
	        &ncv, v, &ldv, iparam, ipntr, workd, workl,
	        &lworkl, &info);
	/*
        printf("info = %d\n", info);
        printf("Arnoldi size = %d\n",iparam[4]);
        printf("iparam[6] = %d\n", iparam[6]);
	*/
        if(iparam[6] == 3)
        {
            //printf("ido = %d\n", ido);
            switch(ido)
            {	
	        case -1:
	            /* Perform y <-- OP(x) = inv[A-sigma*M]*M*x to 
	            force starting vector into the range of OP.
	            User should supply 
	            a matrix vector multiplication routine [input: workd(ipntr(1))] and
	            a linear system solver [output: workd(ipntr(2))]. */
	            
	            Multi_dmatcsr_dvec(M, workd+ipntr[0]-1, rhs);
		    //Linear_solver_amg(amg, level, rhs, workd+ipntr[1]-1, param_local, &resi_norm_amg, &ncycle);
		    //if(resi_norm_amg > param_local.amgsolver_tol) 
			//printf("DNEigenSolverAMG (ido = -1) ncycle = %3d, residual norm = %18.16f\n", ncycle, resi_norm_amg);
		    Linear_solver_pcg_amg(amg, level, rhs, workd+ipntr[1]-1, param_local, NULL, &resi_norm_amg, &niter, &ncycle);
		    if(resi_norm_amg > param_local.pcg_amg_tol) 
			printf("DNEigenSolverPCG (ido = -1) niter= %3d, namgcycle = %3d, rnorm = %18.16f\n", niter, ncycle, resi_norm_amg);
	            break;
	        case 1:
	            /* Perform y <-- OP(x) = inv[A-sigma*M]*M*x where
	            M*x has been saved in workd(ipntr(3)).
	            User only need 
	            a linear system solver [input: workd(ipntr(3))] [output: workd(ipntr(2))]. */
	            
		    //Linear_solver_amg(amg, level, workd+ipntr[2]-1, workd+ipntr[1]-1, param_local, &resi_norm_amg, &ncycle);
		    //if(resi_norm_amg > param_local.amgsolver_tol) 
			//printf("DNEigenSolverAMG (ido =  1) ncycle = %3d, residual norm = %18.16f\n", ncycle, resi_norm_amg);
		    Linear_solver_pcg_amg(amg, level, workd+ipntr[2]-1, workd+ipntr[1]-1, param_local, NULL, &resi_norm_amg, &niter, &ncycle);
	            if(resi_norm_amg > param_local.pcg_amg_tol) 
			printf("DNEigenSolverPCG (ido =  1) niter= %3d, namgcycle = %3d, rnorm = %18.16f\n", niter, ncycle, resi_norm_amg);
                    break;
                case 2:
                    /* Perform y <-- M*x.
	            User only need 
	            matrix vector multiplication routine
	            [input: workd(ipntr(1))] [output: workd(ipntr(2))]. */
	            Multi_dmatcsr_dvec(M, workd+ipntr[0]-1, workd+ipntr[1]-1);
		    //printf("DNEigenSolverPCG (ido =  2)\n");
	            break;
	        default:
	            break;
            }
        }
        else if(iparam[6] == 1)
        {
            switch(ido)
            {	
	        case -1:
	        case  1:
	            /* Perform y <-- A*x.
	            User only need 
	            matrix vector multiplication routine
	            [input: workd(ipntr(1))] [output: workd(ipntr(2))]. */
	            Multi_dmatcsr_dvec(A, workd+ipntr[0]-1, workd+ipntr[1]-1);
	            break;
	        default:
	            break;
	    }
	}
    }while((ido==1) || (ido==-1) || (ido==2));
    
    if(info != 0)
    {
        printf("Error with dnaupd_, info = %d\n",info);
        printf("Converged Ritz value number = %d\n", iparam[4]);
        printf("Check documentation in dnaupd_\n\n");
        assert(1 == 0);
        //exit(-1);
    } 

    //=================================================================
    // compute eigenvector
    //=================================================================
    dneupd_(&rvec, &howmny, select, dr, di, v, &ldv, &sigmar, &sigmai, workev,
	    &bmat, &n, which, &nev_local, &tol, resid, &ncv, v, &ldv, 
	    iparam, ipntr, workd, workl, &lworkl, &ierr);
    if (ierr != 0) 
    {
        printf("Error with zneupd, info = %d\n",ierr);
        printf("Check the documentation of dneupd_.\n\n");
        assert(1 == 0);
        //exit(-1);
    }

    if(nev_local != nev)
    {
	assert(nev_local == nev + 1);
	printf("Warning: there may exist complex eigenvalue(s).\n");
    }
    
    /* for output */
    for (i=0; i<nev; i++) eval[i] = dr[i];

    if(NULL != evec)
    {
        for (i=0; i<nev; i++)
        {
	    for (j=0; j<n; j++)
	        evec[i][j] = v[i*n+j];
        }
    }   

    free(workev);
    free(di);
    free(dr);
    free(select);
    free(rhs);
    free(workl);
    free(workd);
    free(ipntr);
    free(iparam);
    free(v);
    free(resid);
    
    Insertion_ascend_sort_dvec_dvecvec(eval, evec, 0, nev-1);
}

#if 0
//the subrutine for the symmetric matrix of the double precision
void dsaupd_(int *ido, //用来指示如何进行逆通讯来进行Arnoldi迭代
	    const char *bmat, // 表示特征值问题的类型， 当M为单位矩阵是bmat='I'，否则 bmat=‘G’
	    int *n, // 表示所要计算的特征值问题的维数（特征向量的维数）
	    const char *which, //指示所需要计算的特征值： ‘LM’，‘SM’
	    int *nev, //表示所需要计算特征值的个数
	    double *tol, //Arnoldi迭代停止的误差量
	    double *resid, //残量向量，每次Arnoldi迭代的初始向量
	    int *ncv,     //Arnoldi迭代过程中，Krylov子空间中的维数 (有多少个向量)
	    double *v,    /* 用来存Krylov子空间的向量，大小是 n*ncv */
	    int *ldv,  	  /* Krylov子空间中的向量V中每个向量的长度*/
	    int *iparam,  /* 参数表，里面记录了定义具体算法的一些参数 */
	    int *ipntr,   /* WORKD和WORKL的起始位置，在具体的算法实现中可以体会到 */
	   /* Note (from ./PARPACK/SRC/MPI/pdsaupd.f):
	    *     -------------------------------------------------------------
	    *     IPNTR(1): pointer to the current operand vector X in WORKD.
	    *     IPNTR(2): pointer to the current result vector Y in WORKD.
	    *     IPNTR(3): pointer to the vector B * X in WORKD when used in
	    *               the shift-and-invert mode.
	    *     IPNTR(4): pointer to the next available location in WORKL
	    *               that is untouched by the program.
	    *     IPNTR(5): pointer to the NCV by 2 tridiagonal matrix T in WORKL.
	    *     IPNTR(6): pointer to the NCV RITZ values array in WORKL.
	    *     IPNTR(7): pointer to the Ritz estimates in array WORKL associated
	    *               with the Ritz values located in RITZ in WORKL.
	    *     IPNTR(11): pointer to the NP shifts in WORKL. See Remark 6 below.
	    *
	    *     Note: IPNTR(8:10) is only referenced by pdseupd . See Remark 2.
	    *     IPNTR(8): pointer to the NCV RITZ values of the original system.
	    *     IPNTR(9): pointer to the NCV corresponding error bounds.
	    *     IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
	    *                of the tridiagonal matrix T. Only referenced by
	    *                pdseupd  if RVEC = .TRUE. See Remarks. */
	    double *workd,  /* 3*n大小的一个向量，支持算法在其中运行 work array of length 3*N */
	    double *workl,  /* work array of length LWORKL */
	    int *lworkl,  /* sizeof WORKL, at least NCV**2 + 8*NCV */
	    int *info    /* 输入信息 Input:
			 *	0:  a randomly initial residual vector is used
			 *	!0: RESID contains the initial residual vector
			 * Output:
			 *	0: normal exit
			 *	1: maximum number of iterations taken.
			 *	... ... */
	    );
		    
void dseupd_(int *rvec,    /* FALSE: Ritz values only, TRUE: Ritz vectors */
	    char *All,     /* 'A': all, 'S': some (specified by SELECT[]) */
	    int *select,   /* logical array of dimension NCV */
	    double *d,     /* eigenvalues, array of dimension NEV (output) */
	    double *v1,    /* Ritz vectors, N by NEV array (may use V) */
	    int *ldv1,     /* leading dimension of Z */
	    double *sigma, /* represents the shift */
	    /* the following arguments are the same as for pdsaupd */
	    const char *bmat, 
	    int *n,
	    const char *which, 
	    int *nev,
	    double *tol, 
	    double *resid, 
	    int *ncv,
	    double *v,
	    int *ldv,
	    int *iparam, 
	    int *ipntr,
	    double *workd,
	    double *workl,
	    int *lworkl,
	    int *ierr
	    );		  

/** calculate the eigenvalues of the system using ARPACK routines*/
//在双精度下求解对称特征值问题，我们这里假设对于被求解的刚度矩阵和质量矩阵已经
// 处理了边界条件，在这里的计算将不再处理边界条件
void DSEigenSolver(MATRIX *A,  MATRIX *M, int nev, int indictor,double tau,
		double *Evals, double **Evecs, bool Change)
{
/* 
c     %--------------------------------------------------------%
c     | The work array WORKL is used in DSAUPD as         |
c     | workspace.  Its dimension LWORKL is set as         |
c     | illustrated below.  The parameter TOL determines   |
c     | the stopping criterion.  If TOL<=0, machine        |
c     | precision is used.  The variable IDO is used for  |  
c     | reverse communication and is initially set to 0.  |
c     | Setting INFO=0 indicates that a random vector is  |
c     | generated in DSAUPD to start the Arnoldi          |
c     | iteration.                                        | 
c     %-----------------------------------------------------% */
  int i,j,k,N_Eqn, ActiveBound, begin, end;
  int n = A->nr; // 特征值计算的维数
  int ido = 0; 
  const char *bmat; //记录特征值计算的类型，‘G’ 或者'I'
  const char *which; //只是需要求解特征值的范围
  char All[4] = "All"; 
  double tol = 0.0;
  double *resid, *rhs;  //残量向量
  resid = malloc(n*sizeof(double)); 
  rhs = malloc(n*sizeof(double));
  int ncv = 4*nev;   //Krylov子空间的向量的个数，一般ncv > 2*nev (nev表示要求特征值的个数)
  if (ncv>n) 
    ncv = n;
  double *v;  //用来存Krylov子空间中的向量
  int ldv = n;
  //注意下面这个向量定义以及它的大小
  v = malloc(ldv*ncv*sizeof(double));
  //下面定义参数iparam，用到的就给出具体的值，不用到的就不要管，因为再算法执行过程中会赋值
  int *iparam;
  iparam = malloc(11*sizeof(int)); // new int[11];
  iparam[0] = 1;      //表示移动策略：1，表示用了shift
  iparam[2] = n;   //最大Arnoldi迭代次数
  iparam[6] = 3;     //只是计算模式
  int *ipntr;        //指示workd和workl的位置参数，与Arnold迭代之后的子问题的求解有关
  ipntr = malloc(11*sizeof(int)); 
  double *workd; //注意workd的大小
  workd = malloc(3*n*sizeof(double)); 
  double *workl; //workl是用来存储作QR分解的矩阵的
  workl = malloc(ncv*(ncv+8)*sizeof(double)); 
  int lworkl = ncv*(ncv+8);  //记录workl的大小
  int info = 0;
  int rvec = 1;    //（用来得到特征向量的指示子： TRUE表示要得到特征向量，FALSE：表示只需要得到特征值）
  int *select;
  select = malloc(ncv*sizeof(int)); 
  double *d;  //用来存储特征值 （QR迭代的时候也用它来存储所有小矩阵的特征值）
  d = malloc(2*ncv*sizeof(double)); 
  double sigma; 
  sigma = tau;
  int ierr; 
  int *Row, *KCol;
  double *Values, *Valus_A;
/*
c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 3 specified in the      |
c     | documentation of DSAUPD is used (IPARAM(7) = 3).  |
c     | All these options may be changed by the user.     |
c     | For details, see the documentation in DSAUPD.     |
c     %---------------------------------------------------% */
  if((M) && (tau>0.0))
  {
   // MatrixAdd(A, M, -tau);  // 计算矩阵 A=A-tau*M
  }
  if(indictor > 0)
    which = "LM";  //计算模最大的特征值
  else 
    which = "SM";  //计算模最小的特征值
  if(M==NULL)
    bmat = "I";   //（如果没有质量矩阵）表示标准特征值问题
  else 
    bmat = "G";   //（如果有质量矩阵）表示一般的特征值问题
    
  //下面进行特征值求解迭代
    /*
    c     %-----------------------------------------------%
    c     | M A I N   L O O P (Reverse communication) |
    c     %-----------------------------------------------%
    */
  do {
      /*
      c        %---------------------------------------------%
      c        | Repeatedly call the routine DSAUPD and take |
      c        | actions indicated by parameter IDO until    |
      c        | either convergence is indicated or maxitr   |
      c        | has been exceeded.                          |
      c        %---------------------------------------------% 
      */
      dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, 
	      &ncv, v, &ldv, iparam, ipntr, workd, workl,
	      &lworkl, &info); //不是指针的变量前面要加上&以取地址
      
      printf("Ipntr[1]=%d, Ipntr[2]=%d, Ipntr[3]=%d\n",ipntr[0],ipntr[1],ipntr[2]);
      
      switch(ido){
	/*
	c           %------------------------------------------------%
	c           | Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x  |
	c           | to force the starting vector into the      |
	c           | range of OP.  The user should supply       |
	c           | his/her own matrix vector multiplication   |
	c           | routine and a linear system solver here.   |
	c           | The matrix vector multiplication routine   |
	c           | takes workd(ipntr(1)) as the input vector. |
	c           | The final result is returned to            |
	c           | workd(ipntr(2)).                           |
	c           %------------------------------------------------%
	*/
	case -1: 
	  //先计算 y = M*x
	  printf("ido=%d, Matrix-vector and solving\n",ido);
	  MatrixDotVec(M, workd+ipntr[0]-1, rhs);
	  //再计算 inv[A-SIGMA*M]y=w
	  CG(A, rhs, workd+ipntr[1]-1,1e-15,1000,NULL); 
	  break;
	case 1:
	  printf("ido=%d, solving\n",ido);
	  //把 workd(ipntr(3))拷贝到 workd(ipntr(2))
	  //call dcopy ( n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
	  AssignVec(rhs,workd+ipntr[1]-1,n);
	  //计算 inv[A-SIGMA*M]y=w, 结果存在 workd(ipntr(2))
	  CG(A,rhs, workd+ipntr[1]-1,1e-15,1000,NULL); 
	  break;
	case 2:
	  //做举证相乘，输入workd(ipntr(1))，输出workd(ipntr(2))
	  //workd(ipntr(2)) = M*workd(ipntr(1))
	  printf("ido=%d, Matrix multipy\n",ido);
	  //计算 y = M*x
	  MatrixDotVec(M, workd+ipntr[0]-1, workd+ipntr[1]-1);
	  break;	  
      }	
    } while ((ido==1)||(ido==-1)||(ido==2));

  if (info<0) {
         printf("Error with dsaupd, info = %d\n ",info);
         printf("Check documentation in dsaupd\n\n");
  } else {
    dseupd_(&rvec, All, select, d, v, &ldv, &sigma, bmat,
	    &n, which, &nev, &tol, resid, &ncv, v, &ldv,
	    iparam, ipntr, workd, workl, &lworkl, &ierr);    

    if (ierr!=0) {
      printf("Error with dseupd, info = %d\n",ierr); 
      printf("Check the documentation of dseupd.\n\n");
    } else if (info==1) {
       printf("Maximum number of iterations reached.\n\n");
    } else if (info==3) {
      printf("No shifts could be applied during implicit\n");
      printf("Arnoldi update, try increasing NCV.\n\n");
    }
    //printf("sigma= %d\n",sigma); 
    //printf("ldv= %d\n",ldv); 
    //printf("n= %d\n",n);
    int i, j;
    for (i=0; i<nev; i++)
    {
      Evals[i] = d[i];
      printf("The %d-th eigenvalue: %18.15f\n",i,Evals[i]);
    }
    
    if(Evecs)
    {
      Evecs = malloc(nev*sizeof(double)); //new double *[nev];    
      for (i=0; i<nev; i++)
      {
	Evecs[i] = malloc(n*sizeof(double)); // new double [n];
	for (j=0; j<n; j++) Evecs[i][j] = v[i*n+j];
      }
    }
  }
  free(resid);
  free(v);
  free(iparam);
  free(ipntr);
  free(workd);
  free(workl);
  free(select);
  free(d);
}
#endif
