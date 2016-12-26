#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "preprocess.h"
#include "amg_param.h"
#include "matrix.h"
#include "par_matrix_vector.h"
#include "cfsplit.h"
#include "linear_algebra.h"

#define PAR_CFSPLIT 0

#define  DEBUG_PAR_CLJP      0
#define ASSERT_PAR_CLJP      1

//static int ndof = 0;
extern int myrank_id;

static int Generate_par_strong_coupling_set_negtive (par_dmatcsr *A, par_imatcsr *S, amg_param param);
//static int Generate_par_strong_coupling_set_positive(par_dmatcsr *A, par_imatcsr *S, amg_param param);

//static int Get_par_independent_set_D(par_imatcsr *S, double *measure, int *upt_vec, int nUPT, int *dof);

static int Generate_par_strong_coupling_set_negtive (par_dmatcsr *A, par_imatcsr *S, amg_param param)
{
    const double sddt = param.strong_diagonally_dominant_threshold;

    //生成 S->diag
    const int A_diag_nr = A->diag->nr;
    const int A_diag_nc = A->diag->nc;
    const int A_diag_nn = A->diag->nn;
    int      *A_diag_ia = A->diag->ia;
    int      *A_diag_ja = A->diag->ja;
    double   *A_diag_va = A->diag->va;

    const int A_offd_nr = A->offd->nr;
    const int A_offd_nc = A->offd->nc;
    const int A_offd_nn = A->offd->nn;
    int      *A_offd_ia = A->offd->ia;
    int      *A_offd_ja = A->offd->ja;
    double   *A_offd_va = A->offd->va;
    
    assert(A_diag_nr == A_offd_nr);

    int *S_diag_ia = (int *)calloc(A_diag_nr+1, sizeof(int));
    int *S_diag_ja = (int *)calloc(A_diag_nn,   sizeof(int));

    int *S_offd_ia = (int *)calloc(A_offd_nr+1, sizeof(int));
    int *S_offd_ja = (int *)calloc(A_offd_nn,   sizeof(int));
    
    memset(S_diag_ja, -1, A_diag_nn*sizeof(int));
    memset(S_offd_ja, -1, A_offd_nn*sizeof(int));

    int nr = A_diag_nr;
    
    int i, j;
    double row_min;//非对角线最小的负元素
    double row_abs_sum;
    double tol; //强影响阈值
    double aii;
    
    for(i=0; i<nr; i++)
    {
	row_min     = 0.0;
	row_abs_sum = 0.0;
	aii         = 0.0;

	for(j=A_diag_ia[i]; j<A_diag_ia[i+1]; j++)
	{
	    row_abs_sum += MABS(A_diag_va[j]);
	    if(A_diag_ja[j]==i) {aii = A_diag_va[j]; continue;}
	    else row_min = MMIN(row_min, A_diag_va[j]);
	}
	for(j=A_offd_ia[i]; j<A_offd_ia[i+1]; j++)
	{
	    row_abs_sum += MABS(A_offd_va[j]);
	    row_min = MMIN(row_min, A_offd_va[j]);
	}
	tol = row_min * param.strong_connection_threshold;//tol <= 0

	if(row_abs_sum >= (1+sddt)*aii)//考虑对角占优
	{
	    for(j=A_diag_ia[i]; j<A_diag_ia[i+1]; j++)
	    {
		if(A_diag_va[j]<=tol && A_diag_ja[j]!=i)
		    S_diag_ja[j] = A_diag_ja[j];
	    }
	    for(j=A_offd_ia[i]; j<A_offd_ia[i+1]; j++)
	    {
		if(A_offd_va[j] <= tol)
		    S_offd_ja[j] = A_offd_ja[j];
	    }
	}
    }

    int index_diag = 0;
    int index_offd = 0;
    for(i=0; i<nr; i++)
    {
	S_diag_ia[i] = index_diag;
	for(j=A_diag_ia[i]; j<A_diag_ia[i+1]; j++)
	{
	    if(S_diag_ja[j] >= 0)
	    {
		S_diag_ja[index_diag] = S_diag_ja[j];
		index_diag++;
	    }
	}

	S_offd_ia[i] = index_offd;
	for(j=A_offd_ia[i]; j<A_offd_ia[i+1]; j++)
	{
	    if(S_offd_ja[j] >= 0)
	    {
		S_offd_ja[index_offd] = S_offd_ja[j];
		index_offd++;
	    }
	}
    }

    S_diag_ia[A_diag_nr] = index_diag;
    S_diag_ja = (int*)realloc(S_diag_ja, index_diag*sizeof(int));
    int  S_diag_nn = index_diag;
    int *S_diag_va = (int*)calloc(index_diag, sizeof(int));
    for(i=0; i<index_diag; i++) S_diag_va[i] = 1; //为了矩阵打印函数形式的统一

    S_offd_ia[A_offd_nr] = index_offd;
    S_offd_ja = (int*)realloc(S_offd_ja, index_offd*sizeof(int));
    int  S_offd_nn = index_offd;
    int *S_offd_va = (int*)calloc(index_offd, sizeof(int));
    for(i=0; i<index_offd; i++) S_offd_va[i] = 1; //为了矩阵打印函数形式的统一

    imatcsr *S_diag = (imatcsr*)malloc(sizeof(imatcsr));
    S_diag->nr = A_diag_nr;
    S_diag->nc = A_diag_nc;
    S_diag->nn = S_diag_nn;
    S_diag->ia = S_diag_ia;
    S_diag->ja = S_diag_ja;
    S_diag->va = S_diag_va;

    imatcsr *S_offd = (imatcsr*)malloc(sizeof(imatcsr));
    S_offd->nr = A_offd_nr;
    S_offd->nc = A_offd_nc;
    S_offd->nn = S_offd_nn;
    S_offd->ia = S_offd_ia;
    S_offd->ja = S_offd_ja;
    S_offd->va = S_offd_va;

    S->diag             = S_diag;
    S->offd             = S_offd;
    S->nr_global        = -1;
    S->nc_global        = -1;
    S->nn_global        = -1;
    S->map_offd_col_l2g = NULL;
    S->row_start        = NULL;
    S->col_start        = NULL;
    S->comm             = A->comm;
    S->comm_info        = NULL;

    return SUCCESS;
}

/* 1. 与fasp区别在于判断是否强连接时，fasp没有考虑A元素为正的情况，
      而且对于元素大小正好等于tol的时候，fasp默认是不算入强连接 
   2. 加入强对角占优的判断，
      if(row_abs_sum >= (1+max_row_sum)*MABS(aii)) 则不是强对角占优的；
	  若强对角占优，则本行没有强连接关系。
      如果不想判断，就把max_row_sum设为-1.
      fasp的判断条件是
      if(row_sum >= (2 - max_row_sum) * MABS(diag.val[i]) )
      即 (my)max_row_sum 相当于 (fasp)1-max_row_sum
*/
int Generate_par_strong_coupling_set(par_dmatcsr *A, par_imatcsr *S, amg_param param)
{
    if(param.positive_connected == NO)
	return Generate_par_strong_coupling_set_negtive (A, S, param);
    //else if(param.positive_connected == YES)
	//return Generate_par_strong_coupling_set_positive(A, S, param);
    else
	return FAIL;
}


#if PAR_CFSPLIT
static int nCPT = 0;
static int nFPT = 0;
static int nSPT = 0;

int Split_par_CLJP(par_dmatcsr *A, par_imatcsr *S, par_ivec *dof);
{
    MPI_Comm comm = A->comm;

    int  myrank, nproc_global;
    MPI_Comm_size(comm, &nproc_global);
    MPI_Comm_rank(comm, &myrank);
    myrank_id = myrank;

    int i, j, k;
    
    imatcsr *S_diag = S->diag;
    imatcsr *S_offd = S->offd;

    imatcsr *ST_diag = (imatcsr*)malloc(sizeof(imatcsr));
    imatcsr *ST_offd = (imatcsr*)malloc(sizeof(imatcsr));
    Transpose_imatcsr_struct(S_diag, ST_diag);
    Transpose_imatcsr_struct(S_offd, ST_offd);
    
    int ST_diag_nr = ST_diag->nr;
    int ST_offd_nr = ST_offd->nr;

    double *lambda_ST = (double*)malloc((ST_diag_nr+ST_offd_nr)*sizeof(double));
    for(i=0; i<ST_diag_nr; i++) lambda_ST[i]            = ST_diag_ia[i+1] - ST_diag_ia[i];
    for(i=0; i<ST_offd_nr; i++) lambda_ST[i+ST_diag_nr] = ST_offd_ia[i+1] - ST_offd_ia[i];

    
    srand(1);
    int *A_row_start = A->row_start; 
    int tmp;
    for(i=0; i<A_row_start[myrank]; i++) tmp = rand();

    double *lambda_ST = (double*)malloc(ST_nr*sizeof(double));
    for(i=0; i<ST_nr; i++) lambda_ST[i] = ST_ia[i+1]-ST_ia[i] + (double)rand()/RAND_MAX;

    int  A_nc  = A->nc;
    int *S_ia  = S->ia;
    int *S_ja  = S->ja;
    int *ST_ia = ST->ia;
    int   S_nn = S->nn;
    

    int nUPT = 0;
    int *upt_vec = (int*)calloc(ndof, sizeof(int));
    for(i=0; i<A_nc; i++)//for(i=0; i<ST_nr; i++)
    {
	/* fasp判断方法：if(S_ia[i+1] == S_ia[i])，参考fasp, coarsening_rs.c, cfsplitting_cls */
	if(S_ia[i+1] == S_ia[i]) /* i不受别的点影响（SPT或CPT），这个条件等价于A的第i行只有对角线为非零元，说明i可以直接求解出来，则i应为SPT.  */
	{
	    dof[i] = SPT;//SPT: 2
	    nSPT++;
	    lambda_ST[i] = 0.0;
	}
	else if(ST_ia[i+1] == ST_ia[i])/* i不影响别的点，但被别的点影响, 则i应为FPT或SPT, 但上一个if语句中判断是不是SPT，如果来到这里，说明为FPT */
	{
	    dof[i] = FPT;//FPT: 1
	    nFPT++;
	    lambda_ST[i] = 0.0;
	}
	else /* 如果上面两个if语句都不满足，说明i影响一些点，也有一些点影响i，UPT */
	{
	    dof[i] = UPT;
	    upt_vec[nUPT++] = i;
	    //nUPT++;
	}
    }
#if DEBUG_CLJP > 5
    printf("CLJP INI : ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d, nUPT = %d\n", ndof, nCPT, nFPT, nSPT, nUPT);
#endif

    int iUPT, jS, kS;
    int iwhile = 0;
    while(1)
    {
	iwhile++;
	//printf("head   -- iwhile = %d, nUPT = %d\n", iwhile, nUPT);
	//printf("CLJP: ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d, nUPT = %d\n", ndof, nCPT, nFPT, nSPT, nUPT);
	for(iUPT=0; iUPT<nUPT; iUPT++)
	{
	    i = upt_vec[iUPT];
	    //printf("dof[%d] = %d\n", i, dof[i]);
	    if((dof[i]!=CPT) && (lambda_ST[i]<1))
	    {
		dof[i] = FPT;
		//nFPT++;
		//make sure all dependencies have been accounted for
		for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
		{
		    if(S_ja[jS] > -1) 
		    {
			//if(SPT != dof[S_ja[jS]]) dof[i] = UPT;
			dof[i] = UPT;
			//nFPT--;
		    }
		}
	    }

	    if(CPT == dof[i])
	    {
		lambda_ST[i] = 0.0;
		nCPT++;
		nUPT--;
		upt_vec[iUPT] = upt_vec[nUPT];
		upt_vec[nUPT] = i;
		iUPT--;
	    }
	    else if(FPT == dof[i])
	    {
		lambda_ST[i] = 0.0;
		nFPT++;
		nUPT--;
		upt_vec[iUPT] = upt_vec[nUPT];
		upt_vec[nUPT] = i;
		iUPT--;
	    }
	}
	
	//printf("CLJP: ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d, nUPT = %d\n", ndof, nCPT, nFPT, nSPT, nUPT);
	//printf("middle -- iwhile = %d, nUPT = %d\n", iwhile, nUPT);
#if ASSERT_CLJP
	//if(nCPT+nFPT+nSPT+nUPT != ndof)
	//    printf("nCPT = %d, nFPT = %d, nSPT = %d, nUPT = %d, ndof = %d\n\n", nCPT, nFPT, nSPT, nUPT, ndof);
	//else
	//    printf("\n");
	assert(nCPT+nFPT+nSPT+nUPT == ndof);
#endif

	if(nUPT == 0) break;

	if(!Get_independent_set_D(S, ST, lambda_ST, upt_vec, nUPT, dof))
	{
	    printf("ERROR: Cannot find independent set.\n");
	    exit(-1);
	}

	for(iUPT=0; iUPT<nUPT; iUPT++)
	{
	    i = upt_vec[iUPT];
	    //printf("i = %2d, dof[%2d] = %d\n", i, i, dof[i]);
	    if((dof[i]==DPT) || (dof[i]==CPT) || (dof[i]==COMMON_CPT))
	    //if(dof[i] > 0)
	    {
		dof[i] = CPT;
		//nCPT++;
		for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
		{
		    j = S_ja[jS];
		    //printf("    j = %2d\n", j);
		    if(j > -1)
		    {
			S_ja[jS] = -S_ja[jS] - 1;
			if(dof[j] != DPT)
			{
			    lambda_ST[j]--;
			    //printf("        i = %2d, dof[%2d] = %d, lambda_ST[%2d] = %f\n", i, j, dof[j], j, lambda_ST[j]);
			}
		    }
		}
		//printf("\n");
	    }
	    else
	    {
		for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
		{
		    j = S_ja[jS];
		    if(j < 0) j = -j - 1;
		    if((dof[j]==DPT) || (dof[j]==CPT) || (dof[j]==COMMON_CPT))
		    //if(dof[i] > 0)
		    {
			if(S_ja[jS] > -1) S_ja[jS] = -S_ja[jS] - 1;
			dof[j] = COMMON_CPT;
		    }
		    else if(SPT == dof[j])
		    {
			if(S_ja[jS] > -1) S_ja[jS] = -S_ja[jS] - 1;
		    }
		}

		for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
		{
		    j = S_ja[jS];
		    if(j > -1)
		    {
			for(kS=S_ia[j]; kS<S_ia[j+1]; kS++)
			{
			    k = S_ja[kS];
			    if(k < 0) k = -k - 1;
			    if(dof[k] == COMMON_CPT)
			    {
				S_ja[jS] = -S_ja[jS] - 1;
				lambda_ST[j]--;
				break;
			    }
			}
		    }
		}
	    }

	    for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
	    {
		j = S_ja[jS];
		if(j < 0) j = -j - 1;
		if(dof[j] == COMMON_CPT) dof[j] = CPT;
	    }
	}
    }

    for(jS=0; jS<S_nn; jS++)
    {
	if(S_ja[jS] < 0) S_ja[jS] = -S_ja[jS] - 1;
    }

    free(lambda_ST);
    Free_imatcsr(ST_diag);
    Free_imatcsr(ST_offd);
    free(upt_vec);

#if DEBUG_CLJP > 5
    printf("CLJP: ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d\n", ndof, nCPT, nFPT, nSPT);
#endif
#if ASSERT_CLJP
    assert(nCPT+nFPT+nSPT == ndof);
#endif
    return nCPT;
}

//
//
//int Generate_par_sparsity_P_dir(par_imatcsr *S, par_ivec *dof, par_dmatcsr *P);
//int Generate_par_P_dir(par_dmatcsr *A, par_imatcsr *S, par_ivec *dof, par_dmatcsr *P);
//
//void Truncate_par_P(par_dmatcsr *P, amg_param param);


static int Get_independent_set_D(imatcsr *S, imatcsr *ST, double *measure, int *upt_vec, int nUPT, int *dof)
{
    int *S_ia = S->ia;
    int *S_ja = S->ja;
    //int *ST_ia = ST->ia;
    //int *ST_ja = ST->ja;

    int iUPT, i, jS, j;
    //int jST;
    for(iUPT=0; iUPT<nUPT; iUPT++)
    {
	i = upt_vec[iUPT];
	if(measure[i] > 1) dof[i] = DPT;
    }

    for(iUPT=0; iUPT<nUPT; iUPT++)
    {
	i = upt_vec[iUPT];
	for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
	{
	    j = S_ja[jS];
	    if(j < 0) j = -j - 1;

	    if(measure[j] > 1)
	    {
		if(measure[i] > measure[j])
		    dof[j] = NOT_DPT;
		else if(measure[i] < measure[j])
		    dof[i] = NOT_DPT;
	    }
	}

	/*
	for(jST=ST_ia[i]; jST<ST_ia[i+1]; jST++)
	{
	    j = ST_ja[jST];
	    if(j < 0) j = -j - 1;

	    if(measure[j] > 1)
	    {
		if(measure[i] > measure[j])
		    dof[j] = NOT_DPT;
		else if(measure[i] < measure[j])
		    dof[i] = NOT_DPT;
	    }
	}
	*/
    }

    for(iUPT=0; iUPT<nUPT; iUPT++)
    {
	i = upt_vec[iUPT];
	if(DPT == dof[i]) return TRUE;
    }
    return FALSE;
}

void Reset_dof(dmatcsr *A)
{
    nCPT = 0;
    nFPT = 0;
    nSPT = 0;
    ndof = A->nc;
}


static int Generate_strong_coupling_set_positive(dmatcsr *A, imatcsr *S, amg_param param)
{
    const double sddt = param.strong_diagonally_dominant_threshold;

    const int nr = A->nr;
    const int nc = A->nc;
    const int nn = A->nn;
    int      *ia = A->ia;
    int      *ja = A->ja;
    double   *va = A->va;
    
    ndof = nr;
    
    int i, j, index, row_begin, row_end;
    double row_abs_max;//非对角线绝对值最大的元素
    double row_abs_sum;
    double tol; //强影响阈值
    double aii;
    
    S->nr = nr;
    S->nc = nc;
    S->nn = nn;
    S->ia = (int *)calloc(nr+1, sizeof(int));
    S->ja = (int *)calloc(nn,   sizeof(int));
    S->va = NULL;
    
    memset(S->ja, -1, nn*sizeof(int));

    for(i=0; i<nr; i++)
    {
	row_begin = ia[i];
	row_end   = ia[i+1];

	row_abs_max = 0.0;
	row_abs_sum = 0.0;
	aii = 0.0;
	for(j=row_begin; j<row_end; j++)
	{
	    row_abs_sum += MABS(va[j]);
	    if(ja[j]==i) {aii = va[j]; continue;}
	    else         row_abs_max = MMAX(row_abs_max, MABS(va[j]));
	}
	tol = row_abs_max * param.strong_connection_threshold;
	//tol = MMAX(tol, eps);
	
#if DEBUG_GEN_S > 7
	printf("strong connection tolerance[%d] = %f\n", i, tol);
#endif

	if(row_abs_sum >= (1+sddt)*MABS(aii))//考虑对角占优
	{
	    for(j=row_begin; j<row_end; j++)
	    {
		if(MABS(va[j])>=tol && ja[j]!=i)
		    S->ja[j] = ja[j];
	    }
	}
    }

    index = 0;
    for(i=0; i<nr; i++)
    {
	S->ia[i] = index;
	row_begin = ia[i];
	row_end   = ia[i+1];
	for(j=row_begin; j<row_end; j++)
	{
	    if(S->ja[j] >= 0)
	    {
		S->ja[index] = S->ja[j];
		index++;
	    }
	}
    }

    S->ia[nr] = index;
    S->nn = index;
    S->ja = (int*)realloc(S->ja, index*sizeof(int));
    S->va = (int*)calloc(index, sizeof(int));
    for(i=0; i<index; i++) S->va[i] = 1; //为了矩阵打印函数形式的统一
#if 0    
    FILE *file = fopen("../output/S.dat","w");
    fprintf(file, "%d\n", S->nr);
    fprintf(file, "%d\n", S->nc);
    fprintf(file, "%d\n", S->nn);
    fprintf(file, "\n");
    for(i=0; i<S->nr+1; i++) fprintf(file, "%d\n", S->ia[i]);
    fprintf(file, "\n");
    for(i=0; i<S->nn; i++) fprintf(file, "%d\n", S->ja[i]);
    fclose(file);
#endif
    return SUCCESS;
}


int Generate_par_sparsity_P_dir(imatcsr *S, int *dof, dmatcsr *P)
{
    int i, j;
    
    int  nr   = S->nr;
    int *S_ia = S->ia;
    int *S_ja = S->ja;

    P->nr = ndof;
    P->nc = nCPT;
    P->nn = 0;
    P->ia = (int*)calloc(nr+1, sizeof(int));
    
    int  P_nn = 0;
    int *P_ia = P->ia;
    
    for(i=0; i<ndof; i++)
    {
        if(dof[i] == FPT)
        {
            for(j=S_ia[i]; j<S_ia[i+1]; j++)
	        if(dof[S_ja[j]] == CPT) P_nn++;
	}
	else if(dof[i] == CPT)
        {
            P_nn++;
        }
        else if(dof[i] == SPT)
        {
            /* do nothing */
        }
	P_ia[i+1] = P_nn;
    }
    
    P->nn = P_nn;
    P->ja = (int*)   calloc(P_nn, sizeof(int));
    P->va = (double*)calloc(P_nn, sizeof(double));
    
    int    *P_ja = P->ja;
    
    //compute the CPT index map from global to local 
    int *map_F2C = (int*)calloc(ndof, sizeof(int));
    memset(map_F2C, -1, ndof*sizeof(int));

    int CPT_index = 0;
    for(i=0; i<ndof; i++)
    {
        if(dof[i] == CPT) map_F2C[i] = CPT_index++;
    }
    
    int index = 0;
    for(i=0; i<ndof; i++)
    {
        if(dof[i] == FPT)
        {
            for(j=S_ia[i]; j<S_ia[i+1]; j++)
	        if(dof[S_ja[j]] == CPT) P_ja[index++] = map_F2C[S_ja[j]];
	}
	else if(dof[i] == CPT)
        {
            P_ja[index++] = map_F2C[i];
        }
        else if(dof[i] == SPT)
        {
            /* do nothing */
        }
    }
    
    free(map_F2C);
#if ASSERT_GEN_SPA_P
    assert(index == P_nn);
#endif
    return SUCCESS;
}

int Generate_P_dir(dmatcsr *A, imatcsr *S, int *dof, dmatcsr *P)
{
    int i, j, k, m;
    
    int *map_C2F = (int*)calloc(nCPT, sizeof(int));
    for(j=i=0; i<ndof; i++)
    {
        if(dof[i] == CPT) map_C2F[j++] = i;
    }
    
    int    *A_ia = A->ia;
    int    *A_ja = A->ja;
    double *A_va = A->va;
    
    int    *S_ia = S->ia;
    int    *S_ja = S->ja;
    
    int    *P_ia = P->ia;
    int    *P_ja = P->ja;
    double *P_va = P->va;
    
    int *Ci = (int*)malloc(ndof*sizeof(int));
    memset(Ci, -1, ndof*sizeof(int));
    
    double aii = 0.0;
    double spn = 0.0;//sum of positive N
    double snn = 0.0;//sum of negtive  N
    double spp = 0.0;//sum of positive P
    double snp = 0.0;//sum of negtive  P
    double alpha = 0.0;
    double beta  = 0.0;
    int    npc   = 0;//num_positive_coupling
    for(i=0; i<ndof; i++)
    {
        if(dof[i] == FPT)
        {
            aii = spn = snn = spp = snp = alpha = beta = 0.0;
            npc = 0;
            for(j=S_ia[i]; j<S_ia[i+1]; j++)
            {
                k = S_ja[j];
                if(dof[k] == CPT) Ci[k] = i;
            }
            
            for(j=A_ia[i]; j<A_ia[i+1]; j++)
            {
                if(A_ja[j] == i)
                {
                    aii = A_va[j];
                }
                else
                {
                    if(A_va[j] > 0.0)
                    {
                        spn += A_va[j];
                        if(Ci[A_ja[j]] == i)
                        /* whether A_ja[j] \in P_i=C_i^s */
                        {
                            spp += A_va[j];
                            npc++;
                        }
                    }
                    else//if(A_va[j] < 0.0)
                    {
                        snn += A_va[j];
                        if(Ci[A_ja[j]] == i)
                        {
                            snp += A_va[j];
                        }
                    }
                }
            }
            alpha = snn/snp;
            if(npc > 0)
            {
                beta = spn/spp;
            }
            else
            {
                beta = 0.0;
                aii += spn;
            }
            
            for(j=P_ia[i]; j<P_ia[i+1]; j++)
            {
                k = map_C2F[P_ja[j]];
                for(m=A_ia[i]; m<A_ia[i+1]; m++)
                {
                    if(A_ja[m] == k) break;
                }
                if(A_va[m] < 0.0)
                    P_va[j] = -alpha * A_va[m] / aii;
                else
                    P_va[j] = -beta  * A_va[m] / aii;
            }
        }
        else if(dof[i] == CPT)
        {         
            P_va[P_ia[i]] = 1;/* For CPT i, P(i,:) = (0, ..., 0, 1, 0, ..., 0) */
        }
        else if(dof[i] == SPT)
        {
            /* do nothing */
        }
    }

    free(Ci);
    free(map_C2F);
    return SUCCESS;
}

void Truncate_P(dmatcsr *P, amg_param param)
{
    double trunc = param.truncation_threshold;
    if(trunc >= 1.0)
    {
        printf("amg truncation threshold %f >= 1.0!\n", trunc);
        exit(-1);
    }
    
    int nr = P->nr;
    
    int new_nn = 0;
    int new_ja_index = 0;
    int new_va_index = 0;
    
    double minn, maxp;//min negtive/max positive
    double sn, sp; //sum of negtive/positive
    double tsn, tsp;//truncate sum of negtive/postive
    double val;
    double zoomn, zoomp;
    
    int i, j;
    int rb, re;//row begin/end
    
    for(i=0; i<nr; i++)
    {
        minn = maxp = 0.0;
        sn = sp = tsn = tsp = 0.0;
        
        rb = P->ia[i];
        re = P->ia[i+1];
        P->ia[i] = new_nn;
        
        for(j=rb; j<re; j++)
        {
            val = P->va[j];
            
            if(val < 0)
            {
                minn = MMIN(minn, val);
                sn  += val;
            }
            else if(val > 0)
            {
                maxp = MMAX(maxp, val);
                sp  += val;
            }
        }
        
        minn *= trunc;
        maxp *= trunc;
        
        for(j=rb; j<re; j++)
        {
            val = P->va[j];
            
            if(val <= minn)
            {
                new_nn++;
                tsn += val;
                P->ja[new_ja_index++] = P->ja[j];
            }
            else if(val >= maxp)
            {
                new_nn++;
                tsp += val;
                P->ja[new_ja_index++] = P->ja[j];
            }
        }
        
        zoomn = (tsn < -EPS) ? sn/tsn : 1.0;
        zoomp = (tsp >  EPS) ? sp/tsp : 1.0;
        
        for(j=rb; j<re; j++)
        {
            val = P->va[j];
            if     (val <= minn) P->va[new_va_index++] = val * zoomn;
            else if(val >= maxp) P->va[new_va_index++] = val * zoomp;
        }
    }
    
    P->nn     = new_nn;
    P->ia[nr] = new_nn;
    P->ja = (int*)   realloc(P->ja, new_nn*sizeof(int));
    P->va = (double*)realloc(P->va, new_nn*sizeof(double));
}

#endif
