#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "preprocess.h"
#include "matrix.h"
#include "linear_algebra.h"
#include "cfsplit.h"
#include "list.h"
#include "amg_param.h"

#define  DEBUG_GEN_S     0
#define ASSERT_GEN_S     1

#define  DEBUG_PRE       0
#define ASSERT_PRE       1

#define  DEBUG_POST      0 
#define ASSERT_POST      1

#define  DEBUG_GEN_SPA_P 0
#define ASSERT_GEN_SPA_P 1

#define  DEBUG_GEN_P     0
#define ASSERT_GEN_P     1

#define  DEBUG_TRUNC     0
#define ASSERT_TRUNC     1

#define  DEBUG_CLJP      0
#define ASSERT_CLJP      1

static int nCPT = 0;
static int nFPT = 0;
static int nSPT = 0;
static int ndof = 0;

static int Generate_strong_coupling_set_negtive (dmatcsr *A, imatcsr *S, amg_param param);
static int Generate_strong_coupling_set_positive(dmatcsr *A, imatcsr *S, amg_param param);

static int Get_independent_set_D(imatcsr *S, imatcsr *ST, double *measure, int *upt_vec, int nUPT, int *dof);

void Reset_dof(dmatcsr *A)
{
    nCPT = 0;
    nFPT = 0;
    nSPT = 0;
    ndof = A->nc;
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
int Generate_strong_coupling_set(dmatcsr *A, imatcsr *S, amg_param param)
{
    if(param.positive_connected == NO)
	return Generate_strong_coupling_set_negtive (A, S, param);
    else if(param.positive_connected == YES)
	return Generate_strong_coupling_set_positive(A, S, param);
    else
	return FAIL;
}

static int Generate_strong_coupling_set_negtive(dmatcsr *A, imatcsr *S, amg_param param)
{
    const double sddt        = param.strong_diagonally_dominant_threshold;

    const int nr = A->nr;
    const int nc = A->nc;
    const int nn = A->nn;
    int      *ia = A->ia;
    int      *ja = A->ja;
    double   *va = A->va;
    
    ndof = nr;
    
    int i, j, index, row_begin, row_end;
    double row_min;//非对角线最小的负元素
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

	row_min     = 0.0;
	row_abs_sum = 0.0;
	aii         = 0.0;
	for(j=row_begin; j<row_end; j++)
	{
	    row_abs_sum += MABS(va[j]);
	    if(ja[j]==i) {aii = va[j]; continue;}
	    else         row_min = MMIN(row_min, va[j]);
	}
	tol = row_min * param.strong_connection_threshold;//tol <= 0
	
#if DEBUG_GEN_S > 7
	printf("strong connection tolerance[%d] = %f\n", i, tol);
#endif

	if(row_abs_sum >= (1+sddt)*aii)//考虑对角占优
	{
	    for(j=row_begin; j<row_end; j++)
	    {
		if(va[j]<=tol && ja[j]!=i)
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

int Pre_split(dmatcsr *A, imatcsr *S, int *dof)
{
    int i, j, k, t;
    
    imatcsr *ST = (imatcsr*)malloc(sizeof(imatcsr));
    Transpose_imatcsr_struct(S, ST);
    
    int  A_nc  = A->nc;
    int *S_ia  = S->ia;
    int *S_ja  = S->ja;
    int  ST_nr = ST->nr;
    int *ST_ia = ST->ia;
    int *ST_ja = ST->ja;
    
    int *lambda_ST = (int*)malloc(ST_nr*sizeof(int));
    for(i=0; i<ST_nr; i++) lambda_ST[i] = ST_ia[i+1]-ST_ia[i];

    int nUPT = 0;
    for(i=0; i<A_nc; i++)//for(i=0; i<ST_nr; i++)
    {
	/* fasp判断方法：if(S_ia[i+1] == S_ia[i])，参考fasp, coarsening_rs.c, cfsplitting_cls */
	if(S_ia[i+1] == S_ia[i]) /* i不受别的点影响（SPT或CPT），这个条件等价于A的第i行只有对角线为非零元，说明i可以直接求解出来，则i应为SPT.  */
	{
	    dof[i] = SPT;//SPT: 2
	    nSPT++;
	    lambda_ST[i] = 0;
	}
	else if(ST_ia[i+1] == ST_ia[i])/* i不影响别的点, 则i应为FPT或SPT, 但上一个if语句中判断是不是SPT，如果来到这里，说明为FPT */
	{
	    dof[i] = FPT;//FPT: 1
	    nFPT++;
	    lambda_ST[i] = 0;
	    for(k=S_ia[i]; k<S_ia[i+1]; k++)
	        if(dof[S_ja[k]] < 1)//neither SPT nor FPT
		    lambda_ST[S_ja[k]]++;
	}
	else /* 如果上面两个if语句都不满足，说明i影响一些点，也有一些点影响i，UPT */
	{
	    dof[i] = UPT;
	    nUPT++;
	}
    }
#if DEBUG_PRE > 5
    printf("INI : ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d\n", ndof, nCPT, nFPT, nSPT);
#endif

    Node  *node       = (Node *)malloc(ndof   *sizeof(Node));
    Node **local_head = (Node**)calloc(ndof*2, sizeof(Node*));
    Node **local_tail = (Node**)calloc(ndof*2, sizeof(Node*));
    List list;
    List_init(&list, node, local_head, local_tail, lambda_ST, ndof);
    
#if ASSERT_PRE
    assert(list.nlist == nUPT);
#endif

    int max_influence = ndof;
    int max_influence_pos = -1;
    
    while(nUPT > 0)
    {
	max_influence     = list.head->value;
	max_influence_pos = list.head->index;
#if ASSERT_PRE
        assert(nCPT+nFPT+nSPT+nUPT == ndof);
#endif
	if(max_influence == 0) break;
	
	if(UPT == dof[max_influence_pos])
	{
	    dof[max_influence_pos] = CPT;
	    nCPT++;
	    nUPT--;
	    lambda_ST[max_influence_pos] = 0;
	    List_delete(&list, node+max_influence_pos, local_head, local_tail);
	    node[max_influence_pos].value = 0;
#if ASSERT_PRE
	    assert(list.nlist == nUPT);
#endif
	    for(i=ST_ia[max_influence_pos]; i<ST_ia[max_influence_pos+1]; i++)
	    {
		j = ST_ja[i];
		if(dof[j] == UPT)
		{
		    dof[j] = FPT;
		    nFPT++;
		    nUPT--;
		    lambda_ST[j] = 0;
		    List_delete(&list, node+j, local_head, local_tail);
		    node[j].value = 0;
		    for(k=S_ia[j]; k<S_ia[j+1]; k++)
		    {
			if(dof[S_ja[k]] == UPT)
			{
			    lambda_ST[S_ja[k]]++;
			    List_delete(&list, &node[S_ja[k]], local_head, local_tail);
			    node[S_ja[k]].value++;
			    List_insert(&list, &node[S_ja[k]], local_head, local_tail);
			}
		    }
		}
	    }    
	    for(i=S_ia[max_influence_pos]; i<S_ia[max_influence_pos+1]; i++)
	    {
		j = S_ja[i];
		if(dof[j] == UPT)
		{
		    lambda_ST[j]--;
		    List_delete(&list, node+j, local_head, local_tail);
		    node[j].value--;
		    if(lambda_ST[j] > 0)
		    {
		        List_insert(&list, &node[j], local_head, local_tail);
		    }
		    else
		    {
		        /* 很少来到这一步，或者说测试有限元矩阵时候从未来到这一步，对应的 fasp 部分也是这样 */
		        dof[j] = FPT;
		        nFPT++;
		        nUPT--;
		        assert(list.nlist == nUPT);
		        for(t=S_ia[j]; t<S_ia[j+1]; t++)
		        {
		            k = S_ja[t];
		            if(dof[k] == UPT)
		            {
		                lambda_ST[k]++;
		                List_delete(&list, node+k, local_head, local_tail);
		                node[k].value++;
		                List_insert(&list, &node[k], local_head, local_tail);
		            }
		        }
		    }
                }
	    }
	}
    }
    free(local_tail);
    free(local_head);
    free(node);
    free(lambda_ST);
    Free_imatcsr(ST);

#if DEBUG_PRE > 5
    printf("PRE : ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d\n", ndof, nCPT, nFPT, nSPT);
#endif
#if ASSERT_PRE
    assert(nCPT+nFPT+nSPT == ndof);
#endif
    return nCPT;
}

int Pre_split_fasp(dmatcsr *A, imatcsr *S, int *dof)
{
    int i, j, k, t;
    
    imatcsr *ST = (imatcsr*)malloc(sizeof(imatcsr));
    Transpose_imatcsr_struct(S, ST);
    
    int  A_nc  = A->nc;
    int *S_ia  = S->ia;
    int *S_ja  = S->ja;
    int  ST_nr = ST->nr;
    int *ST_ia = ST->ia;
    int *ST_ja = ST->ja;
    
    int *lambda_ST = (int*)malloc(ST_nr*sizeof(int));
    for(i=0; i<ST_nr; i++) lambda_ST[i] = ST_ia[i+1]-ST_ia[i];

    Node  *node       = (Node *)malloc(ndof *sizeof(Node));
    Node **local_head = (Node**)calloc(ndof, sizeof(Node*));
    Node **local_tail = (Node**)calloc(ndof, sizeof(Node*));
    List list;
    list.nlist = 0;
    list.head  = NULL;
    
    int influence;

    int nUPT = 0;
    for(i=0; i<A_nc; i++)
    {
	if(S_ia[i+1]==S_ia[i])
	{
	    dof[i] = SPT;//SPT: 2
	    nSPT++;
	    lambda_ST[i] = 0;
	}
	else
	/* 如果上面两个if语句都不满足，说明i影响一些点，也有一些点影响i，UPT */
	{
	    dof[i] = UPT;
	    nUPT++;
	}
    }
    
    for(i=0; i<A_nc; i++)//for(i=0; i<ST_nr; i++)
    {
	if(dof[i] == SPT) continue;
	influence = lambda_ST[i];
	if(influence > 0)
	{
	    node[i].index = i;
	    node[i].value = lambda_ST[i];
	    List_insert(&list, &node[i], local_head, local_tail);
	}
	else
	{
	    dof[i] = FPT;
	    nFPT++;
	    nUPT--;
	    for(k=S_ia[i]; k<S_ia[i+1]; k++)
	    {
	        j = S_ja[k];
	        if(dof[j] == SPT) continue;
	        
	        if(j < i)
	        {
		    lambda_ST[j]++;
	            List_delete(&list, &node[j], local_head, local_tail);
	            node[j].value++;
	            List_insert(&list, &node[j], local_head, local_tail);
	        }
	        else
	        {
	            lambda_ST[j]++;
	        }
	    }
	}
    }
	
#if DUBUG_PRE > 5
    printf("INI : ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d\n", ndof, nCPT, nFPT, nSPT);
#endif
    
#if ASSERT_PRE
    assert(list.nlist == nUPT);
#endif

    int max_influence = ndof;
    int max_influence_pos = -1;
    
    while(nUPT > 0)
    {
	max_influence     = list.head->value;
	max_influence_pos = list.head->index;
#if ASSERT_PRE
        assert(nCPT+nFPT+nSPT+nUPT == ndof);
#endif
	if(max_influence == 0) break;
	
	if(UPT == dof[max_influence_pos])
	{
	    dof[max_influence_pos] = CPT;
	    nCPT++;
	    nUPT--;
	    lambda_ST[max_influence_pos] = 0;
	    List_delete(&list, node+max_influence_pos, local_head, local_tail);
	    node[max_influence_pos].value = 0;
#if ASSERT_PRE
	    assert(list.nlist == nUPT);
#endif
	    for(i=ST_ia[max_influence_pos]; i<ST_ia[max_influence_pos+1]; i++)
	    {
		j = ST_ja[i];
		if(dof[j] == UPT)
		{
		    dof[j] = FPT;
		    nFPT++;
		    nUPT--;
		    lambda_ST[j] = 0;
		    List_delete(&list, node+j, local_head, local_tail);
		    node[j].value = 0;
		    for(k=S_ia[j]; k<S_ia[j+1]; k++)
		    {
			if(dof[S_ja[k]] == UPT)
			{
			    lambda_ST[S_ja[k]]++;
			    List_delete(&list, &node[S_ja[k]], local_head, local_tail);
			    node[S_ja[k]].value++;
			    List_insert(&list, &node[S_ja[k]], local_head, local_tail);
			}
		    }
		}
	    }    
	    for(i=S_ia[max_influence_pos]; i<S_ia[max_influence_pos+1]; i++)
	    {
		j = S_ja[i];
		if(dof[j] == UPT)
		{
		    lambda_ST[j]--;
		    List_delete(&list, node+j, local_head, local_tail);
		    node[j].value--;
		    if(lambda_ST[j] > 0)
		    {
		        List_insert(&list, &node[j], local_head, local_tail);
		    }
		    else
		    {
		        /* 很少来到这一步，或者说测试有限元矩阵时候从未来到这一步，对应的 fasp 部分也是这样 */
		        dof[j] = FPT;
		        nFPT++;
		        nUPT--;
		        for(t=S_ia[j]; t<S_ia[j+1]; t++)
		        {
		            k = S_ja[t];
		            if(dof[k] == UPT)
		            {
		                lambda_ST[k]++;
		                List_delete(&list, node+k, local_head, local_tail);
		                node[k].value++;
		                List_insert(&list, &node[k], local_head, local_tail);
		            }
		        }
		    }
                }
	    }
	}
    }

    free(local_tail);
    free(local_head);
    free(node);
    free(lambda_ST);
    Free_imatcsr(ST);

#if DUBUG_PRE > 5
    printf("PRE : ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d\n", ndof, nCPT, nFPT, nSPT);
#endif
#if ASSERT_PRE
    assert(nCPT+nFPT+nSPT == ndof);
#endif
    return nCPT;
}

int Post_split(imatcsr *S, int *dof)
{
    int i, j, k, m;

    int *ia = S->ia;
    int *ja = S->ja;

    int *Ci = (int*)calloc(ndof, sizeof(int));
    memset(Ci,  -1, ndof*sizeof(int));
    
    int is_first = TRUE;
    int is_Sk_cap_Ci_empty = FALSE;
    int k_backup = -1;
    
    for(i=0; i<ndof; i++)
    {
        if(dof[i] == FPT)
        {
            for(j=ia[i]; j<ia[i+1]; j++)
                if(dof[ja[j]] == CPT) Ci[ja[j]] = i;
            
            is_first = TRUE;
            
            for(j=ia[i]; j<ia[i+1]; j++)
            {
                k = ja[j];
                is_Sk_cap_Ci_empty = TRUE;
                if(dof[k] == FPT)//k \in D_i^s
                {
                    for(m=ia[k]; m<ia[k+1]; m++)
                    {
                        if(Ci[ja[m]] == i)
                        {
                            is_Sk_cap_Ci_empty = FALSE;
                            break;
                        }
                    }
                    if(TRUE == is_Sk_cap_Ci_empty)
                    {
                        if(TRUE == is_first)
                        {
                            Ci[k] = i;
                            is_first = FALSE;
                            dof[k] = CPT;
                            nCPT++;
                            nFPT--;
                            k_backup = k;
                        }
                        else
                        {
                            dof[i] = CPT;
                            nCPT++;
                            nFPT--;
                            dof[k_backup] = FPT;
                            nFPT++;
                            nCPT--;
                            break;
                        }
                    }
                }
            }
        }
    }
    
#if DEBUG_POST > 5
    printf("POST: ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d\n\n", ndof, nCPT, nFPT, nSPT);
#endif
#if ASSERT_POST
    assert(nCPT+nFPT+nSPT == ndof);
#endif
    
    free(Ci);
    return nCPT;
}


int Generate_sparsity_P_dir(imatcsr *S, int *dof, dmatcsr *P)
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

int Generate_sparsity_P_std(imatcsr *S, int *dof, dmatcsr *P)
{
    int i, j, k, m, n;
    
    int  nr   = S->nr;
    int *S_ia = S->ia;
    int *S_ja = S->ja;

    P->nr = ndof;
    P->nc = nCPT;
    P->nn = 0;
    P->ia = (int*)calloc(nr+1, sizeof(int));
    
    int  P_nn = 0;
    int *P_ia = P->ia;
    
    int *visited = (int*)calloc(nr, sizeof(int));
    memset(visited, -1, nr*sizeof(int));
    
    for(i=0; i<ndof; i++)
    {
        if(dof[i] == FPT)
        {
            for(j=S_ia[i]; j<S_ia[i+1]; j++)
            {
	        k = S_ja[j];
	        if((dof[k]==CPT) && (visited[k]!=i))
	        {
	            visited[k] = i;
	            P_nn++;
	        }
	        else if(dof[k]==FPT)
	        {
	            for(m=S_ia[k]; m<S_ia[k+1]; m++)
	            {
	                n = S_ja[m];
	                if((dof[n]==CPT) && (visited[n]!=i))
	                {
	                    visited[n] = i;
	                    P_nn++;
	                }
	            }
	        }
	    }
        }
        else if(dof[i] == CPT)
        {
            P_nn++;
        }
        else
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
    memset(visited, -1, nr*sizeof(int));
    for(i=0; i<ndof; i++)
    {
        if(dof[i] == FPT)
        {
            for(j=S_ia[i]; j<S_ia[i+1]; j++)
            {
	        k = S_ja[j];
	        if((dof[k]==CPT) && (visited[k]!=i))
	        {
	            visited[k] = i;
	            P_ja[index++] = map_F2C[k];
	        }
	        else if(dof[k]==FPT)
	        {
	            for(m=S_ia[k]; m<S_ia[k+1]; m++)
	            {
	                n = S_ja[m];
	                if((dof[n]==CPT) && (visited[n]!=i))
	                {
	                    visited[n] = i;
	                    P_ja[index++] = map_F2C[n];
	                }
	            }
	        }
	    }
	}
        else if(dof[i] == CPT)
        {
            P_ja[index++] = map_F2C[i];
        }
    }
    
    free(map_F2C);
    free(visited);
#if ASSERT_GEN_SPA_P
    assert(index == P_nn);
#endif
    return SUCCESS;
}

/* 与fasp区别主要有两个：
   1.这里对元素符号做了分开考虑，fasp没有；
   2.fasp（interp_STD 函数，interpolation.c文件）里面针对 RS_C1 条件做了修改，
     #if RS_C1
            alN = psum[i];（默认，不计算SPT的值）
     #else
            alN = nsum[i];
     #endif
     即对 Ni 求和时，考不考虑 SPT 点的值。
     
  如果不对 SPT 剔除（RS_C1 == 0），并且
  对于非对角线元素非正的矩阵（L矩阵），这里的程序和fasp的输出是一样的。
  另见 Generate_strong_coupling_set 的函数说明。
*/
int Generate_P_std(dmatcsr *A, imatcsr *S, int *dof, dmatcsr *P)
{
    int i, j, k, m, n;
    
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
    
    /*strong coupling CPT*/
    int *Ci = (int*)malloc(ndof*sizeof(int));
    memset(Ci, -1, ndof*sizeof(int));
    
    /*all sum of positive N*/
    double *aspn = (double*)calloc(ndof, sizeof(double));
    /*all sum of negtive  N*/
    double *asnn = (double*)calloc(ndof, sizeof(double));
    /*all sum of positive C*/
    double *aspc = (double*)calloc(ndof, sizeof(double));
    /*all sum of negtive  C*/
    double *asnc = (double*)calloc(ndof, sizeof(double));
    /*all diagonal element of A*/
    double *diag = (double*)calloc(ndof, sizeof(double));
    /*all number of positive_couplings*/
    int    *anpc = (int*)   calloc(ndof, sizeof(int));
    
    for(i=0; i<ndof; i++)
    {
        for(j=S_ia[i]; j<S_ia[i+1]; j++)
        {
            k = S_ja[j];
            if(dof[k] == CPT) Ci[k] = i;
        }
        for(j=A_ia[i]; j<A_ia[i+1]; j++)
        {
            k = A_ja[j];
            if(k != i)//not diag entry
            {
                if(A_va[j] < 0.0)
                {
                    asnn[i] += A_va[j];
                    if(Ci[k] == i)
                    {
                        asnc[i] += A_va[j];
                    }
                }
                else if(A_va[j] > 0.0)
                {
                    aspn[i] += A_va[j];
                    if(Ci[k] == i)
                    {
                        aspc[i] += A_va[j];
                        anpc[i] ++;
                    }
                }
            }
            else//diag entry
            {
                diag[i] = A_va[j];
            }
        }
    }
    
    double    *Ahat = (double*)calloc(ndof, sizeof(double));
    
    int *map_i_j2ja = (int*)   calloc(ndof, sizeof(int));
    memset(map_i_j2ja, -1, ndof*sizeof(int));

    int *map_k_j2ja = (int*)   calloc(ndof, sizeof(int));
    memset(map_k_j2ja, -1, ndof*sizeof(int));
    
    double aik, akk, aki, akn, ratio;
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
            alpha = beta = 0.0;
            npc = anpc[i];
            spn = aspn[i];
            snn = asnn[i];
            spp = aspc[i];
            snp = asnc[i];
            
            for(j=A_ia[i]; j<A_ia[i+1]; j++)
                map_i_j2ja[A_ja[j]] = j;
            for(j=P_ia[i]; j<P_ia[i+1]; j++)
                Ahat[map_C2F[P_ja[j]]] = 0.0;
            Ahat[i] = diag[i];
            
            for(j=S_ia[i]; j<S_ia[i+1]; j++)
            {
                k = S_ja[j];
                /*S在第i行的所有列指标是A在第i行所有列指标的子集，
                  所以下面一行对S的列指标做映射可以保证在上面给
                  map_i_j2ja赋值的范围内*/
                aik = A_va[map_i_j2ja[k]];
                
                if(dof[k] == FPT)
                {
                    akk = diag[k];
                    ratio = aik/akk;
                    for(m=A_ia[k]; m<A_ia[k+1]; m++)
                        map_k_j2ja[A_ja[m]] = m;
                    
                    /*对F_i^s中的k点做标准插值*/
                    /*1.先把N_k中等于i的点加入到Ahat[i]*/
                    aki = 0.0;
                    for(m=A_ia[k]; m<A_ia[k+1]; m++)
                    {
                        if(A_ja[m] == i)
                        {
                            aki = A_va[m];
                            Ahat[i] -= ratio*aki;
                        }
                    }
                    /*2.把N_k中的其他点用C_k中的点表示*/
                    for(m=S_ia[k]; m<S_ia[k+1]; m++)
                    {
                        n = S_ja[m];
                        akn = A_va[map_k_j2ja[n]];
                        if(dof[n] == CPT)
                        {
                            Ahat[n] -= ratio*akn;
                        }
                    }
                    
                    if(aik < 0)//ratio < 0, k \in F_i^{s,-}
                    {
                        snn -= ratio*(akk+asnn[k]);
                        if(aki < 0)
                            snn += ratio*aki;
                        spn -= ratio*aspn[k];
                        if(aki > 0)
                            spn += ratio*aki;
                            
                        snp -= ratio*asnc[k];
                        spp -= ratio*aspc[k];
                        
                        if(aspc[k] > 0)
                           npc = 1;
                    }
                    else if(aik > 0)//ratio > 0, k \in F_i^{s,+}
                    {
                        snn -= ratio*aspn[k];
                        if(aki > 0)
                            snn += ratio*aki;
                        spn -= ratio*(akk+asnn[k]);
                        if(aki < 0)
                            spn += ratio*aki;
                            
                        snp -= ratio*aspc[k];
                        spp -= ratio*asnc[k];
                        
                        if(asnc[k] < 0)
                           npc = 1;
                    }
                }
                else if(dof[k] == CPT)
                {
                    Ahat[k] += aik;
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
                Ahat[i] += spn;
            }
            
            for(j=P_ia[i]; j<P_ia[i+1]; j++)
            {
                k = map_C2F[P_ja[j]];
                if(Ahat[k] < 0.0)
                {
                    P_va[j] = -alpha * Ahat[k] / Ahat[i];
                }
                else
                {
                    P_va[j] = -beta  * Ahat[k] / Ahat[i];
                }
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

    free(map_k_j2ja);
    free(map_i_j2ja);
    free(Ahat);
    free(anpc);
    free(diag);
    free(asnc);
    free(aspc);
    free(asnn);
    free(aspn);
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


int CLJP_split(dmatcsr *A, imatcsr *S, int *dof)
{
    //printf("CLJP spliting...\n");

    int i, j, k;
    
    imatcsr *ST = (imatcsr*)malloc(sizeof(imatcsr));
    Transpose_imatcsr_struct(S, ST);
    
    int  A_nc  = A->nc;
    int *S_ia  = S->ia;
    int *S_ja  = S->ja;
    int  ST_nr = ST->nr;
    int *ST_ia = ST->ia;
    int   S_nn = S->nn;
    
    srand(1);
    double *lambda_ST = (double*)malloc(ST_nr*sizeof(double));
    for(i=0; i<ST_nr; i++) lambda_ST[i] = ST_ia[i+1]-ST_ia[i] + (double)rand()/RAND_MAX;

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
    Free_imatcsr(ST);
    free(upt_vec);

#if DEBUG_CLJP > 5
    printf("CLJP: ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d\n", ndof, nCPT, nFPT, nSPT);
#endif
#if ASSERT_CLJP
    assert(nCPT+nFPT+nSPT == ndof);
#endif
    return nCPT;
}

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
