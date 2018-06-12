#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "linear_algebra.h"
#include "preprocess.h"
#include "tool.h"

void Transpose_dmatcsr(dmatcsr *A, dmatcsr *AT)
{
    const int nr = A->nr;
    const int nc = A->nc;
    const int nn = A->nn;

    int    *A_ia = A->ia;
    int    *A_ja = A->ja;
    double *A_va = A->va;


    AT->nr = nc;;
    AT->nc = nr;
    AT->nn = nn;

    AT->ia = (int*)   calloc(AT->nr+1, sizeof(int));
    AT->ja = (int*)   calloc(nn,       sizeof(int));
    AT->va = (double*)calloc(nn,       sizeof(double));

    int    *AT_ia = AT->ia;
    int    *AT_ja = AT->ja;
    double *AT_va = AT->va;

    int i, j, k, m, row_begin, row_end;
    
    for(j=0; j<nn; j++)
    {
	i = A_ja[j];
	AT_ia[i+1]++;
    }
    for(i=2; i<=nc; i++)//i begins at 1 is OK, for AT->ia[0] = 0.
    {
	AT_ia[i] += AT_ia[i-1];
    }

    for(i=0; i<nr; i++)
    {
	row_begin = A_ia[i];
	row_end   = A_ia[i+1];

	for(k=row_begin; k<row_end; k++)
	{
	    j = A_ja[k];
	    m = AT_ia[j];
	    AT_ia[j] = m+1;
	    AT_ja[m] = i;
	    AT_va[m] = A_va[k];
	}
    }

    for(i=nc; i>0; i--)
    {
	AT_ia[i] = AT_ia[i-1];
    }
    AT_ia[0] = 0;
}

void Transpose_imatcsr(imatcsr *A, imatcsr *AT)
{
    const int nr = A->nr;
    const int nc = A->nc;
    const int nn = A->nn;

    int *A_ia = A->ia;
    int *A_ja = A->ja;
    int *A_va = A->va;


    AT->nr = nc;;
    AT->nc = nr;
    AT->nn = nn;
    
    AT->ia = (int*)calloc(AT->nr+1, sizeof(int));
    AT->ja = (int*)calloc(nn,       sizeof(int));
    AT->va = (int*)calloc(nn,       sizeof(int));
    
    int *AT_ia = AT->ia;
    int *AT_ja = AT->ja;
    int *AT_va = AT->va;

    int i, j, k, m, row_begin, row_end;
    for(j=0; j<nn; j++)
    {
	i = A_ja[j];
	AT_ia[i+1]++;
    }
    for(i=2; i<=nc; i++)// i begins at 1 is OK, for AT->ia[0] = 0.
    {
	AT_ia[i] += AT_ia[i-1];
    }
    
    for(i=0; i<nr; i++)
    {
	row_begin = A_ia[i];
	row_end   = A_ia[i+1];

	for(k=row_begin; k<row_end; k++)
	{
	    j = A_ja[k];
	    m = AT_ia[j];
	    AT_ia[j] = m+1;
	    AT_ja[m] = i;
	    AT_va[m] = A_va[k];
	}
    }
    for(i=nc; i>0; i--)
    {
	AT_ia[i] = AT_ia[i-1];
    }
    AT_ia[0] = 0;
}

void Transpose_imatcsr_struct(imatcsr *A, imatcsr *AT)
{
    const int nr = A->nr;
    const int nc = A->nc;
    const int nn = A->nn;

    int *A_ia = A->ia;
    int *A_ja = A->ja;

    AT->nr = nc;;
    AT->nc = nr;
    AT->nn = nn;
    AT->ia = (int*)calloc(AT->nr+1, sizeof(int));
    AT->ja = (int*)calloc(nn,       sizeof(int));
    AT->va = (int*)calloc(nn,       sizeof(int));
    int *AT_ia = AT->ia;
    int *AT_ja = AT->ja;

    int i, j, k, m, row_begin, row_end;
    for(j=0; j<nn; j++)
    {
	i = A_ja[j];
	AT_ia[i+1]++;
    }
    for(i=2; i<=nc; i++)// i begins at 1 is OK, for AT->ia[0] = 0.
    {
	AT_ia[i] += AT_ia[i-1];
    }
    
    for(i=0; i<nr; i++)
    {
	row_begin = A_ia[i];
	row_end   = A_ia[i+1];

	for(k=row_begin; k<row_end; k++)
	{
	    j = A_ja[k];
	    m = AT_ia[j];
	    AT_ia[j] = m+1;
	    AT_ja[m] = i;
	    //AT_va[m] = A_va[k];
	}
    }
   
    for(i=nc; i>0; i--)
    {
	AT_ia[i] = AT_ia[i-1];
    }
    AT_ia[0] = 0;
}

void Remove_zero_dmatcsr(dmatcsr *A)
{
    int nr = A->nr;
    int nn = A->nn;
    
    int    *ia_new = (int*)   calloc(nr+1, sizeof(int));
    int    *ja_new = (int*)   calloc(nn,   sizeof(int));
    double *va_new = (double*)calloc(nn, sizeof(double));
    
    memset(ja_new, -1, nn*sizeof(int));
    
    int i, j, rb, re, ip1;
    int nn_new = 0;
    for(i=0; i<nr; i++)
    {
        
        ip1 = i + 1;
        rb = A->ia[i];
        re = A->ia[ip1];
        
        ia_new[ip1] = ia_new[i];
        
        for(j=rb; j<re; j++)
        {
            if(MABS(A->va[j]) > MYAMGEPS)
            {
                ia_new[ip1]++;
                ja_new[nn_new] = A->ja[j];
                va_new[nn_new] = A->va[j];
                nn_new++;
            }
        }
    }
    
    int    *ja_new_tmp = realloc(ja_new, nn_new*sizeof(int));
    double *va_new_tmp = realloc(va_new, nn_new*sizeof(double));
    
    if((NULL != ja_new_tmp) && (NULL != va_new_tmp))
    {
        free(A->ia);
        free(A->ja);
        free(A->va);
        A->nn = nn_new;
        A->ia = ia_new;
        A->ja = ja_new_tmp;
        A->va = va_new_tmp;
    }
    else
    {
        free(ia_new);
        free(ja_new);
        free(va_new);
    }
}

dmatcsr *Expand_dmatcsr_struct(dmatcsr *A, int expan)
{
    dmatcsr *B = (dmatcsr*)malloc(sizeof(dmatcsr));
    
    int nr = A->nr + expan;
    int nc = A->nc + expan;
    int nn = A->nn + expan*A->nr + expan*A->nc + expan*expan;
    
    double *va = malloc(nn    *sizeof(double));
    int    *ja = malloc(nn    *sizeof(int));
    int    *ia = malloc((nr+1)*sizeof(int));

    int i, j, k;
    int start, end, length;
    ia[0] = 0;
    for(i=0; i<A->nr; i++)
    {
	start   = A->ia[i];
	end     = A->ia[i+1];
	length  = end - start;
	ia[i+1] = ia[i] + length + expan;
	memcpy(&ja[ia[i]], &A->ja[start], length*sizeof(int));
	memcpy(&va[ia[i]], &A->va[start], length*sizeof(double));
	for(j=0; j<expan; j++)
	{
	    k = ia[i] + length;
	    ja[k+j] = A->nc+j;
	    va[k+j] = 0;
	}
    }

    for(i=0; i<expan; i++)
    {
	k = A->nr+i;
	ia[k+1] = ia[k] + A->nc + expan;
	for(j=0; j<A->nc; j++)
	{
	    ja[ia[k]+j] = j;
	    va[ia[k]+j] = 0;
	}
	for(j=0; j<expan; j++)
	{
	    ja[ia[k]+A->nc+j] = A->nc+j;
	    va[ia[k]+A->nc+j] = 0;
	}
    }    
    B->nr = nr;
    B->nc = nc;
    B->nn = nn;
    
    B->ia = ia;
    B->ja = ja;
    B->va = va;

    return B;
}

void Expand_dmatcsr(dmatcsr *A, int expan, double **vec, double **mat)
{
    double *va = A->va;
    int    *ia = A->ia;

    int i, j, k;
    for(i=0; i<A->nr-expan; i++)
    {
	for(j=0; j<expan; j++)
	{
	    k = ia[i+1] - expan;
	    va[k+j] = vec[j][i];
	}
    }

    for(i=0; i<expan; i++)
    {
        k = A->nr - expan + i;
	for(j=0; j<A->nc-expan; j++)
	    va[ia[k]+j] = vec[i][j];
	for(j=0; j<expan; j++)
	    va[ia[k]+A->nc-expan+j] = mat[i][j];
    }
}

void Multi_dmatcsr_dmatcsr(dmatcsr *A, dmatcsr *B, dmatcsr *C)
{
    Multi_dmatcsr_dmatcsr_quick_sorted(A, B, C);
}

/**
 * A: m n
 * B: n p
 * C: m p
 */
void Multi_dmatcsr_dmatcsr_unsorted(dmatcsr *A, dmatcsr *B, dmatcsr *C)
{
    int m = A->nr;
    //int n = A->nc;
    int p = B->nc;
    
    int *A_ia = A->ia;
    int *B_ia = B->ia;
    
    int *A_ja = A->ja;
    int *B_ja = B->ja;
    
    double *A_va = A->va;
    double *B_va = B->va;

    C->nr = m;
    C->nc = p;
    int *C_ia = (int*)calloc(m+1, sizeof(int));

    int *work = (int*)calloc(p, sizeof(int));
    memset(work, -1, p*sizeof(int));
    
    int i, j, k, s, t;
    double ab;
    
    int row_begin;
    int count = 0;
    
    for(i=0; i<m; i++)
    {
        row_begin = count;
        C_ia[i] = row_begin;
        for(s=A_ia[i]; s<A_ia[i+1]; s++)
        {
            k = A_ja[s];
            for(t=B_ia[k]; t<B_ia[k+1]; t++)
            {
		j = B_ja[t];
		if(work[j] < row_begin)
		{
	            work[j] = count;
		    count++;
		}
            }
        }
    }
    C_ia[m] = count;
    C->nn   = count;
    
    int    *C_ja = (int*)   calloc(C->nn, sizeof(int));
    double *C_va = (double*)calloc(C->nn, sizeof(double));
    
    memset(work, -1, p*sizeof(int));
    count = 0;
    for(i=0; i<m; i++)
    {
        row_begin = count;
        for(s=A_ia[i]; s<A_ia[i+1]; s++)
        {
            k = A_ja[s];
            for(t=B_ia[k]; t<B_ia[k+1]; t++)
            {
                j = B_ja[t];
                ab = A_va[s] * B_va[t];
		if(work[j] < row_begin)
		{
		    work[j] = count;
		    C_va[count] = ab;
		    C_ja[count] = j;
		    count++;
		}
		else
		{
		    C_va[work[j]] += ab;
		}
            }
        }
    }
    
    C->ia = C_ia;
    C->ja = C_ja;
    C->va = C_va;
    
    free(work);
}

void Multi_dmatcsr_dmatcsr_quick_sorted(dmatcsr *A, dmatcsr *B, dmatcsr *C)
{
    Multi_dmatcsr_dmatcsr_unsorted(A, B, C);
    int C_nr = C->nr;
    int i;
    for(i=0; i<C_nr; i++)
        Quick_ascend_sort_ivec_dvec(C->ja, C->va, C->ia[i], C->ia[i+1]-1);
}

/*
 * A: m n
 * B: n p
 * C: p q
 * D: m q
 *
 * d_{ij] = \sum_{k=1}*p\sum_{s=1}^n a_{is}b_{sk}c_{kj}
 *
 * remark1: 按照上式遍历时，对第i行来说，要分别对 s, k, j 做循环遍历。
 *          每次更新 s 指标时，k, j, 都会重复循环。但首先确定非零元 d_{ij}
 *          个数时最关键的是 i, j, 因此如果 k 和 j 已经做过遍历，则无需
 *          二次做遍历，这就是 work1 和 work2 的作用
 *
 */
void Multi_dmatcsr_dmatcsr_dmatcsr_unsorted(dmatcsr *A, dmatcsr *B, dmatcsr *C, dmatcsr *D)
{
    int m = A->nr;
    //int n = A->nc;
    int p = B->nc;
    int q = C->nc;
    
    int *A_ia = A->ia;
    int *B_ia = B->ia;
    int *C_ia = C->ia;
    
    int *A_ja = A->ja;
    int *B_ja = B->ja;
    int *C_ja = C->ja;
    
    double *A_va = A->va;
    double *B_va = B->va;
    double *C_va = C->va;
    
    int *work1 = (int*)calloc(p, sizeof(int));
    int *work2 = (int*)calloc(q, sizeof(int));
    memset(work1, -1, p*sizeof(int));//mark B
    memset(work2, -1, q*sizeof(int));//mark C
    D->nr = m;
    D->nc = q;
    int *D_ia = (int*)calloc(m+1, sizeof(int));
    
    int i, j, k, s, i1, i2, i3;
    double ab, abc;
    
    int row_begin;
    int count = 0;
    
    for(i=0; i<m; i++)
    {
        row_begin = count;
        D_ia[i] = row_begin;
        for(i1=A_ia[i]; i1<A_ia[i+1]; i1++)
        {
            s = A_ja[i1];
            for(i2=B_ia[s]; i2<B_ia[s+1]; i2++)
            {
                k = B_ja[i2];
                if(work1[k] != i)//remark1
                {
                    work1[k] = i;
                    for(i3=C_ia[k]; i3<C_ia[k+1]; i3++)
                    {
                        j = C_ja[i3];
                        if(work2[j] < row_begin)
                        {
                            work2[j] = count;
                            count++;
                        }
                    }
                }
            }
        }
    }
    D_ia[m] = count;
    D->nn   = count;
    
    int    *D_ja = (int*)   calloc(D->nn, sizeof(int));
    double *D_va = (double*)calloc(D->nn, sizeof(double));
    
    count = 0;
    memset(work2, -1, q*sizeof(int));//work2[j] denotes j's location in ja
    for(i=0; i<m; i++)
    {
        row_begin = count;
        for(i1=A_ia[i]; i1<A_ia[i+1]; i1++)
        {
            s = A_ja[i1];
            for(i2=B_ia[s]; i2<B_ia[s+1]; i2++)
            {
                k = B_ja[i2];
                ab = A_va[i1] * B_va[i2];
                for(i3=C_ia[k]; i3<C_ia[k+1]; i3++)
                {
                    j = C_ja[i3];
                    abc = ab * C_va[i3];
                    if(work2[j] < row_begin)
                    {
                        work2[j] = count;
                        D_va[count] = abc;
                        D_ja[count] = j;
                        count++;
                    }
                    else
                    {
                        D_va[work2[j]] += abc;
                    }
                }
            }
        }
    }
    
    D->ia = D_ia;
    D->ja = D_ja;
    D->va = D_va;
    
    free(work2);
    free(work1);
}

void Multi_dmatcsr_dmatcsr_dmatcsr_quick_sorted(dmatcsr *A, dmatcsr *B, dmatcsr *C, dmatcsr *D)
{
    Multi_dmatcsr_dmatcsr_dmatcsr_unsorted(A, B, C, D);
    int D_nr = D->nr;
    int i;
    for(i=0; i<D_nr; i++)
        Quick_ascend_sort_ivec_dvec(D->ja, D->va, D->ia[i], D->ia[i+1]-1);
}

void Multi_dmatcsr_dmatcsr_dmatcsr_insertion_sorted(dmatcsr *A, dmatcsr *B, dmatcsr *C, dmatcsr *D)
{
    Multi_dmatcsr_dmatcsr_dmatcsr_unsorted(A, B, C, D);
    int D_nr = D->nr;
    int i;
    for(i=0; i<D_nr; i++)
        Insertion_ascend_sort_ivec_dvec(D->ja, D->va, D->ia[i], D->ia[i+1]-1);
}

#define UNSORTED        0
#define QUICKSORTED     1
#define INSERTIONSORTED 2
#define TOWSTEP         3
#define MULTIVERSION    QUICKSORTED
void Multi_dmatcsr_dmatcsr_dmatcsr(dmatcsr *A, dmatcsr *B, dmatcsr *C, dmatcsr *D)
{
    if(MULTIVERSION == QUICKSORTED)
        Multi_dmatcsr_dmatcsr_dmatcsr_quick_sorted(A, B, C, D);
    else if(MULTIVERSION == INSERTIONSORTED)
        Multi_dmatcsr_dmatcsr_dmatcsr_insertion_sorted(A, B, C, D);
    else if(MULTIVERSION == UNSORTED)
        Multi_dmatcsr_dmatcsr_dmatcsr_unsorted(A, B, C, D);
    else if(MULTIVERSION == TOWSTEP)
    {
	dmatcsr *tmp = (dmatcsr*)malloc(sizeof(dmatcsr));
	Multi_dmatcsr_dmatcsr(A,   B, tmp);
	Multi_dmatcsr_dmatcsr(tmp, C, D);
	Free_dmatcsr(tmp);
    }
}

/**
 * A: nr*nc
 * x: nc*1
 * y: nr*1
 */
void Multi_dmatcsr_dvec(dmatcsr *A, double *x, double *y)
{
    int nr = A->nr;
    int *ia = A->ia;
    int *ja = A->ja;
    double *va = A->va;
    int i, j, row_begin, row_end, row_length;
    
    memset(y, 0, nr*sizeof(double));
    for(i=0; i<nr; i++)
    {
        row_begin  = ia[i];
        row_end    = ia[i+1];
        row_length = row_end - row_begin;
        for(j=0; j<row_length; j++)
            y[i] += va[row_begin+j] * x[ja[row_begin+j]];
    }
}
void Multi_dmatcsr_dvec_hold(dmatcsr *A, double *x, double *y)
{
    int nr = A->nr;
    int *ia = A->ia;
    int *ja = A->ja;
    double *va = A->va;
    int i, j, row_begin, row_end, row_length;
    
    for(i=0; i<nr; i++)
    {
        row_begin  = ia[i];
        row_end    = ia[i+1];
        row_length = row_end - row_begin;
        for(j=0; j<row_length; j++)
            y[i] += va[row_begin+j] * x[ja[row_begin+j]];
    }
}

double Multi_dvec_dvec(double *x, double *y, int length)
{
    double value = 0.0;
    int i;
    for(i=0; i<length; i++) value += x[i] * y[i];
    return value;
}

double Multi_dvec_dmatcsr_dvec(double *x, dmatcsr *A, double *y)
{
  double *tmp = (double*)calloc(A->nr, sizeof(double));
  Multi_dmatcsr_dvec(A, y, tmp);
  double val = Multi_dvec_dvec(x, tmp, A->nr);
  free(tmp); tmp = NULL;
  return val;
}

/* z = a*x + b*y */
void Sum_dvec_axpby    (double *x, double a, double *y, double b, double *z, int length)
{
    int i;
    for(i=0; i<length; i++) z[i] = a*x[i] + b*y[i];
}

/* y = a*x + b*y */
void Sumself_dvec_axpby(double *x, double a, double *y, double b, int length)
{
    int i;
    for(i=0; i<length; i++) y[i] = a*x[i] + b*y[i];
}

int Max_dvec(double *x, int length, double *max, int *pos)
{
    int index = -1;
    int i;
    *max = x[0];
    for(i=1; i<length; i++)
    {
        if(*max<x[i])
        {
            *max = x[i];
            index = i;
        }
    }
    if(NULL != pos) *pos = index;
    return index;
}

int Max_ivec(int *x, int length, int *max, int *pos)
{
    int index = -1;
    int i;
    *max = x[0];
    for(i=1; i<length; i++)
    {
        if(*max<x[i])
        {
            *max = x[i];
            index = i;
        }
    }
    
    if(NULL != pos) *pos = index;
    return index;
}

int Max_abs_dvec(double *x, int length, double *max, int *pos)
{
    int index = 0;
    int i;
    double absv = 0.0;
    *max = fabs(x[0]);
    for(i=1; i<length; i++)
    {
        absv = fabs(x[i]);
        if(*max < absv)
        {
            *max = absv;
            index = i;
        }
    }
    *max = x[index];
    if(NULL != pos) *pos = index;
    return index;
}

void Normalize_dvec(double *x, int length)
{
    double max;
    Max_abs_dvec(x, length, &max, NULL);
    Scale_dvec(x, 1.0/max, length);
    
    if(fabs(max) < MYAMGEPS)
    {
        fprintf(stderr, "Warning in function \"NormalizeVec\": maximum of double vector [%18.15f] maybe equals to 0!\n", fabs(max));
    }
}

void Scale_dvec(double *x, double a, int length)
{
    int i;
    for(i=0; i<length; i++) x[i] *= a;
}

double Get_dvec_2norm(double *x, int length)
{
    return sqrt(Multi_dvec_dvec(x, x, length));
}

int Equal_dmatcsr_struct(dmatcsr *A, dmatcsr *M)
{
    if((A->nr == M->nr) && (A->nc == M->nc) && (A->nn == M->nn))
    {
	int nr = A->nr;
	int nn = A->nn;

	int *A_ia = A->ia;
	int *A_ja = A->ja;
	int *M_ia = M->ia;
	int *M_ja = M->ja;

	int i;

	for(i=0; i<nr+1; i++) if(A_ia[i] != M_ia[i]) return 0;
	for(i=1; i<nn;   i++) if(A_ja[i] != M_ja[i]) return 0;

	return 1;
    }
    return 0;
}

dmatcsr *Sum_dmatcsr_mApnB(double m, dmatcsr *A, double n, dmatcsr *B)
{
    if(!Equal_dmatcsr_struct(A, B))
    {
	printf("The sparsity structure between A and B is not the same!\n");
	printf("Return A instead.\n");
	return Copy_dmatcsr(A);
    }

    dmatcsr *C = Copy_dmatcsr(A);
    Sum_dvec_axpby(A->va, m, B->va, n, C->va, C->nn);
    return C;
}

dmatcsr *Sum_dmatcsr_dmatcsr(dmatcsr *A, dmatcsr *B)
{
    dmatcsr *C = (dmatcsr*)malloc(sizeof(dmatcsr));
    C->nr = A->nr;
    C->nc = A->nc;
    C->nn = 0;
    C->ia = (int*)calloc(C->nr+1, sizeof(int));

    int *A_ia = A->ia;
    int *B_ia = B->ia;
    int *A_ja = A->ja;
    int *B_ja = B->ja;
    double *A_va = A->va;
    double *B_va = B->va;

    int nr = C->nr;
    int nc = C->nc;
	
    int i;
    int repeat = 0;
    int aflag, amax, bflag, bmax;
    int acol, bcol;

    for (i=0; i<nr; ++i)
    {
	aflag = A_ia[i];
	amax  = A_ia[i+1]-1;

	bflag = B_ia[i];
	bmax  = B_ia[i+1]-1;

	while (aflag <= amax && bflag <= bmax)
	{
	    acol = A_ja[aflag];
	    bcol = B_ja[bflag];
	    if(bcol == acol)
	    {
		repeat++;
		aflag++;
		bflag++;
	    }
	    else if(bcol < acol)
		bflag++;
	    else
		aflag++;
	}
	C->ia[i+1] = A_ia[i+1] + B_ia[i+1] -repeat;
    }
    
    C->nn = C->ia[nr];
    C->va = (double*)malloc(C->nn * sizeof(double));
    C->ja = (int*)   malloc(C->nn * sizeof(int));

    int cflag = 0;
    for(i=0; i<nr; ++i)
    {
	aflag = A_ia[i];
	amax  = A_ia[i+1]-1;

	bflag = B_ia[i];
	bmax  = B_ia[i+1]-1;

	while(aflag<=amax || bflag<=bmax)
	{
	    if(aflag <= amax) 
		acol = A_ja[aflag];
	    else 
		acol = nc+1; 

	    if(bflag <= bmax)
		bcol = B_ja[bflag];
	    else
		bcol = nc+1;

	    if(bcol == acol)
	    {
		C->va[cflag] = A_va[aflag] + B_va[bflag];
		aflag++;
		bflag++;
		C->ja[cflag] = acol;
		cflag++;
	    }
	    else if(bcol < acol)
	    {
		C->va[cflag] = B_va[bflag];
		C->ja[cflag] = bcol;
		bflag++;
		cflag++;
	    }
	    else
	    {
		C->va[cflag] = A_va[aflag];
		C->ja[cflag] = acol;
		aflag++;
		cflag++;
	    }
	}
    }

    assert(cflag == C->nn);

    return C;
}

void Sort_dmatcsr_ja_by_insertion(dmatcsr *A)
{
    int nr = A->nr;
    int i;
    for(i=0; i<nr; i++)
        Insertion_ascend_sort_ivec_dvec(A->ja, A->va, A->ia[i], A->ia[i+1]-1);
}
