#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <malloc.h>

#include "tool.h"
#include "preprocess.h"

double Get_time()
{
    struct timeval tv;
    if(gettimeofday(&tv, NULL) == -1) printf("Error in Get_time()!\n");
    return tv.tv_sec + (double)tv.tv_usec/1000000;
}

double Get_memory()
{
    struct mallinfo mi = mallinfo();
    return (double)(mi.usmblks+mi.uordblks)/1024.0/1024.0;
}

void Insertion_ascend_sort_dvec(double *a, int left, int right)
{
    int i, j;
    double tmp;
    for(i=left; i<=right; i++)
    {
        tmp = *(a+i);
        for(j=i; j>left && *(a+j-1)>tmp; j--)
        {
            *(a+j) = *(a+j-1);
        }
        *(a+j) = tmp;
    }
}

void Insertion_ascend_sort_ivec(int *a, int left, int right)
{
    int i, j;
    int tmp;
    for(i=left; i<=right; i++)
    {
        tmp = *(a+i);
        for(j=i; j>left && *(a+j-1)>tmp; j--)
        {
            *(a+j) = *(a+j-1);
        }
        *(a+j) = tmp;
    }
}

void Quick_ascend_sort_ivec(int *a, int left, int right)
{
    if ( left < right )
    {
	int i   = left;
	int j   = right;
	int key = *(a+left);
	int tmp;
	while ( i<j )
	{
	    while(i<j && *(a+j)>key) j--;
	    if(i<j)//swap_indexvec(a+i++, a+j);
	    {
	        tmp    = *(a+i);
	        *(a+i) = *(a+j);
	        *(a+j) = tmp;

	        i++;
	    }
	    while(i<j && *(a+i)<key) i++;
	    if(i<j)//swap_indexvec(a+j--, a+i);
	    {
	        tmp    = *(a+i);
	        *(a+i) = *(a+j);
	        *(a+j) = tmp;

	        j--;
	    }
	}
	*(a+i) = key;//swap_indexvec(a+i, &key);
	Quick_ascend_sort_ivec(a, left, i-1);
	Quick_ascend_sort_ivec(a, i+1, right);
    }
}

void Insertion_ascend_sort_ivec_dvec(int *a, double *b, int left, int right)
{
    int i, j;
    int tmp1;
    double tmp2;
    for(i=left; i<=right; i++)
    {
        tmp1 = *(a+i);
        tmp2 = *(b+i);
        for(j=i; j>left && *(a+j-1)>tmp1; j--)
        {
            *(a+j) = *(a+j-1);
            *(b+j) = *(b+j-1);
        }
        *(a+j) = tmp1;
        *(b+j) = tmp2;
    }
}

void Quick_ascend_sort_ivec_dvec(int *a, double *b, int left, int right)
{
    if ( left < right )
    {
	int i   = left;
	int j   = right;
	int key = *(a+left);
	int tmp1;
	double tmp2;
	while ( i<j )
	{
	    while(i<j && *(a+j)>key) j--;
	    if(i<j)//swap_indexvec(a+i++, a+j);
	    {
	        tmp1   = *(a+i);
	        *(a+i) = *(a+j);
	        *(a+j) = tmp1;
	        
	        tmp2   = *(b+i);
	        *(b+i) = *(b+j);
	        *(b+j) = tmp2;
	        
	        i++;
	    }
	    while(i<j && *(a+i)<key) i++;
	    if(i<j)//swap_indexvec(a+j--, a+i);
	    {
	        tmp1   = *(a+i);
	        *(a+i) = *(a+j);
	        *(a+j) = tmp1;
	        
	        tmp2   = *(b+i);
	        *(b+i) = *(b+j);
	        *(b+j) = tmp2;
	        
	        j--;
	    }
	}
	*(a+i) = key;//swap_indexvec(a+i, &key);
	Quick_ascend_sort_ivec_dvec(a, b, left, i-1);
	Quick_ascend_sort_ivec_dvec(a, b, i+1, right);
    }
}

int Is_ascend_sorted_ivec(int *vec, int length)
{
    int i;
    for(i=0; i<length-1; i++)
    {
        if(vec[i] > vec[i+1]) return FALSE;
    }
    return TRUE;
}



/* assume elements of both A and B are in ascend order */
int A_cap_B_sorted_ascend(int *A, int lengthA, int *B, int lengthB)
{
    int flag = 0;
    int i, j;
    for(i=0; i<lengthA; i++)
    {
        for(j=flag; j<lengthB; j++)
        {
            if(B[j] == A[i])
                return TRUE;
            else if(B[j] < A[i])
                flag++;
            else
                break;
        }
    }

    return FALSE;
}

/* assume elements of both A and B are in ascend order */
int Get_ivec_cap_ivec(int *A, int lengthA, int *B, int lengthB, int *ncap, int **cap, int **indexA, int **indexB)
{
    if((0==lengthA) && (0==lengthB))
    {
	ncap = 0;
	if(NULL != cap) *cap = NULL;
	if(NULL != indexA) *indexA = NULL;
	if(NULL != indexB) *indexB = NULL;
	return TRUE;
    }

    int i, j;

    int *idxA = (int*)malloc(lengthA * sizeof(int));
    int *idxB = (int*)malloc(lengthB * sizeof(int));
    for(i=0; i<lengthA; i++) idxA[i] = -1;
    for(i=0; i<lengthB; i++) idxB[i] = -1;

    int max_num_cap = (lengthA>lengthB) ? lengthB : lengthA;
    int *cap_vec = (int*)calloc(max_num_cap, sizeof(int));

    int num_cap = 0;
    int flag = 0;
    for(i=0; i<lengthA; i++)
    {
        for(j=flag; j<lengthB; j++)
        {
            if(B[j] == A[i])
	    {
		idxA   [num_cap] = i;
		idxB   [num_cap] = j;
		cap_vec[num_cap] = A[i];
		num_cap++;
	    }
            else if(B[j] < A[i])
                flag++;
            else
                break;
        }
    }

    *ncap = num_cap;

    if(NULL != cap) *cap = cap_vec;
    else            {free(cap_vec); cap_vec = NULL;}

    if(NULL != indexA) *indexA = idxA;
    else            {free(idxA); idxA = NULL;}

    if(NULL != indexB) *indexB = idxB;
    else            {free(idxB); idxB = NULL;}

    return TRUE;
}

int A_cap_B(int *A, int lengthA, int *B, int lengthB)
{
    int i, j;
    for(i=0; i<lengthA; i++)
    {
        for(j=0; j<lengthB; j++)
        {
            if(B[j] == A[i])
                return TRUE;
        }
    }

    return FALSE;
}

void Insertion_ascend_sort_dvec_dvecvec(double *a, double **b, int left, int right)
{
    int i, j;
    double  tmp1;
    double *tmp2;
    for(i=left; i<=right; i++)
    {
        tmp1 = *(a+i);
        tmp2 = *(b+i);
        for(j=i; j>left && *(a+j-1)>tmp1; j--)
        {
            *(a+j) = *(a+j-1);
            *(b+j) = *(b+j-1);
        }
        *(a+j) = tmp1;
        *(b+j) = tmp2;
    }
}

int Get_nonnegtive_ivec_all_value_ascend(int *x, int length, int max_val, int *nvalue, int *value)
{
    int i, v, n;
    int *work = (int*)malloc(max_val * sizeof(int));
    for(i=0; i<max_val; i++) work[i] = -1;

    n = 0;
    for(i=0; i<length; i++)
    {
	v = x[i];
	assert(v >= 0);
	if(-1 == work[v])
	{
	    work[v] = v;
	    n++;
	}
    }

    *nvalue = n;

    int index = 0;
    for(i=0; i<max_val; i++)
    {
	if(work[i] > -1)
	    value[index++] = work[i];
    }
    assert(index == n);
    *nvalue = n;

    free(work); work = NULL;

    return n;
}

//假定 a 和 b 中的元素都是从小到大排序的，并且 a 中的全体元素都介于 b 某两个相邻元素之间。
//将 a 插入到 b 中，返回插入后的向量，长度为 length_a+length_b
int *Insert_ascend_ivec_to_ivec(int *a, int length_a, int *b, int length_b, int *position)
{
    //assert((length_a!=0) || (length_b!=0));
    if(length_a==0 && length_b==0)
    {
	if(NULL != position) *position = -1;
	return NULL;
    }

    int i, j;
    int *c = (int*)malloc((length_a+length_b) * sizeof(int));

    if(0 == length_a)
    {
	for(i=0; i<length_b; i++)
	    c[i] = b[i];
	if(NULL != position) *position = -1;
	return c;
    }
    else if(0 == length_b)
    {
	for(i=0; i<length_a; i++)
	    c[i] = a[i];
	if(NULL != position) *position = -1;
	return c;
    }

    int index_c = 0;
    int is_a_insert = 0;
    for(i=0; i<length_b; i++)
    {
	if(0 == is_a_insert)
	{
	    if(a[0] < b[i])
	    {
		assert(a[length_a-1] < b[i]);
		if(NULL != position) *position = index_c;
		for(j=0; j<length_a; j++) 
		{
		    c[index_c++] = a[j];
		}

		is_a_insert = 1;
	    }
	}
	c[index_c++] = b[i];
    }
    if(0 == is_a_insert) 
    {
	for(j=0; j<length_a; j++) 
	    c[index_c++] = a[j];
	if(NULL != position) *position = length_b;
    }

    assert(index_c == length_a + length_b);
    return c;
}

