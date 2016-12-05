#include <stdlib.h>
#include <stdio.h>
#include "constant.h"
#include "matrix.h"
#include "tool.h"


void QuickSort(indexvec *a, int left, int right)
{
    if ( left < right )
    {
	int i = left;
	int j = right;
	indexvec key = *(a+left);
	while ( i<j )
	{
	    while(i<j && (a+j)->value<key.value) j--;
	    if(i<j) swap_indexvec(a+i++, a+j);
	    while(i<j && (a+i)->value>key.value) i++;
	    if(i<j) swap_indexvec(a+j--, a+i);
	}
	//swap_indexvec(a+i, &key);
	*(a+i) = key;
	QuickSort(a, left, i-1);
	QuickSort(a, i+1, right);
    }
}

void Insertion_sort(indexvec *a, int left, int right)
{
    int i, j;
    indexvec tmp;
    for(i=left; i<=right; i++)
    {
        tmp = *(a+i);
        for(j=i; j>left && (a+j-1)->value<tmp.value; j--)
            *(a+j) = *(a+j-1);
        *(a+j) = tmp;
    }
}

int Cmp(indexvec *a, indexvec *b)
{
    if(a->value < b->value)
        return 0;
    else if(a->value==b->value && a->index>b->index)
        return 0;
    else
        return 1;
}

void Insertion_sort2(indexvec *a, int left, int right)
{
    int i, j;
    indexvec tmp;
    for(i=left; i<=right; i++)
    {
        tmp = *(a+i);
        for(j=i; j>left; j--)
        {
            if((a+j-1)->value<tmp.value || ((a+j-1)->value==tmp.value && (a+j-1)->index>tmp.index))
            //if(!Cmp(a+j-1, &tmp))
            {
                *(a+j) = *(a+j-1);
            }
            else
                break;
        }
        *(a+j) = tmp;
    }
}

void swap_indexvec(indexvec *a, indexvec *b)
{
    indexvec tmp;
    tmp = *a;
    *a  = *b;
    *b  = tmp;
}

void swap_indexvec_c(indexvec *a, indexvec *b, int *c)
{
    indexvec tmp;
    tmp = *a;
    *a  = *b;
    *b  = tmp;
    
    int tmp_i;
    tmp_i = c[a->index];
    c[a->index] = c[b->index];
    c[b->index] = tmp_i;
}

int Search_index(indexvec *a, int left, int right, int index)
{
    int i;
    for(i=left; i<=right; i++)
    {
        if((a+i)->index == index)
            return i;
    }
    
    return -1;
}

void bubble_ascend_indexvec(indexvec *a, int where, int *c)
{
    while(where > 0)
    {
        if((a+where)->value > (a+where-1)->value)
        {
            swap_indexvec_c(a+where, a+where-1, c);
            where--;
        }
        else
            break;
    }
}

void bubble_descend_indexvec(indexvec *a, int len, int where, int *c)
{
    while(where < len-1)
    {
        if((a+where)->value < (a+where+1)->value)
        {
            swap_indexvec_c(a+where, a+where+1, c);
            where++;
        }
        else
            break;
    }
}

int test_indexvec(void)
{
#define len 9
    int i;
    int vec[len] = {3, 5, 8, 5, 9, 0, 5, 1, 2};
    indexvec idv[len];
    for(i=0; i<len; i++)
    {
        idv[i].index = i;
        idv[i].value = vec[i];
    }
    printf("before sorted...\n");
    for(i=0; i<len; i++)
    {
        printf("i = %d, index = %d, value = %d\n", i, idv[i].index, idv[i].value);
    }
    
    QuickSort(idv, 0, len-1);
    printf("after quick sorted...\n");
    for(i=0; i<len; i++)
    {
        printf("i = %d, index = %d, value = %d\n", i, idv[i].index, idv[i].value);
    }
    Insertion_sort2(idv, 0, len-1);
    printf("after insertion2 sorted...\n");
    for(i=0; i<len; i++)
    {
        printf("i = %d, index = %d, value = %d\n", i, idv[i].index, idv[i].value);
    }
    
    for(i=0; i<len; i++)
    {
        idv[i].index = i;
        idv[i].value = vec[i];
    }
    
    Insertion_sort(idv, 0, len-1);
    printf("after insertion sorted...\n");
    for(i=0; i<len; i++)
    {
        printf("i = %d, index = %d, value = %d\n", i, idv[i].index, idv[i].value);
    }
    
    int *c = (int*)malloc(len*sizeof(int));
    for(i=0; i<len; i++)
    {
        c[(idv+i)->index] = i;
    }
    
    printf("inherent index...\n");
    for(i=0; i<len; i++)
    {
        printf("i = %d, index = %d\n", i, c[i]);
    }
    
    idv[1].value = 0;
    idv[4].value = 10;
    bubble_ascend_indexvec(idv, 4, c);
    printf("after bubbled...\n");
    for(i=0; i<len; i++)
    {
        printf("i = %d, index = %d, value = %d\n", i, idv[i].index, idv[i].value);
    }
    printf("inherent index...\n");
    for(i=0; i<len; i++)
    {
        printf("i = %d, index = %d\n", i, c[i]);
    }
    
    
    free(c);
    return SUCCESS;
}
