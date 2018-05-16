#ifndef __TOOL__
#define __TOOL__

#include "matrix.h"

typedef struct INDEX_VEC
{
    int index;
    int value;
} indexvec;

void QuickSort(indexvec *a, int left, int right);
void Insertion_sort(indexvec *a, int left, int right);
void Insertion_sort2(indexvec *a, int left, int right);
int Search_index(indexvec *a, int left, int right, int index);

// 1: >=
// 0: <
int Cmp(indexvec *a, indexvec *b);

void swap_indexvec(indexvec *a, indexvec *b);

void bubble_ascend_indexvec (indexvec *a, int where, int *c);
void bubble_descend_indexvec(indexvec *a, int len, int where, int *c);

//int test_indexvec(dmatcsr *A, imatcsr *S, int *dof);

int test_indexvec(void);







#endif
