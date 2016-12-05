#include <string.h>
#include <stdlib.h>
#include "matrix.h"

dmatcsr *Create_dmatcsr(const int nr, const int nc, const int nn)
{
    dmatcsr *A = (dmatcsr*)malloc(sizeof(dmatcsr));

    A->nr = nr;
    A->nc = nc;
    A->nn = nn;

    if(nr > 0) A->ia = (int*)calloc(nr+1, sizeof(int));
    else       A->ia = NULL;

    if(nc > 0){A->ja = (int*)calloc(nn, sizeof(int));
               memset(A->ja, -1, nn*sizeof(int));}
    else       A->ja = NULL;

    if(nn > 0) A->va = (double*)calloc(nn, sizeof(double));
    else       A->va = NULL;

    return A;
}

imatcsr *Create_imatcsr(const int nr, const int nc, const int nn)
{
    imatcsr *A = (imatcsr*)malloc(sizeof(imatcsr));

    A->nr = nr;
    A->nc = nc;
    A->nn = nn;

    if(nr > 0) A->ia = (int*)calloc(nr+1, sizeof(int));
    else       A->ia = NULL;

    if(nc > 0){A->ja = (int*)calloc(nn, sizeof(int));
               memset(A->ja, -1, nn*sizeof(int));}
    else       A->ja = NULL;

    if(nn > 0) A->va = (int*)calloc(nn, sizeof(int));
    else       A->va = NULL;
    
    return A;
}

void Free_dmatcsr(dmatcsr *A)
{
    free(A->ia); A->ia = NULL;
    free(A->ja); A->ja = NULL;
    free(A->va); A->va = NULL;
    
    A->nr = 0;
    A->nc = 0;
    A->nn = 0;
    
    free(A);
    A = NULL;
}

void Free_imatcsr(imatcsr *A)
{
    free(A->ia); A->ia = NULL;
    free(A->ja); A->ja = NULL;
    free(A->va); A->va = NULL;
    
    A->nr = 0;
    A->nc = 0;
    A->nn = 0;
    
    free(A);
    A = NULL;
}

dmatcsr *Copy_dmatcsr(dmatcsr *A)
{
    dmatcsr *B = (dmatcsr*)malloc(sizeof(dmatcsr));
    
    B->nr = A->nr;
    B->nc = A->nc;
    B->nn = A->nn;
    
    B->ia = (int*)   malloc((A->nr+1)*sizeof(int));
    B->ja = (int*)   malloc( A->nn   *sizeof(int));
    B->va = (double*)malloc( A->nn   *sizeof(double));
    
    memcpy(B->ia, A->ia, (A->nr+1)*sizeof(int));
    memcpy(B->ja, A->ja,  A->nn   *sizeof(int));
    memcpy(B->va, A->va,  A->nn   *sizeof(double));
    
    return B;
}
