#ifndef __MATRIX_H__
#define __MATRIX_H__

typedef struct DOUBLE_MATRIX_CSR_
{
    int nr;
    int nc;
    int nn;
    
    int    *ia;                                 
    int    *ja;
    double *va;
} dmatcsr;

typedef struct INT_MATRIX_CSR_
{
    int nr;
    int nc;
    int nn;
    
    int *ia;                                 
    int *ja;
    int *va;
} imatcsr;

dmatcsr *Create_dmatcsr(const int nrow, const int ncol, const int nnz);
imatcsr *Create_imatcsr(const int nrow, const int ncol, const int nnz);

void Free_dmatcsr(dmatcsr *A);
void Free_imatcsr(imatcsr *A);

dmatcsr *Copy_dmatcsr(dmatcsr *A);

#endif
