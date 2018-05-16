#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "preprocess.h"
#include "io.h"
#include "matrix.h"
#include "arpack_interface.h"
#include "tool.h"

int print_rank = 0;
int main()
{
    //char Afile[256] = "../output/AH.dat";
    //char Mfile[256] = "../output/MH.dat";
    char Afile[256] = "../../../FEM_SOFT/dat/gmg_A_refine6.dat";
    char Mfile[256] = "../../../FEM_SOFT/dat/gmg_M_refine6.dat";
    dmatcsr *A = Read_dmatcsr(Afile);
    dmatcsr *M = Read_dmatcsr(Mfile);
    Print_dmatcsr(A);

    int i;
    int nev = 3;
    double  *eval = (double *)calloc(nev,  sizeof(double));
    double **evec = (double**)malloc(nev * sizeof(double*));
    for(i=0; i<nev; i++)
	evec[i] = (double*)calloc(A->nr, sizeof(double));
    
    Eigen_solver_arpack_dn(A, M, nev, eval, evec);
    
    printf("eigenvalues: \n");
    for(i=0; i<nev; i++)
	printf("%18.15f\n", eval[i]);

    for(i=0; i<nev; i++)
    {
	free(evec[i]);
	evec[i] = NULL;
    }
    free(evec);
    evec = NULL;

    free(eval);
    eval = NULL;

    Free_dmatcsr(M);
    Free_dmatcsr(A);
    return 0;
}
