#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "preprocess.h"
#include "matrix.h"
#include "io.h"
#include "linear_algebra.h"
#include "tool.h"

#define file_prefix_41 "fdm2d_9P_bnd_16384x16384"

#define file_prefix_55 "fem_poisson_refine_5"
#define file_prefix_56 "fem_poisson_refine_6"
#define file_prefix_57 "fem_poisson_refine_7"
#define file_prefix_58 "fem_poisson_refine_8"
#define file_prefix_59 "fem_poisson_refine_9"

#define file_prefix    file_prefix_41


int main()
{
#if 0
    printf("=============== test matrix ===============\n");
    /*===================================================================*/
    /*===================================================================*/
    printf("1: Create_dmatcsr...\n");
    dmatcsr *mat_cre = Create_dmatcsr(10, 8, 5);
    Write_dmatcsr_csr(mat_cre, "../output/mat_cre.dmatcsr");
    /*********************************************************************/
    Free_dmatcsr(mat_cre);
    /*===================================================================*/
    /*===================================================================*/
    printf("2: Read, Write, Print, Copy...\n");
    dmatcsr *mat_read1 = Read_dmatcsr("../dat/"file_prefix"_A.dmatcsr");
    //Print_dmatcsr    (mat_read1);
    Write_dmatcsr_csr(mat_read1, "../output/"file_prefix"_read1.dmatcsr");
    /*-------------------------------------------------------------------*/
    dmatcsr *mat_read2 = Read_dmatcsr("../output/"file_prefix"_read1.dmatcsr");
    Write_dmatcsr_csr(mat_read2, "../output/"file_prefix"_read2.dmatcsr");
    /*-------------------------------------------------------------------*/
    dmatcsr *mat_copy = Copy_dmatcsr(mat_read2);
    Write_dmatcsr_csr(mat_copy, "../output/"file_prefix"_copy.dmatcsr");
    /*********************************************************************/
    Free_dmatcsr(mat_copy);
    Free_dmatcsr(mat_read2);
    Free_dmatcsr(mat_read1);
    /*===================================================================*/
    /*===================================================================*/
    printf("3: Transpose...\n");
    dmatcsr *mat_read3 = Read_dmatcsr("../dat/"file_prefix"_A.dmatcsr");
    dmatcsr *mat_trans = malloc(sizeof(dmatcsr));
    Transpose_dmatcsr(mat_read3, mat_trans);
    //Print_dmatcsr    (mat_trans);
    Write_dmatcsr_csr(mat_trans, "../output/"file_prefix"_trans.dmatcsr");
    /*********************************************************************/
    Free_dmatcsr(mat_trans);
    Free_dmatcsr(mat_read3);
    /*===================================================================*/
    /*===================================================================*/
    printf("4: Multiply...\n");
    dmatcsr *A_multi = Read_dmatcsr("../dat/"file_prefix"_A.dmatcsr");
    dmatcsr *B_multi = Read_dmatcsr("../dat/"file_prefix"_A.dmatcsr");
    dmatcsr *C_multi = malloc(sizeof(dmatcsr));
    double tmb = Get_time();
    Multi_dmatcsr_dmatcsr(A_multi, B_multi, C_multi);
    double tme = Get_time();
    Write_dmatcsr_csr(C_multi, "../output/"file_prefix"_multi.dmatcsr");
    printf("multi: %f\n", tme-tmb);
    /*********************************************************************/
    Free_dmatcsr(C_multi);
    Free_dmatcsr(B_multi);
    Free_dmatcsr(A_multi);
    /*===================================================================*/
    dmatcsr *A_ABC_1  = Read_dmatcsr("../dat/"file_prefix"_A.dmatcsr");
    dmatcsr *A_ABC_2  = Read_dmatcsr("../dat/"file_prefix"_A.dmatcsr");
    dmatcsr *A_ABC_3  = Read_dmatcsr("../dat/"file_prefix"_A.dmatcsr");
    
    dmatcsr *M_ABC_1 = malloc(sizeof(dmatcsr));
    double tabcb1 = Get_time();
    Multi_dmatcsr_dmatcsr_dmatcsr(A_ABC_1, A_ABC_2, A_ABC_3, M_ABC_1);
    double tabce1 = Get_time();
    Write_dmatcsr_csr(M_ABC_1, "../output/"file_prefix"_ABC_1.dmatcsr");
    Write_dmatcsr_bmp(M_ABC_1, "../output/ABC_M_0.bmp", 800, 800, ColorMapMatStruct);
    printf("MULTI ABC: %f\n", tabce1-tabcb1);

    dmatcsr *M_ABC_2 = malloc(sizeof(dmatcsr));
    dmatcsr *M_ABC_t = malloc(sizeof(dmatcsr));
    double tabcb2 = Get_time();
    Multi_dmatcsr_dmatcsr(A_ABC_1, A_ABC_2, M_ABC_t);
    double tabctb = Get_time();
    Multi_dmatcsr_dmatcsr(M_ABC_t, A_ABC_3, M_ABC_2);
    double tabce2 = Get_time();
    Write_dmatcsr_csr(M_ABC_2, "../output/"file_prefix"_ABC_2.dmatcsr");
    printf("MULTI ABC2: %f = %f + %f\n", tabce2-tabcb2, tabctb-tabcb2, tabce2-tabctb); 
    Free_dmatcsr(M_ABC_t);
    Free_dmatcsr(M_ABC_2);
    Free_dmatcsr(M_ABC_1);
    Free_dmatcsr(A_ABC_3);
    Free_dmatcsr(A_ABC_2);
    Free_dmatcsr(A_ABC_1);
    printf("===========================================\n");
    /*===================================================================*/
    /*===================================================================*/
    printf("Remove zero elements...\n");
    dmatcsr *A_zeros = Read_dmatcsr("../dat/matrix_data_continuous_coef/gmg_A_refine5.m");
    Print_dmatcsr(A_zeros);
    Remove_zero_dmatcsr(A_zeros);
    Print_dmatcsr(A_zeros);
    Write_dmatcsr_csr(A_zeros, "../output/Azeros.dmatcsr");
    Free_dmatcsr(A_zeros);
#endif

    dmatcsr *A_fem = Read_dmatcsr("../dat/fem2d_poisson_lshape/gmg_A_refine4.m");
    dmatcsr *M_fem = Read_dmatcsr("../dat/fem2d_poisson_lshape/gmg_M_refine4.m");
    dmatcsr *C_fem = Sum_dmatcsr_mApnB(2.0, A_fem, -3.0, M_fem);
    Write_dmatcsr_csr(C_fem, "../output/AminusB.dmatcsr");
    Free_dmatcsr(A_fem);
    Free_dmatcsr(M_fem);
    Free_dmatcsr(C_fem);

    return 0;
}
