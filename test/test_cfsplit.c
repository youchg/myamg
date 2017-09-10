#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "preprocess.h"
#include "matrix.h"
#include "io.h"
#include "linear_algebra.h"
#include "tool.h"
#include "par_matrix_vector.h"
#include "par_cfsplit.h"
#include "par_linear_algebra.h"

int print_rank = 0;
int print_num1 = 0;
int print_num2 = 0;

#define file_prefix_55 "gmg_A_refine5"
#define file_prefix_59 "gmg_A_refine9"
#define file_prefix_510 "gmg_A_refine10"

#define file_prefix    file_prefix_510

int main()
{
    dmatcsr *A    = Read_dmatcsr("../../dat/fem2d_poisson_square/"file_prefix".m");
    imatcsr *S    = (imatcsr*)malloc(sizeof(imatcsr));
    int     *dof  = (int*)    calloc(A->nr, sizeof(int));
    dmatcsr *P    = (dmatcsr*)malloc(sizeof(dmatcsr));
    dmatcsr *R    = (dmatcsr*)malloc(sizeof(dmatcsr));
    dmatcsr *AH   = (dmatcsr*)malloc(sizeof(dmatcsr));
    dmatcsr *AH2  = (dmatcsr*)malloc(sizeof(dmatcsr));
    int      ncpt = 0;;

    amg_param param;
    Init_amg_param(&param);

    double tb_gen_S = Get_time();
    Generate_strong_coupling_set(A, S, param);
    double te_gen_S = Get_time();

    double tb_pre_split = Get_time();
    ncpt = Pre_split(A, S, dof);
    double te_pre_split = Get_time();
    
    double tb_post_split = Get_time();
    ncpt = Post_split(S, dof);
    double te_post_split = Get_time();

    //CLJP_split(AA, SS, dof);

    double tb_gen_P = Get_time();
    Generate_sparsity_P_dir(S, dof, P);
    Generate_P_dir(A, S, dof, P);
    double te_gen_P = Get_time();

    double tb_gen_R = Get_time();
    Transpose_dmatcsr(P, R);
    double te_gen_R = Get_time();

    //Write_dmatcsr_csr(P, "../output/P.dat");
    //Write_dmatcsr_csr(R, "../output/R.dat");

    double tb_AH = Get_time();
    Multi_dmatcsr_dmatcsr_dmatcsr(R, A, P, AH);
    double te_AH = Get_time();
    //Write_dmatcsr_csr(RAP, "../output/AH.dat");

    dmatcsr *RA = malloc(sizeof(dmatcsr));
    double tb_AH2 = Get_time();
    Multi_dmatcsr_dmatcsr(R,  A, RA);
    Multi_dmatcsr_dmatcsr(RA, P, AH2);
    double te_AH2 = Get_time();
    Free_dmatcsr(RA);

    printf("NFPT      : %d\n", A->nr);
    printf("NCPT      : %d\n", ncpt);
    printf("GEN  S    : %f\n", te_gen_S      - tb_gen_S);
    printf("PRE  SPLIT: %f\n", te_pre_split  - tb_pre_split);
    printf("POST SPLIT: %f\n", te_post_split - tb_post_split);
    printf("GEN  P    : %f\n", te_gen_P      - tb_gen_P);
    printf("GEN  R    : %f\n", te_gen_R      - tb_gen_R);
    printf("RAP       : %f\n", te_AH         - tb_AH);
    printf("RAP2      : %f\n", te_AH2        - tb_AH2);

    free(dof);
    Free_dmatcsr(AH);
    Free_dmatcsr(R);
    Free_dmatcsr(P);
    Free_imatcsr(S);
    Free_dmatcsr(A);

    return 0;
}
