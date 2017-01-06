#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "preprocess.h"
#include "matrix.h"
#include "io.h"
#include "linear_algebra.h"
#include "tool.h"
#include "par_matrix_vector.h"
#include "par_cfsplit.h"
#include "par_linear_algebra.h"

#include "mpi.h"

int print_rank = 5;

int my_CLJP_split(imatcsr *S);

int print_num1 = 0;
int print_num2 = 0;

int main(int argc, char *argv[])
{
    amg_param param;
    Init_amg_param(&param);

    MPI_Init(&argc, &argv);

    int  myrank;
    //int nproc_global, myname_len;
    //char myname[MPI_MAX_PROCESSOR_NAME];
    //MPI_Comm_size(MPI_COMM_WORLD, &nproc_global);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    //MPI_Get_processor_name(myname, &myname_len);

    //char file[256] = "../../dat/fem3d/hydrogen-stiff-4913.dat";
    //char file[256] = "../../dat/fem2d_poisson_square/gmg_A_refine5.m";
    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine5.m";
    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine4.m";
    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine3.m";
    //char file[256] = "../../dat/fdm2d9pt/A_fdm9pt_6400x6400.dat";
    char file[256] = "../dat/fdm2d9pt/A_fdm9pt_49x49.dat";
    par_dmatcsr *A = Read_par_dmatcsr(file, MPI_COMM_WORLD);

    par_imatcsr *S = (par_imatcsr*)malloc(sizeof(par_imatcsr));
    Generate_par_strong_coupling_set(A, S, param);
    par_ivec *dof = Init_par_ivec_length_comm(A->diag->nr, A->comm);
    Split_par_CLJP(A, S, dof);
    par_dmatcsr *P = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
    Get_par_interpolation_direct(A, S, dof, P);
    par_dmatcsr *AP = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));

    //Print_par_dmatcsr(P, 0);

    MPI_Barrier(MPI_COMM_WORLD);
    double tb_multi_par = MPI_Wtime();
    Multi_par_dmatcsr_dmatcsr(A, P, AP);
    MPI_Barrier(MPI_COMM_WORLD);
    double te_multi_par = MPI_Wtime();
    if(myrank == print_rank) printf("Multi_par A*P time: %f\n", te_multi_par-tb_multi_par);
    Print_par_dmatcsr(AP, 3);
    //Write_par_dmatcsr_csr(AP, "../output/AP.par.dat", 0);
    //Write_par_dmatcsr_csr(P, "../output/P.par.dat", 0);
    Free_par_dmatcsr(AP);
#if 0
    if(myrank == print_rank)
    {
	printf("P(%d:%d,%d:%d)\n", P->row_start[myrank]+1, P->row_start[myrank+1]-1+1, 
		P->col_start[myrank]+1, P->col_start[myrank+1]-1+1);
	Write_dmatcsr_csr(P->diag, "../output/P_diag.dat");
	Write_dmatcsr_csr(P->offd, "../output/P_offd.dat");
	Write_dmatcsr_csr(A->diag, "../output/A_diag.dat");
	Write_dmatcsr_csr(A->offd, "../output/A_offd.dat");
    }
#endif

    print_num1 = A->col_start[myrank];
    print_num2 = A->col_start[myrank+1];

    Free_par_dmatcsr(P);
    Free_par_ivec(dof);
    //Free_imatcsr(S_ext);
    Free_par_imatcsr(S);
    Free_par_dmatcsr(A);

#if 1
    if(myrank == print_rank)
    {
	fflush(stdout);
	//printf("print from %d to %d\n", print_num1, print_num2-1);

	dmatcsr *AA = Read_dmatcsr(file);
	imatcsr *SS = (imatcsr*)malloc(sizeof(imatcsr));
	int *dof = (int*)calloc(AA->nr, sizeof(int));
	Generate_strong_coupling_set(AA, SS, param);
	CLJP_split(AA, SS, dof);
	dmatcsr *PP = (dmatcsr*)malloc(sizeof(dmatcsr));
	Generate_sparsity_P_dir(SS, dof, PP);
	Generate_P_dir(AA, SS, dof, PP);
	//Write_dmatcsr_csr(AA, "../output/A.dat");
	//Write_dmatcsr_csr(PP, "../output/P.dat");
	dmatcsr *AAPP = (dmatcsr*)malloc(sizeof(dmatcsr));

	double tb_multi = Get_time();
	Multi_dmatcsr_dmatcsr(AA, PP, AAPP);
	double te_multi = Get_time();
	printf("Multi A*P time: %f\n", te_multi-tb_multi);

	//Write_dmatcsr_csr(AAPP, "../output/AAPP.dat");
	//Write_imatcsr_csr(SS, "../output/S.dat");
	//my_CLJP_split(SS);
	free(dof);
	Free_dmatcsr(PP);
	Free_imatcsr(SS);
	Free_dmatcsr(AA);
	Free_dmatcsr(AAPP);
    }
#endif

    MPI_Finalize();

    return 0;
}

#if 0
int my_CLJP_split(imatcsr *S)
{
    //printf("CLJP spliting...\n");

    int i;
    
    imatcsr *ST = (imatcsr*)malloc(sizeof(imatcsr));
    Transpose_imatcsr_struct(S, ST);
    
    int  ST_nr = ST->nr;
    int *ST_ia = ST->ia;
    
    srand(1);
    int *lambda_ST = (int*)malloc(ST_nr*sizeof(int));
    for(i=0; i<ST_nr; i++) lambda_ST[i] = ST_ia[i+1]-ST_ia[i];

    //printf("print from %d to %d\n", print_num1, print_num2-1);
    //Print_ivec(lambda_ST+print_num1, print_num2-print_num1);


    free(lambda_ST);
    Free_imatcsr(ST);

    return 1;
}
#endif
