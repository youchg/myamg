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

#include "mpi.h"

int print_rank = 2;

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

    //char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine5.m";
    char file[256] = "../../dat/fem2d_poisson_lshape/gmg_A_refine3.m";
    //char file[256] = "../../dat/fdm2d9pt/A_fdm9pt_6400x6400.dat";
    //char file[256] = "../dat/fdm2d9pt/A_fdm9pt_49x49.dat";
    par_dmatcsr *A = Read_par_dmatcsr(file, MPI_COMM_WORLD);

    par_imatcsr *S = (par_imatcsr*)malloc(sizeof(par_imatcsr));
    Generate_par_strong_coupling_set(A, S, param);
    par_ivec *dof = Init_par_ivec_length_comm(A->diag->nr, A->comm);
    Split_par_CLJP(A, S, dof);
    //imatcsr *S_ext = Get_S_ext(A, S);

    print_num1 = A->col_start[myrank];
    print_num2 = A->col_start[myrank+1];

    Free_par_ivec(dof);
    //Free_imatcsr(S_ext);
    Free_par_imatcsr(S);
    Free_par_dmatcsr(A);

    if(myrank == print_rank)
    {
	//printf("print from %d to %d\n", print_num1, print_num2-1);

	dmatcsr *AA = Read_dmatcsr(file);
	imatcsr *SS = (imatcsr*)malloc(sizeof(imatcsr));
	Generate_strong_coupling_set(AA, SS, param);
	Write_imatcsr_csr(SS, "../output/S.dat");
	my_CLJP_split(SS);
	Free_imatcsr(SS);
	Free_dmatcsr(AA);
    }


    MPI_Finalize();

    return 0;
}

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
