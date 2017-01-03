#ifndef __PAR_MATRIX_VECTOR_H__
#define __PAR_MATRIX_VECTOR_H__

#define WITH_MPI 1

#include "matrix.h"
#include "mpi.h"

typedef struct PAR_COMMUNICATION_INFO_
{
    int nproc_neighbor; //邻居进程个数
    int *proc_neighbor; //邻居进程编号

    int *nindex_row; //每个邻居进程包含非零行的个数
    int **index_row; //每个邻居进程包含非零行的编号

    int *nindex_col; //每个邻居进程包含非零列的个数
    int **index_col; //每个邻居进程包含非零列的编号
} par_comm_info;

/* Assume that A is symmetric */
typedef struct DOUBLE_PAR_MATRIX_CSR_
{
    dmatcsr *diag;
    dmatcsr *offd;

    int nr_global;
    int nc_global;
    int nn_global;

    int *map_offd_col_l2g;
    int *row_start;
    int *col_start;

    MPI_Comm comm;
    par_comm_info *comm_info;
} par_dmatcsr;

typedef struct INT_PAR_MATRIX_CSR_
{
    imatcsr *diag;
    imatcsr *offd;

    int nr_global;
    int nc_global;
    int nn_global;

    int *map_offd_col_l2g;
    int *row_start;
    int *col_start;

    MPI_Comm comm;
    par_comm_info *comm_info;
} par_imatcsr;

typedef struct PAR_DOUBLE_VECTOR_
{
    int     length_global;
    int     length;
    double *value;

    MPI_Comm comm;
    par_comm_info *comm_info;
    double **send_data;
    double  *recv_data;
    int     *recv_data_start;
} par_dvec;

typedef struct PAR_INT_VECTOR_
{
    int length_global;
    int length;
    int *value;

    MPI_Comm comm;
    par_comm_info *comm_info;
    double **send_data;
    double  *recv_data;
    int     *recv_data_start;
} par_ivec;

par_dmatcsr *Read_par_dmatcsr(const char *filename, MPI_Comm comm);
void Free_par_dmatcsr(par_dmatcsr *A);
void Free_par_imatcsr(par_imatcsr *A);
void Free_par_comm_info(par_comm_info *info);

//for multiply Ax
par_dvec *Init_par_dvec_mv(par_dmatcsr *A);
par_ivec *Init_par_ivec_length_comm(int length, MPI_Comm comm);
void Free_par_dvec(par_dvec *x);
void Free_par_ivec(par_ivec *x);


void Print_par_comm_info(par_comm_info *info);

par_comm_info *Copy_par_comm_info (par_comm_info *info);



void Get_dmatcsr_global_size(const char *filename, int *nr, int *nc, int *nn);
void Get_par_dmatcsr_row_start(int nr_global, int nprocs, int *row_start);
void Separate_dmatcsr_to_diag_offd(dmatcsr *A, int col_idx_min, int col_idx_max, dmatcsr *diag, dmatcsr *offd);
void Get_par_dmatcsr_diag(dmatcsr *diag, int col_idx_min, int col_idx_max);
void Get_par_dmatcsr_offd_and_map_col_offd_l2g(dmatcsr*offd, int *map_offd_col_l2g);
void Get_par_dmatcsr_comm_info(int nproc_global, dmatcsr *offd, int *col_start, int *map_offd_col_l2g, par_comm_info *comm_info);
void Get_neighbor_proc(int nproc_global, dmatcsr *offd, int *col_start, int *map_offd_col_l2g, int *nproc_neighbor, int *proc_neighbor);
void Get_par_dmatcsr_comm_row_info(dmatcsr *offd, int nproc_neighbor, int *proc_neighbor, int *col_start, int *map_offd_col_l2g, int *nidx, int **idx);
void Get_par_dmatcsr_comm_col_info(dmatcsr *offd, int nproc_neighbor, int *proc_neighbor, int *col_start, int *map_offd_col_l2g, int *nidx, int **idx);

void Write_par_dmatcsr_csr(par_dmatcsr *A, const char *filename, int nametype);
#endif
