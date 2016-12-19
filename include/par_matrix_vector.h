#ifndef __PAR_MATRIX_VECTOR_H__
#define __PAR_MATRIX_VECTOR_H__

#define WITH_MPI 1

#include "matrix.h"
#include "mpi.h"

typedef struct PAR_COMMUNICATION_INFO__
{
    int  nproc;
    int *proc;
    int *start;
} par_comm_info;

typedef struct PAR_COMMUNICATION_DATA__
{
    MPI_Datatype  type;
    void         *data;
} par_comm_data;

typedef struct DOUBLE_PAR_MATRIX_CSR_
{
    /* local info */
    //int nr;
    //int nc;
    //int nn;

    dmatcsr *diag;
    dmatcsr *offd;

    /* global info */
    int nr_global;
    int nc_global;
    int nn_global;

    int *map_offd_col_l2g;
    int *row_start;
    int *col_start;

    /* mpi communication */
    MPI_Comm comm;
    par_comm_info *send_info;
    par_comm_data *send_data;
    par_comm_info *recv_info;
    par_comm_data *recv_data;
} par_dmatcsr;



#if 0 //WITH_MPI
typedef struct PAR_DOUBLE_VECTOR_
{
    MPI_Comm comm;

    int length_golbal;

    int length;
    double *value;
    par_comm_info *pcomm_info;
} par_dvec;

double Get_par_dvec_2norm(par_dvec *x);
par_dmatcsr *Read_par_dmatcsr(const char *filename, MPI_Comm comm);
void Free_par_dmatcsr(par_dmatcsr *A);
void Free_par_dvec(par_dvec *x);
void Free_par_comm_info(par_comm_info *pcomm_info);

#if 0
typedef struct INT_PAR_MATRIX_CSR_
{
    int nr;
    int nc;
    int nn;
    
    int *ia;                                 
    int *ja;
    int *va;
} par_imatcsr;

dmatcsr *Create_dmatcsr(const int nrow, const int ncol, const int nnz);
imatcsr *Create_imatcsr(const int nrow, const int ncol, const int nnz);

void Free_dmatcsr(dmatcsr *A);
void Free_imatcsr(imatcsr *A);

dmatcsr *Copy_dmatcsr(dmatcsr *A);
#endif

#endif

//par_dmatcsr *Read_par_dmatcsr(const char *filename, MPI_Comm comm);
void Read_par_dmatcsr(const char *filename, MPI_Comm comm);
void Separate_dmatcsr_to_diag_offd(dmatcsr *A, int col_idx_min, int col_idx_max, dmatcsr *diag, dmatcsr *offd);

void Free_par_dmatcsr(par_dmatcsr *A);
void Free_par_comm_data(par_comm_data *comm_data);
void Free_par_comm_info(par_comm_info *info);

void Print_par_comm_info(par_comm_info *info);
#endif
