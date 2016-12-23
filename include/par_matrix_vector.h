#ifndef __PAR_MATRIX_VECTOR_H__
#define __PAR_MATRIX_VECTOR_H__

#define WITH_MPI 1

#include "matrix.h"
#include "mpi.h"

typedef struct PAR_COMMUNICATION_RECV_INFO_
{
    int  nproc;
    int *proc;
    int *start;
} par_comm_recv_info;

typedef struct PAR_COMMUNICATION_SEND_INFO_
{
    int  nproc;
    int *proc;
    int **index;
    int *nindex;
} par_comm_send_info;

typedef struct PAR_COMMUNICATION_DATA_
{
    MPI_Datatype  type;
    void         *data;
} par_comm_data;

/*
typedef struct PAR_COMMUNICATION_RECV_DATA__
{
    par_comm_recv_info *info;
    void *data;
} par_comm_recv_data;

typedef struct PAR_COMMUNICATION_SEND_DATA__
{
    par_comm_send_info *info;
    void **data;
} par_comm_send_data;
*/

/* Assume that A is symmetric */
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
    par_comm_send_info *send_info;
    par_comm_recv_info *recv_info;

    //par_comm_data *send_data;
    //par_comm_data *recv_data;
} par_dmatcsr;

typedef struct INT_PAR_MATRIX_CSR_
{
    imatcsr *diag;
    imatcsr *offd;

    /* global info */
    int nr_global;
    int nc_global;
    int nn_global;

    int *map_offd_col_l2g;
    int *row_start;
    int *col_start;

    /* mpi communication */
    MPI_Comm comm;
    par_comm_send_info *send_info;
    par_comm_recv_info *recv_info;

    //par_comm_data *send_data;
    //par_comm_data *recv_data;
} par_imatcsr;

typedef struct PAR_DOUBLE_VECTOR_
{
    int     length_global;
    int     length;
    double *value;

    MPI_Comm comm;
    par_comm_send_info *send_info;
    par_comm_recv_info *recv_info;
    double **send_data;
    double  *recv_data;
} par_dvec;

par_dmatcsr *Read_par_dmatcsr(const char *filename, MPI_Comm comm);
void Free_par_dmatcsr(par_dmatcsr *A);
void Free_par_imatcsr(par_imatcsr *A);

par_dvec *Init_par_dvec_from_par_dmatcsr(par_dmatcsr *A);
void Free_par_dvec(par_dvec *x);

void Free_par_comm_data      (par_comm_data      *comm_data);
void Free_par_comm_send_info (par_comm_send_info *info);
void Free_par_comm_recv_info (par_comm_recv_info *info);

void Print_par_comm_recv_info(par_comm_recv_info *info);
void Print_par_comm_send_info(par_comm_send_info *info);

par_comm_send_info *Copy_par_comm_send_info (par_comm_send_info *send_info);
par_comm_recv_info *Copy_par_comm_recv_info (par_comm_recv_info *recv_info);

//void Free_par_comm_send_data(par_comm_send_data *send_data);
//void Free_par_comm_recv_data(par_comm_recv_data *recv_data);

#if 0 //WITH_MPI
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

#endif

#endif
#endif
