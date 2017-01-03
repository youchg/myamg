#include "par_matrix_vector.h"

#if WITH_MPI

#include "par_linear_algebra.h"
#include "linear_algebra.h"
#include "io.h"
#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

extern int print_rank;

/* 必要的时候，检查 A, x, y 的 comm_info, comm_info 是否一致
 * y = A_diag*x_diag + A_offd*x_recv;
 * 
 * 发送：邻居进程 t 上，如果与本进程相关的某列 j (全局编号) 不全为 0，
 *       则说明需要 x 第 j 行的数据.
 *       反映在本进程上，则对应与进程 t 相关的某行 j (全局编号) 不全为0.
 *       因此选择发送数据的时候，需要依据邻居进程的行信息，
 *       把不全为 0 的行 k (局部编号即可) 对应的 x 上第 k 个数据 (相应的局部编号) 发送出去。
 *
 * 接收：直接接收即可。接收的数据按照邻居进程的顺序组成一个数组，即可直接与 offd 相乘。
 */
void Multi_par_dmatcsr_dvec(par_dmatcsr *A, par_dvec *x, par_dvec *y)
{
    MPI_Comm comm = A->comm;

    int  myrank;
    MPI_Comm_rank(comm, &myrank);

    int i, j;

    int nproc_neighbor = A->comm_info->nproc_neighbor;
    int *proc_neighbor = A->comm_info->proc_neighbor;

    int    **index_send = x->comm_info->index_row;
    int    *nindex_send = x->comm_info->nindex_row;
    double **x_send     = x->send_data;
    for(i=0; i<nproc_neighbor; i++)
    {
	for(j=0; j<nindex_send[i]; j++)
	    x_send[i][j] = x->value[index_send[i][j]];
    }

    int    *start_recv = x->recv_data_start;
    double *x_recv     = x->recv_data;

    //MPI_Sendrecv(void *sendbuf,int sendcount,MPI_Datatype sendtype,int dest,  int sendtag,
    //             void *recvbuf,int recvcount,MPI_Datatype recvtype,int source,int recvtag,
    //             MPI_Comm comm,MPI_Status *status)
    for(i=0; i<nproc_neighbor; i++)
    {
	MPI_Sendrecv(x_send[i],            nindex_send[i],                MPI_DOUBLE, proc_neighbor[i], myrank+proc_neighbor[i]*1000, 
		     x_recv+start_recv[i], start_recv[i+1]-start_recv[i], MPI_DOUBLE, proc_neighbor[i], proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
    }
    
    Multi_dmatcsr_dvec     (A->diag, x->value, y->value);
    Multi_dmatcsr_dvec_hold(A->offd, x_recv,   y->value);
}

double Get_par_dvec_2norm(par_dvec *x)
{
    double mynorm_square = Multi_dvec_dvec(x->value, x->value, x->length);
    double   norm_square = 0.0;
    MPI_Allreduce(&mynorm_square, &norm_square, 1, MPI_DOUBLE, MPI_SUM, x->comm);
    return sqrt(norm_square);
}

/*
 * C_diag = A_diag*B_diag + A_offd*B_ext_diag
 * C_offd = A_diag*B_offd + A_offd*B_ext_offd
 */
void Multi_par_dmatcsr_dmatcsr(par_dmatcsr *A, par_dmatcsr *B, par_dmatcsr *C)
{
    assert(A->nc_global == A->nc_global);
    assert(A->diag->nc  == B->diag->nr);
    assert(A->comm == B->comm);

    int  myrank, nproc_global;
    MPI_Comm comm = A->comm;
    MPI_Comm_size(comm, &nproc_global);
    MPI_Comm_rank(comm, &myrank);

    par_comm_info *A_comm_info = A->comm_info;
    int A_nproc_neighbor = A->comm_info->nproc_neighbor;
    int *A_proc_neighbor = A->comm_info->proc_neighbor;
    int *A_map_offd_col_l2g = A->map_offd_col_l2g;

    int *B_map_offd_col_l2g = B->map_offd_col_l2g;

    int *B_diag_ia = B->diag->ia;
    int *B_diag_ja = B->diag->ja;
    double *B_diag_va = B->diag->va;
    int *B_offd_ia = B->offd->ia;
    int *B_offd_ja = B->offd->ja;
    double *B_offd_va = B->offd->va;

    int i, j, k, m;

    int *nrow_need_by_neighbor = (int*) calloc(A_nproc_neighbor, sizeof(int));
    for(i=0; i<A_nproc_neighbor; i++)
    {
	MPI_Sendrecv(&A_comm_info->nindex_col[i], 1, MPI_INT, A_proc_neighbor[i], myrank+A_proc_neighbor[i]*1000, 
		     &nrow_need_by_neighbor[i],   1, MPI_INT, A_proc_neighbor[i], A_proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
    }

    int **row_need_by_neighbor = (int**)malloc(A_nproc_neighbor * sizeof(int*));
    for(i=0; i<A_nproc_neighbor; i++)
    {
	row_need_by_neighbor[i] = (int*)malloc(nrow_need_by_neighbor[i] * sizeof(int));
	int *row_need_by_self = (int*)malloc(A_comm_info->nindex_col[i] * sizeof(int));
	for(j=0; j<A_comm_info->nindex_col[i]; j++)
	    row_need_by_self[j] = A_map_offd_col_l2g[A_comm_info->index_col[i][j]]; 
	MPI_Sendrecv(row_need_by_self,
		     A_comm_info->nindex_col[i], MPI_INT, A_proc_neighbor[i], myrank+A_proc_neighbor[i]*1000, 
		     row_need_by_neighbor[i],   
		     nrow_need_by_neighbor[i],   MPI_INT, A_proc_neighbor[i], A_proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);

	free(row_need_by_self); row_need_by_self = NULL;
    }

    int *B_ext_recv_ia = (int*)calloc(A->offd->nc+1, sizeof(int));
    int *B_ext_send_ia = (int*)calloc(A->diag->nr,   sizeof(int));
    int count_recv = 0;
    for(i=0; i<A_nproc_neighbor; i++)
    {
	for(j=0; j<nrow_need_by_neighbor[i]; j++)
	{
	    k = row_need_by_neighbor[i][j] - A->row_start[myrank];
	    B_ext_send_ia[j] = B_diag_ia[k+1] - B_diag_ia[k] + B_offd_ia[k+1] - B_offd_ia[k];
	}

	MPI_Sendrecv(B_ext_send_ia,
		     nrow_need_by_neighbor[i],   MPI_INT, A_proc_neighbor[i], myrank+A_proc_neighbor[i]*1000, 
		     B_ext_recv_ia + count_recv + 1, 
		     A_comm_info->nindex_col[i], MPI_INT, A_proc_neighbor[i], A_proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);

	count_recv += A_comm_info->nindex_col[i];
    }
    assert(count_recv == A->offd->nc);

    for(i=0; i<A->offd->nc; i++) B_ext_recv_ia[i+1] += B_ext_recv_ia[i];

    int *B_ext_recv_ja = (int*)malloc(B_ext_recv_ia[A->offd->nc] * sizeof(int));
    int *B_ext_send_ja = (int*)calloc(A->diag->nc+A->offd->nc,   sizeof(int));
    int count_send_ja = 0;
    int count_recv_ja = 0;
    for(i=0; i<A_nproc_neighbor; i++)
    {
	count_send_ja = 0;
	for(j=0; j<nrow_need_by_neighbor[i]; j++)
	{
	    k = row_need_by_neighbor[i][j] - A->row_start[myrank];
	    for(m=B_diag_ia[k]; m<B_diag_ia[k+1]; m++)
		B_ext_send_ja[count_send_ja++] = B_diag_ja[m] + B->col_start[myrank];
	    for(m=B_offd_ia[k]; m<B_offd_ia[k+1]; m++)
		B_ext_send_ja[count_send_ja++] = B_map_offd_col_l2g[B_offd_ja[m]];
	}


	int recv_ja_begin = B_ext_recv_ia[count_recv_ja];
	int recv_ja_end   = B_ext_recv_ia[count_recv_ja+A_comm_info->nindex_col[i]];
	MPI_Sendrecv(B_ext_send_ja,
		     count_send_ja, MPI_INT, A_proc_neighbor[i], myrank+A_proc_neighbor[i]*1000, 
		     B_ext_recv_ja+recv_ja_begin,
		     recv_ja_end-recv_ja_begin, MPI_INT, A_proc_neighbor[i], A_proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);

	count_recv_ja += A_comm_info->nindex_col[i];
    }

    double *B_ext_recv_va = (double*)malloc(B_ext_recv_ia[A->offd->nc] * sizeof(double));
    double *B_ext_send_va = (double*)calloc(A->diag->nc+A->offd->nc,     sizeof(double));
    int count_send_va = 0;
    int count_recv_va = 0;
    for(i=0; i<A_nproc_neighbor; i++)
    {
	count_send_va = 0;
	for(j=0; j<nrow_need_by_neighbor[i]; j++)
	{
	    k = row_need_by_neighbor[i][j] - A->row_start[myrank];
	    for(m=B_diag_ia[k]; m<B_diag_ia[k+1]; m++)
		B_ext_send_va[count_send_va++] = B_diag_va[m];
	    for(m=B_offd_ia[k]; m<B_offd_ia[k+1]; m++)
		B_ext_send_va[count_send_va++] = B_offd_va[m];
	}

	int recv_va_begin = B_ext_recv_ia[count_recv_va];
	int recv_va_end   = B_ext_recv_ia[count_recv_va+A_comm_info->nindex_col[i]];
	MPI_Sendrecv(B_ext_send_va,
		     count_send_va, MPI_DOUBLE, A_proc_neighbor[i], myrank+A_proc_neighbor[i]*1000, 
		     B_ext_recv_va+recv_va_begin,
		     recv_va_end-recv_va_begin, MPI_DOUBLE, A_proc_neighbor[i], A_proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
	count_recv_va += A_comm_info->nindex_col[i];

    }

    dmatcsr *B_ext = (dmatcsr*)malloc(sizeof(dmatcsr));
    B_ext->nr = A->offd->nc;
    B_ext->nc = B->nc_global;
    B_ext->nn = B_ext_recv_ia[A->offd->nc];
    B_ext->ia = B_ext_recv_ia;
    B_ext->ja = B_ext_recv_ja;
    B_ext->va = B_ext_recv_va;

    dmatcsr *C_ext = (dmatcsr*)malloc(sizeof(dmatcsr));
    Multi_dmatcsr_dmatcsr(A->offd, B_ext, C_ext);

    dmatcsr *C_diag1 = (dmatcsr*)malloc(sizeof(dmatcsr));
    Multi_dmatcsr_dmatcsr(A->diag, B->diag, C_diag1);

    dmatcsr *C_offd1 = (dmatcsr*)malloc(sizeof(dmatcsr));
    Multi_dmatcsr_dmatcsr(A->diag, B->offd, C_offd1);
    for(i=0; i<C_offd1->nn; i++)
	C_offd1->ja[i] = B_map_offd_col_l2g[C_offd1->ja[i]];
    C_offd1->nc = B->nc_global;

    dmatcsr *C_offd2 = (dmatcsr*)malloc(sizeof(dmatcsr));
    C_offd2 = Sum_dmatcsr_dmatcsr(C_offd1, C_ext);

    dmatcsr *C_diag2 = (dmatcsr*)malloc(sizeof(dmatcsr));
    dmatcsr *C_offd = (dmatcsr*)malloc(sizeof(dmatcsr));
    Separate_dmatcsr_to_diag_offd(C_offd2, B->col_start[myrank], B->col_start[myrank+1]-1, C_diag2, C_offd);

    Get_par_dmatcsr_diag(C_diag2, B->col_start[myrank], B->col_start[myrank+1]-1);
    dmatcsr *C_diag = (dmatcsr*)malloc(sizeof(dmatcsr));
    C_diag = Sum_dmatcsr_dmatcsr(C_diag1, C_diag2);

    int *C_map_offd_col_l2g = (int*)malloc(C_offd->nn * sizeof(int));
    for(i=0; i<C_offd->nn; i++) C_map_offd_col_l2g[i] = -1;
    Get_par_dmatcsr_offd_and_map_col_offd_l2g(C_offd, C_map_offd_col_l2g);
    C_map_offd_col_l2g = (int*)realloc(C_map_offd_col_l2g, C_offd->nc*sizeof(int));

    assert(C_map_offd_col_l2g != NULL);
    assert(C_offd->nr == C_diag2->nr);

    C->diag = C_diag; 
    C->offd = C_offd;
    C->nr_global = A->nr_global;
    C->nc_global = B->nc_global;
    C->nn_global = -1;
    C->map_offd_col_l2g = C_map_offd_col_l2g;
    C->comm = comm;
    C->comm_info = NULL;
    C->row_start = (int*)calloc(nproc_global+1, sizeof(int));
    for(i=0; i<=nproc_global; i++) C->row_start[i] = A->row_start[i];
    C->col_start = (int*)calloc(nproc_global+1, sizeof(int));
    for(i=0; i<=nproc_global; i++) C->col_start[i] = B->col_start[i];

    free(nrow_need_by_neighbor); nrow_need_by_neighbor = NULL;
    for(i=0; i<A_nproc_neighbor; i++)
    {
	free(row_need_by_neighbor[i]);
	row_need_by_neighbor[i] = NULL;
    }
    free(row_need_by_neighbor);
    row_need_by_neighbor = NULL;
    free(B_ext_send_ia);
    free(B_ext_send_ja);
    free(B_ext_send_va);
    Free_dmatcsr(B_ext);
    Free_dmatcsr(C_ext);
    Free_dmatcsr(C_offd2);
    Free_dmatcsr(C_diag1);
    Free_dmatcsr(C_diag2);
    Free_dmatcsr(C_offd1);
}

#endif
