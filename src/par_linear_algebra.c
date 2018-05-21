#ifdef WITH_MPI

#include "par_matrix_vector.h"

#include "par_linear_algebra.h"
#include "linear_algebra.h"
#include "io.h"
#include "mpi.h"
#include "tool.h"
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
    int  myrank;
    MPI_Comm comm = A->comm;
    MPI_Comm_rank(comm, &myrank);
    assert(MPI_COMM_NULL != comm);

    int i, j;

    int nproc_neighbor = A->comm_info->nproc_neighbor;
    int *proc_neighbor = A->comm_info->proc_neighbor;
    int *nindex_recv   = A->comm_info->nindex_col;
    int **index_recv   = A->comm_info->index_col;
    int *map_offd_col_l2g = A->map_offd_col_l2g;

    int *nsend = (int*)malloc(nproc_neighbor * sizeof(int));
    for(i=0; i<nproc_neighbor; i++)
    {
	MPI_Sendrecv(nindex_recv+i, 1, MPI_INT, proc_neighbor[i], myrank+proc_neighbor[i]*1000, 
		     nsend+i,       1, MPI_INT, proc_neighbor[i], proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
	//printf("%d to %d: send = %d, recv = %d\n", myrank, proc_neighbor[i], nindex_recv[i], nsend[i]);
    }

    int **col_send = (int**)malloc(nproc_neighbor * sizeof(int*));
    for(i=0; i<nproc_neighbor; i++) col_send[i] = (int*)malloc(nsend[i] * sizeof(int));

    for(i=0; i<nproc_neighbor; i++)
    {
	int offset;
	if(nindex_recv[i] == 0) offset = 0;
	else                    offset = index_recv[i][0];
	MPI_Sendrecv(map_offd_col_l2g+offset, nindex_recv[i], MPI_INT, proc_neighbor[i], myrank+proc_neighbor[i]*1000, 
		     col_send[i],             nsend[i],       MPI_INT, proc_neighbor[i], proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
	//printf("%d to %d: send = %d, recv = %d\n", myrank, proc_neighbor[i], nindex_recv[i], nsend[i]);
    }
    //MPI_Barrier(comm);

    for(i=0; i<nproc_neighbor; i++)
    {
	for(j=0; j<nsend[i]; j++)
	{
	    assert((col_send[i][j]>=A->col_start[myrank]) && (col_send[i][j]<A->col_start[myrank+1]));
	    col_send[i][j] -= A->col_start[myrank];
	}
    }

    double **x_send = (double**)malloc(nproc_neighbor * sizeof(double*));
    for(i=0; i<nproc_neighbor; i++) x_send[i] = (double*)calloc(nsend[i], sizeof(double));
    for(i=0; i<nproc_neighbor; i++)
	for(j=0; j<nsend[i]; j++)
	    x_send[i][j] = x->value[col_send[i][j]];

    int nrecv_data = 0;
    for(i=0; i<nproc_neighbor; i++) nrecv_data += nindex_recv[i];
    double *x_recv = (double*)calloc(nrecv_data, sizeof(double));

    int *start_recv = (int*)calloc(nproc_neighbor+1, sizeof(int));
    for(i=1; i<=nproc_neighbor; i++) start_recv[i]  = nindex_recv[i-1];
    for(i=1; i<=nproc_neighbor; i++) start_recv[i] +=  start_recv[i-1];

    //MPI_Sendrecv(void *sendbuf,int sendcount,MPI_Datatype sendtype,int dest,  int sendtag,
    //             void *recvbuf,int recvcount,MPI_Datatype recvtype,int source,int recvtag,
    //             MPI_Comm comm,MPI_Status *status)
    for(i=0; i<nproc_neighbor; i++)
    {
	//printf("%d to %d: send = %d, recv = %d\n", myrank, proc_neighbor[i], nsend[i], start_recv[i+1]-start_recv[i]);
	MPI_Sendrecv(x_send[i],            nsend[i],                MPI_DOUBLE, proc_neighbor[i], myrank+proc_neighbor[i]*1000, 
		     x_recv+start_recv[i], start_recv[i+1]-start_recv[i], MPI_DOUBLE, proc_neighbor[i], proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
    }

    if(A->diag->nc > 0)
	Multi_dmatcsr_dvec     (A->diag, x->value, y->value);
    if(A->offd->nc > 0)
	Multi_dmatcsr_dvec_hold(A->offd, x_recv,   y->value);

    for(i=0; i<nproc_neighbor; i++)
    {
	free(col_send[i]);
	col_send[i] = NULL;
    }
    free(col_send); col_send = NULL;
    free(nsend); nsend = NULL;
    free(start_recv); start_recv = NULL;
    free(x_recv); x_recv = NULL;
    for(i=0; i<nproc_neighbor; i++)
    {
	free(x_send[i]);
	x_send[i] = NULL;
    }
    free(x_send); x_send = NULL;
}

double Get_par_dvec_2norm(par_dvec *x)
{
#if 0
    int  myrank, nproc;
    MPI_Comm comm = x->comm;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nproc);
#endif
    double mynorm_square = Multi_dvec_dvec(x->value, x->value, x->length);
#if 0
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    printf("rank %d in %d: mynorm_square = %15.12f\n", myrank, nproc, mynorm_square);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
#endif
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

    //par_comm_info *B_comm_info = B->comm_info;
    int B_nproc_neighbor = B->comm_info->nproc_neighbor;
    int *B_proc_neighbor = B->comm_info->proc_neighbor;
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
    //if(myrank == print_rank) Print_ivec(nrow_need_by_neighbor, A_nproc_neighbor);

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
    int *B_ext_send_ia = (int*)calloc(B->diag->nr,   sizeof(int));
    int count_recv = 0;
    for(i=0; i<A_nproc_neighbor; i++)
    {
	for(j=0; j<nrow_need_by_neighbor[i]; j++)
	{
	    k = row_need_by_neighbor[i][j] - B->row_start[myrank];
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
    //按 nn 申请的空间可能有些浪费
    int *B_ext_send_ja = (int*)calloc(B->diag->nn+B->offd->nn,   sizeof(int));
    int count_send_ja = 0;
    int count_recv_ja = 0;
    for(i=0; i<A_nproc_neighbor; i++)
    {
	count_send_ja = 0;
	for(j=0; j<nrow_need_by_neighbor[i]; j++)
	{
	    k = row_need_by_neighbor[i][j] - B->row_start[myrank];
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
    double *B_ext_send_va = (double*)calloc(B->diag->nn+B->offd->nn,     sizeof(double));
    int count_send_va = 0;
    int count_recv_va = 0;
    for(i=0; i<A_nproc_neighbor; i++)
    {
	count_send_va = 0;
	for(j=0; j<nrow_need_by_neighbor[i]; j++)
	{
	    k = row_need_by_neighbor[i][j] - B->row_start[myrank];
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

    dmatcsr *C_diag2 = Init_empty_dmatcsr(C_offd2->nr);
    dmatcsr *C_offd  = Init_empty_dmatcsr(C_offd2->nr);
    Separate_dmatcsr_to_diag_offd(C_offd2, B->col_start[myrank], B->col_start[myrank+1]-1, C_diag2, C_offd);

    Get_par_dmatcsr_diag(C_diag2, B->col_start[myrank], B->col_start[myrank+1]-1);
    dmatcsr *C_diag = (dmatcsr*)malloc(sizeof(dmatcsr));
    C_diag = Sum_dmatcsr_dmatcsr(C_diag1, C_diag2);

    int *C_map_offd_col_l2g = (int*)malloc(C_offd->nn * sizeof(int));
    for(i=0; i<C_offd->nn; i++) C_map_offd_col_l2g[i] = -1;
    Get_par_dmatcsr_offd_and_map_col_offd_l2g(C_offd, C_map_offd_col_l2g);
    C_map_offd_col_l2g = (int*)realloc(C_map_offd_col_l2g, C_offd->nc*sizeof(int));

    if(0 != C_offd->nc) assert(C_map_offd_col_l2g != NULL);
    assert(C_offd->nr == C_diag2->nr);

    C->diag = C_diag; 
    C->offd = C_offd;
    C->nr_global = A->nr_global;
    C->nc_global = B->nc_global;
    C->nn_global = -1;
    C->map_offd_col_l2g = C_map_offd_col_l2g;
    C->comm = comm;

    C->row_start = (int*)calloc(nproc_global+1, sizeof(int));
    for(i=0; i<=nproc_global; i++) C->row_start[i] = A->row_start[i];
    C->col_start = (int*)calloc(nproc_global+1, sizeof(int));
    for(i=0; i<=nproc_global; i++) C->col_start[i] = B->col_start[i];

    //if(myrank == print_rank) Print_ivec(C->row_start, nproc_global+1);
    //if(myrank == print_rank) Print_ivec(C->col_start, nproc_global+1);


    C->comm_info = Init_par_comm_info();
    par_comm_info *C_comm_info = C->comm_info;
    //C_comm_info->nproc_neighbor = 0;
    
    //if(A_nproc_neighbor > 0)
    //{
    //找到A的邻居进程中B需要多少邻居进程
    int *B_nproc_neighbor_neighbor = (int*)calloc(A_nproc_neighbor, sizeof(int));
    for(i=0; i<A_nproc_neighbor; i++)
    {
	MPI_Sendrecv(&B_nproc_neighbor,           1, MPI_INT, A_proc_neighbor[i], myrank+A_proc_neighbor[i]*1000, 
		     B_nproc_neighbor_neighbor+i, 1, MPI_INT, A_proc_neighbor[i], A_proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
    }
    //找到B的邻居进程中A需要多少邻居进程
    int *A_nproc_neighbor_neighbor = (int*)calloc(B_nproc_neighbor, sizeof(int));
    for(i=0; i<B_nproc_neighbor; i++)
    {
	MPI_Sendrecv(&A_nproc_neighbor,           1, MPI_INT, B_proc_neighbor[i], myrank+B_proc_neighbor[i]*1000, 
		     A_nproc_neighbor_neighbor+i, 1, MPI_INT, B_proc_neighbor[i], B_proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
    }

    //分配空间，找到A的邻居进程中B需要的邻居进程号
    int C_nproc_neighbor = 0;
    for(i=0; i<A_nproc_neighbor; i++)
	C_nproc_neighbor += B_nproc_neighbor_neighbor[i];
    for(i=0; i<B_nproc_neighbor; i++)
	C_nproc_neighbor += A_nproc_neighbor_neighbor[i];
    C_nproc_neighbor += A_nproc_neighbor;
    /*
     * 目前不确定当 0==C_nproc_neighbor 时是否会出现问题，
     * 暂时先防止这一点。
     */
    //assert(0 != C_nproc_neighbor);
    int *C_proc_neighbor  = (int*)malloc(C_nproc_neighbor * sizeof(int));
    for(i=0; i<A_nproc_neighbor; i++) C_proc_neighbor[i] = A_proc_neighbor[i];

    C_nproc_neighbor = A_nproc_neighbor;
    for(i=0; i<A_nproc_neighbor; i++)
    {
	MPI_Sendrecv(B_proc_neighbor,           
		     B_nproc_neighbor,             MPI_INT, A_proc_neighbor[i], myrank+A_proc_neighbor[i]*1000, 
		     C_proc_neighbor + C_nproc_neighbor, 
		     B_nproc_neighbor_neighbor[i], MPI_INT, A_proc_neighbor[i], A_proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
	C_nproc_neighbor += B_nproc_neighbor_neighbor[i];
    }
    for(i=0; i<B_nproc_neighbor; i++)
    {
	MPI_Sendrecv(A_proc_neighbor,           
		     A_nproc_neighbor,             MPI_INT, B_proc_neighbor[i], myrank+B_proc_neighbor[i]*1000, 
		     C_proc_neighbor + C_nproc_neighbor, 
		     A_nproc_neighbor_neighbor[i], MPI_INT, B_proc_neighbor[i], B_proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
	C_nproc_neighbor += A_nproc_neighbor_neighbor[i];
    }

    //将目前 C_proc_neighbor 中重复的或等于本进程号的找出并去掉
    //if(myrank == print_rank) Print_ivec(C_proc_neighbor, C_nproc_neighbor);
    Insertion_ascend_sort_ivec(C_proc_neighbor, 0, C_nproc_neighbor-1);
    int C_nproc_neighbor_actual = 0;
    if(C_nproc_neighbor > 0)
    {
	if(C_proc_neighbor[0] != myrank)
	    C_nproc_neighbor_actual = 1;
	for(i=1; i<C_nproc_neighbor; i++)
	{
	    if((C_proc_neighbor[i] != C_proc_neighbor[i-1]) && (C_proc_neighbor[i] != myrank))
		C_nproc_neighbor_actual++;
	}
	C_comm_info->proc_neighbor  = (int*)malloc(C_nproc_neighbor_actual * sizeof(int));
	C_comm_info->nproc_neighbor = 0;
	
	if(C_proc_neighbor[0] != myrank)
	{
	    C_comm_info->nproc_neighbor = 1;
	    C_comm_info->proc_neighbor[0] = C_proc_neighbor[0];
	}
	for(i=1; i<C_nproc_neighbor; i++)
	{
	    if((C_proc_neighbor[i] != C_proc_neighbor[i-1]) && (C_proc_neighbor[i] != myrank))
		C_comm_info->proc_neighbor[C_comm_info->nproc_neighbor++] = C_proc_neighbor[i];
	}
	
	//if(myrank == print_rank) Print_ivec(C_comm_info->proc_neighbor, C_comm_info->nproc_neighbor);

	C->comm_info->nindex_row = (int*) calloc(C->comm_info->nproc_neighbor,  sizeof(int));
	C->comm_info->index_row  = (int**)malloc(C->comm_info->nproc_neighbor * sizeof(int*));
	for(i=0; i<C->comm_info->nproc_neighbor; i++)
	{
	    C->comm_info->index_row[i] = (int*)malloc(C->offd->nr * sizeof(int));
	    for(j=0; j<C->offd->nr; j++) C->comm_info->index_row[i][j] = -1;
	}
	Get_par_dmatcsr_comm_row_info(C->offd, 
				      C->comm_info->nproc_neighbor, C->comm_info->proc_neighbor, 
				      C->col_start,                 C->map_offd_col_l2g, 
				      C->comm_info->nindex_row,     C->comm_info->index_row);
	for(i=0; i<C->comm_info->nproc_neighbor; i++)
	{
	    C->comm_info->index_row[i] = (int*)realloc(C->comm_info->index_row[i], C->comm_info->nindex_row[i]*sizeof(int));
	    if(0 != C->comm_info->nindex_row[i]) assert(NULL != C->comm_info->index_row[i]);
	}

	C->comm_info->nindex_col = (int*) calloc(C->comm_info->nproc_neighbor,  sizeof(int));
	C->comm_info->index_col  = (int**)malloc(C->comm_info->nproc_neighbor * sizeof(int*));
	for(i=0; i<C->comm_info->nproc_neighbor; i++)
	{
	    C->comm_info->index_col[i] = (int*)malloc(C->offd->nc * sizeof(int));
	    for(j=0; j<C->offd->nc; j++) C->comm_info->index_col[i][j] = -1;
	}
	Get_par_dmatcsr_comm_col_info(C->offd, 
				      C->comm_info->nproc_neighbor, C->comm_info->proc_neighbor, 
				      C->col_start,                 C->map_offd_col_l2g, 
				      C->comm_info->nindex_col,     C->comm_info->index_col);
	for(i=0; i<C->comm_info->nproc_neighbor; i++)
	{
	    C->comm_info->index_col[i] = (int*)realloc(C->comm_info->index_col[i], C->comm_info->nindex_col[i]*sizeof(int));
	    if(0 != C->comm_info->nindex_col[i]) assert(NULL != C->comm_info->index_col[i]);
	}

	//可作为 par_dmatcsr->comm_info 的检测工具之一
	for(i=0; i<C->comm_info->nproc_neighbor; i++)
	{
	    if(0 == C->comm_info->nindex_col[i])
		assert(0 == C->comm_info->nindex_row[i]);
	    if(0 != C->comm_info->nindex_col[i])
		assert(0 != C->comm_info->nindex_row[i]);
	}
    }

    free(C_proc_neighbor); C_proc_neighbor = NULL;

    free(A_nproc_neighbor_neighbor); A_nproc_neighbor_neighbor  = NULL;
    free(B_nproc_neighbor_neighbor); B_nproc_neighbor_neighbor  = NULL;

    //}

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

/*
 * 将 diag 转置
 * 将offd按邻居进程转置，然后发送到邻居进程
 */
void Transpose_par_dmatcsr(par_dmatcsr *A, par_dmatcsr *AT)
{
    int  myrank, nproc_global;
    MPI_Comm comm = A->comm;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nproc_global);

    int i, j, k, m;

    AT->nr_global = A->nc_global;
    AT->nc_global = A->nr_global;
    AT->nn_global = A->nn_global;

    AT->row_start = (int*)calloc(nproc_global+1, sizeof(int));
    for(i=0; i<=nproc_global; i++) AT->row_start[i] = A->col_start[i];
    AT->col_start = (int*)calloc(nproc_global+1, sizeof(int));
    for(i=0; i<=nproc_global; i++) AT->col_start[i] = A->row_start[i];

    dmatcsr *AT_diag = (dmatcsr*)malloc(sizeof(dmatcsr));
    Transpose_dmatcsr(A->diag, AT_diag);
    AT->diag = AT_diag;
    //if(myrank == 1) Write_dmatcsr_csr(AT->diag, "../output/Rdiag.dat");

    dmatcsr *AT_offd = (dmatcsr*)malloc(sizeof(dmatcsr));
    AT_offd->nr = AT->diag->nr;
    AT_offd->ia = (int*)calloc(AT_offd->nr+1, sizeof(int));

    AT->offd = AT_offd;


    /* 将 A->offd 转置为 A_offdT
     * 注意 A_offdT 的列号为 A_offd 的局部行号，行号为A_offd的局部列号
     * 在对A_offdT的数据进行发送时，应当将行号和列号的全局编号信息也发送出去
     */
    dmatcsr *A_offdT = (dmatcsr*)malloc(sizeof(dmatcsr));
    Transpose_dmatcsr(A->offd, A_offdT);
    //将A_offd的列号转化为全局编号，这样就不用再单独发送全局编号了
    //而行号的全局编号需要单独发送一次
    for(i=0; i<A_offdT->nn; i++)
	A_offdT->ja[i] += A->row_start[myrank];

    par_comm_info *A_comm_info = A->comm_info;
    int   nproc_neighbor = A_comm_info->nproc_neighbor;
    int *A_proc_neighbor = A_comm_info->proc_neighbor;
    int *A_map_offd_col_l2g = A->map_offd_col_l2g;

    int send_nr_nc_nn[3];
    int send_nr_mark = 0;
    int *recv_nr_nc_nn = (int*)calloc(nproc_neighbor*3, sizeof(int));
    for(i=0; i<nproc_neighbor; i++)
    {
	send_nr_nc_nn[0] = A_comm_info->nindex_col[i];
	send_nr_nc_nn[1] = A_comm_info->nindex_row[i];
	send_nr_nc_nn[2] = A_offdT->ia[send_nr_mark+send_nr_nc_nn[0]] - A_offdT->ia[send_nr_mark];

	MPI_Sendrecv(send_nr_nc_nn,     3, MPI_INT, A_proc_neighbor[i], myrank+A_proc_neighbor[i]*1000, 
		     recv_nr_nc_nn+3*i, 3, MPI_INT, A_proc_neighbor[i], A_proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);

	send_nr_mark += send_nr_nc_nn[0];
    }

    //if(myrank == print_rank) Print_ivec(recv_nr_nc_nn, nproc_neighbor*3);

    int AT_offd_nc = 0;
    int AT_offd_nn = 0;
    for(i=0; i<nproc_neighbor; i++)
    {
	AT_offd_nc += recv_nr_nc_nn[3*i+1];
	AT_offd_nn += recv_nr_nc_nn[3*i+2];
    }

    int recv_nr_total = 0;
    for(i=0; i<nproc_neighbor; i++)
	recv_nr_total += recv_nr_nc_nn[3*i];

    int *send_row_index_nc = (int*)malloc(A_offdT->nr*2 * sizeof(int));
    int *recv_row_index_nc = (int*)malloc(recv_nr_total*2 * sizeof(int));
    int recv_row_mark = 0;
    for(i=0; i<nproc_neighbor; i++)
    {
	for(j=0; j<A_comm_info->nindex_col[i]; j++)
	{
	    send_row_index_nc[2*j  ] = A_map_offd_col_l2g[A_comm_info->index_col[i][j]];
	    send_row_index_nc[2*j+1] = A_offdT->ia[A_comm_info->index_col[i][j]+1]-A_offdT->ia[A_comm_info->index_col[i][j]];;
	}

	MPI_Sendrecv(send_row_index_nc, 
		     A_comm_info->nindex_col[i]*2, MPI_INT, A_proc_neighbor[i], myrank+A_proc_neighbor[i]*1000, 
		     recv_row_index_nc+recv_row_mark*2, 
		     recv_nr_nc_nn[3*i]*2,         MPI_INT, A_proc_neighbor[i], A_proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);

	recv_row_mark += recv_nr_nc_nn[3*i];
    }
    //if(myrank == print_rank) Print_ivec(recv_row_index_nc, recv_nr_total*2);

    for(i=0; i<recv_nr_total; i++)
    {
	assert((recv_row_index_nc[2*i]>=AT->row_start[myrank]) &&
	       (recv_row_index_nc[2*i]< AT->row_start[myrank+1]));
    }
    for(i=0; i<recv_nr_total; i++)
	recv_row_index_nc[2*i] -= AT->row_start[myrank];

    for(i=0; i<recv_nr_total; i++)
	AT_offd->ia[recv_row_index_nc[2*i]+1] += recv_row_index_nc[2*i+1];

    for(i=0; i<AT_offd->nr; i++) AT_offd->ia[i+1] += AT_offd->ia[i];

    AT_offd->nc = AT_offd_nc;
    AT_offd->nn = AT_offd_nn;
    AT_offd->ja = (int*)   malloc(AT_offd_nn * sizeof(int));
    AT_offd->va = (double*)malloc(AT_offd_nn *  sizeof(double));

    int send_row_mark = 0;
    int *count_ja = (int*)calloc(AT_offd->nr, sizeof(int));
    int *recv_col_index = (int*)malloc(AT_offd->nn * sizeof(int));
    recv_row_mark = 0;
    for(i=0; i<nproc_neighbor; i++)
    {
	MPI_Sendrecv(A_offdT->ja + A_offdT->ia[send_row_mark], 
		     A_offdT->ia[send_row_mark+A_comm_info->nindex_col[i]] - A_offdT->ia[send_row_mark], 
		     MPI_INT, A_proc_neighbor[i], myrank+A_proc_neighbor[i]*1000, 
		     recv_col_index, recv_nr_nc_nn[3*i+2],
		     MPI_INT, A_proc_neighbor[i], A_proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);

	send_row_mark += A_comm_info->nindex_col[i];

	int t = 0;
	for(j=0; j<recv_nr_nc_nn[3*i]; j++)
	{
	    for(k=0; k<recv_row_index_nc[2*recv_row_mark+2*j+1]; k++)
	    {
		m = recv_row_index_nc[2*recv_row_mark+2*j];
		AT_offd->ja[AT_offd->ia[m]+count_ja[m]++] = recv_col_index[t++];
		//if(myrank == print_rank)printf("row %d, %d, %d\n", m, count_ja[m], recv_col_index[t-1]);
	    }
	}
	//if(myrank == print_rank)printf("....................\n");
	assert(t == recv_nr_nc_nn[3*i+2]);
	recv_row_mark += recv_nr_nc_nn[3*i];
    }

    send_row_mark = 0;
    int    *count_va = (int*)   calloc(AT_offd->nr,  sizeof(int));
    double *recv_nnz = (double*)malloc(AT_offd->nn * sizeof(double));
    recv_row_mark = 0;
    for(i=0; i<nproc_neighbor; i++)
    {
	MPI_Sendrecv(A_offdT->va + A_offdT->ia[send_row_mark], 
		     A_offdT->ia[send_row_mark+A_comm_info->nindex_col[i]] - A_offdT->ia[send_row_mark], 
		     MPI_DOUBLE, A_proc_neighbor[i], myrank+A_proc_neighbor[i]*1000, 
		     recv_nnz, recv_nr_nc_nn[3*i+2],
		     MPI_DOUBLE, A_proc_neighbor[i], A_proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);

	send_row_mark += A_comm_info->nindex_col[i];

	int t = 0;
	for(j=0; j<recv_nr_nc_nn[3*i]; j++)
	{
	    for(k=0; k<recv_row_index_nc[2*recv_row_mark+2*j+1]; k++)
	    {
		m = recv_row_index_nc[2*recv_row_mark+2*j];
		AT_offd->va[AT_offd->ia[m]+count_va[m]++] = recv_nnz[t++];
		//if(myrank == print_rank)printf("row %d, %d, %d\n", m, count_ja[m], recv_col_index[t-1]);
	    }
	}
	//if(myrank == print_rank)printf("....................\n");
	assert(t == recv_nr_nc_nn[3*i+2]);
	recv_row_mark += recv_nr_nc_nn[3*i];
    }
    //if(myrank == print_rank) Write_dmatcsr_csr(AT_offd, "../output/AT_offd.dat");

    int *AT_map_offd_col_l2g = (int*)malloc(AT_offd->nn * sizeof(int));
    for(i=0; i<AT_offd->nn; i++) AT_map_offd_col_l2g[i] = -1;
    //下面一行是为了Get_par_dmatcsr_offd_and_map_col_offd_l2g()的需要，后面对该函数进行修改
    AT_offd->nc = AT->nc_global;
    Get_par_dmatcsr_offd_and_map_col_offd_l2g(AT_offd, AT_map_offd_col_l2g);
    //if(AT_offd->nc == 0) printf("rank %d\n", myrank);
    //assert(AT_offd->nc > 0);
    AT_map_offd_col_l2g = (int*)realloc(AT_map_offd_col_l2g, AT_offd->nc*sizeof(int));
    //assert(AT_map_offd_col_l2g != NULL);
    assert(AT_offd->nr == AT_diag->nr);

    AT->map_offd_col_l2g = AT_map_offd_col_l2g;
    //AT->comm_info = AT_comm_info;
    AT->comm = comm;

    AT->comm_info = Init_par_comm_info();
    AT->comm_info->nproc_neighbor = nproc_neighbor;
    if(AT->comm_info->nproc_neighbor > 0)
    {
	AT->comm_info->proc_neighbor = (int*)malloc(nproc_neighbor * sizeof(int));
	for(i=0; i<nproc_neighbor; i++) 
	    AT->comm_info->proc_neighbor[i] = A_proc_neighbor[i];

	AT->comm_info->nindex_row = (int*) calloc(AT->comm_info->nproc_neighbor,  sizeof(int));
	AT->comm_info->index_row  = (int**)malloc(AT->comm_info->nproc_neighbor * sizeof(int*));
	for(i=0; i<AT->comm_info->nproc_neighbor; i++)
	{
	    AT->comm_info->index_row[i] = (int*)malloc(AT->offd->nr * sizeof(int));
	    for(j=0; j<AT->offd->nr; j++) AT->comm_info->index_row[i][j] = -1;
	}
	Get_par_dmatcsr_comm_row_info(AT->offd, 
				      AT->comm_info->nproc_neighbor, AT->comm_info->proc_neighbor, 
				      AT->col_start,                 AT->map_offd_col_l2g, 
				      AT->comm_info->nindex_row,     AT->comm_info->index_row);

	for(i=0; i<AT->comm_info->nproc_neighbor; i++)
	{
	    AT->comm_info->index_row[i] = (int*)realloc(AT->comm_info->index_row[i], AT->comm_info->nindex_row[i]*sizeof(int));
	    if(0 != AT->comm_info->nindex_row[i]) assert(NULL != AT->comm_info->index_row[i]);
	}

	AT->comm_info->nindex_col = (int*) calloc(AT->comm_info->nproc_neighbor,  sizeof(int));
	AT->comm_info->index_col  = (int**)malloc(AT->comm_info->nproc_neighbor * sizeof(int*));
	for(i=0; i<AT->comm_info->nproc_neighbor; i++)
	{
	    AT->comm_info->index_col[i] = (int*)malloc(AT->offd->nc * sizeof(int));
	    for(j=0; j<AT->offd->nc; j++) AT->comm_info->index_col[i][j] = -1;
	}
	Get_par_dmatcsr_comm_col_info(AT->offd, 
				      AT->comm_info->nproc_neighbor, AT->comm_info->proc_neighbor, 
				      AT->col_start,                 AT->map_offd_col_l2g, 
				      AT->comm_info->nindex_col,     AT->comm_info->index_col);
	for(i=0; i<AT->comm_info->nproc_neighbor; i++)
	{
	    AT->comm_info->index_col[i] = (int*)realloc(AT->comm_info->index_col[i], AT->comm_info->nindex_col[i]*sizeof(int));
	    if(0 != AT->comm_info->nindex_col[i]) assert(NULL != AT->comm_info->index_col[i]);
	}

	//ATrint_par_dmatcsr(A);
	//ATrint_par_dmatcsr(AT);

	//可作为 par_dmatcsr->comm_info 的检测工具之一
	for(i=0; i<AT->comm_info->nproc_neighbor; i++)
	{
	    if(0 == AT->comm_info->nindex_col[i])
		assert(0 == AT->comm_info->nindex_row[i]);
	    if(0 != AT->comm_info->nindex_col[i])
	    {

		assert(0 != AT->comm_info->nindex_row[i]);
#if 0
		if(0 == AT->comm_info->nindex_row[i])
		{
		    printf("rank %d index_col %d: ", myrank, i);
		    for(j=0; j<AT->comm_info->nindex_col[i]; j++)
			printf("%d ", AT->map_offd_col_l2g[AT->comm_info->index_col[i][j]]);
		    printf("\n\n");
		}
#endif
	    }
	}
	//根据 AT->comm_info->index_col 进一步确定 AT->comm_info->proc_neighbor
	Remove_par_dmatcsr_extra_proc_neighbor(AT);
#if 0
	int *AT_proc_neighbor_flag_recv = (int*)malloc(AT->comm_info->nproc_neighbor * sizeof(int));
	for(i=0; i<AT->comm_info->nproc_neighbor; i++)
	{
	    MPI_Sendrecv(&AT->comm_info->nindex_col[i],  1, MPI_INT, 
			 AT->comm_info->proc_neighbor[i], myrank+AT->comm_info->proc_neighbor[i]*1000, 
			 &AT_proc_neighbor_flag_recv[i], 1, MPI_INT, 
			 AT->comm_info->proc_neighbor[i], AT->comm_info->proc_neighbor[i]+myrank*1000, 
			 comm, MPI_STATUS_IGNORE);
	}
	int AT_nproc_neighbor = 0;
	for(i=0; i<AT->comm_info->nproc_neighbor; i++)
	{
	    if((0!=AT->comm_info->nindex_col[i]) || (0!=AT_proc_neighbor_flag_recv[i]))
	    {
		AT->comm_info->proc_neighbor[AT_nproc_neighbor] = AT->comm_info->proc_neighbor[i];
		AT->comm_info->nindex_col[AT_nproc_neighbor] = AT->comm_info->nindex_col[i];
		AT->comm_info->nindex_row[AT_nproc_neighbor] = AT->comm_info->nindex_row[i];

		int *tmp;
		tmp = AT->comm_info->index_col[AT_nproc_neighbor];
		AT->comm_info->index_col[AT_nproc_neighbor] = AT->comm_info->index_col[i];
		AT->comm_info->index_col[i] = tmp;

		tmp = AT->comm_info->index_row[AT_nproc_neighbor];
		AT->comm_info->index_row[AT_nproc_neighbor] = AT->comm_info->index_row[i];
		AT->comm_info->index_row[i] = tmp;

		AT_nproc_neighbor++;
	    }
	}

	AT->comm_info->proc_neighbor = (int*)realloc(AT->comm_info->proc_neighbor, AT_nproc_neighbor*sizeof(int));
	if(0 != AT_nproc_neighbor) assert(NULL != AT->comm_info->proc_neighbor);

	if(0 == AT_nproc_neighbor)
	{
	    Free_par_comm_info(AT->comm_info);
	}
	else
	{
	    for(i=AT_nproc_neighbor; i<AT->comm_info->nproc_neighbor; i++)
	    {
		free(AT->comm_info->index_col[i]);
		AT->comm_info->index_col[i] = NULL;

		free(AT->comm_info->index_row[i]);
		AT->comm_info->index_row[i] = NULL;
	    }
	}
	AT->comm_info->nproc_neighbor = AT_nproc_neighbor;
	
	free(AT_proc_neighbor_flag_recv); AT_proc_neighbor_flag_recv = NULL;
#endif
    }


    free(count_ja); count_ja = NULL;
    free(count_va); count_va = NULL;
    free(recv_nnz); recv_nnz = NULL;
    free(recv_col_index); recv_col_index = NULL;
    free(recv_nr_nc_nn); recv_nr_nc_nn = NULL;
    free(recv_row_index_nc); recv_row_index_nc = NULL;
    free(send_row_index_nc); send_row_index_nc = NULL;

    Free_dmatcsr(A_offdT);
}

void Sumself_par_dvec_axpby(par_dvec *x, double a, par_dvec *y, double b)
{
    assert(x->length_global == y->length_global);
    assert(x->length        == y->length);

    int i;
    for(i=0; i<x->length; i++) y->value[i] = a*x->value[i] + b*y->value[i];
}

void Sumself_par_dvec_axpby_length(par_dvec *x, double a, par_dvec *y, double b, int length)
{
    int i;
    for(i=0; i<length; i++) y->value[i] = a*x->value[i] + b*y->value[i];
}

void Copy_par_dvec(par_dvec *x, par_dvec *x_copy)
{
    assert(x->length_global == x_copy->length_global);
    assert(x->length        == x_copy->length);

    int i;
    for(i=0; i<x->length; i++) x_copy->value[i] = x->value[i];
}

double Multi_par_dvec_dvec(par_dvec *x, par_dvec *y)
{
#if 0
    int  myrank, nproc;
    MPI_Comm comm = x->comm;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nproc);
#endif
    assert(x->length_global == y->length_global);
    assert(x->length        == y->length);
    double my_inner_product = Multi_dvec_dvec(x->value, y->value, x->length);
#if 0
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    printf("rank %d in %d: mynorm_square = %15.12f\n", myrank, nproc, mynorm_square);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
#endif
    double inner_product = 0.0;
    MPI_Allreduce(&my_inner_product, &inner_product, 1, MPI_DOUBLE, MPI_SUM, x->comm);
    return inner_product;
}

double Multi_par_dvec_dvec_length(par_dvec *x, par_dvec *y, int length, MPI_Comm comm)
{
    double my_inner_product = Multi_dvec_dvec(x->value, y->value, length);
    double inner_product = 0.0;
    MPI_Allreduce(&my_inner_product, &inner_product, 1, MPI_DOUBLE, MPI_SUM, comm);
    return inner_product;
}

void Scale_par_dvec(par_dvec *x, double a)
{
    int i;
    for(i=0; i<x->length; i++) x->value[i] *= a;
}

void Normalize_par_dvec(par_dvec *x)
{
    double my_max, max;
    Max_abs_dvec(x->value, x->length, &my_max, NULL);
    MPI_Allreduce(&my_max, &max, 1, MPI_DOUBLE, MPI_MAX, x->comm);
    Scale_dvec(x->value, 1.0/max, x->length);
    
    if(fabs(max) < MYAMGEPS)
    {
	int myrank;
	MPI_Comm_rank(x->comm, &myrank);
	if(0 == myrank)
	    fprintf(stderr, "Warning in function \"NormalizeVec\": maximum of double vector [%18.15f] maybe equals to 0!\n", fabs(max));
    }
}

void Normalize_par_dvec_length(par_dvec *x, int length, MPI_Comm comm)
{
    double my_max, max, max1, max2;
    Max_abs_dvec(x->value, length, &my_max, NULL);
#if 0
    int myrank;
    MPI_Comm_rank(comm, &myrank);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    if(myrank == 0) printf("mymax = %18.15f\n", my_max);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    if(myrank == 1) printf("mymax = %18.15f\n", my_max);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    if(myrank == 2) printf("mymax = %18.15f\n", my_max);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    if(myrank == 3) printf("mymax = %18.15f\n", my_max);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
    MPI_Barrier(comm);
#endif
    MPI_Allreduce(&my_max, &max1, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&my_max, &max2, 1, MPI_DOUBLE, MPI_MIN, comm);
    max = (MABS(max1)-MABS(max2)<MYAMGEPS)? max2: max1;
    Scale_dvec(x->value, 1.0/max, length);
    
    if(fabs(max) < MYAMGEPS)
    {
	int myrank;
	MPI_Comm_rank(comm, &myrank);
	if(0 == myrank)
	    fprintf(stderr, "Warning in function \"Normalize_par_dvec_length()\": maximum of double vector [%18.15f] maybe equals to 0!\n", fabs(max));
    }
}
#endif
