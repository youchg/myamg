#include "par_matrix_vector.h"

#if WITH_MPI

#include "io.h"
#include "mpi.h"
//#include "linear_algebra.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#if 0
double Get_par_dvec_2norm(par_dvec *x);
void Free_par_dmatcsr(par_dmatcsr *A);
void Free_par_dvec(par_dvec *x);
void Free_par_comm_info(par_comm_info *pcomm_info);
#endif

static int  myrank_id  = 0;
static int  print_rank = 4;
static void Get_dmatcsr_global_size(const char *filename, int *nr, int *nc, int *nn);
static void Get_par_dmatcsr_row_start(int nr_global, int nprocs, int *row_start);
static void Get_par_dmatcsr_comm_recv_info(int nprocs, int *col_start, dmatcsr *offd, int *map_offd_col_l2g, par_comm_recv_info *recv_info);
static void Get_par_dmatcsr_comm_send_info(int nprocs, int *col_start, dmatcsr *offd, int *map_offd_col_l2g, par_comm_send_info *send_info);
static void Get_par_dmatcsr_offd_and_map_col_offd_l2g(dmatcsr*offd, int *map_offd_col_l2g);
static void Get_par_dmatcsr_diag(dmatcsr *diag, int col_idx_min, int col_idx_max);
static void Get_neighbor_proc(int nprocs, int *col_start, dmatcsr *offd, int *map_offd_col_l2g, int *nproc_neighbor, int *proc_neighbor);
static void Separate_dmatcsr_to_diag_offd(dmatcsr *A, int col_idx_min, int col_idx_max, dmatcsr *diag, dmatcsr *offd);


par_dmatcsr *Read_par_dmatcsr(const char *filename, MPI_Comm comm)
{
    int  myrank, nproc_global, myname_len;
    char myname[MPI_MAX_PROCESSOR_NAME];

    MPI_Comm_size(comm, &nproc_global);
    MPI_Comm_rank(comm, &myrank);
    MPI_Get_processor_name(myname, &myname_len);
    myrank_id = myrank;

    double time_start, time_end;
    time_start = MPI_Wtime();
    int i;

    int nr_global, nc_global, nn_global;
    Get_dmatcsr_global_size(filename, &nr_global, &nc_global, &nn_global);

    int *row_start = (int*)calloc(nproc_global+1, sizeof(int));
    Get_par_dmatcsr_row_start(nr_global, nproc_global, row_start);

    int *col_start = (int*)calloc(nproc_global+1, sizeof(int));
    for(i=0; i<=nproc_global; i++) col_start[i] = row_start[i];

    assert(row_start[nproc_global] == nr_global);
    assert(col_start[nproc_global] == nc_global);

    int nr = row_start[myrank+1] - row_start[myrank];


    dmatcsr *A_local = Read_dmatcsr_part(filename, row_start[myrank], row_start[myrank+1]-1);

    dmatcsr *diag = (dmatcsr*)malloc(sizeof(dmatcsr));
    dmatcsr *offd = (dmatcsr*)malloc(sizeof(dmatcsr));
    Separate_dmatcsr_to_diag_offd(A_local, col_start[myrank], col_start[myrank+1]-1, diag, offd);
    Free_dmatcsr(A_local);

    Get_par_dmatcsr_diag(diag, col_start[myrank], col_start[myrank+1]-1);

    int *map_offd_col_l2g = (int*)malloc(offd->nn * sizeof(int));
    for(i=0; i<offd->nn; i++) map_offd_col_l2g[i] = -1;
    Get_par_dmatcsr_offd_and_map_col_offd_l2g(offd, map_offd_col_l2g);
    map_offd_col_l2g = (int*)realloc(map_offd_col_l2g, offd->nc*sizeof(int));
    assert(map_offd_col_l2g != NULL);
    assert(offd->nr == diag->nr);

    par_comm_recv_info *recv_info = (par_comm_recv_info*)malloc(sizeof(par_comm_recv_info));
    Get_par_dmatcsr_comm_recv_info(nproc_global, col_start, offd, map_offd_col_l2g, recv_info);

    par_comm_send_info *send_info = (par_comm_send_info*)malloc(sizeof(par_comm_send_info));
    Get_par_dmatcsr_comm_send_info(nproc_global, col_start, offd, map_offd_col_l2g, send_info);

    time_end = MPI_Wtime();
    if(myrank == print_rank)
    {
	//Write_dmatcsr_csr(offd, "../output/offd.dat");
	//Write_dmatcsr_csr(diag, "../output/diag.dat");
	printf("\n");
	printf("global info: nr = %d, nc = %d, nn = %d\n\n", nr_global, nc_global, nn_global);
	printf("myrank = %d, nr = %d, row_start = %d, row_end = %d\n\n", myrank, nr, row_start[myrank], row_start[myrank+1]-1);
	printf("read matrix time: %f\n\n", time_end-time_start);
	Print_par_comm_recv_info(recv_info);
	Print_par_comm_send_info(send_info);
	printf("\n");
    }

    par_dmatcsr *A = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
    A->diag = diag;
    A->offd = offd;
    A->nr_global = nr_global;
    A->nc_global = nc_global;
    A->nn_global = nn_global;
    A->map_offd_col_l2g = map_offd_col_l2g;
    A->row_start = row_start;
    A->col_start = col_start;
    A->comm = comm;
    A->send_info = send_info;
    A->recv_info = recv_info;
    //A->send_data = NULL;
    //A->recv_data = NULL;

    return A;
}

//proc_neighbor 应全部初始化为 -1
static void Get_neighbor_proc(int nprocs, int *col_start, dmatcsr *offd, int *map_offd_col_l2g, 
	                      int *nproc_neighbor, int *proc_neighbor)
{
    int np_neighbor = 0;
    //int *proc_neighbor = (int*)malloc(nprocs * sizeof(int));
    //for(i=0; i<nprocs; i++) proc_neighbor[i] = -1;

    assert(map_offd_col_l2g[0] >= col_start[0]);
    assert(0 == col_start[0]);

    int nc = offd->nc;
    int col_global_idx;
    int proc_mark = -1;
    int i, j;
    for(i=0; i<nc; i++)
    {
	col_global_idx = map_offd_col_l2g[i];
	if(col_global_idx >= col_start[proc_mark+1])
	{
	    for(j=proc_mark+1; j<nprocs; j++)
	    {
		if(col_global_idx < col_start[j+1])
		{
		    proc_mark = j;
		    break;
		}
	    }
	    proc_neighbor[np_neighbor++] = proc_mark;
	}
    }

    *nproc_neighbor = np_neighbor;
}

// 假定 A 是对称的
// 则 send_proc 与 recv_proc 应该相同
//
// 矩阵乘向量 A*x 时，y = A_diag*x_diag + A_offd*x_recv;
// 若 offd 对应 t 进程的第 k 行（全局编号）不全为 0， 
// 则说明 t 进程中的第 k 列（全局编号）不为0，
// 因此需要将 x_diag 中的第 k 行（全局编号）发送到 t 进程。
static void Get_par_dmatcsr_comm_send_info(int nprocs, int *col_start, dmatcsr *offd, 
	                                   int *map_offd_col_l2g, par_comm_send_info *send_info)
{
    int i, j, k;

    int  nr = offd->nr;
    int *ia = offd->ia;
    int *ja = offd->ja;

    /* 找到 nsend_proc 与 send_proc */
    int nsend_proc = 0;
    int *send_proc = (int*)malloc(nprocs * sizeof(int));
    for(i=0; i<nprocs; i++) send_proc[i] = -1;
    Get_neighbor_proc(nprocs, col_start, offd, map_offd_col_l2g, &nsend_proc, send_proc);
    send_proc = (int*)realloc(send_proc, nsend_proc*sizeof(int));
    assert(send_proc != NULL);

    int *nidx = (int*) calloc(nsend_proc,  sizeof(int));
    int **idx = (int**)malloc(nsend_proc * sizeof(int*));
    for(i=0; i<nsend_proc; i++)
    {
	idx[i] = (int*)malloc(nr * sizeof(int));
	for(k=0; k<nr; k++) idx[i][k] = -1;
    }

    int **row_status = (int**)malloc(nsend_proc * sizeof(int*));
    for(i=0; i<nsend_proc; i++) row_status[i] = (int*)calloc(nr, sizeof(int));

    int proc_mark = 0;
    int col_global_idx;

    for(i=0; i<nr; i++)
    {
	proc_mark = 0;
	for(k=ia[i]; k<ia[i+1]; k++)
	{
	    col_global_idx = map_offd_col_l2g[ja[k]];

	    if((col_global_idx >= col_start[send_proc[proc_mark]]) && 
	       (col_global_idx <  col_start[send_proc[proc_mark]+1]))
	    {
		if(0 == row_status[proc_mark][i])
		{
		    row_status[proc_mark][i] = 1;
		    idx[proc_mark][nidx[proc_mark]] = i;
		    nidx[proc_mark]++;
		}
	    }
	    else
	    {
		for(j=proc_mark+1; j<nsend_proc; j++)
		{
		    if(col_global_idx < col_start[send_proc[j]+1])
		    {
			proc_mark = j;
			break;
		    }
		}
		if(0 == row_status[proc_mark][i])
		{
		    row_status[proc_mark][i] = 1;
		    idx[proc_mark][nidx[proc_mark]] = i;
		    nidx[proc_mark]++;
		}
	    }
	}
    }

    //idx 是否也要 realloc ?
    send_info->nproc  = nsend_proc;
    send_info->proc   = send_proc;
    send_info->nindex = nidx;
    send_info->index  = idx;

    for(i=0; i<nsend_proc; i++) {free(row_status[i]); row_status[i] = NULL;}
    free(row_status); row_status = NULL;
}

//通过遍历 map_offd_col_l2g 获得 recv_info
//遍历 offd->ja 是不对的
//因为 offd->ja 中的列是按行排的, 不是按邻居进程
//而 map_offd_col_l2g 中的元素则是把每个邻居进程的列从小到大排在一起
//
// 矩阵乘向量 A*x 时，y = A_diag*x_diag + A_offd*x_recv;
// 由于矩阵的 CSR 格式中，不会有全为 0 的列，
// 将 offd 的列号从小到大排序，第 j 列（全局编号）对应全局 x 的第 j 行，
// 需要找出第 j 列（全局编号）对应的进程 t，
// 然后接收进程 t 上 x_diag 的第 j 行（全局编号）。
static void Get_par_dmatcsr_comm_recv_info(int nprocs, int *col_start, dmatcsr *offd, int *map_offd_col_l2g, par_comm_recv_info *recv_info)
{
    int i;

    /* 找到 nrecv_proc 与 recv_proc */
    int nrecv_proc = 0;
    int *recv_proc = (int*)malloc(nprocs * sizeof(int));
    for(i=0; i<nprocs; i++) recv_proc[i] = -1;
    Get_neighbor_proc(nprocs, col_start, offd, map_offd_col_l2g, &nrecv_proc, recv_proc);
    recv_proc = (int*)realloc(recv_proc, nrecv_proc*sizeof(int));
    assert(recv_proc != NULL);

    int *map_recv_start = (int*)calloc(nrecv_proc+1, sizeof(int));

    int nc = offd->nc;
    int proc_mark = 0;
    int count_map_recv_start = 0;
    int col_global_idx;

    for(i=0; i<nc; i++)
    {
	col_global_idx = map_offd_col_l2g[i];
	if((col_global_idx >= col_start[recv_proc[proc_mark]]) && 
	   (col_global_idx <  col_start[recv_proc[proc_mark]+1]))
	{
	    count_map_recv_start++;
	}
	else
	{
	    proc_mark++;
	    map_recv_start[proc_mark] = count_map_recv_start;
	    count_map_recv_start++;
	}
    }
    //printf("nrecv_proc = %d, count = %d\n", nrecv_proc, count_map_recv_start);
    map_recv_start[nrecv_proc] = count_map_recv_start;

    assert(count_map_recv_start == nc);
    assert(proc_mark == nrecv_proc-1);

    recv_info->nproc = nrecv_proc;
    recv_info->proc  = recv_proc;
    recv_info->start = map_recv_start;
}

//map_offd_col_l2g 应全部初始化为 -1
//此时 offd->nc 等于 nc_global
//最后将 offd->nc 从 nc_global 改为正确的列数
static void Get_par_dmatcsr_offd_and_map_col_offd_l2g(dmatcsr*offd, int *map_offd_col_l2g)
{
    int  nn = offd->nn;
    int  nc = offd->nc;
    int *ja = offd->ja;

    int k;

    int *isoffd = (int*)calloc(nc, sizeof(int));
    for(k=0; k<nc; k++) isoffd[k]     = -1;
    for(k=0; k<nn; k++) isoffd[ja[k]] =  1;

    int count_map = 0;
    for(k=0; k<nc; k++)
    {
	if(isoffd[k] > 0)
	{
	    map_offd_col_l2g[count_map] = k;
	    isoffd[k] = count_map;
	    count_map++;
	}
    }

    for(k=0; k<nn; k++) ja[k] = isoffd[ja[k]];

    int ja_max = -1;
    for(k=0; k<nn; k++) if(ja[k] > ja_max) ja_max = ja[k];

    //printf("nc = %d, ja_max = %d, count_map = %d\n", offd->nc, ja_max, count_map);
    offd->nc = count_map;
    assert(count_map == ja_max+1);

    free(isoffd);
} 

static void Get_par_dmatcsr_diag(dmatcsr *diag, int col_idx_min, int col_idx_max)
{
    int  nn = diag->nn;
    int *ja = diag->ja;
    diag->nc = col_idx_max - col_idx_min + 1;

    assert(diag->nc == diag->nr);

    int k;
    for(k=0; k<nn; k++) ja[k] -= col_idx_min;
}

static void Get_dmatcsr_global_size(const char *filename, int *nr, int *nc, int *nn)
{
    FILE *file = My_fopen(filename, "r");
    assert(EOF!=fscanf(file, "%d\n",   nr));
    assert(EOF!=fscanf(file, "%d\n",   nc));
    assert(EOF!=fscanf(file, "%d\n\n", nn));
    fclose(file); file = NULL;
}

static void Get_par_dmatcsr_row_start(int nr_global, int nprocs, int *row_start)
{
    int proc_quotient = nr_global / nprocs;
    int proc_residual = nr_global % nprocs;
    int i;
    for(i=0; i<nprocs; i++)
    {
	if(i < proc_residual) row_start[i+1] = row_start[i] + proc_quotient + 1;
	else                  row_start[i+1] = row_start[i] + proc_quotient;
    }
}

static void Separate_dmatcsr_to_diag_offd(dmatcsr *A, int col_idx_min, int col_idx_max, dmatcsr *diag, dmatcsr *offd)
{
    int     nr = A->nr;
    int     nc = A->nc;
    int     nn = A->nn;
    int    *ia = A->ia;
    int    *ja = A->ja;
    double *va = A->va;

    int  i, j, k;
    int  nn_diag = 0;
    int  nn_offd = 0;
    int *ia_diag = (int*)calloc(nr+1, sizeof(int));
    int *ia_offd = (int*)calloc(nr+1, sizeof(int));
    int *isdiag  = (int*)calloc(nn, sizeof(int));
    for(i=0; i<nr; i++)
    {
	for(k=ia[i]; k<ia[i+1]; k++)
	{
	    j = ja[k];
	    if(j>=col_idx_min && j<=col_idx_max) {nn_diag++; isdiag[k] = 1;}
	    else                                 {nn_offd++;}
	}
	ia_diag[i+1] = nn_diag;
	ia_offd[i+1] = nn_offd;
    }

    assert(nn_diag + nn_offd == nn);

    int    *ja_diag = (int*)   malloc(nn_diag * sizeof(int));
    int    *ja_offd = (int*)   malloc(nn_offd * sizeof(int));
    double *va_diag = (double*)malloc(nn_diag * sizeof(double));
    double *va_offd = (double*)malloc(nn_offd * sizeof(double));

    int count_diag = 0;
    int count_offd = 0;
    for(k=0; k<nn; k++)
    {
	if(isdiag[k]) ja_diag[count_diag++] = ja[k];
	else           ja_offd[count_offd++] = ja[k];
    }
    assert(count_diag == nn_diag);
    assert(count_offd == nn_offd);

    count_diag = 0;
    count_offd = 0;
    for(k=0; k<nn; k++)
    {
	if(isdiag[k]) va_diag[count_diag++] = va[k];
	else          va_offd[count_offd++] = va[k];
    }
    assert(count_diag == nn_diag);
    assert(count_offd == nn_offd);

    diag->nr = nr;
    diag->nc = nc;
    diag->nn = nn_diag;
    diag->ia = ia_diag;
    diag->ja = ja_diag;
    diag->va = va_diag;

    offd->nr = nr;
    offd->nc = nc;
    offd->nn = nn_offd;
    offd->ia = ia_offd;
    offd->ja = ja_offd;
    offd->va = va_offd;

    free(isdiag);
}

void Free_par_dmatcsr(par_dmatcsr *A)
{
    Free_dmatcsr(A->diag);
    Free_dmatcsr(A->offd);

    free(A->map_offd_col_l2g); A->map_offd_col_l2g = NULL;
    free(A->row_start); A->row_start = NULL;
    free(A->col_start); A->col_start = NULL;

    Free_par_comm_send_info(A->send_info);
    Free_par_comm_recv_info(A->recv_info);
    //Free_par_comm_data(A->send_data);
    //Free_par_comm_data(A->recv_data);

    free(A); A = NULL;
}

void Free_par_imatcsr(par_imatcsr *A)
{
    Free_imatcsr(A->diag);
    Free_imatcsr(A->offd);

    free(A->map_offd_col_l2g); A->map_offd_col_l2g = NULL;
    free(A->row_start); A->row_start = NULL;
    free(A->col_start); A->col_start = NULL;

    Free_par_comm_send_info(A->send_info);
    Free_par_comm_recv_info(A->recv_info);
    //Free_par_comm_data(A->send_data);
    //Free_par_comm_data(A->recv_data);

    free(A); A = NULL;
}

void Free_par_comm_recv_info(par_comm_recv_info *info)
{
    if(NULL != info)
    {
	free(info->proc);  info->proc  = NULL;
	free(info->start); info->start = NULL;
	free(info);        info        = NULL;
    }
}

void Free_par_comm_send_info(par_comm_send_info *info)
{
    if(NULL != info)
    {
	int i;
	free(info->proc);   info->proc  = NULL;
	free(info->nindex); info->nindex = NULL;
	for(i=0; i<info->nproc; i++) {free(info->index[i]); info->index[i] = NULL;}
	free(info->index);  info->index = NULL;
	free(info);         info        = NULL;
    }
}

void Free_par_comm_data(par_comm_data *comm_data)
{
    if(NULL != comm_data)
    {
	free(comm_data->data); comm_data->data = NULL;
	free(comm_data);       comm_data       = NULL;
    }
}

void Print_par_comm_recv_info(par_comm_recv_info *info)
{
    int  nproc = info->nproc;
    int *proc  = info->proc;
    int *start = info->start;
    int i;
    printf("\n----------------------------------\n");
    printf("par_comm_recv_info:\n");
    printf("nproc = %d\n", nproc);
    printf("proc  = ");
    for(i=0; i< nproc; i++) printf("%d ", proc[i]);
    printf("\n");
    printf("start = ");
    for(i=0; i<=nproc; i++) printf("%d ", start[i]);
    printf("\n");
    printf("----------------------------------\n\n");
}

void Print_par_comm_send_info(par_comm_send_info *info)
{
    int  nproc  = info->nproc;
    int *proc   = info->proc;
    int *nindex = info->nindex;
    int **index = info->index;
    int i, j;
    printf("\n----------------------------------\n");
    printf("par_comm_send_info:\n");
    printf("nproc  = %d\n", nproc);
    printf("proc   = ");
    for(i=0; i<nproc; i++) printf("%d ", proc[i]);
    printf("\n");
    printf("nindex = ");
    for(i=0; i<nproc; i++) printf("%d ", nindex[i]);
    printf("\n");
    printf("index  = \n");
    for(i=0; i<nproc; i++) 
    {
	printf("         ");
	for(j=0; j<nindex[i]; j++) 
	    printf("%02d ", index[i][j]); 
	printf("\n");
    }
    printf("----------------------------------\n\n");
}

par_dvec *Init_par_dvec_from_par_dmatcsr(par_dmatcsr *A)
{
    par_dvec *x = (par_dvec*)malloc(sizeof(par_dvec));
    x->length_global = A->nr_global;
    x->length        = A->diag->nr;
    x->value         = (double*)calloc(x->length, sizeof(double));
    x->comm          = A->comm;
    x->send_info     = Copy_par_comm_send_info(A->send_info);
    x->recv_info     = Copy_par_comm_recv_info(A->recv_info);

    int i;
    int nproc_send   = x->send_info->nproc;
    x->send_data     = (double**)malloc(nproc_send * sizeof(double*));
    for(i=0; i<nproc_send; i++)
    {
	x->send_data[i] = (double*)calloc(x->send_info->nindex[i], sizeof(double));
    }

    int nproc_recv   = x->recv_info->nproc;
    x->recv_data     = (double*)calloc(x->recv_info->start[nproc_recv], sizeof(double));

    return x;
}

void Free_par_dvec(par_dvec *x)
{
    if(NULL != x)
    {
	free(x->value); x->value = NULL;

	int i;
	int nproc_send = x->send_info->nproc;
	for(i=0; i<nproc_send; i++)
	{
	    free(x->send_data[i]);
	    x->send_data[i] = NULL;
	}
	free(x->send_data); x->send_data = NULL;

	free(x->recv_data); x->recv_data = NULL;

	Free_par_comm_send_info(x->send_info);
	Free_par_comm_recv_info(x->recv_info);

	free(x); x = NULL;
    }
}

par_comm_send_info *Copy_par_comm_send_info (par_comm_send_info *send_info)
{
    int i, j;

    int  nproc  = send_info->nproc;
    int  *proc  = send_info->proc;
    int **index = send_info->index;
    int *nindex = send_info->nindex;

    int *proc_copy = (int*)malloc(nproc * sizeof(int));
    for(i=0; i<nproc; i++) proc_copy[i] = proc[i];

    int *nindex_copy = (int*)malloc(nproc * sizeof(int));
    for(i=0; i<nproc; i++) nindex_copy[i] = nindex[i];

    int **index_copy = (int**)malloc(nproc * sizeof(int*));
    for(i=0; i<nproc; i++)
    {
	index_copy[i] = (int*)malloc(nindex[i] * sizeof(int));
	for(j=0; j<nindex[i]; j++) index_copy[i][j] = index[i][j];
    }

    par_comm_send_info *send_info_copy = (par_comm_send_info*)malloc(sizeof(par_comm_send_info));
    send_info_copy->nproc  =  nproc;
    send_info_copy->proc   =   proc_copy;
    send_info_copy->index  =  index_copy;
    send_info_copy->nindex = nindex_copy;

    return send_info_copy;
}

par_comm_recv_info *Copy_par_comm_recv_info (par_comm_recv_info *recv_info)
{
    int  nproc  = recv_info->nproc;
    int  *proc  = recv_info->proc;
    int  *start = recv_info->start;

    int i;

    int *proc_copy = (int*)malloc(nproc * sizeof(int));
    for(i=0; i<nproc; i++) proc_copy[i] = proc[i];

    int *start_copy = (int*)malloc((nproc+1) * sizeof(int));
    for(i=0; i<=nproc; i++) start_copy[i] = start[i];

    par_comm_recv_info *recv_info_copy = (par_comm_recv_info*)malloc(sizeof(par_comm_recv_info));
    recv_info_copy->nproc  = nproc;
    recv_info_copy->proc   =  proc_copy;
    recv_info_copy->start  = start_copy;

    return recv_info_copy;
}

#if 0
void Free_par_comm_send_data(par_comm_send_data *send_data)
{
    if(NULL != send_data)
    {
	int i;
	for(i=0; i<info->nproc; i++) {free(send_data->data[i]); send_data->data[i] = NULL;}
	free(send_data->data);  send_data->data = NULL;

	//free(send_data->ndata); ndata = NULL;

	if(NULL != send_data->info) 
	    Free_par_comm_send_info(send_data->info);

	free(send_data); send_data = NULL;
    }
}

void Free_par_comm_recv_data(par_comm_recv_data *recv_data)
{
    if(NULL != recv_data)
    {
	free(recv_data->data);  recv_data->data = NULL;

	if(NULL != recv_data->info) 
	    Free_par_comm_recv_info(recv_data->info);

	free(recv_data); recv_data = NULL;
    }
}

par_comm_send_data *Create_par_comm_send_data(par_comm_send_info *info)
{
    int i, j;
    int **data = (int**)malloc(info->nproc * sizeof(int*));
    for(i=0; i<info->nproc; i++)
    {
	data[i] = (int*)calloc(info->nindex[i], sizeof(int));
    }

    par_comm_send_data *send_data = (par_comm_send_data*)malloc(sizeof(par_comm_send_data));
    send_data->info = Copy_par_comm_send_info(info);
    send_data->data = data;

    return send_data;
}

par_comm_recv_data *Create_par_comm_recv_data(par_comm_recv_info *info)
{
    par_comm_recv_data *recv_data = (par_comm_recv_data*)malloc(sizeof(par_comm_recv_data));
    recv_data->info = Copy_par_comm_recv_info(info);
    recv_data->data = malloc;

    return recv_data;
}
#endif


#endif
