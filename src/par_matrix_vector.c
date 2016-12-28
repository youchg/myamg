#include "par_matrix_vector.h"

#if WITH_MPI

#include "io.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

extern int  print_rank;

static void Get_dmatcsr_global_size(const char *filename, int *nr, int *nc, int *nn);
static void Get_par_dmatcsr_row_start(int nr_global, int nprocs, int *row_start);
static void Separate_dmatcsr_to_diag_offd(dmatcsr *A, int col_idx_min, int col_idx_max, dmatcsr *diag, dmatcsr *offd);
static void Get_par_dmatcsr_diag(dmatcsr *diag, int col_idx_min, int col_idx_max);
static void Get_par_dmatcsr_offd_and_map_col_offd_l2g(dmatcsr*offd, int *map_offd_col_l2g);

static void Get_par_dmatcsr_comm_info(int nproc_global, dmatcsr *offd, int *col_start, int *map_offd_col_l2g, par_comm_info *comm_info);
static void Get_neighbor_proc(int nproc_global, dmatcsr *offd, int *col_start, int *map_offd_col_l2g, int *nproc_neighbor, int *proc_neighbor);
static void Get_par_dmatcsr_comm_row_info(dmatcsr *offd, int nproc_neighbor, int *proc_neighbor, int *col_start, int *map_offd_col_l2g, int *nidx, int **idx);
static void Get_par_dmatcsr_comm_col_info(dmatcsr *offd, int nproc_neighbor, int *proc_neighbor, int *col_start, int *map_offd_col_l2g, int *nidx, int **idx);

par_dmatcsr *Read_par_dmatcsr(const char *filename, MPI_Comm comm)
{
    int  myrank, nproc_global;
    MPI_Comm_size(comm, &nproc_global);
    MPI_Comm_rank(comm, &myrank);

    double time_start, time_end;
    time_start = MPI_Wtime();

    int i;

    //-----------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------
    int nr_global, nc_global, nn_global;
    Get_dmatcsr_global_size(filename, &nr_global, &nc_global, &nn_global);

    int *row_start = (int*)calloc(nproc_global+1, sizeof(int));
    Get_par_dmatcsr_row_start(nr_global, nproc_global, row_start);

    int *col_start = (int*)calloc(nproc_global+1, sizeof(int));
    for(i=0; i<=nproc_global; i++) col_start[i] = row_start[i];

    assert(row_start[nproc_global] == nr_global);
    assert(col_start[nproc_global] == nc_global);

    //-----------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------
    dmatcsr *A_local = Read_dmatcsr_part(filename, row_start[myrank], row_start[myrank+1]-1);

    dmatcsr *diag = (dmatcsr*)malloc(sizeof(dmatcsr));
    dmatcsr *offd = (dmatcsr*)malloc(sizeof(dmatcsr));
    //将 A_local 分裂成 diag 和 offd
    //此时 diag 和 offd 中的列还是全局编号
    Separate_dmatcsr_to_diag_offd(A_local, col_start[myrank], col_start[myrank+1]-1, diag, offd);

    Free_dmatcsr(A_local);

    //-----------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------
    //将 diag 中的列号转化为局部独立编号
    Get_par_dmatcsr_diag(diag, col_start[myrank], col_start[myrank+1]-1);

    //-----------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------
    //map_offd_col_l2g 要先初始化为 -1
    int *map_offd_col_l2g = (int*)malloc(offd->nn * sizeof(int));
    for(i=0; i<offd->nn; i++) map_offd_col_l2g[i] = -1;
    //将 offd 中的列号转化为局部独立编号，同时生成局部编号和全局编号的映射
    Get_par_dmatcsr_offd_and_map_col_offd_l2g(offd, map_offd_col_l2g);
    map_offd_col_l2g = (int*)realloc(map_offd_col_l2g, offd->nc*sizeof(int));

    assert(map_offd_col_l2g != NULL);
    assert(offd->nr == diag->nr);

    //-----------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------
    par_comm_info *comm_info = (par_comm_info*)malloc(sizeof(par_comm_info));
    Get_par_dmatcsr_comm_info(nproc_global, offd, col_start, map_offd_col_l2g, comm_info);

    time_end = MPI_Wtime();
    if(myrank == print_rank)
    {
	//Write_dmatcsr_csr(offd, "../output/offd.dat");
	//Write_dmatcsr_csr(diag, "../output/diag.dat");
	printf("\n");
	printf("global info: nr = %d, nc = %d, nn = %d\n\n", nr_global, nc_global, nn_global);
	printf("myrank = %d, nr = %d, row_start = %d, row_end = %d\n\n", 
		myrank, row_start[myrank+1]-row_start[myrank], row_start[myrank], row_start[myrank+1]-1);
	printf("read matrix time: %f\n\n", time_end-time_start);
	//Print_dmatcsr(diag);
	//Print_dmatcsr(offd);
	Print_par_comm_info(comm_info);
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
    A->comm_info = comm_info;

    return A;
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

static void Get_par_dmatcsr_diag(dmatcsr *diag, int col_idx_min, int col_idx_max)
{
    int  nn = diag->nn;
    int *ja = diag->ja;
    diag->nc = col_idx_max - col_idx_min + 1;

    assert(diag->nc == diag->nr);

    int k;
    for(k=0; k<nn; k++) ja[k] -= col_idx_min;
}

/*
* map_offd_col_l2g 应全部初始化为 -1
* 此时 offd->nc 等于 nc_global
* 最后将 offd->nc 从 nc_global 改为正确的列数
*/
static void Get_par_dmatcsr_offd_and_map_col_offd_l2g(dmatcsr*offd, int *map_offd_col_l2g)
{
    int  nn = offd->nn;
    int  nc = offd->nc;
    int *ja = offd->ja;

    int k;

    //好的做法应该是将 offd->ja 排序
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

    offd->nc = count_map;
    assert(count_map == ja_max+1);

    free(isoffd);
} 

static void Get_par_dmatcsr_comm_info(int nproc_global, dmatcsr *offd, int *col_start, int *map_offd_col_l2g, par_comm_info *comm_info)
{
    int i, k;

    int  nr = offd->nr;
    int  nc = offd->nc;

    int nproc_neighbor = 0;
    int *proc_neighbor = (int*)malloc(nproc_global * sizeof(int));
    for(i=0; i<nproc_global; i++) proc_neighbor[i] = -1;
    Get_neighbor_proc(nproc_global, offd, col_start, map_offd_col_l2g, &nproc_neighbor, proc_neighbor);
    proc_neighbor = (int*)realloc(proc_neighbor, nproc_neighbor*sizeof(int));
    assert(proc_neighbor != NULL);

    int *nidx_row = (int*) calloc(nproc_neighbor,  sizeof(int));
    int **idx_row = (int**)malloc(nproc_neighbor * sizeof(int*));
    for(i=0; i<nproc_neighbor; i++)
    {
	idx_row[i] = (int*)malloc(nr * sizeof(int));
	for(k=0; k<nr; k++) idx_row[i][k] = -1;
    }
    Get_par_dmatcsr_comm_row_info(offd, nproc_neighbor, proc_neighbor, col_start, map_offd_col_l2g, nidx_row, idx_row);
    for(i=0; i<nproc_neighbor; i++)
    {
	idx_row[i] = (int*)realloc(idx_row[i], nidx_row[i]*sizeof(int));
	assert(idx_row[i] != NULL);
    }

    int *nidx_col = (int*) calloc(nproc_neighbor,  sizeof(int));
    int **idx_col = (int**)malloc(nproc_neighbor * sizeof(int*));
    for(i=0; i<nproc_neighbor; i++)
    {
	idx_col[i] = (int*)malloc(nc * sizeof(int));
	for(k=0; k<nc; k++) idx_col[i][k] = -1;
    }
    Get_par_dmatcsr_comm_col_info(offd, nproc_neighbor, proc_neighbor, col_start, map_offd_col_l2g, nidx_col, idx_col);
    for(i=0; i<nproc_neighbor; i++)
    {
	idx_col[i] = (int*)realloc(idx_col[i], nidx_col[i]*sizeof(int));
	assert(idx_col[i] != NULL);
    }
    
    comm_info->nproc_neighbor = nproc_neighbor;
    comm_info->proc_neighbor  = proc_neighbor;
    comm_info->nindex_row     = nidx_row;
    comm_info->index_row      = idx_row;
    comm_info->nindex_col     = nidx_col;
    comm_info->index_col      = idx_col;
}


/*
* proc_neighbor 应全部初始化为 -1
*/
static void Get_neighbor_proc(int nproc_global, dmatcsr *offd, int *col_start, int *map_offd_col_l2g, 
	                      int *nproc_neighbor, int *proc_neighbor)
{
    int np_neighbor = 0;

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
	    for(j=proc_mark+1; j<nproc_global; j++)
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

/*
* 假定 A 是对称的
* 则 send_proc 与 recv_proc 应该相同
*
* 矩阵乘向量 A*x 时，y = A_diag*x_diag + A_offd*x_recv;
* 若 offd 对应 t 进程的第 k 行（全局编号）不全为 0， 
* 则说明 t 进程中的第 k 列（全局编号）不为0，
* 因此需要将 x_diag 中的第 k 行（全局编号）发送到 t 进程。
*/
static void Get_par_dmatcsr_comm_row_info(dmatcsr *offd, 
	                                  int nproc_neighbor, int  *proc_neighbor, 
	                                  int *col_start,     int  *map_offd_col_l2g, 
					  int *nidx,          int **idx)
{
    int i, j, k;

    int  nr = offd->nr;
    int *ia = offd->ia;
    int *ja = offd->ja;

    int **row_status = (int**)malloc(nproc_neighbor * sizeof(int*));
    for(i=0; i<nproc_neighbor; i++) row_status[i] = (int*)calloc(nr, sizeof(int));

    int proc_mark = 0;
    int col_global_idx;

    for(i=0; i<nr; i++)
    {
	proc_mark = 0;
	for(k=ia[i]; k<ia[i+1]; k++)
	{
	    col_global_idx = map_offd_col_l2g[ja[k]];

	    assert(col_global_idx >= col_start[proc_neighbor[proc_mark]]);

	    if(col_global_idx >= col_start[proc_neighbor[proc_mark]+1])
	    {
		for(j=proc_mark+1; j<nproc_neighbor; j++)
		{
		    if(col_global_idx < col_start[proc_neighbor[j]+1])
		    {
			proc_mark = j;
			break;
		    }
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

    for(i=0; i<nproc_neighbor; i++) {free(row_status[i]); row_status[i] = NULL;}
    free(row_status); row_status = NULL;
}

/*
* 通过遍历 map_offd_col_l2g 获得 recv_info
* 遍历 offd->ja 是不对的
* 因为 offd->ja 中的列是按行排的, 不是按邻居进程
* 而 map_offd_col_l2g 中的元素则是把每个邻居进程的列从小到大排在一起

* 矩阵乘向量 A*x 时，y = A_diag*x_diag + A_offd*x_recv;
* 由于矩阵的 CSR 格式中，不会有全为 0 的列，
* 将 offd 的列号从小到大排序，第 j 列（全局编号）对应全局 x 的第 j 行，
* 需要找出第 j 列（全局编号）对应的进程 t，
* 然后接收进程 t 上 x_diag 的第 j 行（全局编号）。
*/
static void Get_par_dmatcsr_comm_col_info(dmatcsr *offd, 
	                                  int nproc_neighbor, int *proc_neighbor,
					  int *col_start,     int *map_offd_col_l2g, 
					  int *nidx,          int **idx)
{
    int i;

    int nc = offd->nc;
    int proc_mark = 0;
    int col_global_idx;

    for(i=0; i<nc; i++)
    {
	col_global_idx = map_offd_col_l2g[i];
	assert(col_global_idx >= col_start[proc_neighbor[proc_mark]]);

	if(col_global_idx >= col_start[proc_neighbor[proc_mark]+1]) proc_mark++;

	idx[proc_mark][nidx[proc_mark]] = i;
	nidx[proc_mark]++;
    }

    assert(proc_mark == nproc_neighbor-1);
}

void Free_par_dmatcsr(par_dmatcsr *A)
{
    Free_dmatcsr(A->diag);
    Free_dmatcsr(A->offd);

    if(NULL != A->map_offd_col_l2g)
    {
	free(A->map_offd_col_l2g);
	A->map_offd_col_l2g = NULL;
    }
    if(NULL != A->row_start)
    {
	free(A->row_start);
	A->row_start = NULL;
    }
    if(NULL != A->col_start)
    {
	free(A->col_start);
	A->col_start = NULL;
    }

    Free_par_comm_info(A->comm_info);

    free(A); A = NULL;
}

void Free_par_imatcsr(par_imatcsr *A)
{
    Free_imatcsr(A->diag);
    Free_imatcsr(A->offd);

    if(NULL != A->map_offd_col_l2g)
    {
	free(A->map_offd_col_l2g);
	A->map_offd_col_l2g = NULL;
    }
    if(NULL != A->row_start)
    {
	free(A->row_start);
	A->row_start = NULL;
    }
    if(NULL != A->col_start)
    {
	free(A->col_start);
	A->col_start = NULL;
    }

    Free_par_comm_info(A->comm_info);

    free(A); A = NULL;
}

void Free_par_comm_info(par_comm_info *info)
{
    if(NULL != info)
    {
	int i;
	free(info->proc_neighbor);   info->proc_neighbor = NULL;

	free(info->nindex_row); info->nindex_row = NULL;
	for(i=0; i<info->nproc_neighbor; i++)
	{
	    free(info->index_row[i]); 
	    info->index_row[i] = NULL;
	}
	free(info->index_row);  info->index_row = NULL;

	free(info->nindex_col); info->nindex_col = NULL;
	for(i=0; i<info->nproc_neighbor; i++)
	{
	    free(info->index_col[i]); 
	    info->index_col[i] = NULL;
	}
	free(info->index_col);  info->index_col = NULL;

	free(info);         info        = NULL;
    }
}

void Print_par_comm_info(par_comm_info *info)
{
    int  nproc_neighbor = info->nproc_neighbor;
    int *proc_neighbor  = info->proc_neighbor;
    int *nindex_row     = info->nindex_row;
    int **index_row     = info->index_row;
    int *nindex_col     = info->nindex_col;
    int **index_col     = info->index_col;

    int i, j;
    printf("\n----------------------------------\n");
    printf("par_comm_info:\n");
    printf("nproc_neighbor  = %d\n", nproc_neighbor);
    printf("proc_neighbor   = ");
    for(i=0; i<nproc_neighbor; i++) printf("%d ", proc_neighbor[i]);
    printf("\n");

    printf("nindex_row = ");
    for(i=0; i<nproc_neighbor; i++) printf("%d ", nindex_row[i]);
    printf("\n");
    printf("index_row  = \n");
    for(i=0; i<nproc_neighbor; i++) 
    {
	printf("             ");
	for(j=0; j<nindex_row[i]; j++) 
	    printf("%02d ", index_row[i][j]); 
	printf("\n");
    }
    printf("\n");

    printf("nindex_col = ");
    for(i=0; i<nproc_neighbor; i++) printf("%d ", nindex_col[i]);
    printf("\n");
    printf("index_col  = \n");
    for(i=0; i<nproc_neighbor; i++) 
    {
	printf("             ");
	for(j=0; j<nindex_col[i]; j++) 
	    printf("%02d ", index_col[i][j]); 
	printf("\n");
    }
    printf("----------------------------------\n\n");
}

par_dvec *Init_par_dvec_mv(par_dmatcsr *A)
{
    par_dvec *x = (par_dvec*)malloc(sizeof(par_dvec));
    x->length_global = A->nr_global;
    x->length        = A->diag->nr;
    x->value         = (double*)calloc(x->length, sizeof(double));
    x->comm          = A->comm;
    x->comm_info     = Copy_par_comm_info(A->comm_info);

    int nproc_neighbor   = x->comm_info->nproc_neighbor;

    int i;
    x->send_data = (double**)malloc(nproc_neighbor * sizeof(double*));
    for(i=0; i<nproc_neighbor; i++)
    {
	x->send_data[i] = (double*)calloc(x->comm_info->nindex_row[i], sizeof(double));
    }

    int nrecv_data = 0;
    for(i=0; i<nproc_neighbor; i++) nrecv_data += x->comm_info->nindex_col[i];
    x->recv_data = (double*)calloc(nrecv_data, sizeof(double));

    x->recv_data_start = (int*)calloc(nproc_neighbor+1, sizeof(int));
    for(i=1; i<=nproc_neighbor; i++) x->recv_data_start[i] = x->comm_info->nindex_col[i-1];
    for(i=1; i<=nproc_neighbor; i++) x->recv_data_start[i] += x->recv_data_start[i-1];

    return x;
}

par_ivec *Init_par_ivec_length_comm(int length, MPI_Comm comm)
{
    par_ivec *x = (par_ivec*)malloc(sizeof(par_ivec));
    x->length_global   = -1;
    x->length          = length;
    x->value           = (int*)calloc(x->length, sizeof(int));
    x->comm            = comm;
    x->comm_info       = NULL;
    x->send_data       = NULL;
    x->recv_data       = NULL;
    x->recv_data_start = NULL;

    return x;
}

void Free_par_dvec(par_dvec *x)
{
    if(NULL != x)
    {
	free(x->value); x->value = NULL;

	int i;
	if(NULL != x->send_data)
	{
	    assert(NULL != x->comm_info);
	    int nproc_neighbor = x->comm_info->nproc_neighbor;

	    for(i=0; i<nproc_neighbor; i++)
	    {
		free(x->send_data[i]);
		x->send_data[i] = NULL;
	    }
	    free(x->send_data); x->send_data = NULL;
	}

	if(NULL != x->recv_data)
	{
	    free(x->recv_data);
	    x->recv_data = NULL;
	}

	Free_par_comm_info(x->comm_info);

	free(x->recv_data_start); x->recv_data_start = NULL;

	free(x); x = NULL;
    }
}

void Free_par_ivec(par_ivec *x)
{
    if(NULL != x)
    {
	free(x->value); x->value = NULL;

	int i;
	if(NULL != x->send_data)
	{
	    assert(NULL != x->comm_info);
	    int nproc_neighbor = x->comm_info->nproc_neighbor;

	    for(i=0; i<nproc_neighbor; i++)
	    {
		free(x->send_data[i]);
		x->send_data[i] = NULL;
	    }
	    free(x->send_data); x->send_data = NULL;
	}

	if(NULL != x->recv_data)
	{
	    free(x->recv_data);
	    x->recv_data = NULL;
	}

	Free_par_comm_info(x->comm_info);

	free(x->recv_data_start); x->recv_data_start = NULL;

	free(x); x = NULL;
    }
}

par_comm_info *Copy_par_comm_info (par_comm_info *info)
{
    int i, j;

    int  nproc_neighbor = info->nproc_neighbor;
    int  *proc_neighbor = info->proc_neighbor;

    int *proc_neighbor_copy = (int*)malloc(nproc_neighbor * sizeof(int));
    for(i=0; i<nproc_neighbor; i++) proc_neighbor_copy[i] = proc_neighbor[i];

    int *nindex_row = info->nindex_row;
    int *nindex_row_copy = (int*)malloc(nproc_neighbor * sizeof(int));
    for(i=0; i<nproc_neighbor; i++) nindex_row_copy[i] = nindex_row[i];

    int **index_row = info->index_row;
    int **index_row_copy = (int**)malloc(nproc_neighbor * sizeof(int*));
    for(i=0; i<nproc_neighbor; i++)
    {
	index_row_copy[i] = (int*)malloc(nindex_row[i] * sizeof(int));
	for(j=0; j<nindex_row[i]; j++) index_row_copy[i][j] = index_row[i][j];
    }

    int *nindex_col = info->nindex_col;
    int *nindex_col_copy = (int*)malloc(nproc_neighbor * sizeof(int));
    for(i=0; i<nproc_neighbor; i++) nindex_col_copy[i] = nindex_col[i];

    int **index_col = info->index_col;
    int **index_col_copy = (int**)malloc(nproc_neighbor * sizeof(int*));
    for(i=0; i<nproc_neighbor; i++)
    {
	index_col_copy[i] = (int*)malloc(nindex_col[i] * sizeof(int));
	for(j=0; j<nindex_col[i]; j++) index_col_copy[i][j] = index_col[i][j];
    }

    par_comm_info *info_copy = (par_comm_info*)malloc(sizeof(par_comm_info));
    info_copy->nproc_neighbor = nproc_neighbor;
    info_copy->proc_neighbor  =  proc_neighbor_copy;
    info_copy->index_row      =  index_row_copy;
    info_copy->nindex_row     = nindex_row_copy;
    info_copy->index_col      =  index_col_copy;
    info_copy->nindex_col     = nindex_col_copy;

    return info_copy;
}
#endif
