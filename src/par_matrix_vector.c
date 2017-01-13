#ifdef WITH_MPI

#include "par_matrix_vector.h"

#include "io.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

extern int  print_rank;

/*
 * 假定所读取矩阵为对称方阵
 */
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

    //dmatcsr *diag = (dmatcsr*)malloc(sizeof(dmatcsr));
    //dmatcsr *offd = (dmatcsr*)malloc(sizeof(dmatcsr));
    dmatcsr *diag = Init_empty_dmatcsr(col_start[myrank+1]-col_start[myrank]);
    dmatcsr *offd = Init_empty_dmatcsr(col_start[myrank+1]-col_start[myrank]);
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

    if(0 != offd->nn) assert(map_offd_col_l2g != NULL);
    assert(offd->nr == diag->nr);

    //-----------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------
    par_comm_info *comm_info = Init_par_comm_info();
    Get_par_dmatcsr_comm_info(nproc_global, offd, col_start, map_offd_col_l2g, comm_info);

    time_end = MPI_Wtime();
#if 1
    if(myrank == print_rank)
    {
	//Write_dmatcsr_csr(offd, "../output/offd.dat");
	//Write_dmatcsr_csr(diag, "../output/diag.dat");
	printf("\n");
	printf("global info: nr = %d, nc = %d, nn = %d\n", nr_global, nc_global, nn_global);
	//printf("myrank = %d, nr = %d, row_start = %d, row_end = %d\n\n", myrank, row_start[myrank+1]-row_start[myrank], row_start[myrank], row_start[myrank+1]-1);
	printf("read matrix time: %f\n\n", time_end-time_start);
	//Print_dmatcsr(diag);
	//Print_dmatcsr(offd);
	//Print_par_comm_info(comm_info, 5);
	printf("\n");
    }
#endif

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

void Get_dmatcsr_global_size(const char *filename, int *nr, int *nc, int *nn)
{
    FILE *file = My_fopen(filename, "r");
    assert(EOF!=fscanf(file, "%d\n",   nr));
    assert(EOF!=fscanf(file, "%d\n",   nc));
    assert(EOF!=fscanf(file, "%d\n\n", nn));
    fclose(file); file = NULL;
}

void Get_par_dmatcsr_row_start(int nr_global, int nprocs, int *row_start)
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

void Separate_dmatcsr_to_diag_offd(dmatcsr *A, int col_idx_min, int col_idx_max, dmatcsr *diag, dmatcsr *offd)
{
    if(0 == A->nn) return;

    int     nr = A->nr;
    int     nc = A->nc;
    int     nn = A->nn;
    int    *ia = A->ia;
    int    *ja = A->ja;
    double *va = A->va;

    int  i, j, k;
    int  nn_diag = 0;
    int  nn_offd = 0;
    //int *ia_diag = (int*)calloc(nr+1, sizeof(int));
    //int *ia_offd = (int*)calloc(nr+1, sizeof(int));
    int *ia_diag = diag->ia;
    int *ia_offd = offd->ia;
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

    //printf("nn_offd = %d\n", nn_offd);

    int *ja_diag;
    int *ja_offd;
    double *va_diag;
    double *va_offd;

    if(nn_diag > 0)
    {
	ja_diag = (int*)   malloc(nn_diag * sizeof(int));
	va_diag = (double*)malloc(nn_diag * sizeof(double));
    }
    else
    {
	ja_diag = NULL;
	va_diag = NULL;
    }
    if(nn_offd > 0)
    {
	ja_offd = (int*)   malloc(nn_offd * sizeof(int));
	va_offd = (double*)malloc(nn_offd * sizeof(double));
    }
    else
    {
	ja_offd = NULL;
	va_offd = NULL;
    }

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
    //diag->ia = ia_diag;
    diag->ja = ja_diag;
    diag->va = va_diag;

    offd->nr = nr;
    offd->nc = nc;
    offd->nn = nn_offd;
    //offd->ia = ia_offd;
    offd->ja = ja_offd;
    offd->va = va_offd;

    free(isdiag);
}

void Get_par_dmatcsr_diag(dmatcsr *diag, int col_idx_min, int col_idx_max)
{
    int  nn = diag->nn;
    int *ja = diag->ja;
    int  k;
    diag->nc = col_idx_max - col_idx_min + 1;
    for(k=0; k<nn; k++) ja[k] -= col_idx_min;
}

/*
* map_offd_col_l2g 应全部初始化为 -1
* 此时 offd->nc 等于 nc_global
* 最后将 offd->nc 从 nc_global 改为正确的列数
*/
void Get_par_dmatcsr_offd_and_map_col_offd_l2g(dmatcsr*offd, int *map_offd_col_l2g)
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

void Get_par_dmatcsr_comm_info(int nproc_global, dmatcsr *offd, int *col_start, int *map_offd_col_l2g, par_comm_info *comm_info)
{
    if(0 == offd->nn) return;

    int i, k;

    int  nr = offd->nr;
    int  nc = offd->nc;

    int nproc_neighbor = 0;
    int *proc_neighbor = (int*)malloc(nproc_global * sizeof(int));
    for(i=0; i<nproc_global; i++) proc_neighbor[i] = -1;
    Get_neighbor_proc(nproc_global, offd, col_start, map_offd_col_l2g, &nproc_neighbor, proc_neighbor);
    proc_neighbor = (int*)realloc(proc_neighbor, nproc_neighbor*sizeof(int));

    if(nproc_neighbor > 0) 
    {
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
}

/*
* proc_neighbor 应全部初始化为 -1
*/
void Get_neighbor_proc(int nproc_global, dmatcsr *offd, int *col_start, int *map_offd_col_l2g, 
	                      int *nproc_neighbor, int *proc_neighbor)
{
    if(0 == offd->nn)
    {
	*nproc_neighbor = 0;
	return;
    }

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
void Get_par_dmatcsr_comm_row_info(dmatcsr *offd, 
				   int nproc_neighbor, int  *proc_neighbor, 
				   int *col_start,     int  *map_offd_col_l2g, 
				   int *nidx,          int **idx)
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if(0 == offd->nn)
    {
	nidx = 0;
	return;
    }

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
void Get_par_dmatcsr_comm_col_info(dmatcsr *offd, 
	                                  int nproc_neighbor, int *proc_neighbor,
					  int *col_start,     int *map_offd_col_l2g, 
					  int *nidx,          int **idx)
{
    if(0 == offd->nn)
    {
	*nidx = 0;
	return;
    }
    int i, j;

    int nc = offd->nc;
    int proc_mark = 0;
    int col_global_idx;

    for(i=0; i<nc; i++)
    {
	col_global_idx = map_offd_col_l2g[i];
	assert(col_global_idx >= col_start[proc_neighbor[proc_mark]]);

	//if(col_global_idx >= col_start[proc_neighbor[proc_mark]+1]) proc_mark++;
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

	idx[proc_mark][nidx[proc_mark]] = i;
	nidx[proc_mark]++;
    }


    /*
     * 一般情况下，该函数执行到这里有 proc_mark == nproc_neighbor-1;
     * 因为默认矩阵 offd 的列是“压缩”过的，即不含全为0的列。
     *
     * 但有时候矩阵的 proc_neighbor 比真正的要多，
     * 比如生成插值矩阵时 P 的 proc_neighbor 被设置为 与 A 的相同，
     * 就会出现 proc_mark < nproc_neighbor-1
     */

    assert(proc_mark <= nproc_neighbor-1);
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
	if(NULL != info->proc_neighbor)
	{
	    free(info->proc_neighbor);
	    info->proc_neighbor = NULL;
	}

	if(NULL != info->nindex_row)
	{
	    free(info->nindex_row); 
	    info->nindex_row = NULL;
	}

	if(NULL != info->index_row)
	{
	    for(i=0; i<info->nproc_neighbor; i++)
	    {
		free(info->index_row[i]); 
		info->index_row[i] = NULL;
	    }
	    free(info->index_row);
	    info->index_row = NULL;
	}
	if(NULL != info->nindex_col)
	{
	    free(info->nindex_col);
	    info->nindex_col = NULL;
	}

	if(NULL != info->index_col)
	{
	    for(i=0; i<info->nproc_neighbor; i++)
	    {
		free(info->index_col[i]); 
		info->index_col[i] = NULL;
	    }
	    free(info->index_col);
	    info->index_col = NULL;
	}

	info->nproc_neighbor = 0;
	free(info); 
	info = NULL;
    }
}

void Print_par_comm_info(par_comm_info *info, int print_level)
{
    int  nproc_neighbor = info->nproc_neighbor;
    int *proc_neighbor  = info->proc_neighbor;
    int *nindex_row     = info->nindex_row;
    int **index_row     = info->index_row;
    int *nindex_col     = info->nindex_col;
    int **index_col     = info->index_col;

    int i, j;
    if(print_level > 2) printf("----------------------------------\n\n");
    if(print_level > 2) printf("par_comm_info:\n");
    printf("nproc_neighbor  = %d\n", nproc_neighbor);
    printf("proc_neighbor   = ");
    for(i=0; i<nproc_neighbor; i++) printf("%d ", proc_neighbor[i]);
    printf("\n");

    if(print_level > 2)
    {
	if(nindex_row != NULL)
	{
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
	}

	if(nindex_row != NULL)
	{
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
	}
    }
    //if(print_level > 2) printf("----------------------------------\n\n");
}

par_dvec *Init_par_dvec_mv(par_dmatcsr *A)
{
    par_dvec *x = (par_dvec*)malloc(sizeof(par_dvec));
    x->length_global = A->nr_global;
    x->length        = A->diag->nr;
    x->value         = (double*)calloc(x->length, sizeof(double));
    x->comm          = A->comm;

#if 0
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
#endif

    return x;
}

par_ivec *Init_par_ivec_length_comm(int length, MPI_Comm comm)
{
    par_ivec *x = (par_ivec*)malloc(sizeof(par_ivec));
    x->length_global   = -1;
    x->length          = length;
    x->value           = (int*)calloc(x->length, sizeof(int));
    x->comm            = comm;
#if 0
    x->comm_info       = NULL;
    x->send_data       = NULL;
    x->recv_data       = NULL;
    x->recv_data_start = NULL;
#endif

    return x;
}

void Free_par_dvec(par_dvec *x)
{
    if(NULL != x)
    {
	free(x->value); x->value = NULL;

#if 0
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
#endif
	free(x); x = NULL;
    }
}

void Free_par_ivec(par_ivec *x)
{
    if(NULL != x)
    {
	free(x->value); x->value = NULL;

#if 0
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
#endif
	free(x); x = NULL;
    }
}

par_comm_info *Copy_par_comm_info (par_comm_info *info)
{
    par_comm_info *info_copy = Init_par_comm_info();
    if(0 == info->nproc_neighbor) return info_copy;

    int  nproc_neighbor = info->nproc_neighbor;

    int  *proc_neighbor = info->proc_neighbor;

    int i, j;
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

    info_copy->nproc_neighbor = nproc_neighbor;
    info_copy->proc_neighbor  =  proc_neighbor_copy;
    info_copy->index_row      =  index_row_copy;
    info_copy->nindex_row     = nindex_row_copy;
    info_copy->index_col      =  index_col_copy;
    info_copy->nindex_col     = nindex_col_copy;

    return info_copy;
}

par_dmatcsr *Copy_par_dmatcsr(par_dmatcsr *A)
{
    par_dmatcsr *A_copy = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
    A_copy->diag = Copy_dmatcsr(A->diag);
    A_copy->offd = Copy_dmatcsr(A->offd);
    A_copy->nr_global = A->nr_global;
    A_copy->nc_global = A->nc_global;
    A_copy->nn_global = A->nn_global;
    A_copy->comm      = A->comm;
    A_copy->comm_info = Copy_par_comm_info(A->comm_info);

    A_copy->map_offd_col_l2g = NULL;

    int i;
    if(A->offd->nc > 0)
    {
	A_copy->map_offd_col_l2g = (int*)malloc(A->offd->nc * sizeof(int));
	for(i=0; i<A->offd->nc; i++)
	    A_copy->map_offd_col_l2g[i] = A->map_offd_col_l2g[i];
    }

    int  nproc_global;
    MPI_Comm_size(A->comm, &nproc_global);

    A_copy->row_start = (int*)malloc((nproc_global+1) * sizeof(int));
    for(i=0; i<=nproc_global; i++)
	A_copy->row_start[i] = A->row_start[i];
    A_copy->col_start = (int*)malloc((nproc_global+1) * sizeof(int));
    for(i=0; i<=nproc_global; i++)
	A_copy->col_start[i] = A->col_start[i];
    
    return A_copy;
}


void Write_par_dmatcsr_csr(par_dmatcsr *A, const char *filename, int nametype)
{
    int  myrank;
    MPI_Comm_rank(A->comm, &myrank);
    char crank[128];
    if(nametype == 1)
	sprintf(crank, ".rank%02d-(%d,%d)", myrank, A->row_start[myrank], A->row_start[myrank+1]-1);
    else
	sprintf(crank, ".rank%02d", myrank);

    char filename_par[128] = "";
    strcat(filename_par, filename);
    strcat(filename_par, crank);

    FILE *file = fopen(filename_par,"w");
    if(!file)
    {
        printf("\nError: Cannot open %s!\n", filename_par);
	exit(0);
    }
    
    int i, j;
    
    dmatcsr *A_diag = A->diag;
    dmatcsr *A_offd = A->offd;
    
    fprintf(file, "%d\n", A_diag->nr);
    fprintf(file, "%d\n", A->nc_global);
    fprintf(file, "%d\n", A_diag->nn+A_offd->nn);
    fprintf(file, "\n");
    
    int nr = A_diag->nr;
    
    for(i=0; i<=nr; i++) fprintf(file, "%d\n",      A_diag->ia[i] + A_offd->ia[i]);
    fprintf(file, "\n");

    for(i=0; i<nr; i++)
    {
	for(j=A_diag->ia[i]; j<A_diag->ia[i+1]; j++)
	    fprintf(file, "%d\n", A_diag->ja[j] + A->col_start[myrank]);
	for(j=A_offd->ia[i]; j<A_offd->ia[i+1]; j++)
	    fprintf(file, "%d\n", A->map_offd_col_l2g[A_offd->ja[j]]);
    }
    fprintf(file, "\n");
    for(i=0; i<nr; i++)
    {
	for(j=A_diag->ia[i]; j<A_diag->ia[i+1]; j++)
	    fprintf(file, "%15.12f\n", A_diag->va[j]);
	for(j=A_offd->ia[i]; j<A_offd->ia[i+1]; j++)
	    fprintf(file, "%15.12f\n", A_offd->va[j]);
    }
    
    fclose(file);
}

void Write_par_imatcsr_csr(par_imatcsr *A, const char *filename, int nametype)
{
    int  myrank;
    MPI_Comm_rank(A->comm, &myrank);
    char crank[128];
    if(nametype == 1)
	sprintf(crank, ".rank%02d-(%d,%d)", myrank, A->row_start[myrank], A->row_start[myrank+1]-1);
    else
	sprintf(crank, ".rank%02d", myrank);

    char filename_par[128] = "";
    strcat(filename_par, filename);
    strcat(filename_par, crank);

    FILE *file = fopen(filename_par,"w");
    if(!file)
    {
        printf("\nError: Cannot open %s!\n", filename_par);
	exit(0);
    }
    
    int i, j;
    
    imatcsr *A_diag = A->diag;
    imatcsr *A_offd = A->offd;
    
    fprintf(file, "%d\n", A_diag->nr);
    fprintf(file, "%d\n", A->nc_global);
    fprintf(file, "%d\n", A_diag->nn+A_offd->nn);
    fprintf(file, "\n");
    
    int nr = A_diag->nr;
    
    for(i=0; i<=nr; i++) fprintf(file, "%d\n",      A_diag->ia[i] + A_offd->ia[i]);
    fprintf(file, "\n");

    for(i=0; i<nr; i++)
    {
	for(j=A_diag->ia[i]; j<A_diag->ia[i+1]; j++)
	    fprintf(file, "%d\n", A_diag->ja[j] + A->col_start[myrank]);
	for(j=A_offd->ia[i]; j<A_offd->ia[i+1]; j++)
	    fprintf(file, "%d\n", A->map_offd_col_l2g[A_offd->ja[j]]);
    }
    fprintf(file, "\n");
    for(i=0; i<nr; i++)
    {
	for(j=A_diag->ia[i]; j<A_diag->ia[i+1]; j++)
	    fprintf(file, "%d\n", A_diag->va[j]);
	for(j=A_offd->ia[i]; j<A_offd->ia[i+1]; j++)
	    fprintf(file, "%d\n", A_offd->va[j]);
    }
    
    fclose(file);
}

void Print_par_dmatcsr(par_dmatcsr *A, int print_level)
{
    MPI_Comm comm = A->comm;
    int  myrank, nproc_global;
    MPI_Comm_size(comm, &nproc_global);
    MPI_Comm_rank(comm, &myrank);

    if(myrank == 0)
	printf("global info: nr = %d, nc = %d, nn = %d\n\n", A->nr_global, A->nc_global, A->nn_global);

    int i;
    for(i=0; i<nproc_global; i++)
    {
	MPI_Barrier(comm);
	if(myrank == i)
	{
	    printf("---------------------------------------------------------------------------\n");
	    printf("rank = %d, nr = %d from %d to %d, nc = %d from %d to %d\n", 
		    myrank, A->row_start[myrank+1]-A->row_start[myrank], A->row_start[myrank], A->row_start[myrank+1]-1, 
			    A->col_start[myrank+1]-A->col_start[myrank], A->col_start[myrank], A->col_start[myrank+1]-1);
	    Print_par_comm_info(A->comm_info, print_level);
	    //printf("---------------------------------------------------------------------------\n");
	}
	MPI_Barrier(comm);
    }

    MPI_Barrier(comm);
    if(myrank == 0)
	printf("===========================================================================\n");
    MPI_Barrier(comm);
}

dmatcsr *Init_empty_dmatcsr(int nr)
{
    dmatcsr *A = (dmatcsr*)malloc(sizeof(dmatcsr));

    A->nr = nr;
    A->nc = 0;
    A->nn = 0;

    A->ia = (int*)calloc(nr+1, sizeof(int));
    A->ja = NULL;
    A->va = NULL;

    return A;
}
imatcsr *Init_empty_imatcsr(int nr)
{
    imatcsr *A = (imatcsr*)malloc(sizeof(imatcsr));

    A->nr = nr;
    A->nc = 0;
    A->nn = 0;

    A->ia = (int*)calloc(nr+1, sizeof(int));
    A->ja = NULL;
    A->va = NULL;

    return A;
}

par_comm_info *Init_par_comm_info(void)
{
    par_comm_info *comm_info = (par_comm_info*)malloc(sizeof(par_comm_info));

    comm_info->nproc_neighbor = 0;
    comm_info->proc_neighbor  = NULL;

    comm_info->nindex_row = NULL;
    comm_info->index_row  = NULL;

    comm_info->nindex_col = NULL;
    comm_info->index_col  = NULL;

    return comm_info;
}

void Remove_par_dmatcsr_extra_proc_neighbor(par_dmatcsr *A)
{
    if(A->comm_info->nproc_neighbor == 0) return;

    int  myrank;
    MPI_Comm comm = A->comm;
    MPI_Comm_rank(comm, &myrank);

    par_comm_info *comm_info = A->comm_info;
    int nproc_neighbor = comm_info->nproc_neighbor;
    int *proc_neighbor = comm_info->proc_neighbor;

    int i;

    int *proc_neighbor_flag = (int*)malloc(nproc_neighbor * sizeof(int));
    for(i=0; i<nproc_neighbor; i++) proc_neighbor_flag[i] = 1;

    for(i=0; i<nproc_neighbor; i++)
    {
	if(0 == comm_info->nindex_row[i])
	    assert(0 == comm_info->nindex_col[i]);
	if(0 < comm_info->nindex_row[i])
	    assert(0 < comm_info->nindex_col[i]);

	MPI_Sendrecv(&comm_info->nindex_row[i], 1, MPI_INT, proc_neighbor[i], myrank+proc_neighbor[i]*1000, 
		     proc_neighbor_flag+i,      1, MPI_INT, proc_neighbor[i], proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
    }

    
    comm_info->nproc_neighbor = 0;
    for(i=0; i<nproc_neighbor; i++)
    {
	if((comm_info->nindex_row[i]>0) || (proc_neighbor_flag[i]>0))
	{
	    comm_info->proc_neighbor[comm_info->nproc_neighbor] = comm_info->proc_neighbor[i];
	    comm_info->nindex_col[   comm_info->nproc_neighbor] = comm_info->nindex_col[i];
	    comm_info->nindex_row[   comm_info->nproc_neighbor] = comm_info->nindex_row[i];

	    int *tmp;
	    tmp = comm_info->index_col[comm_info->nproc_neighbor];
	    comm_info->index_col[comm_info->nproc_neighbor] = comm_info->index_col[i];
	    comm_info->index_col[i] = tmp;

	    tmp = comm_info->index_row[comm_info->nproc_neighbor];
	    comm_info->index_row[comm_info->nproc_neighbor] = comm_info->index_row[i];
	    comm_info->index_row[i] = tmp;

	    comm_info->nproc_neighbor++;
	}
    }
    

    comm_info->proc_neighbor = (int*)realloc(comm_info->proc_neighbor, comm_info->nproc_neighbor*sizeof(int));
    if(0 != comm_info->nproc_neighbor) assert(NULL != comm_info->proc_neighbor);

    if(0 == comm_info->nproc_neighbor)
    {
	Free_par_comm_info(A->comm_info);
	A->comm_info = Init_par_comm_info();
    }
    else
    {
	for(i=comm_info->nproc_neighbor; i<nproc_neighbor; i++)
	{
	    free(comm_info->index_col[i]);
	    comm_info->index_col[i] = NULL;

	    free(comm_info->index_row[i]);
	    comm_info->index_row[i] = NULL;
	}
    }

    free(proc_neighbor_flag); proc_neighbor_flag = NULL;
}


#endif
