#include "par_matrix_vector.h"

#if WITH_MPI

#include "io.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#if 0
double Get_par_dvec_2norm(par_dvec *x);
void Free_par_dmatcsr(par_dmatcsr *A);
void Free_par_dvec(par_dvec *x);
void Free_par_comm_info(par_comm_info *pcomm_info);
#endif

static void Get_dmatcsr_global_size(const char *filename, int *nr, int *nc, int *nn);
static void Get_par_dmatcsr_row_start(int nr_global, int nprocs, int *row_start);
static void Get_par_dmatcsr_comm_send_info(int nprocs, int *col_start, dmatcsr *offd, int *map_offd_col_l2g, par_comm_info *send_info);
static void Get_par_dmatcsr_offd_and_map_col_offd_l2g(dmatcsr*offd, int *map_offd_col_l2g);
static void Get_par_dmatcsr_diag(dmatcsr *diag, int col_idx_min);


//par_dmatcsr *Read_par_dmatcsr(const char *filename, MPI_Comm comm)
void Read_par_dmatcsr(const char *filename, MPI_Comm comm)
{
    int  myrank, nprocs, myname_len;
    char myname[MPI_MAX_PROCESSOR_NAME];
    //MPI_Status mpi_status;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myrank);
    MPI_Get_processor_name(myname, &myname_len);

    //double time_start, time_end;
    //time_start = MPI_Wtime();
    int i;

    int nr_global, nc_global, nn_global;
    Get_dmatcsr_global_size(filename, &nr_global, &nc_global, &nn_global);

    int *row_start = (int*)calloc(nprocs+1, sizeof(int));
    Get_par_dmatcsr_row_start(nr_global, nprocs, row_start);

    int *col_start = (int*)calloc(nprocs+1, sizeof(int));
    for(i=0; i<=nprocs; i++) col_start[i] = row_start[i];

    assert(row_start[nprocs] == nr_global);
    assert(col_start[nprocs] == nc_global);

    int nr = row_start[myrank+1] - row_start[myrank];


    dmatcsr *A_local = Read_dmatcsr_part(filename, row_start[myrank], row_start[myrank+1]-1);

    dmatcsr *diag = (dmatcsr*)malloc(sizeof(dmatcsr));
    dmatcsr *offd = (dmatcsr*)malloc(sizeof(dmatcsr));
    Separate_dmatcsr_to_diag_offd(A_local, col_start[myrank], col_start[myrank+1]-1, diag, offd);

    Get_par_dmatcsr_diag(diag, col_start[myrank]);

    int *map_offd_col_l2g = (int*)malloc(offd->nn * sizeof(int));
    Get_par_dmatcsr_offd_and_map_col_offd_l2g(offd, map_offd_col_l2g);
    map_offd_col_l2g = (int*)realloc(map_offd_col_l2g, offd->nc*sizeof(int));
    assert(map_offd_col_l2g != NULL);
    //if(myrank == 0) Print_ivec(map_offd_col_l2g, offd->nc);

    Write_dmatcsr_csr(offd, "../output/offd.dat");
    Write_dmatcsr_csr(diag, "../output/diag.dat");
    if(myrank == 0) printf("offd->nc = %d\n", offd->nc);

    if(myrank == 0) printf("global info: nr = %d, nc = %d, nn = %d\n", nr_global, nc_global, nn_global);
    if(myrank == 0) printf("myrank = %d, nr = %d, row_start = %d, row_end = %d\n", myrank, nr, row_start[myrank], row_start[myrank+1]-1);

    par_comm_info *send_info = (par_comm_info*)malloc(sizeof(par_comm_info));
    Get_par_dmatcsr_comm_send_info(nprocs, col_start, offd, map_offd_col_l2g, send_info);
    if(myrank == 0) Print_par_comm_info(send_info);
    Free_par_comm_info(send_info);

    Free_dmatcsr(A_local);

    free(map_offd_col_l2g);
    Free_dmatcsr(diag);
    Free_dmatcsr(offd);

    free(row_start);
    free(col_start);

    //par_dmatcsr *A = (par_dmatcsr*)malloc(sizeof(par_dmatcsr));
    //A->diag = diag;
    //A->offd = offd;
    //A->nr_global = nr_global;
    //A->nc_global = nc_global;
    //A->nn_global = nn_global;
    //A->map_offd_col_l2g = map_offd_col_l2g;
    //A->row_start = row_start;
    //A->col_start = col_start;
    //A->comm = comm;
    //A->send_info = send_info;
    //A->send_data = NULL;
    //A->recv_info = NULL;
    //A->recv_data = NULL;

    //return A;
}

static void Get_par_dmatcsr_comm_send_info(int nprocs, int *col_start, dmatcsr *offd, int *map_offd_col_l2g, par_comm_info *send_info)
{
    int i, j;

    int nsend_proc = 0;

    int *send_proc = (int*)malloc(nprocs * sizeof(int));
    for(i=0; i<nprocs; i++) send_proc[i] = -1;

    int *map_send_start = (int*)calloc(nprocs+1, sizeof(int));

    int  nn = offd->nn;
    int *ja = offd->ja;
    int proc_mark = 0;
    int col_global_idx;
    int count_map_send_start = 0;

    //找到第一个邻居进程号
    col_global_idx = map_offd_col_l2g[ja[0]];
    for(j=proc_mark; j<nprocs; j++)
    {
	if(col_global_idx < col_start[j+1])
	{
	    proc_mark = j;
	    send_proc[nsend_proc] = proc_mark;
	    nsend_proc++;
	    break;
	}
    }

    for(i=0; i<nn; i++)
    {
	col_global_idx = map_offd_col_l2g[ja[i]];
	if(col_global_idx>=col_start[proc_mark] && col_global_idx<col_start[proc_mark+1])
	{
	    count_map_send_start++;
	}
	else
	{
	    for(j=proc_mark+1; j<nprocs; j++)
	    {
		if(col_global_idx < col_start[j+1])
		{
		    proc_mark = j;
		    break;
		}
	    }
	    count_map_send_start++;

	    send_proc[nsend_proc] = proc_mark;
	    map_send_start[nsend_proc] = count_map_send_start;
	    nsend_proc++;
	}
    }
    map_send_start[nsend_proc+1] = count_map_send_start;

    assert(count_map_send_start == nn);

    send_info->nproc = nsend_proc;
    send_info->proc  = send_proc;
    send_info->start = map_send_start;
}

static void Get_par_dmatcsr_offd_and_map_col_offd_l2g(dmatcsr*offd, int *map_offd_col_l2g)
{
    int  nn = offd->nn;
    int  nc = offd->nc;
    int *ja = offd->ja;

    int k;
    for(k=0; k<nn; k++) map_offd_col_l2g[k] = -1;

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

static void Get_par_dmatcsr_diag(dmatcsr *diag, int col_idx_min)
{
    int  nn = diag->nn;
    int *ja = diag->ja;

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

void Separate_dmatcsr_to_diag_offd(dmatcsr *A, int col_idx_min, int col_idx_max, dmatcsr *diag, dmatcsr *offd)
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

    Free_par_comm_info(A->send_info);
    Free_par_comm_data(A->send_data);
    Free_par_comm_info(A->recv_info);
    Free_par_comm_data(A->recv_data);

    free(A); A = NULL;
}

void Free_par_comm_info(par_comm_info *info)
{
    if(NULL != info)
    {
	free(info->proc);  info->proc  = NULL;
	free(info->start); info->start = NULL;
	free(info);        info        = NULL;
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

void Print_par_comm_info(par_comm_info *info)
{
    int  nproc = info->nproc;
    int *proc  = info->proc;
    int *start = info->start;
    int i;
    printf("\n----------------------------------\n");
    printf("par_comm_info:\n");
    printf("nproc = %d\n", info->nproc);
    printf("proc =\n");
    for(i=0; i< nproc; i++) printf("       %d\n", proc[i]);
    printf("start =\n");
    for(i=0; i<=nproc; i++) printf("        %d\n", start[i]);
    printf("----------------------------------\n\n");
}
#endif
