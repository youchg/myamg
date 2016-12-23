#include "par_matrix_vector.h"

#if WITH_MPI

#include "par_linear_algebra.h"
#include "linear_algebra.h"
#include "mpi.h"
#include <math.h>

//必要的时候，检查 A, x, y 的 send_info, recv_info 是否一致
// y = A_diag*x_diag + A_offd*x_recv;
void Multi_par_dmatcsr_dvec(par_dmatcsr *A, par_dvec *x, par_dvec *y)
{
    MPI_Comm comm = A->comm;

    int  myrank;
    MPI_Comm_rank(comm, &myrank);

    int i, j;

    int nproc = A->send_info->nproc;
    int *proc = A->send_info->proc;

    int    **index_send = x->send_info->index;
    int    *nindex_send = x->send_info->nindex;
    double **x_send     = x->send_data;
    for(i=0; i<nproc; i++)
    {
	for(j=0; j<nindex_send[i]; j++)
	    x_send[i][j] = x->value[index_send[i][j]];
    }

    int    *start_recv = x->recv_info->start;
    double *x_recv     = x->recv_data;

    //MPI_Sendrecv(void *sendbuf,int sendcount,MPI_Datatype sendtype,int dest,  int sendtag,
    //             void *recvbuf,int recvcount,MPI_Datatype recvtype,int source,int recvtag,
    //             MPI_Comm comm,MPI_Status *status)
    for(i=0; i<nproc; i++)
    {
	MPI_Sendrecv(x_send[i],            nindex_send[i],                MPI_DOUBLE, proc[i], myrank+proc[i]*1000, 
		     x_recv+start_recv[i], start_recv[i+1]-start_recv[i], MPI_DOUBLE, proc[i], proc[i]+myrank*1000, 
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

#endif
