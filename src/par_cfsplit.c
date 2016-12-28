#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "par_cfsplit.h"
#include "preprocess.h"
#include "amg_param.h"
#include "matrix.h"
#include "par_matrix_vector.h"
#include "cfsplit.h"
#include "linear_algebra.h"
#include "io.h"
#include "tool.h"

#define PAR_CFSPLIT 0

#define  DEBUG_PAR_CLJP      9
#define ASSERT_PAR_CLJP      1

//static int ndof = 0;
extern int  print_rank;

static int nCPT = 0;
static int nFPT = 0;
static int nSPT = 0;

static int Generate_par_strong_coupling_set_negtive (par_dmatcsr *A, par_imatcsr *S, amg_param param);
//static int Generate_par_strong_coupling_set_positive(par_dmatcsr *A, par_imatcsr *S, amg_param param);

//static int Get_par_independent_set_D(par_imatcsr *S, double *measure, int *upt_vec, int nUPT, int *dof);

/* 从全局来说，
 * 假定 A 的稀疏结构对称。生成的邻接矩阵 S 不对称，稀疏结构也不对称。
 * 这是由于强连接不是对称的，就会导致 j 强影响 i (S_ij = 1) 但 i 不强影响 j (S_ji = 0)，
 * 因此也就导致了 S 的稀疏结构不对称。
 *
 * 从编程的角度，本来 S 的稀疏结构与 A 是完全一致的，即将上面所说的 S_ji=0 看做值为 0 的非零元。
 * 这样造成至少不方便：
 * 用 S 中的元素判断时候强连接时，需要考虑 S 中元素的值，即需要判断值是否为 0，才能确定是否强连接。
 * 如果将这些 0 去掉 (对 S 进行压缩)，
 * 就无需考虑 S 元素的值，S->ja 向量中存储的列号 (假定位于第 i 行) 都是强影响 i 的点的编号。
 * 注：Hypre 中 par_strength.c 说也有考虑不压缩 S 的做法, "There are several pros and cons to discuss."
 *
 *  S->diag, S->offd 中列的编号采用与 A->diag, A->offd 的分别相同，
 *  而没有采用类似 A->diag, A->offd 中独立编号的方式，
 *  这是为了在后面考虑影响某点的点的编号时，方便查找。
 *  S->diag->nr = A->diag->nr;
 *  S->diag->nc = A->diag->nc;
 *  S->offd->nr = A->offd->nr;
 *  S->offd->nc = A->offd->nc;
 */
static int Generate_par_strong_coupling_set_negtive (par_dmatcsr *A, par_imatcsr *S, amg_param param)
{
    const double sddt = param.strong_diagonally_dominant_threshold;

    //生成 S->diag
    const int A_diag_nr = A->diag->nr;
    const int A_diag_nc = A->diag->nc;
    const int A_diag_nn = A->diag->nn;
    int      *A_diag_ia = A->diag->ia;
    int      *A_diag_ja = A->diag->ja;
    double   *A_diag_va = A->diag->va;

    const int A_offd_nr = A->offd->nr;
    const int A_offd_nc = A->offd->nc;
    const int A_offd_nn = A->offd->nn;
    int      *A_offd_ia = A->offd->ia;
    int      *A_offd_ja = A->offd->ja;
    double   *A_offd_va = A->offd->va;
    
    assert(A_diag_nr == A_offd_nr);

    int *S_diag_ia = (int *)calloc(A_diag_nr+1, sizeof(int));
    int *S_diag_ja = (int *)calloc(A_diag_nn,   sizeof(int));

    int *S_offd_ia = (int *)calloc(A_offd_nr+1, sizeof(int));
    int *S_offd_ja = (int *)calloc(A_offd_nn,   sizeof(int));
    
    memset(S_diag_ja, -1, A_diag_nn*sizeof(int));
    memset(S_offd_ja, -1, A_offd_nn*sizeof(int));

    int nr = A_diag_nr;
    
    int i, j;
    double row_min;//非对角线最小的负元素
    double row_abs_sum;
    double tol; //强影响阈值
    double aii;
    
    for(i=0; i<nr; i++)
    {
	row_min     = 0.0;
	row_abs_sum = 0.0;
	aii         = 0.0;

	for(j=A_diag_ia[i]; j<A_diag_ia[i+1]; j++)
	{
	    row_abs_sum += MABS(A_diag_va[j]);
	    if(A_diag_ja[j]==i) {aii = A_diag_va[j]; continue;}
	    else row_min = MMIN(row_min, A_diag_va[j]);
	}
	for(j=A_offd_ia[i]; j<A_offd_ia[i+1]; j++)
	{
	    row_abs_sum += MABS(A_offd_va[j]);
	    row_min = MMIN(row_min, A_offd_va[j]);
	}
	tol = row_min * param.strong_connection_threshold;//tol <= 0

	if(row_abs_sum >= (1+sddt)*aii)//考虑对角占优
	{
	    for(j=A_diag_ia[i]; j<A_diag_ia[i+1]; j++)
	    {
		if(A_diag_va[j]<=tol && A_diag_ja[j]!=i)
		    S_diag_ja[j] = A_diag_ja[j];
	    }
	    for(j=A_offd_ia[i]; j<A_offd_ia[i+1]; j++)
	    {
		if(A_offd_va[j] <= tol)
		    S_offd_ja[j] = A_offd_ja[j];
	    }
	}
    }

    int index_diag = 0;
    int index_offd = 0;
    for(i=0; i<nr; i++)
    {
	S_diag_ia[i] = index_diag;
	for(j=A_diag_ia[i]; j<A_diag_ia[i+1]; j++)
	{
	    if(S_diag_ja[j] >= 0)
	    {
		S_diag_ja[index_diag] = S_diag_ja[j];
		index_diag++;
	    }
	}

	S_offd_ia[i] = index_offd;
	for(j=A_offd_ia[i]; j<A_offd_ia[i+1]; j++)
	{
	    if(S_offd_ja[j] >= 0)
	    {
		S_offd_ja[index_offd] = S_offd_ja[j];
		index_offd++;
	    }
	}
    }

    S_diag_ia[A_diag_nr] = index_diag;
    S_diag_ja = (int*)realloc(S_diag_ja, index_diag*sizeof(int));
    int  S_diag_nn = index_diag;
    int *S_diag_va = (int*)calloc(index_diag, sizeof(int));
    for(i=0; i<index_diag; i++) S_diag_va[i] = 1; //为了矩阵打印函数形式的统一

    S_offd_ia[A_offd_nr] = index_offd;
    S_offd_ja = (int*)realloc(S_offd_ja, index_offd*sizeof(int));
    int  S_offd_nn = index_offd;
    int *S_offd_va = (int*)calloc(index_offd, sizeof(int));
    for(i=0; i<index_offd; i++) S_offd_va[i] = 1; //为了矩阵打印函数形式的统一

    imatcsr *S_diag = (imatcsr*)malloc(sizeof(imatcsr));
    S_diag->nr = A_diag_nr;
    S_diag->nc = A_diag_nc;
    S_diag->nn = S_diag_nn;
    S_diag->ia = S_diag_ia;
    S_diag->ja = S_diag_ja;
    S_diag->va = S_diag_va;

    imatcsr *S_offd = (imatcsr*)malloc(sizeof(imatcsr));
    S_offd->nr = A_offd_nr;
    S_offd->nc = A_offd_nc;
    S_offd->nn = S_offd_nn;
    S_offd->ia = S_offd_ia;
    S_offd->ja = S_offd_ja;
    S_offd->va = S_offd_va;

    S->diag             = S_diag;
    S->offd             = S_offd;
    S->nr_global        = -1;
    S->nc_global        = -1;
    S->nn_global        = -1;
    S->map_offd_col_l2g = NULL;
    S->row_start        = NULL;
    S->col_start        = NULL;
    S->comm             = A->comm;
    S->comm_info        = NULL;

    return SUCCESS;
}

/* 1. 与fasp区别在于判断是否强连接时，fasp没有考虑A元素为正的情况，
      而且对于元素大小正好等于tol的时候，fasp默认是不算入强连接 
   2. 加入强对角占优的判断，
      if(row_abs_sum >= (1+max_row_sum)*MABS(aii)) 则不是强对角占优的；
	  若强对角占优，则本行没有强连接关系。
      如果不想判断，就把max_row_sum设为-1.
      fasp的判断条件是
      if(row_sum >= (2 - max_row_sum) * MABS(diag.val[i]) )
      即 (my)max_row_sum 相当于 (fasp)1-max_row_sum
*/
int Generate_par_strong_coupling_set(par_dmatcsr *A, par_imatcsr *S, amg_param param)
{
    if(param.positive_connected == NO)
	return Generate_par_strong_coupling_set_negtive (A, S, param);
    //else if(param.positive_connected == YES)
	//return Generate_par_strong_coupling_set_positive(A, S, param);
    else
	return FAIL;
}

int Split_par_CLJP(par_dmatcsr *A, par_imatcsr *S, par_ivec *dof)
{
    MPI_Comm comm = A->comm;

    int  myrank, nproc_global;
    MPI_Comm_size(comm, &nproc_global);
    MPI_Comm_rank(comm, &myrank);

    int i, j;
    
    imatcsr *S_diag = S->diag;
    imatcsr *S_offd = S->offd;

    /*
     * 某点 i (全局编号) 的影响力，即 ST 第 i 行的行和，亦即 S 第 i 列的列和
     * 因此不考虑将 S 做全局转置的情况下，需要将 S 的所有列求和。
     * 反应到每个进程上，需要将 S 的所有列求和，
     * 然后将 offd 中列和的信息发送，也接收需要的列和信息
     */
    int *S_diag_ja = S_diag->ja;
    int *S_offd_ja = S_offd->ja;
    int  S_diag_nc = S_diag->nc;
    int  S_offd_nc = S_offd->nc;
    int  S_diag_nn = S_diag->nn;
    int  S_offd_nn = S_offd->nn;

    int  A_diag_nc = A->diag->nc;
    par_comm_info *comm_info = A->comm_info;

    /*
     * 排列方式：先存放 S_diag 的列和，再存放 S_offd 的列和
     * S->diag->nc = A->diag->nc;
     * S->offd->nc = A->offd->nc;
     */
    double *measure = (double*)calloc(S_diag_nc+S_offd_nc, sizeof(double));
    for(j=0; j<S_diag_nn; j++) measure[          S_diag_ja[j]] += 1.0;
    for(j=0; j<S_offd_nn; j++) measure[S_diag_nc+S_offd_ja[j]] += 1.0;

    /*
     * 发送：将 S_offd 中的列和发送出去。
     *       需要 A 的进程邻接信息———
     *           全局来看，A 至少要求稀疏结构是对称的，
     *           但一般来说 S 不但不对称，稀疏结构也不对称。
     *       具体地说，设本进程与某进程 t 相邻，需要将 A->offd 在进程 t 上的列和全部发出去
     *
     * 接收：A_offd 关于进程 t 的部分中某行不为 0, 说明全局来看 AT 的该列不为 0.
     *       按照发送的规则，这也说明 ST 的该列不为 0，需要将该列的值加到 ST 在本地该列的列和上.
     *       注意这里列号的判定是基于 A_offd 的行号.
     */
    int nproc_neighbor = comm_info->nproc_neighbor;
    int *nindex_send   = comm_info->nindex_col;
    int **index_send   = comm_info->index_col;
    double **measure_send = (double**)malloc(nproc_neighbor * sizeof(double*));
    for(i=0; i<nproc_neighbor; i++) measure_send[i] = (double*)calloc(nindex_send[i], sizeof(double));
    for(i=0; i<nproc_neighbor; i++)
    {
	for(j=0; j<nindex_send[i]; j++)
	    measure_send[i][j] = measure[S_diag_nc+index_send[i][j]];
    }

    int *nindex_recv   = comm_info->nindex_row;
    int **index_recv   = comm_info->index_row;
    double **measure_recv = (double**)malloc(nproc_neighbor * sizeof(double*));
    for(i=0; i<nproc_neighbor; i++) measure_recv[i] = (double*)calloc(nindex_recv[i], sizeof(double));

    //MPI_Sendrecv(void *sendbuf,int sendcount,MPI_Datatype sendtype,int dest,  int sendtag,
    //             void *recvbuf,int recvcount,MPI_Datatype recvtype,int source,int recvtag,
    //             MPI_Comm comm,MPI_Status *status)
    int *proc_neighbor = comm_info->proc_neighbor;
    for(i=0; i<nproc_neighbor; i++)
    {
	MPI_Sendrecv(measure_send[i], nindex_send[i], MPI_DOUBLE, proc_neighbor[i], myrank+proc_neighbor[i]*1000, 
		     measure_recv[i], nindex_recv[i], MPI_DOUBLE, proc_neighbor[i], proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
    }

#if 0
    if(myrank == print_rank)
    {
	printf("send: \n");
	for(i=0; i<nproc_neighbor; i++)
	{
	    for(j=0; j<nindex_send[i]; j++)
		printf("%02d ", (int)measure_send[i][j]);
	    printf("\n");
	}
	printf("\n");

	printf("recv: \n");
	for(i=0; i<nproc_neighbor; i++)
	{
	    for(j=0; j<nindex_recv[i]; j++)
		printf("%02d ", (int)measure_recv[i][j]);
	    printf("\n");
	}
	printf("\n");
    }
#endif

    for(i=0; i<nproc_neighbor; i++)
    {
	for(j=0; j<nindex_recv[i]; j++)
	    measure[index_recv[i][j]] += measure_recv[i][j];
    }

    for(i=S_diag_nc; i<S_offd_nc; i++) measure[i] = 0.0;

    //通讯完毕之后，再加rand
    srand(1);
    int *A_row_start = A->row_start; 
    for(i=0; i<A_row_start[myrank]; i++) rand(); //保证 measure 的值不随进程数的改变而改变
    for(i=0; i<S_diag_nc; i++) measure[i] += (double)rand()/RAND_MAX;

    //if(myrank == print_rank) Print_dvec(measure, S_diag_nc);

    // measure 初始化完毕，开始进行 cfsplit
    
    //int  A_diag_nc  = A->diag->nc;
    int *S_diag_ia  = S->diag->ia;
    int *S_offd_ia  = S->offd->ia;
    int *cfmarker   = dof->value;
    int  ndof       = dof->length;

    int nUPT = 0;
    int *upt_vec = (int*)calloc(ndof, sizeof(int));
    for(i=0; i<A_diag_nc; i++)
    {
	if((S_diag_ia[i+1]==S_diag_ia[i]) && (S_offd_ia[i+1]==S_offd_ia[i])) 
	    // i 不受别的点影响，即第 i 行非对角线全为 0
	{
	    cfmarker[i] = SPT;
	    measure[i]  = 0.0;
	    nSPT++;
	}
	else if(measure[i] < 1)// i 不影响别的点
	{
	    cfmarker[i] = FPT;
	    measure[i]  = 0.0;
	    nFPT++;
	}
	else // 如果上面两个if语句都不满足，说明i影响一些点，也有一些点影响i，UPT
	{
	    measure[i] = UPT;
	    upt_vec[nUPT++] = i;
	}
    }

#if DEBUG_PAR_CLJP > 5
    int ndof_global = 0;
    int nCPT_global = 0;
    int nFPT_global = 0;
    int nSPT_global = 0;
    int nUPT_global = 0;

    MPI_Allreduce(&ndof, &ndof_global, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&nCPT, &nCPT_global, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&nFPT, &nFPT_global, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&nSPT, &nSPT_global, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&nUPT, &nUPT_global, 1, MPI_INT, MPI_SUM, comm);

    if(myrank == print_rank)
    {
	printf("CLJP INI GLOBAL: ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d, nUPT = %d\n", 
		ndof_global, nCPT_global, nFPT_global, nSPT_global, nUPT_global);
    }
#endif

#if ASSERT_PAR_CLJP
    assert(nCPT+nFPT+nSPT+nUPT == ndof);
#endif

    //???? 有没有可能 S_ext 中的元素全为 0?
    imatcsr *S_ext = Get_S_ext(A, S);

#if 1
    /*
     * A->offd 的所有列中，每个邻居进程的开始位置
     * 待改进：只需要考虑 S->offd 的所有列中，每个邻居进程的开始位置
     *         但是需要向邻居进程通讯三次，第一次为告诉邻居进程需要多少列的信息，
     *         第二次告诉邻居进程需要哪些列的信息，第三次真正传输这些信息。
     *         生成 S_ext 的时候，已经考虑过这些，但这里目前先不采用这种较复杂的方式。
     */
    int *col_start_neighbor = (int*)calloc(nproc_neighbor+1, sizeof(int));
    for(i=0; i<nproc_neighbor; i++) col_start_neighbor[i+1] = col_start_neighbor[i] + comm_info->nindex_col[i];

    int *cfmarker_offd      = (int*)malloc(A->offd->nc * sizeof(int));
    int *cfmarker_diag_send = (int*)malloc(A->diag->nc * sizeof(int));
    int ncfmarker_diag_send = 0;
    int ncfmarker_offd_recv = 0;

    //int iUPT, jS, kS;
    int iwhile = 0;
    while(1)
    {
	if(iwhile >= 1) break;
	iwhile++;

	//进行通信，获得 A->offd 上的所有的 dof 类型
	for(i=0; i<A->offd->nc; i++) cfmarker_offd[i] = EXCEPTION_PT;
	for(i=0; i<nproc_neighbor; i++)
	{

	    ncfmarker_diag_send = comm_info->nindex_row[i];
	    for(j=0; j<ncfmarker_diag_send; j++) cfmarker_diag_send[j] = cfmarker[comm_info->index_row[i][j]];

	    ncfmarker_offd_recv = comm_info->nindex_col[i];
	    MPI_Sendrecv(cfmarker_diag_send,                  ncfmarker_diag_send, MPI_INT, proc_neighbor[i], myrank+proc_neighbor[i]*1000, 
			 cfmarker_offd+col_start_neighbor[i], ncfmarker_offd_recv, MPI_INT, proc_neighbor[i], proc_neighbor[i]+myrank*1000, 
			 comm, MPI_STATUS_IGNORE);
	}
	if(myrank == print_rank)
	{
	    printf("shit cfmarker_offd: ");
	    for(i=0; i<A->offd->nc; i++) printf("%02d ", cfmarker_offd[i]);
	    printf("\n");
	}

	//printf("head   -- iwhile = %d, nUPT = %d\n", iwhile, nUPT);
	//printf("CLJP: ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d, nUPT = %d\n", ndof, nCPT, nFPT, nSPT, nUPT);
	for(iUPT=0; iUPT<nUPT; iUPT++)
	{
	    i = upt_vec[iUPT];
	    //printf("dof[%d] = %d\n", i, dof[i]);
	    if((dof[i]!=CPT) && (lambda_ST[i]<1))
	    {
		dof[i] = FPT;
		//nFPT++;
		//make sure all dependencies have been accounted for
		for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
		{
		    if(S_ja[jS] > -1) 
		    {
			//if(SPT != dof[S_ja[jS]]) dof[i] = UPT;
			dof[i] = UPT;
			//nFPT--;
		    }
		}
	    }

	    if(CPT == dof[i])
	    {
		lambda_ST[i] = 0.0;
		nCPT++;
		nUPT--;
		upt_vec[iUPT] = upt_vec[nUPT];
		upt_vec[nUPT] = i;
		iUPT--;
	    }
	    else if(FPT == dof[i])
	    {
		lambda_ST[i] = 0.0;
		nFPT++;
		nUPT--;
		upt_vec[iUPT] = upt_vec[nUPT];
		upt_vec[nUPT] = i;
		iUPT--;
	    }
	}
    }
#endif

#if 0
    int iUPT, jS, kS;
    int iwhile = 0;
    while(1)
    {
	iwhile++;
	//printf("head   -- iwhile = %d, nUPT = %d\n", iwhile, nUPT);
	//printf("CLJP: ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d, nUPT = %d\n", ndof, nCPT, nFPT, nSPT, nUPT);
	for(iUPT=0; iUPT<nUPT; iUPT++)
	{
	    i = upt_vec[iUPT];
	    //printf("dof[%d] = %d\n", i, dof[i]);
	    if((dof[i]!=CPT) && (lambda_ST[i]<1))
	    {
		dof[i] = FPT;
		//nFPT++;
		//make sure all dependencies have been accounted for
		for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
		{
		    if(S_ja[jS] > -1) 
		    {
			//if(SPT != dof[S_ja[jS]]) dof[i] = UPT;
			dof[i] = UPT;
			//nFPT--;
		    }
		}
	    }

	    if(CPT == dof[i])
	    {
		lambda_ST[i] = 0.0;
		nCPT++;
		nUPT--;
		upt_vec[iUPT] = upt_vec[nUPT];
		upt_vec[nUPT] = i;
		iUPT--;
	    }
	    else if(FPT == dof[i])
	    {
		lambda_ST[i] = 0.0;
		nFPT++;
		nUPT--;
		upt_vec[iUPT] = upt_vec[nUPT];
		upt_vec[nUPT] = i;
		iUPT--;
	    }
	}
	
	//printf("CLJP: ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d, nUPT = %d\n", ndof, nCPT, nFPT, nSPT, nUPT);
	//printf("middle -- iwhile = %d, nUPT = %d\n", iwhile, nUPT);
#if ASSERT_CLJP
	//if(nCPT+nFPT+nSPT+nUPT != ndof)
	//    printf("nCPT = %d, nFPT = %d, nSPT = %d, nUPT = %d, ndof = %d\n\n", nCPT, nFPT, nSPT, nUPT, ndof);
	//else
	//    printf("\n");
	assert(nCPT+nFPT+nSPT+nUPT == ndof);
#endif

	if(nUPT == 0) break;

	if(!Get_independent_set_D(S, ST, lambda_ST, upt_vec, nUPT, dof))
	{
	    printf("ERROR: Cannot find independent set.\n");
	    exit(-1);
	}

	for(iUPT=0; iUPT<nUPT; iUPT++)
	{
	    i = upt_vec[iUPT];
	    //printf("i = %2d, dof[%2d] = %d\n", i, i, dof[i]);
	    if((dof[i]==DPT) || (dof[i]==CPT) || (dof[i]==COMMON_CPT))
	    //if(dof[i] > 0)
	    {
		dof[i] = CPT;
		//nCPT++;
		for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
		{
		    j = S_ja[jS];
		    //printf("    j = %2d\n", j);
		    if(j > -1)
		    {
			S_ja[jS] = -S_ja[jS] - 1;
			if(dof[j] != DPT)
			{
			    lambda_ST[j]--;
			    //printf("        i = %2d, dof[%2d] = %d, lambda_ST[%2d] = %f\n", i, j, dof[j], j, lambda_ST[j]);
			}
		    }
		}
		//printf("\n");
	    }
	    else
	    {
		for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
		{
		    j = S_ja[jS];
		    if(j < 0) j = -j - 1;
		    if((dof[j]==DPT) || (dof[j]==CPT) || (dof[j]==COMMON_CPT))
		    //if(dof[i] > 0)
		    {
			if(S_ja[jS] > -1) S_ja[jS] = -S_ja[jS] - 1;
			dof[j] = COMMON_CPT;
		    }
		    else if(SPT == dof[j])
		    {
			if(S_ja[jS] > -1) S_ja[jS] = -S_ja[jS] - 1;
		    }
		}

		for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
		{
		    j = S_ja[jS];
		    if(j > -1)
		    {
			for(kS=S_ia[j]; kS<S_ia[j+1]; kS++)
			{
			    k = S_ja[kS];
			    if(k < 0) k = -k - 1;
			    if(dof[k] == COMMON_CPT)
			    {
				S_ja[jS] = -S_ja[jS] - 1;
				lambda_ST[j]--;
				break;
			    }
			}
		    }
		}
	    }

	    for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
	    {
		j = S_ja[jS];
		if(j < 0) j = -j - 1;
		if(dof[j] == COMMON_CPT) dof[j] = CPT;
	    }
	}
    }
#endif
    //============================ free ========================================
    free(col_start_neighbor); col_start_neighbor = NULL;
    free(cfmarker_offd); cfmarker_offd = NULL;
    free(cfmarker_diag_send); cfmarker_diag_send = NULL;
    Free_imatcsr(S_ext);

    free(upt_vec); upt_vec = NULL;

    for(i=0; i<nproc_neighbor; i++)
    {
	free(measure_recv[i]);
	measure_recv[i] = NULL;
    }
    free(measure_recv);
    measure_recv = NULL;
    for(i=0; i<nproc_neighbor; i++)
    {
	free(measure_send[i]);
	measure_send[i] = NULL;
    }
    free(measure_send);
    measure_send = NULL;
    free(measure);

    return SUCCESS;
}

imatcsr *Get_S_ext(par_dmatcsr *A, par_imatcsr *S)
{
    /*
     * 1. 通信：自己需要信息的行数，以及这些信息位于进程号的个数(包含本进程)
     * 2. 通信：需要那些行，哪些进程号
     * 3. 计算：遍历需要发送的行，将信息汇总起来
     * 4. 通信：将信息发送出去
     * 5. 计算：将接收的信息汇总，生成 S_ext
     *
     * 在数值算例中发现，有时候本进程与某进程 t 相邻 (假设为第 s 个邻居进程)，
     * 但是并不需要 t 中的数据，即下面 nneed_neighbor[2*s] = 0.
     * 这种情况下，仍会建立 MPI_Sendrecv 通信，通信量为 0, 但未见发生报错的现象。
     *
     * 修改：不将自己的邻居进程号发送出去，而是将 S->offd 的所有列的列号（全局编号）
     *       发送出去。这是因为，自己需要的在邻居进程上的某行，实际上并不需要该行
     *       的全部数据，而只需要 S->diag和S->offd在该行中包含的列 与 
     *       本进程上S->diag和S->offd的全部列的交集。也就是说，只需要那些与本进程相关
     *       的列的数据。对于本进程上S->diag和S->offd的所有列的列号汇总，然后发出。
     * 
     */

    int myrank;
    MPI_Comm comm = A->comm;
    MPI_Comm_rank(comm, &myrank);

    int i, j, k;
    imatcsr *S_offd = S->offd;
    int *S_offd_ja = S_offd->ja;
    int  S_offd_nc = S_offd->nc; // 等于 A_offd_nc
    int  S_offd_nn = S_offd->nn;
    imatcsr *S_diag = S->diag;
    int *S_diag_ja = S_diag->ja;
    int  S_diag_nc = S_diag->nc; // 等于 A_diag_nc
    int  S_diag_nn = S_diag->nn;

    par_comm_info *comm_info = A->comm_info;
    int nproc_neighbor = A->comm_info->nproc_neighbor;
    int *proc_neighbor = A->comm_info->proc_neighbor;
    int *map_offd_col_l2g = A->map_offd_col_l2g;

    // 找到 S_diag 中所有的列号，并按照从小到大排序
    int nc_diag = 0; // actual S->diag->nc
    int *S_diag_col = (int*)malloc(S_diag_nc * sizeof(int));
    for(i=0; i<S_diag_nc; i++) S_diag_col[i] = -1;
    Get_nonnegtive_ivec_all_value_ascend(S_diag_ja, S_diag_nn, S_diag_nc, &nc_diag, S_diag_col);
    // 找到 S_offd 中所有的列号，并按照从小到大排序
    int nc_offd = 0; // actual S->offd->nc
    int *S_offd_col = (int*)malloc(S_offd_nc * sizeof(int));
    for(i=0; i<S_offd_nc; i++) S_offd_col[i] = -1;
    Get_nonnegtive_ivec_all_value_ascend(S_offd_ja, S_offd_nn, S_offd_nc, &nc_offd, S_offd_col);
    // 将 S_diag_col 与 S_offd_col 按全局编号合并
    int  nc = nc_diag + nc_offd;
    int *S_diag_col_global_index = (int*)malloc(S_diag_nc * sizeof(int));
    int col_diag_offset = A->col_start[myrank];
    for(i=0; i<nc_diag; i++) S_diag_col_global_index[i] = S_diag_col[i] + col_diag_offset;
    int *S_offd_col_global_index = (int*)malloc(S_offd_nc * sizeof(int));
    for(i=0; i<nc_offd; i++) S_offd_col_global_index[i] = map_offd_col_l2g[S_offd_col[i]];
    int  position_diag = -1;
    int *S_col_global_index = Insert_ascend_ivec_to_ivec(S_diag_col_global_index, nc_diag, S_offd_col_global_index, nc_offd, &position_diag);

#if 0 //引入了 Insert_ascend_ivec_to_ivec, 下面的代码不需要了
    int *S_col = (int*)malloc(nc * sizeof(int));
    int index_S_col = 0;
    int is_diag_insert = 0;
    int col_diag_offset = A->col_start[myrank];
    for(i=0; i<nc_offd; i++)
    {
	if(0 == is_diag_insert)
	{
	    if(S_diag_col[0]+col_diag_offset < map_offd_col_l2g[S_offd_col[i]])
	    {
		assert(S_diag_col[nc_diag-1]+col_diag_offset < map_offd_col_l2g[S_offd_col[i]]);
		for(j=0; j<nc_diag; j++)
		    S_col[index_S_col++] = S_diag_col[j]+col_diag_offset;
		is_diag_insert = 1;
	    }
	}
	S_col[index_S_col++] = map_offd_col_l2g[S_offd_col[i]];
    }
    if(0 == is_diag_insert)
    {
	for(j=0; j<nc_diag; j++)
	    S_col[index_S_col++] = S_diag_col[j]+col_diag_offset;
    }
    assert(index_S_col == nc);
#endif

#if 1
    if(myrank == print_rank)
    {
	printf("G: S_diag_col = \n");
	printf("             ");
	for(i=0; i<nc_diag; i++)
	    printf("%02d ", i);
	printf("\n");
	printf("             ");
	for(i=0; i<nc_diag; i++)
	    printf("%02d ", S_diag_col[i]+A->col_start[myrank]);
	printf("\n");

	printf("G: S_offd_col = \n");
	printf("             ");
	for(i=0; i<nc_offd; i++)
	    printf("%02d ", i);
	printf("\n");
	printf("             ");
	for(i=0; i<nc_offd; i++)
	    printf("%02d ", map_offd_col_l2g[S_offd_col[i]]);
	printf("\n");
	
	printf("G: S_col      = \n");
	printf("             ");
	for(i=0; i<nc; i++)
	    printf("%02d ", i);
	printf("\n");
	printf("             ");
	for(i=0; i<nc; i++)
	    printf("%02d ", S_col_global_index[i]);
	printf("\n");

	printf("Insert diag position: %d\n", position_diag);
    }
#endif
#if 1
    int sum_S_diag_col = 0;
    int sum_S_offd_col = 0;
    int sum_S_col = 0;
    for(i=0; i<nc_diag; i++) sum_S_diag_col += S_diag_col_global_index[i];
    for(i=0; i<nc_offd; i++) sum_S_offd_col += S_offd_col_global_index[i];
    for(i=0; i<nc;      i++) sum_S_col += S_col_global_index[i];
    //if(myrank == print_rank)
    //	printf("sum_diag = %d, sum_offd = %d, sum_diag_offd = %d, sum_total = %d\n",
    //		sum_S_diag_col, sum_S_offd_col, sum_S_diag_col+sum_S_offd_col, sum_S_col);
    assert(sum_S_col == sum_S_diag_col + sum_S_offd_col);
#endif
    
    // 本进程需要每个邻居进程的行数，列数和进程数
    // 即每个邻居进程上的列数，S的总列数，和进程数+1
    int *nneed_self     = (int*) calloc(nproc_neighbor*3, sizeof(int));
    // 自己需要每个邻居进程的行号（全局编号）
    int **row_need_self = (int**)malloc(nproc_neighbor * sizeof(int*));
    for(i=0; i<nproc_neighbor; i++)
    {
	row_need_self[i] = (int*)malloc(nc_diag+nc_offd * sizeof(int));
	for(j=0; j<nc_offd; j++) row_need_self[i][j] = -1;
    }
    int index = 0;
    int index_row_need_self = 0;
    int index_col_last = -1; // 邻居进程的最后一列列号（局部编号）
    for(i=0; i<nproc_neighbor; i++)
    {
	index_row_need_self = 0;
	// 计算 S->offd 列在邻居进程的数目
	for(j=index; j<nc_offd; j++)
	{
	    index_col_last = comm_info->index_col[i][comm_info->nindex_col[i]-1];
	    if(S_offd_col[j] <= index_col_last)
	    {
		nneed_self[3*i]++;
		row_need_self[i][index_row_need_self++] = map_offd_col_l2g[S_offd_col[j]];
		index++;
	    }
	    else
	    {
		break;
	    }
	}
	// 需要在邻居进程上查找的列数
	nneed_self[3*i+1] = nc;
	// 需要在邻居进程上查找本进程以及本进程邻居上的数据
	nneed_self[3*i+2] = nproc_neighbor + 1;

	assert(index_row_need_self == nneed_self[3*i]);
    }
    int count = 0;
    for(i=0; i<nproc_neighbor; i++) count += nneed_self[3*i];
    assert(count == nc_offd);

    // 接收    ：邻居进程需要的行数，列数和需要的进程数
    // 排列规则：[nrow_need_neighbor, ncol_need_neighbor, nproc_need_neighbor]
    int *nneed_neighbor = (int*)calloc(nproc_neighbor*3, sizeof(int));

    // 将需要信息的行数，以及邻居进程的个数发送出去
    // 同时接收别的进程需要本进程的行数，以及邻居进程的个数
    for(i=0; i<nproc_neighbor; i++)
    {
	MPI_Sendrecv(nneed_self+3*i,     3, MPI_INT, proc_neighbor[i], myrank+proc_neighbor[i]*1000, 
		     nneed_neighbor+3*i, 3, MPI_INT, proc_neighbor[i], proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
    }

    // 接收：邻居需要哪些行
    int **row_need_neighbor = (int**)malloc(nproc_neighbor * sizeof(int*));
    for(i=0; i<nproc_neighbor; i++)
    {
	row_need_neighbor[i] = (int*)malloc(nneed_neighbor[3*i] * sizeof(int));
    }

    /* 
     * 发送自己需要那些行
     * 接收邻居需要哪些行
     */
    for(i=0; i<nproc_neighbor; i++)
    {
	MPI_Sendrecv(row_need_self[i],     nneed_self[3*i],     MPI_INT, proc_neighbor[i], myrank+proc_neighbor[i]*1000, 
		     row_need_neighbor[i], nneed_neighbor[3*i], MPI_INT, proc_neighbor[i], proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
    }

    // 接收：邻居需要哪些列
    int **col_need_neighbor = (int**)malloc(nproc_neighbor * sizeof(int*));
    for(i=0; i<nproc_neighbor; i++)
    {
	col_need_neighbor[i] = (int*)malloc(nneed_neighbor[3*i+1] * sizeof(int));
    }

    /* 
     * 发送自己需要那些列
     * 接收邻居需要哪些列
     */
    for(i=0; i<nproc_neighbor; i++)
    {
	MPI_Sendrecv(S_col_global_index,   nneed_self[3*i+1],     MPI_INT, proc_neighbor[i], myrank+proc_neighbor[i]*1000, 
		     col_need_neighbor[i], nneed_neighbor[3*i+1], MPI_INT, proc_neighbor[i], proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
    }

    // 发送：自己需要哪些进程号（自己进程号，以及邻居进程号）
    int *proc_need_self = Insert_ascend_ivec_to_ivec(&myrank, 1, proc_neighbor, nproc_neighbor, NULL);

#if 0 //引入了 Insert_ascend_ivec_to_ivec, 下面的代码不需要了
    int *proc_need_self = (int*)malloc((nproc_neighbor+1) * sizeof(int));
    is_diag_insert = 0;
    int index_proc_need_self = 0;
    for(i=0; i<nproc_neighbor; i++) 
    {
	if(0 == is_diag_insert)
	{
	    if(myrank < proc_neighbor[i])
	    {
		proc_need_self[index_proc_need_self++] = myrank;
		is_diag_insert = 1;
	    }
	}
	proc_need_self[index_proc_need_self++] = proc_neighbor[i];
    }
    if(0 == is_diag_insert)
    {
	proc_need_self[index_proc_need_self++] = myrank;
    }
    assert(index_proc_need_self == nproc_neighbor + 1);
#endif

    //邻居需要哪些进程号（邻居进程本身，以及其邻居进程号）
    int **proc_need_neighbor = (int**)malloc(nproc_neighbor * sizeof(int*));
    for(i=0; i<nproc_neighbor; i++)
    {
	proc_need_neighbor[i] = (int*)malloc(nneed_neighbor[3*i+2] * sizeof(int));
    }

    /* 
     * 发送自己需要哪些进程号
     * 接收邻居需要哪些进程号
     */
    for(i=0; i<nproc_neighbor; i++)
    {
	MPI_Sendrecv(proc_need_self,        nneed_self[3*i+2],     MPI_INT, proc_neighbor[i], myrank+proc_neighbor[i]*1000, 
		     proc_need_neighbor[i], nneed_neighbor[3*i+2], MPI_INT, proc_neighbor[i], proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
    }

    /*
     * 现在初始化工作已经做好。接下来需要做的是，对每个邻居进程，
     * 将本进程 S->offd, S->diag 中的相关行 (通过 row_need_neighbor, 再将其全局编号转化为本地编号) 找出
     * 然后提取这些行在相关进程 (通过 col_need_neighbor 找出) 的部分。
     * 最后将这个小矩阵的 ia 和 ja 生成，然后发送到邻居进程。
     */
    int *row_start = A->row_start;
    int *col_start = A->col_start;
    int *S_diag_ia = S_diag->ia;
    int *S_offd_ia = S_offd->ia;

    int *nS_ext_col_send = (int*) calloc(nproc_neighbor,  sizeof(int));
    int **S_ext_col_send = (int**)malloc(nproc_neighbor * sizeof(int*));
    for(i=0; i<nproc_neighbor; i++)
    {
	// need_neighbor 排列规则：[nrow_need_neighbor, ncol_need_neighbor, nproc_need_neighbor]
	// 将 row_need_neighbor 的全局编号转化为局部编号
	int nrow = nneed_neighbor[3*i];
	int ncol = nneed_neighbor[3*i+1];
	S_ext_col_send[i] = (int*)calloc(nrow*ncol+nrow, sizeof(int));

	int *row_need_neighbor_local_index = (int*)malloc(nrow * sizeof(int));
	for(j=0; j<nrow; j++)
	    row_need_neighbor_local_index[j] = row_need_neighbor[i][j] - row_start[myrank];

	/*
	 * 不能 col_need_neighbor 的全局编号转化为局部编号，因为不确定本进程上包含所有这些列。
	 * 因此要将本进程 S->offd, S->diag 中的相关行（设共有 nrow 行）提取出来，
	 * 将列转化为全局编号，然后再和 col_need_neighbor 找交集。
	 *
	 * 例如，找S的第i行S_i和 col_need_neighbor 的交集cap时，以 col_need_neighbor 为主，
	 * 即 cap 的长度和 col_need_neighbor 相同（设为 ncol），
	 * 如果第j列（全局编号）位于col_need_neighbor的第k个元素，则 cap[k] = j,
	 * cap 的其余元素都置为 -1。
	 *
	 * 从结果来看，最后输出的是一个稠密矩阵，大小 nrow * ncol
	 * 考虑到矩阵的稀疏性，每行的交集 cap 中很可能有大量元素为 -1,
	 * 这是由于每行 cap 中的元素很少造成的。
	 *
	 * 具体传输信息时，将所有 cap 中的元素的位置记录下来，然后传递。
	 * 此时传递的向量长度为 ncap * 2.
	 *
	 * 但愿每行 cap 中的元素真的很少。
	 */
	for(j=0; j<nrow; j++)
	{
	    int row = row_need_neighbor_local_index[j];

	    int nc_diag_need = S_diag_ia[row+1] - S_diag_ia[row];
	    int *ja_diag_need_global_index = (int*)malloc(nc_diag_need * sizeof(int));
	    for(k=0; k<nc_diag_need; k++)
		ja_diag_need_global_index[k] = S_diag_ja[k+S_diag_ia[row]] + col_start[myrank];

	    int nc_offd_need = S_offd_ia[row+1] - S_offd_ia[row];
	    int *ja_offd_need_global_index = (int*)malloc(nc_offd_need * sizeof(int));
	    for(k=0; k<nc_offd_need; k++)
		ja_offd_need_global_index[k] = map_offd_col_l2g[S_offd_ja[k+S_offd_ia[row]]];

	    int *ja_need_global_index = Insert_ascend_ivec_to_ivec(ja_diag_need_global_index, nc_diag_need, ja_offd_need_global_index, nc_offd_need, NULL);

	    int *cap;
	    int *index_cap;
	    int  ncap;
	    Get_ivec_cap_ivec(ja_need_global_index, nc_diag_need+nc_offd_need, col_need_neighbor[i], ncol, &ncap, &cap, NULL, &index_cap);

	    S_ext_col_send[i][nS_ext_col_send[i]] = ncap;
	    nS_ext_col_send[i] += 1;
	    for(k=0; k<ncap; k++) 
	    {
		S_ext_col_send[i][nS_ext_col_send[i]+k] = index_cap[k];
	    }
	    nS_ext_col_send[i] += ncap;

	    free(index_cap); index_cap = NULL;
	    free(cap); cap = NULL;
	    free(ja_need_global_index); ja_need_global_index = NULL;
	    free(ja_offd_need_global_index); ja_offd_need_global_index = NULL;
	    free(ja_diag_need_global_index); ja_diag_need_global_index = NULL;
	}
	
	free(row_need_neighbor_local_index); row_need_neighbor_local_index = NULL;

	S_ext_col_send[i] = (int*)realloc(S_ext_col_send[i], nS_ext_col_send[i]*sizeof(int));
	if(0 != nS_ext_col_send[i]) assert(NULL != S_ext_col_send[i]);
    }

    /*
     * 发送：邻居进程需要的S_ext每行的列数，以及这些列在对应进程上S_col_global_index的编号
     *
     * 接收：  本进程需要的S_ext每行的列数，以及这些列在本进程上的S_col_global_index的编号
     * 
     * 注意接收的是列在本进程上S_col_global_index的编号，而不是列的全局编号。
     * 这样的好处在于容易得到局部编号。
     * 具体做法是，比较接收到的编号和生成S_col_global_index时diag插入的位置position_diag,
     * 然后就可以知道是否属于diag, 然后就可以计算该编号在S_diag_col或S_offd_col的位置，进而得到局部编号.
     * 
     */
    int *nS_ext_col_recv = (int*)calloc(nproc_neighbor,  sizeof(int));
    for(i=0; i<nproc_neighbor; i++)
    {
	MPI_Sendrecv(&nS_ext_col_send[i], 1, MPI_INT, proc_neighbor[i], myrank+proc_neighbor[i]*1000, 
		     &nS_ext_col_recv[i], 1, MPI_INT, proc_neighbor[i], proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
    }

    //开始对 S_ext 的数据进行通信
    int **S_ext_col_recv = (int**)malloc(nproc_neighbor * sizeof(int*));
    for(i=0; i<nproc_neighbor; i++)
    {
	S_ext_col_recv[i] = (int*)malloc(nS_ext_col_recv[i] * sizeof(int));
	for(j=0; j<nS_ext_col_recv[i]; j++)
	    S_ext_col_recv[i][j] = -1;
    }
    for(i=0; i<nproc_neighbor; i++)
    {
	MPI_Sendrecv(S_ext_col_send[i], nS_ext_col_send[i], MPI_INT, proc_neighbor[i], myrank+proc_neighbor[i]*1000, 
		     S_ext_col_recv[i], nS_ext_col_recv[i], MPI_INT, proc_neighbor[i], proc_neighbor[i]+myrank*1000, 
		     comm, MPI_STATUS_IGNORE);
    }

#if 1
    if(myrank == print_rank)
    {
	Print_imatcsr(S->diag);
	Print_imatcsr(S->offd);
	printf("S->offd->nc_actual = %d", nc_offd);

	printf("\n");
	printf("======================== S info ========================\n");
	printf("S_offd_col: ");
	for(i=0; i<nc_offd; i++) printf("%d ", S_offd_col[i]);
	printf("\n");
	printf("........................................................\n");

	printf("nneed_self: \n");
	for(i=0; i<nproc_neighbor; i++)
	    printf("            %02d: %02d  %02d  %02d\n", proc_neighbor[i], nneed_self[3*i], nneed_self[3*i+1], nneed_self[3*i+2]);

	printf("row_need_self (global index): \n");
	for(i=0; i<nproc_neighbor; i++)
	{
	    printf("            %02d: ", proc_neighbor[i]);
	    for(j=0; j<nneed_self[3*i]; j++)
		printf("%02d ", row_need_self[i][j]);
	    printf("\n");
	}

	printf("col_need_self (global index): \n            ");
	for(i=0; i<nc; i++) printf("%02d ", S_col_global_index[i]);
	printf("\n");

	printf("proc_need_self : \n            ");
	for(i=0; i<nproc_neighbor+1; i++) printf("%02d ", proc_need_self[i]);
	printf("\n");
	printf("........................................................\n");

	printf("nneed_neighbor: \n");
	for(i=0; i<nproc_neighbor; i++)
	    printf("                %02d: %02d  %02d  %02d\n", proc_neighbor[i], nneed_neighbor[3*i], nneed_neighbor[3*i+1], nneed_neighbor[3*i+2]);

	printf("row_need_neighbor (global index): \n");
	for(i=0; i<nproc_neighbor; i++)
	{
	    printf("                %02d: ", proc_neighbor[i]);
	    for(j=0; j<nneed_neighbor[3*i]; j++)
		printf("%02d ", row_need_neighbor[i][j]);
	    printf("\n");
	}

	printf("col_need_neighbor (global index): \n");
	for(i=0; i<nproc_neighbor; i++)
	{
	    printf("                %02d: ", proc_neighbor[i]);
	    for(j=0; j<nneed_neighbor[3*i+1]; j++)
		printf("%02d ", col_need_neighbor[i][j]);
	    printf("\n");
	}

	printf("proc_need_neighbor: \n");
	for(i=0; i<nproc_neighbor; i++)
	{
	    printf("                %02d: ", proc_neighbor[i]);
	    for(j=0; j<nneed_neighbor[3*i+2]; j++)
		printf("%02d ", proc_need_neighbor[i][j]);
	    printf("\n");
	}
	printf("........................................................\n");

	printf("S_ext_col_recv: \n");
	for(i=0; i<nproc_neighbor; i++)
	{
	    printf("                %02d: ", proc_neighbor[i]);
	    for(j=0; j<nS_ext_col_recv[i]; j++)
		printf("%02d ", S_ext_col_recv[i][j]);
	    printf("\n");
	}
	printf("========================================================\n");
	printf("\n");
    }
#endif

    /* 通信完毕。接下来对接收到的数据进行最后的处理，生成本进程上的 S_ext.
     * S_ext 是一个 nc_offd行 nc列 的矩阵
     *
     * 先组装 S_ext ，将列化为局部编号，方法前面已经说明。
     * 然后保持 diag 的部分，将属于 offd 的列号去相反数，以示区分。
     */
    imatcsr *S_ext = (imatcsr*)malloc(sizeof(imatcsr));
    S_ext->nr = nc_offd;
    S_ext->nc = nc;

    int S_ext_nn = 0;
    for(i=0; i<nproc_neighbor; i++) S_ext_nn += nS_ext_col_recv[i];

    S_ext->nn = S_ext_nn;
    S_ext->ia = (int*)calloc(S_ext->nr+1, sizeof(int));
    S_ext->ja = (int*)calloc(S_ext->nn,   sizeof(int));

    int index_S_ext_ja = 0;
    int index_row = 0;
    int index_col_recv = 0;
    int ncol_S_ext = 0;
    for(i=0; i<nproc_neighbor; i++)
    {
	index_col_recv = 0;
	for(j=0; j<nneed_self[3*i]; j++)
	{
	    ncol_S_ext = S_ext_col_recv[i][index_col_recv++];
	    S_ext->ia[index_row+1] = S_ext->ia[index_row] + ncol_S_ext;
	    index_row++;
	    for(k=0; k<ncol_S_ext; k++)
		S_ext->ja[index_S_ext_ja++] = S_ext_col_recv[i][index_col_recv++];
	}
    }

    S_ext->va = (int*)calloc(S_ext->nn, sizeof(int));

    if(myrank == print_rank)
    {
	printf("Print S_ext...\n");
	Print_imatcsr(S_ext);
	Write_imatcsr_csr(S_ext, "../output/S_ext.dat");
    }

    //对 S_ext 的列进行处理得到局部编号，并区分属于diag还是offd
    //得到局部编号的方法：在 map_offd_col_l2g 中查找列号，如果找不到，则说明属于 diag.
    //由于 S_ext 每一行中的列号是从小到大排列的，因此每一行遍历一次 map_offd_col_l2g 即可.
    int  S_ext_nr = S_ext->nr;
    int *S_ext_ia = S_ext->ia;
    int *S_ext_ja = S_ext->ja;
    int  jS_ext;
    if(myrank == print_rank)
    {
	printf("\nmap_offd_col_l2g: \n");
	for(i=0; i<S_offd_nc; i++)
	    printf("%02d ", i);
	printf("\n");
	for(i=0; i<S_offd_nc; i++)
	    printf("%02d ", map_offd_col_l2g[i]);
	printf("\n");
	printf("\nS_ext ja: \n");

	for(i=0; i<S_ext_nr; i++)
	{
	    printf("row %02d:\n          ", i);
	    for(jS_ext=S_ext_ia[i]; jS_ext<S_ext_ia[i+1]; jS_ext++)
		printf("%02d ", jS_ext);
	    printf("\n          ");
	    for(jS_ext=S_ext_ia[i]; jS_ext<S_ext_ia[i+1]; jS_ext++)
		printf("%02d ", S_ext_ja[jS_ext]);
	    printf("\n");
	}
    }
    int insert_diag_left  = position_diag;
    int insert_diag_right = position_diag + nc_diag -1;
    for(k=0; k<S_ext_nn; k++)
    {
	j = S_ext_ja[k];
	if(j < insert_diag_left)
	{
	    assert(j >= 0);
	    S_ext_ja[k] = -j - 1;
	}
	else if(j > insert_diag_right)
	{
	    assert(j < nc);
	    S_ext_ja[k] = -(j-nc_diag) - 1;
	}
	else
	{
	    S_ext_ja[k] = j-insert_diag_left;
	}
    }
    if(myrank == print_rank)
    {
	printf("******************************************************\n");
	for(i=0; i<S_ext_nr; i++)
	{
	    printf("row %02d:\n          ", i);
	    for(jS_ext=S_ext_ia[i]; jS_ext<S_ext_ia[i+1]; jS_ext++)
		printf("%02d ", jS_ext);
	    printf("\n          ");
	    for(jS_ext=S_ext_ia[i]; jS_ext<S_ext_ia[i+1]; jS_ext++)
		printf("%02d ", S_ext_ja[jS_ext]);
	    printf("\n");
	}
    }
#if 0  // 错误的代码，待删除
    int  index_map = 0;
    int  j_global_index;
    int  is_find = 0;
    for(i=0; i<S_ext_nr; i++)
    {
	index_map = 0;
	for(jS_ext=S_ext_ia[i]; jS_ext<S_ext_ia[i+1]; jS_ext++)
	{
	    j_global_index = S_ext_ja[jS_ext];
	    is_find = 0;
	    for(k=index_map; k<S_offd_nc; k++) // S_offd_nc 等于A_offd_nc 等于map_offd_col_l2g的长度
	    {
		if((j_global_index>=col_start[myrank]) && (j_global_index<col_start[myrank+1]))
		{
		    S_ext_ja[jS_ext] -= col_start[myrank];
		    is_find = 1;
		    break;
		}
		else if(j_global_index == map_offd_col_l2g[k])
		{
		    S_ext_ja[jS_ext] = -k;
		    is_find = 1;
		    break;
		}
		index_map++;
	    }
	    //assert(1 == is_find);
	}
    }
#endif

    if(myrank == print_rank)
    {
	Write_imatcsr_csr(S_ext, "../output/S_ext2.dat");
    }


    for(i=0; i<nproc_neighbor; i++)
    {
	free(S_ext_col_recv[i]);     S_ext_col_recv[i]     = NULL;
	free(S_ext_col_send[i]);     S_ext_col_send[i]     = NULL;

	free(proc_need_neighbor[i]); proc_need_neighbor[i] = NULL;
	free(row_need_neighbor[i]);  row_need_neighbor[i]  = NULL;
	free(col_need_neighbor[i]);  col_need_neighbor[i]  = NULL;
	free(row_need_self[i]);      row_need_self[i]      = NULL;
    }

    free(S_ext_col_recv);          S_ext_col_recv          = NULL;
    free(S_ext_col_send);          S_ext_col_send          = NULL;

    free(proc_need_neighbor);      proc_need_neighbor      = NULL;
    free(row_need_neighbor);       row_need_neighbor       = NULL;
    free(col_need_neighbor);       col_need_neighbor       = NULL;
    free(row_need_self);           row_need_self           = NULL;

    free(nS_ext_col_recv);         nS_ext_col_recv         = NULL;
    free(nS_ext_col_send);         nS_ext_col_send         = NULL;
    free(proc_need_self);          proc_need_self          = NULL;
    free(nneed_self);              nneed_self              = NULL;
    free(nneed_neighbor);          nneed_neighbor          = NULL;
    free(S_offd_col);              S_offd_col              = NULL;
    free(S_diag_col);              S_diag_col              = NULL;
    free(S_col_global_index);      S_col_global_index      = NULL;
    free(S_offd_col_global_index); S_offd_col_global_index = NULL;
    free(S_diag_col_global_index); S_diag_col_global_index = NULL;

    return S_ext;
}

int Get_nonnegtive_ivec_all_value_ascend(int *x, int length, int max_val, int *nvalue, int *value)
{
    int i, v, n;
    int *work = (int*)malloc(max_val * sizeof(int));
    for(i=0; i<max_val; i++) work[i] = -1;

    n = 0;
    for(i=0; i<length; i++)
    {
	v = x[i];
	assert(v >= 0);
	if(-1 == work[v])
	{
	    work[v] = v;
	    n++;
	}
    }

    *nvalue = n;

    int index = 0;
    for(i=0; i<max_val; i++)
    {
	if(work[i] > -1)
	    value[index++] = work[i];
    }
    assert(index == n);
    *nvalue = n;

    free(work); work = NULL;

    return n;
}

//假定 a 和 b 中的元素都是从小到大排序的，并且 a 中的全体元素都介于 b 某两个相邻元素之间。
//将 a 插入到 b 中，返回插入后的向量，长度为 length_a+length_b
int *Insert_ascend_ivec_to_ivec(int *a, int length_a, int *b, int length_b, int *position)
{
    assert((length_a!=0) || (length_b!=0));

    int i, j;
    int *c = (int*)malloc((length_a+length_b) * sizeof(int));

    if(0 == length_a)
    {
	for(i=0; i<length_b; i++)
	    c[i] = b[i];
	if(NULL != position) *position = -1;
	return c;
    }
    else if(0 == length_b)
    {
	for(i=0; i<length_a; i++)
	    c[i] = a[i];
	if(NULL != position) *position = -1;
	return c;
    }

    int index_c = 0;
    int is_a_insert = 0;
    for(i=0; i<length_b; i++)
    {
	if(0 == is_a_insert)
	{
	    if(a[0] < b[i])
	    {
		assert(a[length_a-1] < b[i]);
		if(NULL != position) *position = index_c;
		for(j=0; j<length_a; j++) 
		{
		    c[index_c++] = a[j];
		}

		is_a_insert = 1;
	    }
	}
	c[index_c++] = b[i];
    }
    if(0 == is_a_insert) 
    {
	for(j=0; j<length_a; j++) 
	    c[index_c++] = a[j];
	if(NULL != position) *position = length_b;
    }

    assert(index_c == length_a + length_b);
    return c;
}


#if PAR_CFSPLIT

int Split_par_CLJP(par_dmatcsr *A, par_imatcsr *S, par_ivec *dof)
{
    MPI_Comm comm = A->comm;

    int  myrank, nproc_global;
    MPI_Comm_size(comm, &nproc_global);
    MPI_Comm_rank(comm, &myrank);

    int i, j, k;
    
    imatcsr *S_diag = S->diag;
    imatcsr *S_offd = S->offd;

    imatcsr *ST_diag = (imatcsr*)malloc(sizeof(imatcsr));
    imatcsr *ST_offd = (imatcsr*)malloc(sizeof(imatcsr));
    Transpose_imatcsr_struct(S_diag, ST_diag);
    Transpose_imatcsr_struct(S_offd, ST_offd);
    
    int ST_diag_nr = ST_diag->nr;
    int ST_offd_nr = ST_offd->nr;

    //是不是可以考虑 lambda_ST 指向的内存类型改为原本的 int
    //这样通信的长度从 double 变为 int 减少了一半
    //后面再将其指向的内存类型转换为 double 的
    //因为 double 类型只是为了 CLJP 算法中加上一个随机的小数
    double *lambda_ST = (double*)malloc((ST_diag_nr+ST_offd_nr)*sizeof(double));
    for(i=0; i<ST_diag_nr; i++) lambda_ST[i]            = ST_diag_ia[i+1] - ST_diag_ia[i];
    for(i=0; i<ST_offd_nr; i++) lambda_ST[i+ST_diag_nr] = ST_offd_ia[i+1] - ST_offd_ia[i];

    
    //通讯完毕之后，再加rand
    srand(1);
    int *A_row_start = A->row_start; 
    int tmp;
    for(i=0; i<A_row_start[myrank]; i++) tmp = rand();

    double *lambda_ST = (double*)malloc(ST_nr*sizeof(double));
    for(i=0; i<ST_nr; i++) lambda_ST[i] = ST_ia[i+1]-ST_ia[i] + (double)rand()/RAND_MAX;

    int  A_nc  = A->nc;
    int *S_ia  = S->ia;
    int *S_ja  = S->ja;
    int *ST_ia = ST->ia;
    int   S_nn = S->nn;
    

    int nUPT = 0;
    int *upt_vec = (int*)calloc(ndof, sizeof(int));
    for(i=0; i<A_nc; i++)//for(i=0; i<ST_nr; i++)
    {
	/* fasp判断方法：if(S_ia[i+1] == S_ia[i])，参考fasp, coarsening_rs.c, cfsplitting_cls */
	if(S_ia[i+1] == S_ia[i]) /* i不受别的点影响（SPT或CPT），这个条件等价于A的第i行只有对角线为非零元，说明i可以直接求解出来，则i应为SPT.  */
	{
	    dof[i] = SPT;//SPT: 2
	    nSPT++;
	    lambda_ST[i] = 0.0;
	}
	else if(ST_ia[i+1] == ST_ia[i])/* i不影响别的点，但被别的点影响, 则i应为FPT或SPT, 但上一个if语句中判断是不是SPT，如果来到这里，说明为FPT */
	{
	    dof[i] = FPT;//FPT: 1
	    nFPT++;
	    lambda_ST[i] = 0.0;
	}
	else /* 如果上面两个if语句都不满足，说明i影响一些点，也有一些点影响i，UPT */
	{
	    dof[i] = UPT;
	    upt_vec[nUPT++] = i;
	    //nUPT++;
	}
    }
#if DEBUG_CLJP > 5
    printf("CLJP INI : ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d, nUPT = %d\n", ndof, nCPT, nFPT, nSPT, nUPT);
#endif

    int iUPT, jS, kS;
    int iwhile = 0;
    while(1)
    {
	iwhile++;
	//printf("head   -- iwhile = %d, nUPT = %d\n", iwhile, nUPT);
	//printf("CLJP: ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d, nUPT = %d\n", ndof, nCPT, nFPT, nSPT, nUPT);
	for(iUPT=0; iUPT<nUPT; iUPT++)
	{
	    i = upt_vec[iUPT];
	    //printf("dof[%d] = %d\n", i, dof[i]);
	    if((dof[i]!=CPT) && (lambda_ST[i]<1))
	    {
		dof[i] = FPT;
		//nFPT++;
		//make sure all dependencies have been accounted for
		for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
		{
		    if(S_ja[jS] > -1) 
		    {
			//if(SPT != dof[S_ja[jS]]) dof[i] = UPT;
			dof[i] = UPT;
			//nFPT--;
		    }
		}
	    }

	    if(CPT == dof[i])
	    {
		lambda_ST[i] = 0.0;
		nCPT++;
		nUPT--;
		upt_vec[iUPT] = upt_vec[nUPT];
		upt_vec[nUPT] = i;
		iUPT--;
	    }
	    else if(FPT == dof[i])
	    {
		lambda_ST[i] = 0.0;
		nFPT++;
		nUPT--;
		upt_vec[iUPT] = upt_vec[nUPT];
		upt_vec[nUPT] = i;
		iUPT--;
	    }
	}
	
	//printf("CLJP: ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d, nUPT = %d\n", ndof, nCPT, nFPT, nSPT, nUPT);
	//printf("middle -- iwhile = %d, nUPT = %d\n", iwhile, nUPT);
#if ASSERT_CLJP
	//if(nCPT+nFPT+nSPT+nUPT != ndof)
	//    printf("nCPT = %d, nFPT = %d, nSPT = %d, nUPT = %d, ndof = %d\n\n", nCPT, nFPT, nSPT, nUPT, ndof);
	//else
	//    printf("\n");
	assert(nCPT+nFPT+nSPT+nUPT == ndof);
#endif

	if(nUPT == 0) break;

	if(!Get_independent_set_D(S, ST, lambda_ST, upt_vec, nUPT, dof))
	{
	    printf("ERROR: Cannot find independent set.\n");
	    exit(-1);
	}

	for(iUPT=0; iUPT<nUPT; iUPT++)
	{
	    i = upt_vec[iUPT];
	    //printf("i = %2d, dof[%2d] = %d\n", i, i, dof[i]);
	    if((dof[i]==DPT) || (dof[i]==CPT) || (dof[i]==COMMON_CPT))
	    //if(dof[i] > 0)
	    {
		dof[i] = CPT;
		//nCPT++;
		for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
		{
		    j = S_ja[jS];
		    //printf("    j = %2d\n", j);
		    if(j > -1)
		    {
			S_ja[jS] = -S_ja[jS] - 1;
			if(dof[j] != DPT)
			{
			    lambda_ST[j]--;
			    //printf("        i = %2d, dof[%2d] = %d, lambda_ST[%2d] = %f\n", i, j, dof[j], j, lambda_ST[j]);
			}
		    }
		}
		//printf("\n");
	    }
	    else
	    {
		for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
		{
		    j = S_ja[jS];
		    if(j < 0) j = -j - 1;
		    if((dof[j]==DPT) || (dof[j]==CPT) || (dof[j]==COMMON_CPT))
		    //if(dof[i] > 0)
		    {
			if(S_ja[jS] > -1) S_ja[jS] = -S_ja[jS] - 1;
			dof[j] = COMMON_CPT;
		    }
		    else if(SPT == dof[j])
		    {
			if(S_ja[jS] > -1) S_ja[jS] = -S_ja[jS] - 1;
		    }
		}

		for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
		{
		    j = S_ja[jS];
		    if(j > -1)
		    {
			for(kS=S_ia[j]; kS<S_ia[j+1]; kS++)
			{
			    k = S_ja[kS];
			    if(k < 0) k = -k - 1;
			    if(dof[k] == COMMON_CPT)
			    {
				S_ja[jS] = -S_ja[jS] - 1;
				lambda_ST[j]--;
				break;
			    }
			}
		    }
		}
	    }

	    for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
	    {
		j = S_ja[jS];
		if(j < 0) j = -j - 1;
		if(dof[j] == COMMON_CPT) dof[j] = CPT;
	    }
	}
    }

    for(jS=0; jS<S_nn; jS++)
    {
	if(S_ja[jS] < 0) S_ja[jS] = -S_ja[jS] - 1;
    }

    free(lambda_ST);
    Free_imatcsr(ST_diag);
    Free_imatcsr(ST_offd);
    free(upt_vec);

#if DEBUG_CLJP > 5
    printf("CLJP: ndof = %d, nCPT = %d, nFPT = %d, nSPT = %d\n", ndof, nCPT, nFPT, nSPT);
#endif
#if ASSERT_CLJP
    assert(nCPT+nFPT+nSPT == ndof);
#endif
    return nCPT;
}

//
//
//int Generate_par_sparsity_P_dir(par_imatcsr *S, par_ivec *dof, par_dmatcsr *P);
//int Generate_par_P_dir(par_dmatcsr *A, par_imatcsr *S, par_ivec *dof, par_dmatcsr *P);
//
//void Truncate_par_P(par_dmatcsr *P, amg_param param);


static int Get_independent_set_D(imatcsr *S, imatcsr *ST, double *measure, int *upt_vec, int nUPT, int *dof)
{
    int *S_ia = S->ia;
    int *S_ja = S->ja;
    //int *ST_ia = ST->ia;
    //int *ST_ja = ST->ja;

    int iUPT, i, jS, j;
    //int jST;
    for(iUPT=0; iUPT<nUPT; iUPT++)
    {
	i = upt_vec[iUPT];
	if(measure[i] > 1) dof[i] = DPT;
    }

    for(iUPT=0; iUPT<nUPT; iUPT++)
    {
	i = upt_vec[iUPT];
	for(jS=S_ia[i]; jS<S_ia[i+1]; jS++)
	{
	    j = S_ja[jS];
	    if(j < 0) j = -j - 1;

	    if(measure[j] > 1)
	    {
		if(measure[i] > measure[j])
		    dof[j] = NOT_DPT;
		else if(measure[i] < measure[j])
		    dof[i] = NOT_DPT;
	    }
	}

	/*
	for(jST=ST_ia[i]; jST<ST_ia[i+1]; jST++)
	{
	    j = ST_ja[jST];
	    if(j < 0) j = -j - 1;

	    if(measure[j] > 1)
	    {
		if(measure[i] > measure[j])
		    dof[j] = NOT_DPT;
		else if(measure[i] < measure[j])
		    dof[i] = NOT_DPT;
	    }
	}
	*/
    }

    for(iUPT=0; iUPT<nUPT; iUPT++)
    {
	i = upt_vec[iUPT];
	if(DPT == dof[i]) return TRUE;
    }
    return FALSE;
}

void Reset_dof(dmatcsr *A)
{
    nCPT = 0;
    nFPT = 0;
    nSPT = 0;
    ndof = A->nc;
}


static int Generate_strong_coupling_set_positive(dmatcsr *A, imatcsr *S, amg_param param)
{
    const double sddt = param.strong_diagonally_dominant_threshold;

    const int nr = A->nr;
    const int nc = A->nc;
    const int nn = A->nn;
    int      *ia = A->ia;
    int      *ja = A->ja;
    double   *va = A->va;
    
    ndof = nr;
    
    int i, j, index, row_begin, row_end;
    double row_abs_max;//非对角线绝对值最大的元素
    double row_abs_sum;
    double tol; //强影响阈值
    double aii;
    
    S->nr = nr;
    S->nc = nc;
    S->nn = nn;
    S->ia = (int *)calloc(nr+1, sizeof(int));
    S->ja = (int *)calloc(nn,   sizeof(int));
    S->va = NULL;
    
    memset(S->ja, -1, nn*sizeof(int));

    for(i=0; i<nr; i++)
    {
	row_begin = ia[i];
	row_end   = ia[i+1];

	row_abs_max = 0.0;
	row_abs_sum = 0.0;
	aii = 0.0;
	for(j=row_begin; j<row_end; j++)
	{
	    row_abs_sum += MABS(va[j]);
	    if(ja[j]==i) {aii = va[j]; continue;}
	    else         row_abs_max = MMAX(row_abs_max, MABS(va[j]));
	}
	tol = row_abs_max * param.strong_connection_threshold;
	//tol = MMAX(tol, eps);
	
#if DEBUG_GEN_S > 7
	printf("strong connection tolerance[%d] = %f\n", i, tol);
#endif

	if(row_abs_sum >= (1+sddt)*MABS(aii))//考虑对角占优
	{
	    for(j=row_begin; j<row_end; j++)
	    {
		if(MABS(va[j])>=tol && ja[j]!=i)
		    S->ja[j] = ja[j];
	    }
	}
    }

    index = 0;
    for(i=0; i<nr; i++)
    {
	S->ia[i] = index;
	row_begin = ia[i];
	row_end   = ia[i+1];
	for(j=row_begin; j<row_end; j++)
	{
	    if(S->ja[j] >= 0)
	    {
		S->ja[index] = S->ja[j];
		index++;
	    }
	}
    }

    S->ia[nr] = index;
    S->nn = index;
    S->ja = (int*)realloc(S->ja, index*sizeof(int));
    S->va = (int*)calloc(index, sizeof(int));
    for(i=0; i<index; i++) S->va[i] = 1; //为了矩阵打印函数形式的统一
#if 0    
    FILE *file = fopen("../output/S.dat","w");
    fprintf(file, "%d\n", S->nr);
    fprintf(file, "%d\n", S->nc);
    fprintf(file, "%d\n", S->nn);
    fprintf(file, "\n");
    for(i=0; i<S->nr+1; i++) fprintf(file, "%d\n", S->ia[i]);
    fprintf(file, "\n");
    for(i=0; i<S->nn; i++) fprintf(file, "%d\n", S->ja[i]);
    fclose(file);
#endif
    return SUCCESS;
}


int Generate_par_sparsity_P_dir(imatcsr *S, int *dof, dmatcsr *P)
{
    int i, j;
    
    int  nr   = S->nr;
    int *S_ia = S->ia;
    int *S_ja = S->ja;

    P->nr = ndof;
    P->nc = nCPT;
    P->nn = 0;
    P->ia = (int*)calloc(nr+1, sizeof(int));
    
    int  P_nn = 0;
    int *P_ia = P->ia;
    
    for(i=0; i<ndof; i++)
    {
        if(dof[i] == FPT)
        {
            for(j=S_ia[i]; j<S_ia[i+1]; j++)
	        if(dof[S_ja[j]] == CPT) P_nn++;
	}
	else if(dof[i] == CPT)
        {
            P_nn++;
        }
        else if(dof[i] == SPT)
        {
            /* do nothing */
        }
	P_ia[i+1] = P_nn;
    }
    
    P->nn = P_nn;
    P->ja = (int*)   calloc(P_nn, sizeof(int));
    P->va = (double*)calloc(P_nn, sizeof(double));
    
    int    *P_ja = P->ja;
    
    //compute the CPT index map from global to local 
    int *map_F2C = (int*)calloc(ndof, sizeof(int));
    memset(map_F2C, -1, ndof*sizeof(int));

    int CPT_index = 0;
    for(i=0; i<ndof; i++)
    {
        if(dof[i] == CPT) map_F2C[i] = CPT_index++;
    }
    
    int index = 0;
    for(i=0; i<ndof; i++)
    {
        if(dof[i] == FPT)
        {
            for(j=S_ia[i]; j<S_ia[i+1]; j++)
	        if(dof[S_ja[j]] == CPT) P_ja[index++] = map_F2C[S_ja[j]];
	}
	else if(dof[i] == CPT)
        {
            P_ja[index++] = map_F2C[i];
        }
        else if(dof[i] == SPT)
        {
            /* do nothing */
        }
    }
    
    free(map_F2C);
#if ASSERT_GEN_SPA_P
    assert(index == P_nn);
#endif
    return SUCCESS;
}

int Generate_P_dir(dmatcsr *A, imatcsr *S, int *dof, dmatcsr *P)
{
    int i, j, k, m;
    
    int *map_C2F = (int*)calloc(nCPT, sizeof(int));
    for(j=i=0; i<ndof; i++)
    {
        if(dof[i] == CPT) map_C2F[j++] = i;
    }
    
    int    *A_ia = A->ia;
    int    *A_ja = A->ja;
    double *A_va = A->va;
    
    int    *S_ia = S->ia;
    int    *S_ja = S->ja;
    
    int    *P_ia = P->ia;
    int    *P_ja = P->ja;
    double *P_va = P->va;
    
    int *Ci = (int*)malloc(ndof*sizeof(int));
    memset(Ci, -1, ndof*sizeof(int));
    
    double aii = 0.0;
    double spn = 0.0;//sum of positive N
    double snn = 0.0;//sum of negtive  N
    double spp = 0.0;//sum of positive P
    double snp = 0.0;//sum of negtive  P
    double alpha = 0.0;
    double beta  = 0.0;
    int    npc   = 0;//num_positive_coupling
    for(i=0; i<ndof; i++)
    {
        if(dof[i] == FPT)
        {
            aii = spn = snn = spp = snp = alpha = beta = 0.0;
            npc = 0;
            for(j=S_ia[i]; j<S_ia[i+1]; j++)
            {
                k = S_ja[j];
                if(dof[k] == CPT) Ci[k] = i;
            }
            
            for(j=A_ia[i]; j<A_ia[i+1]; j++)
            {
                if(A_ja[j] == i)
                {
                    aii = A_va[j];
                }
                else
                {
                    if(A_va[j] > 0.0)
                    {
                        spn += A_va[j];
                        if(Ci[A_ja[j]] == i)
                        /* whether A_ja[j] \in P_i=C_i^s */
                        {
                            spp += A_va[j];
                            npc++;
                        }
                    }
                    else//if(A_va[j] < 0.0)
                    {
                        snn += A_va[j];
                        if(Ci[A_ja[j]] == i)
                        {
                            snp += A_va[j];
                        }
                    }
                }
            }
            alpha = snn/snp;
            if(npc > 0)
            {
                beta = spn/spp;
            }
            else
            {
                beta = 0.0;
                aii += spn;
            }
            
            for(j=P_ia[i]; j<P_ia[i+1]; j++)
            {
                k = map_C2F[P_ja[j]];
                for(m=A_ia[i]; m<A_ia[i+1]; m++)
                {
                    if(A_ja[m] == k) break;
                }
                if(A_va[m] < 0.0)
                    P_va[j] = -alpha * A_va[m] / aii;
                else
                    P_va[j] = -beta  * A_va[m] / aii;
            }
        }
        else if(dof[i] == CPT)
        {         
            P_va[P_ia[i]] = 1;/* For CPT i, P(i,:) = (0, ..., 0, 1, 0, ..., 0) */
        }
        else if(dof[i] == SPT)
        {
            /* do nothing */
        }
    }

    free(Ci);
    free(map_C2F);
    return SUCCESS;
}

void Truncate_P(dmatcsr *P, amg_param param)
{
    double trunc = param.truncation_threshold;
    if(trunc >= 1.0)
    {
        printf("amg truncation threshold %f >= 1.0!\n", trunc);
        exit(-1);
    }
    
    int nr = P->nr;
    
    int new_nn = 0;
    int new_ja_index = 0;
    int new_va_index = 0;
    
    double minn, maxp;//min negtive/max positive
    double sn, sp; //sum of negtive/positive
    double tsn, tsp;//truncate sum of negtive/postive
    double val;
    double zoomn, zoomp;
    
    int i, j;
    int rb, re;//row begin/end
    
    for(i=0; i<nr; i++)
    {
        minn = maxp = 0.0;
        sn = sp = tsn = tsp = 0.0;
        
        rb = P->ia[i];
        re = P->ia[i+1];
        P->ia[i] = new_nn;
        
        for(j=rb; j<re; j++)
        {
            val = P->va[j];
            
            if(val < 0)
            {
                minn = MMIN(minn, val);
                sn  += val;
            }
            else if(val > 0)
            {
                maxp = MMAX(maxp, val);
                sp  += val;
            }
        }
        
        minn *= trunc;
        maxp *= trunc;
        
        for(j=rb; j<re; j++)
        {
            val = P->va[j];
            
            if(val <= minn)
            {
                new_nn++;
                tsn += val;
                P->ja[new_ja_index++] = P->ja[j];
            }
            else if(val >= maxp)
            {
                new_nn++;
                tsp += val;
                P->ja[new_ja_index++] = P->ja[j];
            }
        }
        
        zoomn = (tsn < -EPS) ? sn/tsn : 1.0;
        zoomp = (tsp >  EPS) ? sp/tsp : 1.0;
        
        for(j=rb; j<re; j++)
        {
            val = P->va[j];
            if     (val <= minn) P->va[new_va_index++] = val * zoomn;
            else if(val >= maxp) P->va[new_va_index++] = val * zoomp;
        }
    }
    
    P->nn     = new_nn;
    P->ia[nr] = new_nn;
    P->ja = (int*)   realloc(P->ja, new_nn*sizeof(int));
    P->va = (double*)realloc(P->va, new_nn*sizeof(double));
}

#endif
