#CC    = gcc
#CC    = mpicc
CC = /opt/mpi/mvapich2-1.6_gcc-4.6/bin/mpicc
CFLAG = -g -W -Wall -Wextra -Wfloat-equal -O3 
CFLAG += -DWITH_MPI

SUITE_SPARSE = 
BLAS         = 
FASP         = 

#ifeq ($(CC), mpicc)
#    CFLAG += '-cc=gcc' -DWITH_MPI
#endif

ifdef SUITE_SPARSE
    CFLAG += -DWITH_UMFPACK
endif

ifdef FASP
    CFLAG += -DWITH_FASP
endif

CINC = $(HOME)/include
CLIB = $(HOME)/lib
CBIN = $(HOME)/bin
CDAT = $(HOME)/dat

version = 0.3

######################## include ############################
ifdef SUITE_SPARSE
    INC_UMFPACK = -I$(SuiteSparse)/UMFPACK/Include/ \
                  -I$(SuiteSparse)/SuiteSparse_config/ \
		  -IS(SuiteSparse)/AMD/Include
endif
ifdef FASP 
    INC_FASP    = -I$(HOME)/../../faspsolver/base/include
endif

INCEXTERN   = $(INC_UMFPACK) $(INC_FASP)
######################## extern lib #########################
LIB_ARPACK  = $(HOME)/../../ARPACK/libarpack_CentOS.a -lgfortran

ifdef SUITE_SPARSE
    LIB_UMFPACK = $(SuiteSparse)/UMFPACK/Lib/libumfpack.a \
                  $(SuiteSparse)/AMD/Lib/libamd.a \
		  $(SuiteSparse)/SuiteSparse_config/libsuitesparseconfig.a \
		  $(SuiteSparse)/CHOLMOD/Lib/libcholmod.a \
		  $(SuiteSparse)/COLAMD/Lib/libcolamd.a \
		  $(BLAS)/blas_LINUX.a -lgfortran
endif
ifdef FASP 
    LIB_FASP    = $(FASP)/lib/libfasp.a -lgfortran
endif

LIBEXTERN   = -lm $(LIB_UMFPACK) $(LIB_ARPACK) $(LIB_FASP)
#LIBEXTERN   = -lm $(LIB_UMFPACK) $(LIB_ARPACK) $(LIB_FASP) -L/opt/intel/Compiler/11.1/064/lib/intel64 -lirc
#############################################################
LIBAMG = $(CLIB)/libamg-$(version).a

RM = /bin/rm
RMFLAG = -f

AR = ar
ARFLAG = curs
#RANLIB = ranlib

VALGRIND = valgrind
VALGRINDFLAG = #--tool=memcheck --leak-check=full -v --track-origins=yes --show-reachable=yes --read-var-info=yes
