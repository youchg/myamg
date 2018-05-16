#CC    = gcc
CC    = mpicc
CFLAG = -g -W -Wall -Wextra -Wfloat-equal -O3 
CFLAG += -DWITH_MPI

SUITE_SPARSE = /home/ycg/Software/SuiteSparse
BLAS         = /home/ycg/Software/BLAS
FASP         = /home/ycg/Software/faspsolver

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
    INC_UMFPACK = -I$(SUITE_SPARSE)/include 
endif

ifdef FASP 
    INC_FASP    = -I$(HOME)/../../faspsolver/base/include
endif

INCEXTERN   = $(INC_UMFPACK) $(INC_FASP)
######################## extern lib #########################
LIB_ARPACK  = $(HOME)/../../ARPACK/libarpack_CentOS.a -lgfortran

ifdef SUITE_SPARSE
    LIB_UMFPACK = -lrt -Wl,-rpath=$(SUITE_SPARSE)/lib -L$(SUITE_SPARSE)/lib -lumfpack -lamd -lsuitesparseconfig -lccolamd -lcamd -lgfortran
    LIB_UMFPACK2 = $(SUITE_SPARSE)/UMFPACK/Lib/libumfpack.a \
                  $(SUITE_SPARSE)/AMD/Lib/libamd.a \
		  $(SUITE_SPARSE)/SuiteSparse_config/libsuitesparseconfig.a \
		  $(SUITE_SPARSE)/CAMD/Lib/libcamd.a \
		  $(SUITE_SPARSE)/CHOLMOD/Lib/libcholmod.a \
		  $(SUITE_SPARSE)/COLAMD/Lib/libcolamd.a \
		  $(SUITE_SPARSE)/CCOLAMD/Lib/libccolamd.a \
		  $(BLAS)/../lapack/liblapack.a \
		  $(BLAS)/blas_LINUX.a \
		  -lgfortran
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

