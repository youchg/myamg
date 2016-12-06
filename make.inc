CC    = gcc
#CC    = mpicc
CFLAG = -g -W -Wall -Wextra -Wfloat-equal -O3 

CINC = $(HOME)/include
CLIB = $(HOME)/lib
CBIN = $(HOME)/bin
CDAT = $(HOME)/dat

version = 0.3

######################## include ############################
INC_UMFPACK = -I$(HOME)/../../SuiteSparse/UMFPACK/Include/ \
              -I$(HOME)/../../SuiteSparse/SuiteSparse_config/ \
              -I$(HOME)/../../SuiteSparse/AMD/Include

INC_FASP    = -I$(HOME)/../../faspsolver/base/include

INCEXTERN   = $(INC_UMFPACK) $(INC_FASP)

######################## extern lib #########################
LIB_ARPACK  = $(HOME)/../../ARPACK/libarpack_CentOS.a -lgfortran

LIB_UMFPACK = $(HOME)/../../SuiteSparse/UMFPACK/Lib/libumfpack.a \
              $(HOME)/../../SuiteSparse/AMD/Lib/libamd.a \
              $(HOME)/../../SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a \
              $(HOME)/../../SuiteSparse/CHOLMOD/Lib/libcholmod.a \
              $(HOME)/../../SuiteSparse/COLAMD/Lib/libcolamd.a \
              $(HOME)/../../BLAS-3.5.0/blas_LINUX.a -lgfortran

LIB_FASP    = $(HOME)/../../faspsolver/lib/libfasp.a -lgfortran

LIBEXTERN   = -lm $(LIB_UMFPACK) $(LIB_ARPACK) $(LIB_FASP)
            

#############################################################
LIBAMG = $(CLIB)/libamg-$(version).a

RM = /bin/rm
RMFLAG = -f

AR = ar
ARFLAG = curs
#RANLIB = ranlib

VALGRIND = valgrind
VALGRINDFLAG = #--tool=memcheck --leak-check=full -v --track-origins=yes --show-reachable=yes --read-var-info=yes
