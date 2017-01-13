#ifndef __PREPROCESS_H__
#define __PREPROCESS_H__

/** 
 * \brief dependency relation
 */
#define WITH_BMP     1
//#define WITH_UMFPACK 1

/** 
 * \brief debugging macros
 */
#define DEBUG_LEVEL 0
#define DEBUGk(k, param) { if(DEBUG_LEVEL >= (k)) {printf(param);} }
#define DEBUGm1(param) DEBUGk(-1, param)
#define DEBUG0(param) DEBUGk(0, param)
#define DEBUG1(param) DEBUGk(1, param)
#define DEBUG2(param) DEBUGk(2, param)
#define DEBUG3(param) DEBUGk(3, param)
#define DEBUG4(param) DEBUGk(4, param)
#define DEBUG5(param) DEBUGk(5, param)
#define DEBUG6(param) DEBUGk(6, param)
#define DEBUG7(param) DEBUGk(7, param)
#define DEBUG8(param) DEBUGk(8, param)
#define DEBUG9(param) DEBUGk(9, param)
/* turn off debugging macros */
/* #define DEBUG0(param) */



/** 
 * \brief debugging of printing depth
 */
#define NONE  0
#define BASIC 1 
#define FEW   2
#define MOST  5
#define ALL   9



/** 
 * \brief switchs
 */
#define SUCCESS 1
#define FAIL    0

#define TRUE    1
#define FALSE   0

#define YES     1
#define NO      0



/** 
 * \brief formats for matrix io
 */
#define BMP    1
#define PS     2
#define CSR    3



/** 
 * \brief type of points (dofs) for C/F splitting
 */
#define CPT  1 /**< coarse    grid point */
#define UPT  0 /**< undecided grid point */
#define FPT -2 /**< fine      grid point */
#define SPT -3 /**< solitary  grid point */

#define        DPT  5
#define    NOT_DPT  UPT
#define COMMON_CPT  7

#define EXCEPTION_PT -99

/**
 * \brief coarsening type
 */
#define STD_RS 1
#define CLJP   2

/**
 * \brief interpolation type
 */
#define DIR 1
#define STD 2


/**
 * \brief constant numbers
 */
#define LARGE  9999
#define SMALL -9999
#define EPS    1e-15


/**
 * \brief math definition macros
 */
#define MABS(x)   (((x)>=0) ?(x):-(x))
#define MMAX(x,y) (((x)>(y))?(x): (y))
#define MMIN(x,y) (((x)<(y))?(x): (y))

/**
 * \brief linear solver type
 */
#ifdef WITH_UMFPACK
#define EXTERN_UMFPACK 0
#endif
#define CG 1
#define GS 2




#endif
