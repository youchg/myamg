#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "io.h"
#include "matrix.h"
#include "multigrid.h"
#include "setup_phase.h"
#include "tool.h"

int main(int argc, char* argv[])
{
    dmatcsr *A;
    amg_param param;
    Init_argv(argc, argv, &param, &A, NULL, NULL);
    //Print_amg_param(param);
    multigrid *amg = Build_amg(A, NULL, 15);
    Setup_phase(amg, param);
    Print_amg(amg);
    Free_multigrid(amg);
    Free_dmatcsr(A);
    return 0;
}
