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
    Init_amg_param_argv(argc, argv, &param, &A, NULL, NULL);
    //Print_amg_param(param);
    multigrid *amg = Build_amg(A, NULL, 15);
    double tb_setup = Get_time();
    Setup_phase(amg, param);
    double te_setup = Get_time();
    Print_amg(amg);
    printf("Setup phase time: %f\n", te_setup-tb_setup);
    Free_multigrid(amg);
    Free_dmatcsr(A);
    return 0;
}
