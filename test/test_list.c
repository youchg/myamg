#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "preprocess.h"
#include "list.h"

void print_dvec(double *dvec, int length);

int main()
{
    int i;
    int length     = 10;

    //srand(1);
    srand((unsigned)time(NULL));
    double *dvec = (double*)calloc(length, sizeof(double));
    for(i=0; i<length; i++) dvec[i] = (double)(rand()%length) + (double)rand()/RAND_MAX;
    print_dvec(dvec, length);

    DNode  *node       = (DNode *)malloc(length * sizeof(DNode));
    DNode **local_head = (DNode**)malloc(length * 2 * sizeof(DNode*));
    DNode **local_tail = (DNode**)malloc(length * 2 * sizeof(DNode*));
    DList list;
    DList_init(&list, node, local_head, local_tail, dvec, length);
    DList_print(&list, local_head, local_tail);

    free(local_tail);
    free(local_head);
    free(node);
    free(dvec);

    return 0;
}

void print_dvec(double *dvec, int length)
{
    int i;
    printf("\n");
    for(i=0; i<length; i++) printf("[%02d]  %f\n", i, dvec[i]);
    printf("\n");
}

