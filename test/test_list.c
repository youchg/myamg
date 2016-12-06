#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "preprocess.h"
#include "list.h"
#include "tool.h"

void print_dvec(double *dvec, int length);

int main()
{
    int i;
    int length     = 1000000;
    int max_val    = 100;

    //srand(1);
    srand((unsigned)time(NULL));
    double *dvec = (double*)calloc(length, sizeof(double));
    double tb_rand = Get_time();
    for(i=0; i<length; i++) dvec[i] = (double)(rand()%max_val) + (double)rand()/RAND_MAX;
    //print_dvec(dvec, length);
    double te_rand = Get_time();
    printf("rand time: %f\n", te_rand - tb_rand);

    double tb_list = Get_time();
    DNode  *node       = (DNode *)malloc(length  * sizeof(DNode));
    DNode **local_head = (DNode**)malloc(max_val * sizeof(DNode*));
    DNode **local_tail = (DNode**)malloc(max_val * sizeof(DNode*));
    DList list;
    DList_init(&list, node, local_head, local_tail, dvec, length);
    //DList_print(&list, local_head, local_tail);
    
    int delete_index = rand() % length;
    DList_delete(&list, node+delete_index, local_head, local_tail);
    //DList_print(&list, local_head, local_tail);
    double te_list = Get_time();
    printf("list time: %f\n", te_list - tb_list);

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

