#ifndef __DOUBLE_LINK_LIST_H__
#define __DOUBLE_LINK_LIST_H__

typedef struct NODE
{
    int index;
    int value;
    struct NODE *prev;
    struct NODE *next;
} Node;

typedef struct LIST
{
    int nlist;
    Node *head;
} List;

void List_init(List *list, Node *n, 
	       Node **local_head, Node **local_tail, 
	       int *vec, int len);

void List_insert(List *list, Node *n, 
	         Node **local_head, Node **local_tail);

void List_delete(List *list, Node *n, 
	         Node **local_head, Node **local_tail);

void List_print(List *list, Node **local_head, Node **local_tail);

typedef struct NODE_DOUBLE
{
    int index;
    double value;
    struct NODE_DOUBLE *prev;
    struct NODE_DOUBLE *next;
} DNode;

typedef struct LIST_DOUBLE
{
    int nlist;
    DNode *head;
} DList;

void DList_init(DList *list, DNode *n, 
	        DNode **local_head, DNode **local_tail, 
	        double *vec, int len);

void DList_insert(DList *list, DNode *n, 
	          DNode **local_head, DNode **local_tail);

void DList_delete(DList *list, DNode *n, 
	          DNode **local_head, DNode **local_tail);

void DList_print(DList *list, DNode **local_head, DNode **local_tail);

#endif
