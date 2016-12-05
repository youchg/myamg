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

#endif
