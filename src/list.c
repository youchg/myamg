#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "list.h"

#define ASSERT_LIST 1

//allocate memory outside
void List_init(List *list, Node *n, 
	       Node **local_head, Node **local_tail, 
	       int *vec, int len)
{
    list->nlist = 0;
    list->head  = NULL;
    int i;
    for(i=0; i<len; i++)
    {
	(n+i)->index = i;
	(n+i)->value = vec[i];
	if(vec[i] > 0)
	{
	    List_insert(list, n+i, local_head, local_tail);
	}
    }
}

void List_insert(List *list, Node *n, 
	         Node **local_head, Node **local_tail)
{
    list->nlist++;
    if(NULL == list->head)
    {
	list->head = n;
	n->prev = NULL;
	n->next = NULL;
	local_head[n->value] = n;
	local_tail[n->value] = n;
	return;
    }
    
    Node *q = list->head;
    Node *t;
    while(NULL != q)
    {
	if(n->value > q->value)
	{
	    local_head[n->value] = n;
	    local_tail[n->value] = n;
	    if(NULL != q->prev)
	    {
		q->prev->next = n;
		n->prev = q->prev;
		n->next = q;
		q->prev = n;
		return;
	    }
	    else
	    {
		n->next = q;
		q->prev = n;
		n->prev = NULL;
		list->head = n;
		return;
	    }
	}
	else if(n->value == q->value)
	{
	    Node *tail = local_tail[n->value];
	    local_tail[n->value] = n;
	    if(NULL != tail->next)
	    {
		tail->next->prev = n;
		n->next = tail->next;
		n->prev = tail;
		tail->next = n;
		return;
	    }
	    else
	    {
		n->next = NULL;
		n->prev = tail;
		tail->next = n;
		return;
	    }
	}
	t = local_tail[q->value];
	q = local_tail[q->value]->next;
    }
    local_head[n->value] = n;
    local_tail[n->value] = n;
    t->next = n;
    n->prev = t;
    n->next = NULL;
    return;
}

void List_delete(List *list, Node *n, 
	         Node **local_head, Node **local_tail)
{
    list->nlist--;
    int va = n->value;
    if(n!=local_head[va] && n!=local_tail[va])
    {
#if ASSERT_LIST
	assert(NULL!=n->prev && NULL!=n->next);
#endif
	n->next->prev = n->prev;
	n->prev->next = n->next;
	return;
    }
    else if(n!=local_head[va] && n==local_tail[va])
    {
#if ASSERT_LIST
	assert(NULL != n->prev);
#endif
	local_tail[va] = n->prev;
	if(NULL != n->next)
	{
	    n->next->prev = n->prev;
	    n->prev->next = n->next;
	}
	else
	{
	    n->prev->next = NULL;
	}
	return;
    }
    else if(n==local_head[va] && n!=local_tail[va])
    {
#if ASSERT_LIST
	assert(NULL != n->next);
#endif
	local_head[va] = n->next;
	if(NULL != n->prev)
	{
	    n->next->prev = n->prev;
	    n->prev->next = n->next;
	}
	else
	{
	    n->next->prev = NULL;
	    list->head = n->next;
	}
	return;
    }
    else
    {
#if ASSERT_LIST
	assert(n==local_head[va] && n==local_tail[va]);
#endif
	//pointed out by hongqc, 20160420
	//local_head = NULL;
	//local_tail = NULL;
	if(NULL!=n->prev && NULL!=n->next)
	{
	    n->next->prev = n->prev;
	    n->prev->next = n->next;
	}
	else if(NULL!=n->prev && NULL==n->next)
	{
	    n->prev->next = NULL;
	}
	else if(NULL==n->prev && NULL!=n->next)
	{
	    n->next->prev = NULL;
	    list->head = n->next;
	}
	else
	{
#if ASSERT_LIST
	    assert(NULL==n->prev && NULL==n->next);
#endif
	    list->head = NULL;
	}
	return;
    }
#if ASSERT_LIST
    assert(1 == 0);
#endif
}

void List_print(List *list, Node **local_head, Node **local_tail)
{
    int i=0;
    if(list->head == NULL)
    {
        printf("empty list.\n");
    }
    else
    {
        Node *node = list->head;
        if(node->next == NULL)// single element list
        {
            printf("i = 0, value = %d, index = %d, prev = NULL, next = NULL", node->value, node->index);
            if(local_head[node->value] == node)
                printf(", local head");
            if(local_tail[node->value] == node)
                printf(", local tail");
            printf("\n");
            i++;
        }
        else
        {
            printf("i = 0, value = %d, index = %d, prev = NULL, next = %d",
                   node->value, node->index, node->next->index);
            if(local_head[node->value] == node)
                printf(", local head");
            if(local_tail[node->value] == node)
                printf(", local tail");
            printf("\n");
#if ASSERT_LIST
            assert(node->index == node->next->prev->index);
#endif
            i++;
            node = node->next;
            while(node->next != NULL)// && i<list->nlist)
            {
                printf("i = %d, value = %d, index = %d, prev = %d, next = %d",
                       i++, node->value, node->index, node->prev->index, node->next->index);
                if(local_head[node->value] == node)
                    printf(", local head");
                if(local_tail[node->value] == node)
                    printf(", local tail");
                printf("\n");
#if ASSERT_LIST
                assert(node->index == node->next->prev->index);
                assert(node->index == node->prev->next->index);
#endif
                node = node->next;
            }
            printf("i = %d, value = %d, index = %d, prev = %d, next = NULL",
                   i++, node->value, node->index, node->prev->index);
            if(local_head[node->value] == node)
                printf(", local head");
            if(local_tail[node->value] == node)
                printf(", local tail");
            printf("\n");
#if ASSERT_LIST
            assert(node->index == node->prev->next->index);
#endif
        }
    }
#if ASSERT_LIST
    assert(i==list->nlist);
#endif
}


void DList_init(DList *list, DNode *n, 
	        DNode **local_head, DNode **local_tail, 
	        double *vec, int len)
{
    list->nlist = 0;
    list->head  = NULL;
    int i;
    for(i=0; i<len; i++)
    {
	(n+i)->index = i;
	(n+i)->value = vec[i];
	if(vec[i] > 0)
	{
	    DList_insert(list, n+i, local_head, local_tail);
	}
    }
}

void DList_insert(DList *list, DNode *n, 
	          DNode **local_head, DNode **local_tail)
{
    list->nlist++;
    int ivalue_n = (int)n->value;

    if(NULL == list->head)
    {
	list->head = n;
	n->prev = NULL;
	n->next = NULL;
	local_head[ivalue_n] = n;
	local_tail[ivalue_n] = n;
	return;
    }
    
    DNode *q = list->head;
    int ivalue_q = (int)q->value;
    DNode *t;
    while(NULL != q)
    {
	if(ivalue_n > ivalue_q)
	{
	    local_head[ivalue_n] = n;
	    local_tail[ivalue_n] = n;
	    if(NULL != q->prev)
	    {
		q->prev->next = n;
		n->prev = q->prev;
		n->next = q;
		q->prev = n;
		return;
	    }
	    else
	    {
		n->next = q;
		q->prev = n;
		n->prev = NULL;
		list->head = n;
		return;
	    }
	}
	else if(ivalue_n == ivalue_q)
	{
	    DNode *tail = local_tail[ivalue_n];
	    local_tail[ivalue_n] = n;
	    if(NULL != tail->next)
	    {
		tail->next->prev = n;
		n->next = tail->next;
		n->prev = tail;
		tail->next = n;
		return;
	    }
	    else
	    {
		n->next = NULL;
		n->prev = tail;
		tail->next = n;
		return;
	    }
	}
	t = local_tail[ivalue_q];
	q = local_tail[ivalue_q]->next;
    }
    local_head[ivalue_n] = n;
    local_tail[ivalue_n] = n;
    t->next = n;
    n->prev = t;
    n->next = NULL;
    return;
}

void DList_delete(DList *list, DNode *n, 
	          DNode **local_head, DNode **local_tail)
{
    list->nlist--;
    int ivalue_n = (int)n->value;
    if(n!=local_head[ivalue_n] && n!=local_tail[ivalue_n])
    {
#if ASSERT_LIST
	assert(NULL!=n->prev && NULL!=n->next);
#endif
	n->next->prev = n->prev;
	n->prev->next = n->next;
	return;
    }
    else if(n!=local_head[ivalue_n] && n==local_tail[ivalue_n])
    {
#if ASSERT_LIST
	assert(NULL != n->prev);
#endif
	local_tail[ivalue_n] = n->prev;
	if(NULL != n->next)
	{
	    n->next->prev = n->prev;
	    n->prev->next = n->next;
	}
	else
	{
	    n->prev->next = NULL;
	}
	return;
    }
    else if(n==local_head[ivalue_n] && n!=local_tail[ivalue_n])
    {
#if ASSERT_LIST
	assert(NULL != n->next);
#endif
	local_head[ivalue_n] = n->next;
	if(NULL != n->prev)
	{
	    n->next->prev = n->prev;
	    n->prev->next = n->next;
	}
	else
	{
	    n->next->prev = NULL;
	    list->head = n->next;
	}
	return;
    }
    else
    {
#if ASSERT_LIST
	assert(n==local_head[ivalue_n] && n==local_tail[ivalue_n]);
#endif
	//pointed out by hongqc, 20160420
	//local_head = NULL;
	//local_tail = NULL;
	if(NULL!=n->prev && NULL!=n->next)
	{
	    n->next->prev = n->prev;
	    n->prev->next = n->next;
	}
	else if(NULL!=n->prev && NULL==n->next)
	{
	    n->prev->next = NULL;
	}
	else if(NULL==n->prev && NULL!=n->next)
	{
	    n->next->prev = NULL;
	    list->head = n->next;
	}
	else
	{
#if ASSERT_LIST
	    assert(NULL==n->prev && NULL==n->next);
#endif
	    list->head = NULL;
	}
	return;
    }
#if ASSERT_LIST
    assert(1 == 0);
#endif
}

void DList_print(DList *list, DNode **local_head, DNode **local_tail)
{
    int i=0;
    if(list->head == NULL)
    {
        printf("empty list.\n");
    }
    else
    {
        DNode *node = list->head;
	int ivalue_node = (int)node->value;
        if(node->next == NULL)// single element list
        {
            printf("i = 0, value = %f, index = %d, prev = NULL, next = NULL", node->value, node->index);
            if(local_head[ivalue_node] == node)
                printf(", local head");
            if(local_tail[ivalue_node] == node)
                printf(", local tail");
            printf("\n");
            i++;
        }
        else
        {
            printf("i = 0, value = %f, index = %d, prev = NULL, next = %d",
                   node->value, node->index, node->next->index);
            if(local_head[ivalue_node] == node) printf(", local head");
            if(local_tail[ivalue_node] == node) printf(", local tail");
            printf("\n");
#if ASSERT_LIST
            assert(node->index == node->next->prev->index);
#endif
            i++;
            node = node->next;
            while(node->next != NULL)// && i<list->nlist)
            {
                printf("i = %d, value = %f, index = %d, prev = %d, next = %d",
                       i++, node->value, node->index, node->prev->index, node->next->index);
                if(local_head[ivalue_node] == node) printf(", local head");
                if(local_tail[ivalue_node] == node) printf(", local tail");
                printf("\n");
#if ASSERT_LIST
                assert(node->index == node->next->prev->index);
                assert(node->index == node->prev->next->index);
#endif
                node = node->next;
            }
            printf("i = %d, value = %f, index = %d, prev = %d, next = NULL",
                   i++, node->value, node->index, node->prev->index);
            if(local_head[ivalue_node] == node) printf(", local head");
            if(local_tail[ivalue_node] == node) printf(", local tail");
            printf("\n");
#if ASSERT_LIST
            assert(node->index == node->prev->next->index);
#endif
        }
    }
#if ASSERT_LIST
    assert(i==list->nlist);
#endif
}
