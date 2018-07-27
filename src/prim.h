#ifndef PRIM_H
#define PRIM_H

#include "c_inc.h"

typedef struct heap_node{
	int key;
	int value;
} heap_node;

typedef struct min_heap{
	int size;
	int capacity;
	int * pos;
	heap_node ** heap;
} min_heap;

extern int prim(INC_GRP * meta, MST_GRP * mst);

#endif