#ifndef PRIM_H
#define PRIM_H

typedef struct adj_list_node{
	int node_id;
	node * next;
} adj_node;

// Support O(1) insertion
typedef struct adj_list{
	int size;
	node * start;
	node * end;
} adj_list;

typedef struct mst_tree{
	int n;
	adj_list * list;
} tree;

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

extern int prim(float ** adj_mat, int num_sequence, int * max_mst_weight, adj_list ** mst, int * prim_ordering);

#endif