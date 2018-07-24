#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prim.h"
#include "utilities.h"
#include "options.h"

#define root_idx 			0
#define last_node_index 	(heap->size - 1)
#define left_child(i) 		((i << 1) + 1)
#define right_child(i)		((i << 1) + 2)
#define get_parent(i) 		((i - 1) >> 1)
#define is_root(i) 			(i == root_idx)
#define get_node(i)			(heap->heap[i])
#define get_pos(i) 			(heap->pos[i])
#define get_value(i)		(get_node(i)->value)
#define get_key(i) 			(get_node(i)->key)
#define in_heap(i)			(heap->pos[i] < heap->size)

heap_node * heap_node_constructor(int value, int key){
	heap_node * node = malloc(sizeof(heap_node));

	// Check malloc error

	node->value 	= value;
	node->key 		= key;

	return node;
}

min_heap * heap_constructor(int c){
	int i;

	min_heap * heap = malloc(sizeof(min_heap));

	heap->size 			= 0;
	heap->capacity 		= c;
	heap->pos 			= malloc(c * sizeof(int));
	heap->heap 			= malloc(c * sizeof(heap_node *));

	for(i = 0; i < c; i++){
		heap->heap[i] 	= malloc(sizeof(heap_node));
	}

	return heap;
}

void swap_node(int a, int b, min_heap * pos){
	heap_node * tmp = get_node(a);
	get_node(a) = get_node(b);
	get_node(b) = tmp;

	get_pos(get_key(a)) = b;
	get_pos(get_key(b)) = a;
}

// Heapify down implementation
void heapify(min_heap * heap, int i){
	// Find the key with the smallest value between i and its 2 child (if exists)
	int min_index = i;

	if(left_child(i) < heap->size)
		if(get_value(left_child(i)) < get_value(min_index))
			min_index = left_child(i);

	if(right_child(i) < heap->size)
		if(get_value(right_child(i)) < get_value(min_index))
			min_index = right_child(i);

	// Check if heap property is violated
	if(min_index != i){
		// Swap node and recurse
		swap_node(i, min_index, heap);
		heapify(heap, min_index);
	}
}


heap_node * pop_heap(min_heap * heap){
	if(!heap->size) return NULL;

	// Swap the first and the last node, decrease the size of the heap and heapify the new root if nec.
	swap_node(root_idx, last_node_index, heap);
	heap->size--;
	heapify(heap, root_idx);

	return get_node(last_node_index); // this is the old node
}

void update(min_heap * heap, int key, int new_value){
	int idx = get_pos(key);
	get_value(idx) = new_value;

	// Maintain heap property by explicit heapify up
	while(!is_root(idx) && get_value(idx) < get_value(get_parent(idx))){
		swap_node(idx, get_parent(idx));
		idx = get_parent(idx);
	}
}

adj_list ** adj_list_graph_constructor(int n){
	int i;
	adj_list ** graph;

	graph = malloc(sizeof(adj_list *) * n);
	for(i = 0; i < n; i++){
		graph[i] = malloc(sizeof(adj_list));

		graph[i]->size 	= 0;
		graph[i]->start = NULL;
		graph[i]->end 	= NULL;
	}

	return graph;
}

int * prim_odering_constructor(int n){
	int i;
	int * order = malloc(n * sizeof(int));

	for(i = 0; i < n; i++){
		order[i] = -1; 
	}
	return order;
}

int prim(float ** adj_mat, int num_sequence, float * max_mst_weight, int * prim_ordering, int * adj_in_mst){
	min_heap * heap;
	heap_node * head;
	int i; 
	int j;
	int to_be_adjacent;
	int prim_ordering_counter;
	int cur_key;
	float cur_weight;

	int value_array[num_sequence];
	int parent_array[num_sequence];

	*max_mst_weight = 0.0;
	prim_ordering = prim_odering_constructor(num_sequence);
	adj_in_mst = prim_odering_constructor(num_sequence);

	prim_ordering_counter = 0;

	heap = heap_constructor(num_sequence);
	for(i = 0; i < num_sequence; i++){
		value_array[i] = 1000000; // large number
		get_node(i) 	= heap_node_constructor(value_array[i], i);
		get_pos(i)		= i;
	}

	// Initializing the heap
	value_array[0] = 0;
	heap->size = num_sequence;

	// Prim's MST
	while(heap->size){
		head = pop_heap(heap);
		cur_key = head->key;

		for(i = 0; i < num_sequence; i++){
			if(in_heap(i) && adj_mat[cur_key][i] < value_array[i]){
				value_array[i] = adj_mat[cur_key][i];
				parent_array[i] = cur_key;
				update(heap, i, value_array[i]);
				j = i;
			}
		}
		*max_mst_weight = (*max_mst_weight) < adj_mat[cur_key][j] ? adj_mat[cur_key][j] : (*max_mst_weight);
		prim_ordering[prim_ordering_counter] = j;
		adj_in_mst[prim_ordering_counter] = parent[j];

		prim_ordering_counter++; 
	}
}