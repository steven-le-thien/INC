#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prim.h"
#include "utilities.h"
#include "options.h"

// Helper macro (for this file only)
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

// Internal functions

heap_node * 	heap_node_constructor(int value, int key);
min_heap * 		heap_constructor(int c);
heap_node * 	pop_heap(min_heap * heap);
int 			update(min_heap * heap, int key, int new_value);
int 			swap_node(int a, int b, min_heap * heap);
int 			heapify(min_heap * heap, int i);

/* Construct a MST using heap implementation of Prim's algorithm
 * Input:       meta        meta variables, including the distance matrix (adjacency matrix) 
 *              mst         mst variables, to store the return value
 * Output: 0 on success, ERROR otherwise
 * Effect: allocate some memories, build some trees, open some files, init 2 arrays in map and 1 in meta 
 */
int prim(INC_GRP * meta, MST_GRP * mst){
	// Heap variables
	min_heap * heap;
	heap_node * head;

	// Loop variables
	int i, j; 
	int prim_ordering_counter;
	int cur_key;

	// Temporary stack array
	int v_arr[meta->n_taxa];

	// Mst initialization
	mst->n_taxa 			= meta->n_taxa; 
	mst->max_w 				= 0.0;
	mst->prim_ord 			= malloc(mst->n_taxa * sizeof(int));
	mst->prim_par 			= malloc(mst->n_taxa * sizeof(int)); 

	if(!prim_ord || !prim_par) 		PRINT_AND_RETURN("malloc failed to init mst in prim\n", MALLOC ERROR);
	for(i = 0; i < n_taxa; i ++){
		mst->prim_ord = -1;
		mst->prim_par = -1;
	}

	// Heap initialization
	heap = heap_constructor(mst->n_taxa)  	
	if(!heap) 						PRINT_AND_RETURN("heap initialization failed in prim\n", MALLOC_ERROR);
	for(i = 0; i < mst->n_taxa; i++){
		v_arr[i] = INT_MAX; // large number
		get_node(i) 	= heap_node_constructor(v_arr[i], i);
		if(!get_node(i)) 			PRINT_AND_RETURN("node initialization failed in prim\n", MALLOC_ERROR);
		get_pos(i)		= i;
	}

	// Loop initialization
	prim_ordering_counter = 0;
	value_array[0] = 0;
	heap->size = mst->n_taxa;

	// Prim's MST
	while(heap->size){
		head = pop_heap(heap);
		if(!head) 					PRINT_AND_RETURN("heap popping failed in prim\n", GENERAL_ERROR);
		cur_key = head->key;

		for(i = 0; i < mst->n_taxa; i++){
			if(in_heap(i) && meta->d[cur_key][i] < v_arr[i]){
				v_arr[i] = meta->d[cur_key][i];
				mst_prim_par[i] = cur_key;
				if(update(heap, i, value_array[i]) != SUCCESS) 
									PRINT_AND_RETURN("updating heap value failed in prim\n", GENERAL_ERROR);
				j = i;
			}
		}
		mst->max_w = max(mst->max_w, meta->d[cur_key][j]); 
		mst->prim_ord[prim_ordering_counter++]; 
	}

	return 0;
}

// INTERNAL IMPLEMENTATION
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

int swap_node(int a, int b, min_heap * heap){
	if(!heap) PRINT_AND_RETURN("heap is NULL in swap node\n", GENERAL_ERROR);
	heap_node * tmp = get_node(a);
	get_node(a) = get_node(b);
	get_node(b) = tmp;

	get_pos(get_key(a)) = b;
	get_pos(get_key(b)) = a;

	return 0;
}

// Heapify down implementation
int heapify(min_heap * heap, int i){
	if(!heap) PRINT_AND_RETURN("heap is NULL in swap node\n", GENERAL_ERROR);

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
		if(swap_node(i, min_index, heap) != SUCCESS) PRINT_AND_RETURN("swap node failed in heapify\n", GENERAL_ERROR);
		heapify(heap, min_index);
	}

	return 0;
}


heap_node * pop_heap(min_heap * heap){
	if(!heap->size) return NULL;

	// Swap the first and the last node, decrease the size of the heap and heapify the new root if nec.
	if(swap_node(root_idx, last_node_index, heap) 	!= SUCCESS) PRINT_AND_RETURN("swap node failed in pop_head\n", NULL);
	heap->size--;
	if(heapify(heap, root_idx) 						!= SUCCESS) PRINT_AND_RETURN("heapify faield in pop_head\n", NULL);

	return get_node(last_node_index); // this is the old node
}

int update(min_heap * heap, int key, int new_value){
	int idx = get_pos(key);
	get_value(idx) = new_value;

	// Maintain heap property by explicit heapify up
	while(!is_root(idx) && get_value(idx) < get_value(get_parent(idx))){
		if(swap_node(idx, get_parent(idx)) != SUCCESS) PRINT_AND_RETURN("swap node failed in update in prim\n", GENERAL_ERROR);
		idx = get_parent(idx);
	}

	return 0;
}
