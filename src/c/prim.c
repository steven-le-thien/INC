// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prim.h"
#include "utilities.h"
#include "options.h"
#include "dist.h"

// Helper macro (for this file only)
#define root_idx            0
#define last_node_index     (heap->size - 1)
#define left_child(i)       ((i << 1) + 1)
#define right_child(i)      ((i << 1) + 2)
#define get_parent(i)       ((i - 1) >> 1)
#define is_root(i)          (i == root_idx)
#define get_node(i)         (heap->heap[i])
#define get_pos(i)          (heap->pos[i])
#define get_value(i)        (get_node(i)->value)
#define get_key(i)          (get_node(i)->key)
#define in_heap(i)          (heap->pos[i] < heap->size)

// Internal functions

heap_node *     heap_node_constructor(int value, int key);
min_heap *      heap_constructor(int c);
heap_node *     pop_heap(min_heap * heap);
int             update(min_heap * heap, int key, int new_value);
int             swap_node(int a, int b, min_heap * heap);
int             heapify(min_heap * heap, int i);

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
    float v_arr[meta->n_taxa];

    // Mst initialization
    mst->n_taxa             = meta->n_taxa; 
    mst->max_w              = 0.0;
    mst->prim_ord           = malloc(mst->n_taxa * sizeof(int));
    mst->prim_par           = malloc(mst->n_taxa * sizeof(int)); 

    if(!mst->prim_ord || !mst->prim_par)        PRINT_AND_RETURN("malloc failed to init mst in prim\n", MALLOC_ERROR);
    for(i = 0; i < meta->n_taxa; i ++){
        mst->prim_ord[i] = -1;
        mst->prim_par[i] = -1;
    }

    // Heap initialization
    heap = heap_constructor(mst->n_taxa);   
    if(!heap)                       PRINT_AND_RETURN("heap initialization failed in prim\n", MALLOC_ERROR);
    for(i = 0; i < mst->n_taxa; i++){
        v_arr[i] = INT_MAX; // large number
        get_node(i)     = heap_node_constructor(i ? v_arr[i] : 0, i);
        if(!get_node(i))            PRINT_AND_RETURN("node initialization failed in prim\n", MALLOC_ERROR);
        get_pos(i)      = i;
    }

    // Loop initialization
    prim_ordering_counter = 0;
    v_arr[0] = 0;
    heap->size = mst->n_taxa;
                                                                                            #if DEBUG 
                                                                                                                                                                                                 int dflag = 0;

                                                                                            #endif
    // Prim's MST
    while(heap->size){
                                                                                            #if DEBUG 
                                                                                                //checking that heap property is always maintained at the root at each iteration\n"); 
                                                                                                for(i = 1; i < heap->size; i++){
                                                                                                    if(heap->heap[0]->value > heap->heap[i]->value)
                                                                                                        dflag++;

                                                                                                }

                                                                                            #endif
        head = pop_heap(heap);
        if(!head)                   PRINT_AND_RETURN("heap popping failed in prim\n", GENERAL_ERROR);
        cur_key = head->key;

        for(i = 0; i < mst->n_taxa; i++){
            if(in_heap(i) && meta->d[cur_key][i] - v_arr[i] < EPS && i != cur_key){
                v_arr[i] = meta->d[cur_key][i];
                mst->prim_par[i] = cur_key;
                if(update(heap, i, v_arr[i]) != SUCCESS) 
                                    PRINT_AND_RETURN("updating heap value failed in prim\n", GENERAL_ERROR);
                j = i;
            }  
        }
        mst->max_w = MAX(mst->max_w, meta->d[cur_key][j]); 
        mst->prim_ord[prim_ordering_counter++] = cur_key; 
    }
                                                                                            #if DEBUG
                                                                                                if(!dflag) printf("debug: heap property is good throughout at the root\n");
                                                                                                else printf("debug: heap value is wrong, flag is %d\n", dflag);

                                                                                                // printf("debug: the following print out prim's ordering\n");
                                                                                                // for(i = 0; i < meta->n_taxa; i++)
                                                                                                //     printf("%d ", mst->prim_ord[i]);
                                                                                                // printf("\n");

                                                                                                // printf("debug: the following print out prim's tree (parenting)\n");
                                                                                                // for(i = 0; i < meta->n_taxa; i++)
                                                                                                //     printf("%d ", mst->prim_par[i]);
                                                                                                // printf("\n");

                                                                                                printf("debug: test: the prim ordering is unique\n");
                                                                                                dflag = 0;
                                                                                                for(i = 0; i < meta->n_taxa; i++)
                                                                                                    for(j = 0; j < meta->n_taxa; j++)
                                                                                                        if(mst->prim_ord[i] == mst->prim_ord[j] && i != j)
                                                                                                            dflag ++;
                                                                                                if(!dflag) printf("pass test\n");
                                                                                                else printf("failed test with flag %d\n", dflag);

                                                                                                printf("debug: small test: if someone is A's parent's, it must be after A \n");
                                                                                                dflag = 0;
                                                                                                for(i = 1; i < meta->n_taxa; i++){
                                                                                                    int parent = mst->prim_par[i];
                                                                                                    int iflag = 0;

                                                                                                    for(j = 0; j < meta->n_taxa; j++){
                                                                                                        if(mst->prim_ord[j] == i){
                                                                                                            if(iflag){
                                                                                                                iflag++; // this is fine
                                                                                                            } else dflag++;  // this is not
                                                                                                        }
                                                                                                        if(mst->prim_ord[j] == parent){
                                                                                                            if(iflag) dflag++; 
                                                                                                            else iflag++;
                                                                                                        }
                                                                                                    }
                                                                                                }
                                                                                                if(!dflag) printf("pass test\n");
                                                                                                else printf("failed test with flag %d\n", dflag);
                                                                                                // while(1);
                                                                                            #endif
    return 0;
}

int prim_on_small_graph(int n, GRAPH * graph, MST_GRP * mst, char * distance_model, char ** data){
    // Heap variables
    min_heap * heap;
    heap_node * head;

    char * tmp_dist_data[2];

    // Loop variables
    int i; 
    int prim_ordering_counter;
    int cur_key;

    ADJ_LIST * cur_edge;
    double cur_dist;
    double sto_dist;
    int num_site = strlen(data[0]);

    // Temporary stack array
    float v_arr[n];

    // Mst initialization
    mst->n_taxa             = n; 
    mst->max_w              = 0.0;
    mst->prim_ord           = malloc(n * sizeof(int));
    mst->prim_par           = malloc(n * sizeof(int)); 

    if(!mst->prim_ord || !mst->prim_par)        PRINT_AND_RETURN("malloc failed to init mst in prim\n", MALLOC_ERROR);
    for(i = 0; i < n; i ++){
        mst->prim_ord[i] = -1;
        mst->prim_par[i] = -1;
    }

    // Heap initialization
    heap = heap_constructor(mst->n_taxa);   
    if(!heap)                       PRINT_AND_RETURN("heap initialization failed in prim\n", MALLOC_ERROR);
    for(i = 0; i < mst->n_taxa; i++){
        v_arr[i] = INT_MAX; // large number
        get_node(i)     = heap_node_constructor(i ? v_arr[i] : 0, i);
        if(!get_node(i))            PRINT_AND_RETURN("node initialization failed in prim\n", MALLOC_ERROR);
        get_pos(i)      = i;
    }

    // Loop initialization
    prim_ordering_counter = 0;
    v_arr[0] = 0;
    heap->size = mst->n_taxa;

    // Prim's MST
    while(heap->size){
        printf("we\n");
        head = pop_heap(heap);
        if(!head)                   PRINT_AND_RETURN("heap popping failed in prim\n", GENERAL_ERROR);
        cur_key = head->key;

        cur_edge = graph->adjacency_list[cur_key];
        while(cur_edge){
            i = cur_edge->dest;

            tmp_dist_data[0] = data[cur_key];
            tmp_dist_data[1] = data[i];
            cur_dist = (strcmp(distance_model, "JC") == 0) ? compute_jc_distance(tmp_dist_data, num_site) : compute_logdet_distance(tmp_dist_data, num_site);

            if(in_heap(i) && cur_dist - v_arr[i] < EPS && i != cur_key){
                v_arr[i] = cur_dist;
                mst->prim_par[i] = cur_key;
                if(update(heap, i, v_arr[i]) != SUCCESS) 
                                    PRINT_AND_RETURN("updating heap value failed in prim\n", GENERAL_ERROR);
                // j = i;
                sto_dist = cur_dist; 
            }  
            cur_edge = cur_edge->next;
        }
        mst->max_w = MAX(mst->max_w, sto_dist); 
        mst->prim_ord[prim_ordering_counter++] = cur_key; 
    }
    return 0;
}

// INTERNAL IMPLEMENTATION
heap_node * heap_node_constructor(int value, int key){
    heap_node * node = malloc(sizeof(heap_node));

    // Check malloc error

    node->value     = value;
    node->key       = key;

    return node;
}

min_heap * heap_constructor(int c){
    int i;

    min_heap * heap = malloc(sizeof(min_heap));

    heap->size          = 0;
    heap->capacity      = c;
    heap->pos           = malloc(c * sizeof(int));
    heap->heap          = malloc(c * sizeof(heap_node *));

    for(i = 0; i < c; i++){
        heap->heap[i]   = malloc(sizeof(heap_node));
    }

    return heap;
}

int swap_node(int a, int b, min_heap * heap){
    if(!heap) PRINT_AND_RETURN("heap is NULL in swap node\n", GENERAL_ERROR);
    get_pos(get_key(a)) = b;
    get_pos(get_key(b)) = a;

    heap_node * tmp = get_node(a);
    get_node(a) = get_node(b);
    get_node(b) = tmp;



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

    heap_node * tmp = get_node(0);

    // Swap the first and the last node, decrease the size of the heap and heapify the new root if nec.
    if(swap_node(root_idx, last_node_index, heap)   != SUCCESS) PRINT_AND_RETURN("swap node failed in pop_head\n", NULL);
    heap->size--;
    if(heapify(heap, root_idx)                      != SUCCESS) PRINT_AND_RETURN("heapify faield in pop_head\n", NULL);
    return tmp; // this is the old node
}

int update(min_heap * heap, int key, int new_value){
    int idx = get_pos(key);
    get_value(idx) = new_value;

    // Maintain heap property by explicit heapify up
    while(!is_root(idx) && get_value(idx) < get_value(get_parent(idx))){
        if(swap_node(idx, get_parent(idx), heap) != SUCCESS) PRINT_AND_RETURN("swap node failed in update in prim\n", GENERAL_ERROR);
        idx = get_parent(idx);
    }

    return 0;
}
