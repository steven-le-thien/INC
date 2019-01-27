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
  int i; 
  int prim_ordering_counter;
  int cur_key;

  // Temporary stack array
  double v_arr[meta->n_taxa];

  // Twice the maximum defined distance
  // double m_def;
  double tmp_dist;
  double sto_dist;

  // char * tmp_dist_data[2];
  // int num_site;

  // if(!meta->master_ml_options->use_distance_matrix){
  //     // Compute m_def
  //     m_def = (double) -1e9;
  //     num_site = meta->msa->N;
  //     tmp_dist_data[0] = malloc(num_site);
  //     tmp_dist_data[1] = malloc(num_site);
  //     for(i = 0; i < meta->n_taxa; i++){
  //          printf("i is %d\n", i);
  //         for(j = i + 1; j < meta->n_taxa; j++){
         

  //             strcpy(tmp_dist_data[0], meta->msa->msa[i]);
  //             strcpy(tmp_dist_data[1], meta->msa->msa[j]);

  //             tmp_dist = meta->master_ml_options->distance_model == D_JC ? compute_jc_distance(tmp_dist_data, num_site) : compute_logdet_distance(tmp_dist_data, num_site);
  //             if(tmp_dist >= 0.0 && tmp_dist > m_def) 
  //                 m_def = tmp_dist;

  //             meta->correction = m_def;
  //         }
  //     }
  //     m_def *= 2.0;
 
  // }
  // Mst initialization
  mst->n_taxa             = meta->n_taxa; 
  mst->max_w              = 0.0;
  mst->prim_ord           = SAFE_MALLOC(mst->n_taxa * sizeof(int));
  mst->prim_par           = SAFE_MALLOC(mst->n_taxa * sizeof(int)); 

  for(i = 0; i < meta->n_taxa; i ++){
    mst->prim_ord[i] = -1;
    mst->prim_par[i] = -1;
  }

  // Heap initialization
  heap = heap_constructor(mst->n_taxa);   

  for(i = 0; i < mst->n_taxa; i++){
    v_arr[i] = INT_MAX; // large number
    get_node(i)     = heap_node_constructor(i ? v_arr[i] : 0, i);
    get_pos(i)      = i;
  }

  // Loop initialization
  prim_ordering_counter = 0;
  v_arr[0] = 0;
  heap->size = mst->n_taxa;

  // Prim's MST
  while(heap->size){
    head = pop_heap(heap);
    ASSERT(
      GENERAL_ERROR,
      F_HEAD_POP_IN_PRIM, 
      head
    );
    cur_key = head->key;

    for(i = 0; i < mst->n_taxa; i++){
      if(meta->master_ml_options->use_distance_matrix)
        tmp_dist = meta->d[cur_key][i];
      else
        tmp_dist = dist_from_msa(
            meta->msa, 
            meta->master_ml_options->distance_model, 
            cur_key, 
            i, 
            meta->correction
        );
  
      if(in_heap(i) && tmp_dist - v_arr[i] < EPS && i != cur_key){
        v_arr[i] = tmp_dist;
        mst->prim_par[i] = cur_key;
        FCAL(
            GENERAL_ERROR, 
            F_UPDATE_IN_PRIM,
            update(heap, i, v_arr[i])
        );
        sto_dist = tmp_dist;
      }  
    }
    mst->max_w = MAX(mst->max_w, sto_dist); 
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

int parse_initial_tree_as_mst(INC_GRP * meta, MST_GRP * mst){

  // meta->master_ml_options->init_tree_name; 
  return 0;
}

int prim_on_small_graph(int n, GRAPH * graph, MST_GRP * mst, DIST_MOD distance_model, char ** data){
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
  mst->prim_ord           = SAFE_MALLOC(n * sizeof(int));
  mst->prim_par           = SAFE_MALLOC(n * sizeof(int)); 

  for(i = 0; i < n; i ++){
    mst->prim_ord[i] = -1;
    mst->prim_par[i] = -1;
  }

  // Heap initialization
  heap = heap_constructor(mst->n_taxa);   
  for(i = 0; i < mst->n_taxa; i++){
    v_arr[i] = INT_MAX; // large number
    get_node(i)     = heap_node_constructor(i ? v_arr[i] : 0, i);
    get_pos(i)      = i;
  }

  // Loop initialization
  prim_ordering_counter = 0;
  v_arr[0] = 0;
  heap->size = mst->n_taxa;

  // Prim's MST
  while(heap->size){
    head = pop_heap(heap);
    ASSERT(
      GENERAL_ERROR,
      F_HEAD_POP_IN_PRIM, 
      head
    );
    cur_key = head->key;

    cur_edge = graph->adjacency_list[cur_key];
    while(cur_edge){
      i = cur_edge->dest;

      tmp_dist_data[0] = data[cur_key];
      tmp_dist_data[1] = data[i];
      cur_dist = distance_model == D_JC ? 
          compute_jc_distance(tmp_dist_data, num_site) : 
          compute_logdet_distance(tmp_dist_data, num_site);

      if(in_heap(i) && cur_dist - v_arr[i] < EPS && i != cur_key){
        v_arr[i] = cur_dist;
        mst->prim_par[i] = cur_key;
        FCAL(
            GENERAL_ERROR, 
            F_UPDATE_IN_PRIM,
            update(heap, i, v_arr[i])
        );
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
  heap_node * node = SAFE_MALLOC(sizeof(heap_node));

  node->value     = value;
  node->key       = key;

  return node;
}

min_heap * heap_constructor(int c){
  int i;

  min_heap * heap = SAFE_MALLOC(sizeof(min_heap));

  heap->size          = 0;
  heap->capacity      = c;
  heap->pos           = SAFE_MALLOC(c * sizeof(int));
  heap->heap          = SAFE_MALLOC(c * sizeof(heap_node *));

  for(i = 0; i < c; i++){
    heap->heap[i]   = SAFE_MALLOC(sizeof(heap_node));
  }

  return heap;
}

int swap_node(int a, int b, min_heap * heap){
  heap_node * tmp;
  get_pos(get_key(a)) = b;
  get_pos(get_key(b)) = a;

  tmp = get_node(a);
  get_node(a) = get_node(b);
  get_node(b) = tmp;

  return 0;
}

// Heapify down implementation
int heapify(min_heap * heap, int i){
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
    FCAL(
        GENERAL_ERROR,
        F_SWAP_NODE_IN_HEAPIFY,
        swap_node(i, min_index, heap)
    );
    FCAL(
        GENERAL_ERROR,
        F_SWAP_NODE_IN_HEAPIFY,
        swap_node(i, min_index, heap)
    );
    FCAL(
        GENERAL_ERROR,
        F_REC_IN_HEAPIFY,
        heapify(heap, min_index)
    );
    
  }

  return 0;
}


heap_node * pop_heap(min_heap * heap){
  if(!heap->size) return NULL;

  heap_node * tmp = get_node(0);

  // Swap the first and the last node, decrease the size of the heap 
  //    and heapify the new root if nec.
  FCAL(
      NULL,
      F_SWAP_NODE_IN_POP_HEAD,
      swap_node(root_idx, last_node_index, heap)
  );
  heap->size--;
  FCAL(
      NULL,
      F_HEAPIFY_IN_POP_HEAD,
      heapify(heap, root_idx)
  );
  return tmp; // this is the old node
}

int update(min_heap * heap, int key, int new_value){
  int idx = get_pos(key);
  get_value(idx) = new_value;

  // Maintain heap property by explicit heapify up
  while(!is_root(idx) && get_value(idx) < get_value(get_parent(idx))){
    FCAL(
        GENERAL_ERROR,
        F_SWAP_NODE_IN_UPDATE,
        swap_node(idx, get_parent(idx), heap)
    );
    idx = get_parent(idx);
  }
  return 0;
}
