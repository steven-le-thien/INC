#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "fast_mst.h"
#include "dist.h"
#include "prim.h"

int max_subset_size_mult = 5;
int p = 10;

int int_swap(int* a, int x, int y){
  int z;

  z = a[x];
  a[x] = a[y];
  a[y] = z;

  return 0;
}

// Internal functions 
int init_graph(GRAPH * graph, int n){
  int i;

  // graph = malloc(sizeof(graph));
  // printf("%d\n", n);
  graph->n = n;
  graph->adjacency_list = malloc(n * sizeof(ADJ_LIST*));
    // printf("%d\n",  graph->adjacency_list);

  for(i = 0; i < n; i++){
    // printf("%d\n", i);
    graph->adjacency_list[i] = NULL; // no edges at the beginning
  }
  // printf("awd\n");

  graph->last_node = malloc(n * sizeof(ADJ_LIST*));
  for(i = 0; i < n; i++)
    graph->last_node[i] = NULL;
    // printf("awd\n");

  return 0;
}

int add_edge(GRAPH * graph, int u, int v){
  if(!graph->adjacency_list[u]){
    graph->adjacency_list[u] = malloc(sizeof(ADJ_LIST));
    graph->adjacency_list[u]->dest = v;
    graph->adjacency_list[u]->next = NULL;
    graph->last_node[u] = graph->adjacency_list[u]; 
  } else {
    graph->last_node[u]->next = malloc(sizeof(ADJ_LIST));
    graph->last_node[u]->next->dest = v;
    graph->last_node[u]->next->next = NULL;
    graph->last_node[u] = graph->last_node[u]->next;
  }

  if(!graph->adjacency_list[v]){
    graph->adjacency_list[v] = malloc(sizeof(ADJ_LIST));
    graph->adjacency_list[v]->dest = u;
    graph->adjacency_list[v]->next = NULL;
    graph->last_node[v] = graph->adjacency_list[v]; 
  } else {
    // printf("halg\n");
    graph->last_node[v]->next = malloc(sizeof(ADJ_LIST));
    graph->last_node[v]->next->dest = u;
    graph->last_node[v]->next->next = NULL;
    graph->last_node[v] = graph->last_node[v]->next;
    // printf("end\n");
  }
  return 0;
}

int random_centroids(
    char ** data, 
    DIST_MOD distance_model, 
    int n, 
    int seed,
    int *** disjoint_subset, 
    int ** centroid_arr, 
    double *** small_dist_mat, 
    int *** overlapping_subset)
{
  int * permu; 
  int * centroid_size_count;
  int * fast_disjoint_subset;
  int condition_flag;
  int sqrt_n; // useful precalculation
  int i, j;
  int num_site;
  char * tmp_dist_data[2];

  sqrt_n = (int) sqrt(1.0 * n);
  num_site = strlen(data[0]);

  centroid_size_count = SAFE_MALLOC(sqrt_n * sizeof(int));
  
  // we don't actually construct the subset but maintain an array of indices
  (*disjoint_subset) = SAFE_MALLOC(n * sizeof(int*)); 

  for(i = 0; i < n; i++)
    (*disjoint_subset)[i] = SAFE_MALLOC(sqrt_n * sizeof(int));

  fast_disjoint_subset = SAFE_MALLOC(n * sizeof(int));
  

  permu = SAFE_MALLOC(n * sizeof(int));
  for(i = 0; i < n; i++)
    permu[i] = i;

  (*small_dist_mat) = SAFE_MALLOC(n * sizeof(double*));
  for(i = 0; i < n; i++)
    (*small_dist_mat)[i] = SAFE_MALLOC(sqrt_n * sizeof(double)); 

  srand(seed); // use the random seed
  while(1){
    // Select the centroid at random by shuffling the permutation array. 
    for(i = n - 1; i > 0; i--)
      int_swap(permu, i, rand() % (i + 1)); // Fisher-Yates 

    for(i = 0; i < n; i++)
      for(j = 0; j < sqrt_n; j++)
        (*disjoint_subset)[i][j] = 0; 

    for(i = 0; i < sqrt_n; i++)
      centroid_size_count[i] = 0;

    for(i = 0; i < n; i++)
      fast_disjoint_subset[i] = -1;

    for(j = 0; j < n; j++){
      for(i = 0; i < sqrt_n; i++)
        if(j == permu[i]) 
          (*small_dist_mat)[j][i] = 0.0;
        else{
          if(fast_disjoint_subset[j] == -1){
            fast_disjoint_subset[j] = 0;
            tmp_dist_data[0] = data[permu[0]];
            tmp_dist_data[1] = data[j];
            (*small_dist_mat)[j][0] = distance_model == D_JC ? 
                compute_jc_distance(tmp_dist_data, num_site) : 
                compute_logdet_distance(tmp_dist_data, num_site);
          }
          else {
            tmp_dist_data[0] = data[permu[j]];
            tmp_dist_data[1] = data[permu[i]];
            (*small_dist_mat)[j][i] = distance_model == D_JC ? 
                compute_jc_distance(tmp_dist_data, num_site) : 
                compute_logdet_distance(tmp_dist_data, num_site);

            if((*small_dist_mat)[j][i] < 
                  (*small_dist_mat)[j][fast_disjoint_subset[j]]){
              fast_disjoint_subset[j] = i;
            }
          }
        }
      centroid_size_count[fast_disjoint_subset[j]]++;
    }

    // Check condition
    condition_flag = 1;

    for(i = 0; i < sqrt_n; i++){
      if(centroid_size_count[i] > max_subset_size_mult * sqrt_n){
        condition_flag = 0;
        break;
      }
    }

    if(condition_flag) break;
  }
  for(i = 0; i < n; i++)
    (*disjoint_subset)[i][fast_disjoint_subset[i]] = 1;

  (*centroid_arr) = SAFE_MALLOC(sqrt_n * sizeof(int));
  for(i = 0; i < sqrt_n; i++)
    (*centroid_arr)[i] = permu[i];

  // Centroids are in their own clusters
  for(i = 0; i < sqrt_n; i++)
    for(j = 0; j < sqrt_n; j++)
      if(i == j)
        (*disjoint_subset)[(*centroid_arr)[i]][(*centroid_arr)[j]] = 1;
      else 
        (*disjoint_subset)[(*centroid_arr)[i]][(*centroid_arr)[j]] = 0; 

  // we don't actually construct the subset but maintain an array of indices
  (*overlapping_subset) = SAFE_MALLOC(n * sizeof(int*)); 
  for(i = 0; i < n; i++){
    (*overlapping_subset)[i] = SAFE_MALLOC(sqrt_n * sizeof(int));
    for(j = 0; j < sqrt_n; j++)
      (*overlapping_subset)[i][j] = (*disjoint_subset)[i][j]; 
  }

  return 0;
}

int extend_subset(
    int n, 
    int ** disjoint_subset,
    int p, 
    double ** small_dist_mat)
{ 
  int i, j, k; 
  int sqrt_n;
  int min_idx;

  sqrt_n = (int) sqrt(1.0 * n);

  for(i = 0; i < n; i++){ // the first closest point is already chosen, we only
                          // need the remaining p - 1 point
    for(j = 0; j < p - 1; j++){
      min_idx = -1;
      for(k = 0; k < sqrt_n; k++)
        if(disjoint_subset[i][k]) continue;
        else if(min_idx == -1 || 
              small_dist_mat[i][k] < small_dist_mat[i][min_idx])
          min_idx = k;
      
      disjoint_subset[i][min_idx] = 1; 
    }
  }

  return 0; 
}

int make_cluster_cliques(
    GRAPH * graph, 
    int n, 
    int ** disjoint_subset, 
    int * centroid_arr)
{
  int * clique_v;
  int i, j, k;
  int sqrt_n;
  int clique_count;

  sqrt_n = (int) sqrt(1.0 * n);
  clique_v = SAFE_MALLOC((max_subset_size_mult + p) * sqrt_n * sizeof(int) + 
      sizeof(int)); // one more for the centroid



  for(i = 0; i < sqrt_n; i++){
    clique_count = 0;
    for(j = 0; j < n; j++){
      if(disjoint_subset[j][i]) 
        clique_v[clique_count++] = j;
    }

    clique_v[clique_count++] = centroid_arr[i];

    for(j = 0; j < clique_count; j++)
      for(k = j + 1; k < clique_count; k++){
        // need to check when there are overlaps
        add_edge(graph, clique_v[j], clique_v[k]); 
      }
  }

  free(clique_v);

  return 0;
}

int make_centroid_clique(GRAPH * graph, int n, int * centroid_arr){
  int i, j;
  int sqrt_n = (int) sqrt(1.0 * n);

  for(i = 0; i < sqrt_n; i++)
    for(j = i+1; j < sqrt_n; j++)
      add_edge(graph, centroid_arr[i], centroid_arr[j]);

  return 0; 
}

int make_complete_bipartite(GRAPH * graph, int n, int * centroid_arr){
  int * centroid_map;
  int sqrt_n;
  int i, j;

  sqrt_n = (int) sqrt(1.0 * n);
  centroid_map = malloc(n * sizeof(int)); 

  for(i = 0; i < n; i++)
    centroid_map[i] = 0; 

  for(i = 0; i < sqrt_n; i++) 
    centroid_map[centroid_arr[i]] = 1;

  for(i = 0; i < sqrt_n; i++)
    for(j = 0; j < n; j++)
      if(j != centroid_arr[i] && !centroid_map[j])
        add_edge(graph, j, centroid_arr[i]);

  return 0;
}

int nn_from_sequences(
    GRAPH * graph, 
    int *** disjoint_subset, 
    char ** data, 
    DIST_MOD distance_model, 
    int n, 
    int seed)
{
  int ** overlapping_subset = NULL;
  int * centroid_arr = NULL;
  double ** small_dist_mat = NULL;

  // Randomly select the centroid and create the X_i's, this also copies
  //    disjoint subset to overlapping one

  FCAL(
      GENERAL_ERROR,
      F_RANDOM_CENTROID_IN_NN_FROM_SEQ,
      random_centroids(
          data, 
          distance_model, 
          n, 
          seed, 
          disjoint_subset, 
          &centroid_arr, 
          &small_dist_mat, 
          &overlapping_subset
      )
  );
  printf(STATE_DONE_CENTROID);

  // Create the Y_i's by adding additional close centroids
  FCAL(
      GENERAL_ERROR,
      F_EXT_SUBSET_IN_NN_FROM_SEQ,
      extend_subset(n, overlapping_subset, p, small_dist_mat)
  );
  printf(STATE_DONE_EXTEND_SS);

  // Make cliques from the Y_i's
  FCAL(
      GENERAL_ERROR,
      F_MAKE_CLUSTER_CLIQUES_IN_NN_FROM_SEQ,
      make_cluster_cliques(graph, n, overlapping_subset, centroid_arr)
  );
  printf(STATE_DONE_MAKE_CLUSTER_CLIQUES);

  // Make centroid clique
  FCAL(
      GENERAL_ERROR,
      F_MAKE_CENTROID_CLIQUES_IN_NN_FROM_SEQ,
      make_centroid_clique(graph, n, centroid_arr)
  );

  // Make complete bipartite and centroid clique
  FCAL(
      GENERAL_ERROR,
      F_MAKE_COMPLETE_BP_IN_NN_FROM_SEQ,
      make_complete_bipartite(graph, n, centroid_arr)
  );

  return 0;
}

// This is an approximate mst that does not depend on the distance matrix and 
// runs in time O(n^1.5). Input is a set of sequences in a char array 
// Output is an mst as defined in prim.c (since the last step of this function 
// is always to call the prim algorithm)
int fast_mst(
    char ** data, 
    int n, 
    DIST_MOD distance_model, 
    int seed, 
    MST_GRP * mst, 
    int *** disjoint_subset)
{
  GRAPH small_graph; 

  // Compute NN graph
  FCAL(
      GENERAL_ERROR,
      F_INIT_GRAPH_IN_FAST_MST,
      init_graph(&small_graph, n)
  );
  FCAL(
      GENERAL_ERROR,
      F_NN_FROM_SEQ_IN_FAST_MST,
      nn_from_sequences(
          &small_graph, 
          disjoint_subset, 
          data, 
          distance_model, 
          n, 
          seed
      )
  );
  printf(STATE_DONE_NN);

  FCAL(
      GENERAL_ERROR,
      F_PRIM_ON_SMALL_GRAPH, 
      prim_on_small_graph(n, &small_graph, mst, distance_model, data)
  );
  printf(STATE_DONE_PRIM);
  return 0; 
}