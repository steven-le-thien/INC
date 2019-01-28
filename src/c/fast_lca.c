#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utilities.h"
#include "fast_rmq.h"
#include "fast_lca.h"

// Internal functions
int euler_dfs(
  int cur_node, 
  int * visited, 
  LCA_T * LCA, 
  int * level_counter, 
  int * node_counter);

int fast_lca(LCA_T * LCA, int st, int nd, int * ret);
int llsp(int st, int nd, double * ret, LCA_T * LCA);

// Given a tree, and 2 nodes, we want to find the internode distance in O(1) 
//      time with some preprocessing
int fast_lca_init(LCA_T * LCA){
  int * visited;
  int level_counter = 0;
  int node_counter = 0;

  // 0. Set up the variables 
  LCA->euler_tour = SAFE_MALLOC(3 * LCA->tree->n_node * sizeof(int));
  LCA->level = SAFE_MALLOC(3 * LCA->tree->n_node * sizeof(int));
  LCA->first_occurence = SAFE_MALLOC(LCA->tree->n_node * sizeof(int)); 
  visited = SAFE_MALLOC(LCA->tree->n_node * sizeof(int)); 
  LCA->RMQ = SAFE_MALLOC(sizeof(RMQ_T));
  LCA->d_from_root = SAFE_MALLOC(LCA->tree->n_node * sizeof(double));
  memset(LCA->level, -1, LCA->tree->n_node);
  memset(LCA->first_occurence, -1, LCA->tree->n_node);
  memset(visited, 0, LCA->tree->n_node);

  // 1. Build an Euler Tour via DFS, keeping track of the depth and 
  //    first occurence
  euler_dfs(0, visited, LCA, &level_counter, &node_counter);

  // 2. Build RMQ structure
  fast_rmq_init(node_counter, LCA->level, LCA->RMQ); 

  return 0;
}

// Assuming that the LCA has been initialized, perform fpm in constant time
// The indices passed in are indexed from the guide tree
int fpm_on_tree(LCA_T * LCA, int * u, int * ret){
  int p = u[0], u1 = u[1], u2 = u[2], x = u[3];
  double pu1, pu2, px, u1u2, u1x, u2x;
  double sump, sum1, sum2;

  llsp(p, u1, &pu1, LCA);
  llsp(p, u2, &pu2, LCA);
  llsp(p, x, &px, LCA);
  llsp(u1, u2, &u1u2, LCA);
  llsp(u1, x, &u1x, LCA);
  llsp(u2, x, &u2x, LCA);

  // *weight = MAX(MAX(MAX(MAX(MAX(pu1, pu2), px), u1u2), u1x), u2x);

  sump = px + u1u2;
  sum1 = pu2 + u1x;
  sum2 = pu1 + u2x;

  if(sump - sum1 < -0.000000000001 && sump - sum2 < -0.000000000001)
    *ret = 0;
  else{
    if(sum1 - sum2 < -0.000000000001)
      *ret = 1;
    else 
      *ret = 2;
  }
  return 0;
} 


// Internal implementation
int euler_dfs(
    int cur_node, 
    int * visited, 
    LCA_T * LCA, 
    int * level_counter, 
    int * node_counter)
{
  int i;
  int child_node;

  LCA->euler_tour[*node_counter] = cur_node;
  LCA->level[*node_counter] = *level_counter;
  if(!visited[cur_node]){
    LCA->d_from_root[cur_node] = *level_counter;
    LCA->first_occurence[cur_node] = *node_counter;
    visited[cur_node] = 1;
  }

  for(i = 0; i < LCA->tree->degree[cur_node]; i++){
    child_node = LCA->tree->adj_list[cur_node][i].dest;
    if(!visited[child_node]){
      (*level_counter)++;
      (*node_counter)++;
      euler_dfs(child_node, visited, LCA, level_counter, node_counter); 
      (*level_counter)--;
      (*node_counter)++;
      LCA->euler_tour[*node_counter] = cur_node;
      LCA->level[*node_counter] = *level_counter; 
    }
  }

  return 0;
}

// Given a initialized LCA structure and 2 queries in node index, 
//    return the LCA
int fast_lca(LCA_T * LCA, int st, int nd, int * ret){
  int ret_idx;
  int qs, qe;

  qs = LCA->first_occurence[st] > LCA->first_occurence[nd] ? 
      LCA->first_occurence[nd] : 
      LCA->first_occurence[st];

  qe = LCA->first_occurence[st] > LCA->first_occurence[nd] ? 
      LCA->first_occurence[st] : 
      LCA->first_occurence[nd];

  // Ask for the min level
  fast_rmq(qs, qe, &ret_idx, LCA->RMQ);

  // Return the node that is the min level
  *ret = LCA->euler_tour[ret_idx];

  return 0;
}

// Assuming that the LCA has been initialized, find the shortest distance
// from st to nd
int llsp(int st, int nd, double * ret, LCA_T * LCA){
  int lca;

  fast_lca(LCA, st, nd, &lca);
  *ret = 
      LCA->d_from_root[st] + LCA->d_from_root[nd] - 2.0 * LCA->d_from_root[lca];

  return 0;
}


