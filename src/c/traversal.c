// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "traversal.h"
#include "utilities.h"
#include "options.h"
#include "quartet.h"
#include "fast_lca.h"
#include "dist.h"
#include "tree.h"

int skip_counter =0;

// TODO: figure out ith BFS used master_to_midx

// All trees should be BTs 

// All trees should be BTs 
int dfs_lca_implementation(
    int node, 
    int parent, 
    BT * tree, 
    int * dp, 
    int * in_building, 
    int * lca_parent, 
    int * lca_child, 
    int mode);

int dfs_backtrack(
    int node, 
    int parent, 
    int findme, 
    BT * tree, 
    int * next_to_start, 
    int * found);

int dfs_preorder(int node, int parent, BT * tree, int * in_building, int flag);

int bfs_vote_implementation(
    INC_GRP * meta,
    MAP_GRP * map,
    VOTE_GRP * vote,
    int x,
    double q0,
    int power,
    int ** revote_map,
    int all_quartets,
    int revoting);

int less_than_3_shared_taxa(INC_GRP * meta);

void all_valid(INC_GRP * meta, VOTE_GRP * vote);

int init_bfs(VOTE_GRP * vote, double * vote_a, int * parent_map, int * queue);

int init_revote(int revoting, int ** revote_map, int n);

int update_edge_count(
    BT * tree, 
    int cur, 
    int * parent_map, 
    int * edge_count, 
    int * adj_idx_a);

int run_qmethod(
    INC_GRP * meta, 
    int * master_to_midx,
    char ** name_map, 
    int * u, 
    int * quartet_result,
    double * M);

int do_quartet(
    INC_GRP * meta, 
    MAP_GRP * map, 
    VOTE_GRP * vote, 
    int cur, 
    int * adj_idx_a, 
    double q0, 
    double * itrnl_vote, 
    double * max_vote,
    int * revote_map,
    int revote_power,
    int x,
    int all_quartets);


/* Function to initialize variables relevant to voting to a blank state so 
 *      that a new round of voting can be done 
 * Input:   meta,   including the array of visited taxa in the growing tree,
 *              the growing tree as well as the constraint tree
 *          map     including bookkeeping arrays for fast lookups
 *          mst     including the prim ordering to get the next taxon
 *          i       indexing into the prim ordering
 * Output:  0 on success, ERROR otherwise
 * Effet:   run memset to clear fields, reset meta's visited array, 
 *              identify the correct constraint tree
 */
int init_vote(
    INC_GRP * meta, 
    MAP_GRP * map, 
    MST_GRP * mst, 
    VOTE_GRP * vote, 
    int i)
{
  int j;
  // Set all other special edges to 0 at one go. 
  // WARNING: redo this when the struct changes   
  memset(&(vote->valid_st), -1, 10 * sizeof(int));

  // Clear visited array from bipartition (but still keep the visted flags)
  for(j = 0; j < meta->n_taxa; j++)
    if(meta->visited[j] > 0) 
      meta->visited[j] = -1;

  // Identify the constraint tree
  vote->ctree_idx = map->master_to_ctree[mst->prim_ord[i]];

  return 0;
}

/* Assuming that the correct constraint tree is found and meta's visited array has  
 *      the correct record for leaves in the growing tree, this function 
 *      1. identify the intersection A between the leaf set of the growing 
 *          tree and the correct constraint tree
 *      2. determine the bipartition of the intersection set within the constraint 
 *          tree by removing the component containing a certain taxon in t_c^A
 * Exact implementation detail is as followed:
 *      1. Given L, the set of leaves in the growing tree, root the tree at 
 *          the query taxon and run a tree dp (/a post-order dfs) recursion that, 
 *          at each internal node n, count the number of leaves in the subtree 
 *          rooted at n contained in L, this is the dp value
 *      2. As we move up the tree, we also keep track of the last time dp changes
 *          (by storing the parent and child of the edge on which the value change). 
 *          This gives the LCA of taxa in L within the constraint tree 
 *      3. The LCA, when removed, divides the constraint tree into 3 subtrees, 
 *          1 containing the query taxon and no leaves in L; another containing 
 *          a subset of L_1 leaves in L and the last contain the rest of L.
 *          Note that we want to mark L_1 and L/L_1. 
 *      4. 2 simple preorder dfs, each starting at 2 children of the LCA found 
 *          previously, suffice to identify and mark out this bipartition of L
 * Input:   all meta groups
 * Output:  0 on success, ERROR otherwise
 * Effect:  set some fields in vote
 */ 
int find_bipartition(
    INC_GRP * meta, 
    MAP_GRP * map, 
    MST_GRP * mst, 
    VOTE_GRP * vote, 
    int i)
{
  int placeholder = 0;    // extra variable for dp
  int bp_counter = 0;     // bipartition counter

  BT * lca_tree;
  int lca_node;          
  int lca_parent;
  int cur_child;
  int j; 

  
  if(vote->ctree_idx < 0 
      || meta->ctree[vote->ctree_idx]->n_node < 3) 
    return 0; // trivial tree, don't have any bipartition (from lca)

  // Specifying local variables
  lca_tree = meta->ctree[vote->ctree_idx];

  FCAL(
      GENERAL_ERROR,
      F_DFS_LCA_IMP_IN_FIND_BIPART,
      dfs_lca_implementation(
          map->master_to_cidx[mst->prim_ord[i]],
          -1,
          lca_tree,
          &placeholder,                           
          meta->visited,                          
          &(vote->c_lca.p),                      
          &(vote->c_lca.c),
          0
      )
  );

  // If there are less than 3 shared common taxa then don't bother
  if(placeholder < 3) return 0;

  // Identify the 2 children of the LCA and mark the correct bipartition
  lca_node = vote->c_lca.c;
  lca_parent = vote->c_lca.p;
  for(j = 0; j < get_degree(lca_tree, lca_node); j++){
    cur_child = get_adj(lca_tree, lca_node, j);
    if(cur_child != lca_parent)
      FCAL(
          GENERAL_ERROR,
          F_DFS_PRE_IN_FIND_BIPART,
          dfs_preorder(
              cur_child,
              lca_node,
              lca_tree,
              meta->visited,
              ++bp_counter
          )
      );
  }

  return 0;
}

/* Assuming the correct bipartition of A (LeafSet(t_c) \cap LeafSet(t)) where 
 *    t_c is the constraint tree and t is the growing tree
 *    has been identified, this function find the valid subtree (the 'component 
 *    in t_c') in which the query taxon may be added. 
 * The exact implementation detail is as followed:
 *      1. Let the biparition be {A1, A2}. There guarantees to be a node in A1
 *          (i.e, A1 is nonempty). Find that node, root the tree there and find 
 *          the LCA of A2, call this l2
 *      2. Likewise, find some node in A2, root the tree there and find the LCA 
 *          of A1, call this l1
 *      3. Run any graph traversal (here, a dfs equiped with backtracking) from 
 *          l1 and stops when reaching l2. This attempts to find the correct 
 *          orientation of l1 and l2 in the valid subtree
 *        3.1 Here we denote the parent of l1 is the neighbor of l1 that leads
 *            to l2. This is unique since there's a unqiue path from any 2 
 *            nodes in a tree 
 * Input: all meta variables
 * Output: 0 on success, ERROR otherwise
 * Effect: set fields in vote
 */
int find_valid_subtree(
    INC_GRP * meta, 
    MAP_GRP * map, 
    MST_GRP * mst, 
    VOTE_GRP * vote)
{
  int i;
  int placeholder;

  BT * gtree = meta->gtree;

  if(less_than_3_shared_taxa(meta)){ // if less than 3 shared taxa
    all_valid(meta, vote);
    return 0;
  }

  // Run from the first leaf. We root the tree at the first taxon. 
  // If the taxon is in the constraint tree, we find the other bipartition; 
  //    else we find the 1st biparition
  for(i = 0; i < meta->n_taxa; ++i){
    if(meta->visited[i] == 2){
      break;
    }
  }
  FCAL(
      GENERAL_ERROR,
      F_DFS_LCA_IMP_IN_FIND_VALID,
      dfs_lca_implementation(
          map->master_to_gidx[i],
          -1,
          gtree,
          &placeholder,
          meta->visited,
          &(vote->st_lca.p),  // parent isn't important here                
          &(vote->st_lca.c),
          1
      )
  );

  // Find lca of the other bipartition that was not found in the first phase, 
  //    by starting from the first phase lca, avoiding its own parent ()
  for(i = 0; i < meta->n_taxa; ++i){
    if(meta->visited[i] == 1){
      break;
    }
  }
  // for(i = 0; i < meta->n_taxa; i = meta->visited[i] == 1 ? meta->n_taxa : i + 1);  
  //   ASSERT(
  //     GENERAL_ERROR,
  //     MISSING_BIP,
  //     i == meta->n_taxa
  // );  
  FCAL(
      GENERAL_ERROR,
      F_DFS_LCA_IMP_IN_FIND_VALID,
      dfs_lca_implementation(
          map->master_to_gidx[i],
          -1,
          gtree,
          &placeholder,
          meta->visited,
          &(vote->nd_lca.p),  // parent isn't important here                
          &(vote->nd_lca.c),
          2
      )
  );

  // Find correct orientation
  placeholder = 0;
  FCAL(
      GENERAL_ERROR,
      F_DFS_BACKTRACK_IN_FIND_VALID,
      dfs_backtrack(
          vote->st_lca.c,
          -1,
          vote->nd_lca.c,
          gtree,
          &(vote->st_lca.p),
          &placeholder
      )
  );
  return 0;
}

/* Assuming that the correct valid subtree is identified and the orientation of
 *    its endpoint is correct, this runs a bfs from one of the endpoints 
 *    to update vote counting relative to the starting edge. This is 
 *    described in Theorem 7 of the original paper. Implementation detail is as 
 *    followed
 *  1. Recall that we determine the orientation of the 2 endpoints l1, 
 *      l2 by the neighbor of l1 that leads to l2. Call this vertex lp, the 
 *      algo in Theorem 7 starts with (l1, lp)
 *  2. At each vertex visited in the bfs order, we are only interested in the 
 *      internal nodes (since we want to update values on the edge). Thus 
 *      seeing a leaf or seeing l2 ends the run.
 *  3. At each internal node, we identify the `parent' in the bfs order. 
 *      Theorem 7 gives a way to update the 2 children edge based on the 
 *      parents' value; the validity of the voting quartet and results from 4 
 *      point method. 
 *  4. Node with higher vote than its parent is immediate checked whether it is
 *      the node with the most vote encountered so far.
 * Input: all meta variables
 * Output: 0 on success, ERROR otherwise
 * Effect: set fields in vote
 */
int bfs_vote(
    INC_GRP * meta, 
    MAP_GRP * map, 
    MST_GRP * mst, 
    VOTE_GRP * vote, 
    int i)
{
  int * revote_map;
  if(vote->st_lca.p == vote->nd_lca.c && 
      vote->st_lca.c == vote->nd_lca.p){ // only 1 edge is valid
    vote->ins.c = vote->st_lca.c;
    vote->ins.p = vote->st_lca.p;
    return 0;
  }
  vote->tree = meta->gtree;
  FCAL(
      GENERAL_ERROR,
      F_BFS_VOTE_IMPL_IN_BFS_VOTE,
      bfs_vote_implementation(
          meta,
          map, 
          vote,
          mst->prim_ord[i],
          (double) mst->max_w,
          INV_SQR,
          &revote_map,
          ALL_QUARTET,
          NO_REVOTING
      )
  );

  free(revote_map);
  return 0;
}



// All trees should be BTs 
int dfs_lca_implementation(
    int node, 
    int parent, 
    BT * tree, 
    int * dp, 
    int * in_building, 
    int * lca_parent, 
    int * lca_child, 
    int mode)
{
  int i;
  int child_dp; 
  int child_lca_parent;
  int child_lca_child;
  int flag;

  if(tree->degree[node] == 1 && parent != -1){
    if((!mode && in_building[tree->master_idx_map[node]]) ||
      (mode && in_building[tree->master_idx_map[node]] == mode)){
      *dp = 1;
      *lca_parent = parent;
      *lca_child = node;
      return 0;
    } else {
      *dp = 0;
      *lca_child = -1;
      *lca_parent = -1;
      return 0;
    }
  }

  // Initialization
  *dp = 0;
  flag = 0;
  child_dp = -1;
  child_lca_parent    = -1;
  child_lca_child     = -1;

  for(i = 0; i < tree->degree[node]; i++){
    if(tree->adj_list[node][i].dest == parent) continue;
    FCAL(
        GENERAL_ERROR,
        F_REC_DFS_LCA,
        dfs_lca_implementation(
            tree->adj_list[node][i].dest,
            node,
            tree,
            &child_dp,
            in_building,
            &child_lca_parent,
            &child_lca_child,
            mode
        )
    ); 

    if(child_dp){
      *dp += child_dp;
      if(!flag){
        *lca_child = child_lca_child; 
        *lca_parent = child_lca_parent; flag = 1;
      } else {
        *lca_child = node;  
        *lca_parent = parent;
      }
    }
  }
  return 0;
}

int dfs_backtrack(
    int node, 
    int parent, 
    int findme, 
    BT * tree, 
    int * next_to_start, 
    int * found)
{
  int i;
  int child_found = 0;

  if(*found) return 0; // someone found it first (but nothing should be here)
  if(node == findme){ // gottem 
    *found = 1;
    *next_to_start = node;
  } else { // search next
    for(i = 0; i < tree->degree[node]; i++){
      if(tree->adj_list[node][i].dest == parent) continue;

      FCAL(
          GENERAL_ERROR,
          F_REC_DFS_BACK,
          dfs_backtrack(
              tree->adj_list[node][i].dest,
              node,
              findme,
              tree,
              next_to_start,
              &child_found
          )
      );
      if(child_found && !(*found)){
        *found = 1;
        *next_to_start = tree->adj_list[node][i].dest; 
        break;
      }
    }
  }
  return 0;
}

int dfs_preorder(int node, int parent, BT * tree, int * in_building, int flag){
  int i;
  if(tree->degree[node] == 1)
    if(in_building[tree->master_idx_map[node]])
      in_building[tree->master_idx_map[node]] = flag;

  for(i = 0; i < tree->degree[node]; i++){
    if(tree->adj_list[node][i].dest == parent) continue;
    FCAL(
        GENERAL_ERROR,
        F_REC_DFS_PRR,
        dfs_preorder(
            tree->adj_list[node][i].dest,
            node, 
            tree,
            in_building, 
            flag
        )
    );
  }
  return 0;
}

int bfs_vote_implementation(
    INC_GRP * meta,
    MAP_GRP * map,
    VOTE_GRP * vote,
    int x,
    double q0,
    int power,
    int ** revote_map,
    int all_quartets,
    int revoting)
{
  // Tree statistics
  BT * tree = meta->gtree;
  int n = tree->n_node;
  int e = vote->nd_lca.c;
  int s_p = vote->st_lca.c; 

  // Initialize the first edge to some arbitrary array
  int qs = 0, qe = 0, queue[n], edge_count = 1;; // mock queue
  int parent_map[n];
  double itrnl_vote[n];

  // State
  int cur, i;

  // Quartet method
  int adj_idx_a[3];

  // Trackers
  double max_vote = revoting ? (revote_map[0] ? 0.0 : -10000.0) : 0.0;
  // Init sequence 
  FCAL(
      GENERAL_ERROR,
      F_INIT_BFS_IN_BFS,
      init_bfs(
          vote,
          itrnl_vote,
          parent_map,
          queue
      )
  );

  FCAL(
      GENERAL_ERROR,
      F_INIT_IN_REVOTE_IN_BFS,
      init_revote(
          revoting,
          revote_map, 
          n
      )
  );

  // First vertex in the BFS
  queue[qe++] = vote->st_lca.p; 
  while(qs != qe){
    cur = queue[qs++];

    // We only deal with inner edges
    if(get_degree(tree, cur) == 1 || cur == e || cur == s_p) continue;

    // Updating edge count 
    FCAL(
        GENERAL_ERROR,
        F_UPDATE_EDGE_IN_BFS,
        update_edge_count(
            tree, 
            cur, 
            parent_map, 
            &edge_count,
            adj_idx_a
        )
    );

    // 4 point
    FCAL(
        GENERAL_ERROR,
        F_COMPUTE_M_IN_BFS,
        do_quartet(
            meta,
            map,
            vote,
            cur, 
            adj_idx_a,
            q0,
            itrnl_vote,
            &max_vote,
            (*revote_map),
            power,
            x,
            all_quartets
        )
    );

    for(i = 0; i < get_degree(tree, cur); i++){
      // Stopping condition
      if(get_adj(tree, cur, i) == parent_map[cur]) continue;
      else 
        queue[qe++] = tree->adj_list[cur][i].dest;
    }
  }

  if(!revoting) // first round, write down revote
    for(i = 0; i < n - 1; i++)
      (*revote_map)[i] = 
          (itrnl_vote[i] - max_vote < EPS) && 
          (itrnl_vote[i] - max_vote > -EPS);
  else  // others
    for(i = 0; i < n - 1; i++)
      (*revote_map)[i] = 
          ((*revote_map)[i]) && 
          (itrnl_vote[i] - max_vote < EPS) && 
          (itrnl_vote[i] - max_vote > -EPS);

  return 0;
}

int less_than_3_shared_taxa(INC_GRP * meta){
  int cap_counter = 0;
  int j;

  cap_counter = 0;
  for(j = 0; j < meta->n_taxa; j++)
    cap_counter += meta->visited[j] > 0;

  return cap_counter < 3; 
}

void all_valid(INC_GRP * meta, VOTE_GRP * vote){
  vote->st_lca.p = meta->gtree->adj_list[0][0].dest; 
  vote->st_lca.c = 0;
}

int init_bfs(VOTE_GRP * vote, double * vote_a, int * parent_map, int * queue){
  BT * tree = vote->tree;
  int n = tree->n_node;
  int s = vote->st_lca.p;
  int s_p = vote->st_lca.c; 
  int i;

  memset(vote_a, 0, n* sizeof(double));
  memset(parent_map, -1, n * sizeof(int));
  memset(queue, -1, n * sizeof(int));

  vote->ins.c = s;
  vote->ins.p = s_p;

  parent_map[s] = s_p;
  for(i = 0; i < get_degree(tree, s); i++)
    if(get_adj(tree, s, i) == s_p)
      set_edge_master_idx(tree, s, i, 0);

  for(i = 0; i < get_degree(tree, s_p); i++)
    if(get_adj(tree, s_p, i) == s)
      set_edge_master_idx(tree, s_p, i, 0);

  return 0;
}

int init_revote(int revoting, int ** revote_map, int n){
  int i; 
  if(!revoting){
    *revote_map = SAFE_MALLOC(n * sizeof(int)); 
    for(i = 0; i < n; (*revote_map)[i++] = 1); 
  }
  return 0;
}

int update_edge_count(
    BT * tree, 
    int cur, 
    int * parent_map, 
    int * edge_count, 
    int * adj_idx_a)
{
  const int DEGREE = 3;
  int i, j;
  int adj_counter;
  int nex;

  // Find parent order
  for(i = 0; i < DEGREE; i++)
    if(get_adj(tree, cur, i) == parent_map[cur])
      adj_idx_a[0] = i;

  adj_counter = 1;
  for(i = 0; i < DEGREE; i++)
    if(i != adj_idx_a[0]){
      adj_idx_a[adj_counter++] = i;

      nex = get_adj(tree, cur, i);
      for(j = 0; j < get_degree(tree, nex); j++)
        if(get_adj(tree, nex, j) == cur)
          break; 

      set_edge_master_idx(tree, cur, i, *edge_count);
      set_edge_master_idx(tree, nex, j, *edge_count);
      (*edge_count)++;

      parent_map[nex] = cur;
    }
  return 0;
}

int run_qmethod(
    INC_GRP * meta, 
    int * master_to_midx,
    char ** name_map, 
    int * u, 
    int * quartet_result,
    double * M)
{
  int i;
  float ** d = meta->d, ** dm = meta->dm;

  ml_options * master_ml_options = meta->master_ml_options;
  QTREE_METHOD qtree_method = master_ml_options->qtree_method;

  int use_distance_matrix = master_ml_options->use_distance_matrix;
  LCA_T * LCA = master_ml_options->LCA; 

  switch(qtree_method){
    case Q_FPM:
      FCAL(
        GENERAL_ERROR,
        F_FPM_IN_BFS,
        four_point_method(d, u, quartet_result)
      );
      break;

    case Q_SUBTREE:
      if(use_distance_matrix){
        for(i = 0; i < 4; i++)
          u[i] = master_to_midx[u[i]];
        
        FCAL(
            GENERAL_ERROR,
            F_FPM_MAT_IN_BFS,
            four_point_method(dm, u, quartet_result)
        );
      }
      else
        FCAL(
            GENERAL_ERROR,
            F_FPM_MAT_IN_BFS,
            fpm_on_tree(LCA, u, quartet_result)
        );
      break;

    case Q_RAXML:
      FCAL(
          GENERAL_ERROR,
          F_RAXML_Q_IN_BFS,
          new_quartets_raxml(name_map, u, quartet_result, master_ml_options)
      );
      break;

    case Q_ML:
      FCAL(
          GENERAL_ERROR,
          F_ML_Q_IN_BFS,
          ml_quartet(
              name_map, 
              u, 
              quartet_result, 
              master_ml_options, 
              M          
          )
      );
      break;   
  }
  return 0;
}

int do_quartet(
    INC_GRP * meta, 
    MAP_GRP * map, 
    VOTE_GRP * vote, 
    int cur, 
    int * adj_idx_a, 
    double q0, 
    double * itrnl_vote, 
    double * max_vote,
    int * revote_map,
    int revote_power,
    int x,
    int all_quartets)
{
  const int DEGREE = 3;
  const int NUM_CHILD = 2;
  const int CHILD_OFFSET = 1;
  const int PAR_OFFSET = 0;

  int i, j;
  int nex_idx;
  int u[4]; // 0 is parent, 1 is first child, 2 is second child, 
            // 3 is the insertion taxon
  double M = -1e9;
  BT * tree = meta->gtree;

  int quartet_result;

  int * gidx_to_master = tree->master_idx_map;
  DIST_MOD distance_model = meta->master_ml_options->distance_model;
  int use_distance_matrix = meta->master_ml_options->use_distance_matrix; 

  u[3] = x; 
  for(i = 0; i < DEGREE; i++)
    u[i] = gidx_to_master[get_edge_sample(tree, cur, adj_idx_a[i])];

  for(i = 0; i < 4; ++i)
    for(j = i + 1; j < 4; ++j)
      M = MAX(M, 
          use_distance_matrix ?
              meta->d[u[i]][u[j]] : 
              dist_from_msa(
                  meta->msa, 
                  distance_model,
                  u[i],
                  u[j], 
                  meta->correction
              )
      );

  // Invalid quartet case
  for(i = 0; i < NUM_CHILD; i++)
    itrnl_vote[get_edge_master_idx(tree, cur, adj_idx_a[CHILD_OFFSET + i])]
        = itrnl_vote[get_edge_master_idx(tree, cur, adj_idx_a[PAR_OFFSET])];

  if(all_quartets || M - 8 * q0 <= EPS) {
    FCAL(
        GENERAL_ERROR,
        F_RUN_QMETHOD_IN_DO_Q,
        run_qmethod(
            meta,
            map->master_to_midx,
            map->master_to_name,
            u,
            &quartet_result,
            &M
        )
    );

    if(quartet_result == 0)  // parent wins
      for(i = 0; i < 2; i++)
        itrnl_vote[get_edge_master_idx(tree, cur, adj_idx_a[PAR_OFFSET])] 
            -= POWER(1.0 / M, revote_power);
    for(i = 0; i < NUM_CHILD; i++)
      if(quartet_result == i + CHILD_OFFSET)
        itrnl_vote[get_edge_master_idx(tree, cur, adj_idx_a[i + CHILD_OFFSET])]
            -= POWER(1.0 / M, revote_power);
  }

  for(i = 0; i < 2; i++){
    nex_idx = get_edge_master_idx(tree, cur, adj_idx_a[i + CHILD_OFFSET]);
    if(revote_map[nex_idx] && itrnl_vote[0] - (*max_vote) > EPS){
        (*max_vote) = itrnl_vote[nex_idx];
        vote->ins.c = get_adj(tree, cur, adj_idx_a[i + CHILD_OFFSET]);
        vote->ins.p = cur;
    }
  }

  return 0;
}
