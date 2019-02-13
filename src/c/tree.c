// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tree.h"
#include "utilities.h"
#include "options.h"
#include "tools.h"

// Global variable for quickly intializing a tree
BT_edge adj_list[MAXN][3];
int degree[MAXN];
int parent_map[MAXN];
int master_idx_map[MAXN];

// Private functions
int attach_leaf_to_edge_impl(
    BT * growing_tree, 
    int x, 
    int addition_edge_parent,
    int additional_edge_child, 
    int adjacent_in_mst);

BT * tree_constructor(
    int n, 
    BT_edge ** adj_list, 
    int * degree, 
    int * master_idx_map);

BT * read_newick(
    MAP_GRP * map, 
    char * filename, 
    int tree_idx);

int init();

int check(int node_a, BT * tree);

int swap_adajcency(int nod_a, int node_b, int node_c, BT* tree);

int make_adjacent(int node_a, int node_b, BT * tree);

int make_parent(int node_p, int node_c, BT * tree);

int save_name(int cur_node, char * name, MAP_GRP * map, int tree_idx);

void petite_dfs(
    int node, 
    int parent, 
    char ** name_map, 
    char * builder, 
    BT * tree);

/* Read all constraint trees at once into memory
 * Input:       meta        meta variables, including the pointers to 
 *                          constraint trees in memory
 *              map         map variables, including the indexing map 
 *                          and name map
 *              options     user input options, including the name of 
 *                          constraint trees
 * Output: 0 on success, ERROR otherwise
 * Effect: allocate some memories, build some trees, open some files, 
 *         init 2 arrays in map and 1 in meta 
 */
int parse_tree(INC_GRP * meta, MAP_GRP * map, ml_options * options){
  int i;  // loop variable

  // Init meta's field
  if(meta->master_ml_options->qtree_method == Q_SUBTREE){
    meta->n_ctree = options->num_trees - 1;
    map->master_to_midx    = SAFE_MALLOC(meta->n_taxa  * sizeof(int));

    if(meta->master_ml_options->use_distance_matrix){
      meta->dm = SAFE_MALLOC(meta->n_taxa * sizeof(float*));
      for(i = 0; i < meta->n_taxa; i++){
        meta->dm[i] = SAFE_MALLOC(meta->n_taxa * sizeof(float));
      }

      FCAL(
          GENERAL_ERROR,
          F_MK_UNW_MAT_IN_PARSE_TREE,
          make_unweighted_matrix(
              options->tree_names[0], 
              options->output_prefix, 
              meta->dm, 
              map->master_to_name, 
              map->master_to_midx
          )
      );
    }
  } else meta->n_ctree = options->num_trees;
 
  // printf("wda %d\n", options->num_trees);

  // Malloc sequence
  meta->ctree             = SAFE_MALLOC(meta->n_ctree * sizeof(BT *));
  map->master_to_ctree    = SAFE_MALLOC(meta->n_taxa  * sizeof(int));
  map->master_to_cidx     = SAFE_MALLOC(meta->n_taxa  * sizeof(int));
  map->master_to_gidx     = SAFE_MALLOC(meta->n_taxa  * sizeof(int));

  // Init
  for(i = 0; i < meta->n_taxa; i++){
    map->master_to_ctree[i]     = -1;
    map->master_to_cidx[i]      = -1;
    map->master_to_gidx[i]      = -1;        
  }

  // Sanitize 
  for(i = 0; i < meta->n_ctree; i++){
    if(meta->master_ml_options->qtree_method == Q_SUBTREE)
      FCAL(
          GENERAL_ERROR, 
          F_RM_LBL_IN_MK_FT_CONSTRAINT,
          rm_label_job(options->tree_names[i + 1], options->tree_names[i + 1])
      );
    else 
      FCAL(
          GENERAL_ERROR, 
          F_RM_LBL_IN_MK_FT_CONSTRAINT,
          rm_label_job(options->tree_names[i], options->tree_names[i])
      );    
  }

  // Call the newick reader
  for(i = 0; i < meta->n_ctree; i++){
    if(meta->master_ml_options->qtree_method == Q_SUBTREE)
      meta->ctree[i]      = read_newick(map, options->tree_names[i + 1], i);
    else 
      meta->ctree[i]      = read_newick(map, options->tree_names[i], i);
  }                                               
  return 0;
}


/* Creating a 3-leaf star from the firt 3 nodes visited in the Prim's ordering
 *    of the mst
 * Input:       meta        meta variables, including the pointers to growing 
 *                          tree in memory
 *              mst         mst variables, including the prim ordering and the 
 *                          mst itself
 * Output: 0 on success, ERROR otherwise
 * Effect: allocate some memories, build some trees, init the growing tree and 
 *    visited array in meta
 */
int init_growing_tree(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst){ 
  const int N_START_TAXON = 3;
  int i;

  // Update master->idx mapping
  for(i = 0; i < N_START_TAXON; i++)
    map->master_to_gidx[mst->prim_ord[i]] = i;

  // Check that there must be at least 3 leaves
  ASSERT(
      GENERAL_ERROR,
      LESS_THAN_THREE_TAX,
      meta->n_taxa >= 3
  );

  // Meta initialization
  meta->visited = SAFE_MALLOC(meta->n_taxa * sizeof(int));
  memset(meta->visited, 0, meta->n_taxa * sizeof(int));

  // Growing tree initialization
  meta->gtree = tree_constructor(4 * meta->n_taxa,  NULL, NULL, NULL); 
  ASSERT(
      GENERAL_ERROR,
      N_GTREE,
      meta->gtree
  );

  // Override tree logistics since there are only 4 real nodes in the tree, 
  //    the rest can be viewed as singletons
  meta->gtree->n_node = 4;
  for(i = 0; i < N_START_TAXON; i++)
    FCAL(
        GENERAL_ERROR,
        F_MK_ADJ_IN_INIT_GTREE,
        make_adjacent(i, 3, meta->gtree)
    );

  // Init tree variables  
  for(i = 0; i < 3; i++){
    meta->gtree->master_idx_map[i] = mst->prim_ord[i];
    meta->visited[mst->prim_ord[i]] = -1;
    meta->gtree->adj_list[3][i].sample = i;
  }

  // This is a bit tricky, we need to find adjacent nodes in the mst
  meta->gtree->adj_list[0][0].sample = 1; // fist vertex in prim is obviously 
                                          //    connected to second
  meta->gtree->adj_list[1][0].sample = 0; // vice versa

  // check if the parent of 3rd vertex in mst is adjacent to 0 or 1
  meta->gtree->adj_list[2][0].sample = 
      map->master_to_gidx[mst->prim_par[mst->prim_ord[2]]]; 

  return 0;
}

int attach_leaf_to_edge(
    INC_GRP * meta,  
    MAP_GRP * map,
    MST_GRP * mst, 
    VOTE_GRP * vote, 
    int i)
{
  // Set the new edge to be present in the growing tree
  meta->visited[mst->prim_ord[i]] = -1;

  // Actual attachment
  FCAL(
      GENERAL_ERROR,
      F_ATTACH_LEAF,
      attach_leaf_to_edge_impl(
          meta->gtree,
          mst->prim_ord[i],
          vote->ins.p,
          vote->ins.c,
          map->master_to_gidx[mst->prim_par[mst->prim_ord[i]]]
      ) 
  );

  // Set the correct mapping of the index of the new edge
  //    (the internal tree mapping is set in the function above)
  map->master_to_gidx[mst->prim_ord[i]] = meta->gtree->n_node - 1;
  return 0;

}


int write_newick(BT * tree, char * filename, char ** name_map){
  char buffer[MAX_BUFFER_SIZE * 100];
  FILE * f;

  // Subdivide an edge and make it the root
  tree->n_node++;

  swap_adajcency(0, tree->n_node - 1, tree->adj_list[0][0].dest, tree);
  STR_CLR(buffer);
  petite_dfs(tree->n_node - 1, -1, name_map, buffer, tree);
  strcat(buffer, ";");

  f = fopen(filename, "w");
  fprintf(f, "%s", buffer);
  fclose(f);

  return 0;
}

int get_degree(BT * tree, int idx){
  return tree->degree[idx];
}

void set_degree(BT * tree, int idx, int val){
  tree->degree[idx] = val;
}

void incr_deg(BT * tree, int idx){
  set_degree(tree, idx, get_degree(tree, idx) + 1);
}

int get_adj(BT * tree, int idx, int order){
  return tree->adj_list[idx][order].dest;
}

void set_adj(BT * tree, int idx, int order, int val){
  tree->adj_list[idx][order].dest = val;
}

int get_edge_sample(BT * tree, int idx, int order){
  return tree->adj_list[idx][order].sample;
}

void set_edge_sample(BT * tree, int idx, int order, int val){
  tree->adj_list[idx][order].sample = val;
}

void set_edge_master_idx(BT * tree, int ini, int dest, int val){
  tree->adj_list[ini][dest].master_idx = val;
}

int get_edge_master_idx(BT * tree, int ini, int dest){
  return tree->adj_list[ini][dest].master_idx;
}

// INTERNAL FUNCTION IMPLEMENTATIONS
int attach_leaf_to_edge_impl(
    BT * growing_tree, 
    int x, 
    int additional_edge_parent, 
    int additional_edge_child, 
    int adjacent_in_mst)
{
  int i, n;
  // Changing tree logistics (if this fails or terminate in the middle of some 
  //    process then we are screwed)
  // Because of the way we mallocate the growing tree, there are already memory 
  //    for the extra 2 nodes 
  n = growing_tree->n_node += 2;

  // Subdivide the edge
  FCAL(
      GENERAL_ERROR,
      F_SWAP_ADJ_IN_ATTACH_EDGE,
      swap_adajcency(
          additional_edge_parent,
          n - 2,
          additional_edge_child,
          growing_tree
      )
  );

  FCAL(
      GENERAL_ERROR,
      F_MK_ADJ_IN_ATTACH_EDGE,
      make_adjacent(
          n - 1, 
          n - 2, 
          growing_tree
      )
  );

  for(i = 0; i < 3; i++)
    if(get_adj(growing_tree, n - 2, i) == n - 1)
      set_edge_sample(growing_tree, n - 2, i, n - 1);

  set_edge_sample(growing_tree, n - 1, 0, adjacent_in_mst);

  growing_tree->master_idx_map[growing_tree->n_node - 1] = x;

  return 0;
}

BT * tree_constructor(
    int n, 
    BT_edge ** adj_list, 
    int * degree, 
    int * master_idx_map)
{
  // The arrays mayeb null (in case of null construction)
  int i, j;

  BT * tree = SAFE_MALLOC(sizeof(BT)); 

  tree->n_node = n;

  tree->adj_list              = SAFE_MALLOC(n * sizeof(BT_edge *));
  tree->degree                = SAFE_MALLOC(n * sizeof(int));
  tree->master_idx_map        = SAFE_MALLOC(n * sizeof(int));

  for(i = 0; i < n; i++){
    tree->adj_list[i] = SAFE_MALLOC(3 * sizeof(BT_edge));
  }
  // Initialization 
  for(i = 0; i < n; i++){
    if(adj_list && adj_list[i]) 
      memcpy(tree->adj_list[i], adj_list[i], 3 * sizeof(BT_edge));
    else 
      for(j = 0; j < 3; j++){
        tree->adj_list[i][j].dest = -1;
        tree->adj_list[i][j].sample = -1;
        tree->adj_list[i][j].master_idx = -1;
      }

    tree->master_idx_map[i] = -1;
    tree->degree[i] = 0; 
  }

  if(degree)
    memcpy(tree->degree, degree, n * sizeof(int));

  if(master_idx_map)
    memcpy(tree->master_idx_map, master_idx_map, n * sizeof(int));

  return tree;
}

#define READ_NAME_STATE 0
#define OTHER_STATE     1

#define incr_node(node, max_node) \
    node = ++max_node

#define incr_level(node_p, node_c, max_node)\
    do{\
      node_p = node_c;\
      incr_node(node_c, max_node);\
    } while(0)

#define decr_level(node_p, node_c)\
    do{\
      node_c = node_p;\
      node_p = parent_map[node_c];\
    } while(0)


BT * read_newick(MAP_GRP * map, char * filename, int tree_idx){
  FILE *  f;
  int i, j;

  BT_edge ** mock;

  // FSM variables
  int     cur_state;

  int     cur_node;
  int     max_node;
  int     parent_node;

  char    cur_char;
  char    cur_name[MAX_NAME_SIZE];

  // File opening
  f = SAFE_FOPEN_RD(filename);

  // Intialization
  cur_state       = READ_NAME_STATE;
  cur_node        = 0;
  max_node        = 0;
  parent_node     = -1;

  cur_char        = 0;
  STR_CLR(cur_name);

  // Initialize global variables
  init();

  while(fscanf(f, "%c", &cur_char) == 1){
    switch(cur_char){
      case '(': //start of a new level and start of a new node
        incr_level(parent_node, cur_node, max_node);
        FCAL(
            NULL,
            F_MK_PAR_IN_RD_NW,
            make_parent(parent_node, cur_node, NULL) 
        );
        cur_state = READ_NAME_STATE;
        break;
      case ',': // start of new node end of old node
        FCAL(
            NULL,
            F_SV_NAME_IN_RD_NW,
            save_name(cur_node, cur_name, map, tree_idx)
        ); 
        incr_node(cur_node, max_node);
        FCAL(
            NULL,
            F_MK_PAR_IN_RD_NW,
            make_parent(parent_node, cur_node, NULL) 
        );
        cur_state = READ_NAME_STATE;
        break;
      case ')': // end of level
        FCAL(
            NULL,
            F_SV_NAME_IN_RD_NW,
            save_name(cur_node, cur_name, map, tree_idx)
        ); 
        decr_level(parent_node, cur_node); 
        cur_state = OTHER_STATE;
        break;
      case ':': 
        cur_state = OTHER_STATE;
        break;
      default:
        if(cur_state == READ_NAME_STATE)
          cat_c_to_a(cur_char, cur_name);
    }
  }

  fclose(f);
  // Deal with the case where the root has degree 2
  if(degree[0] == 2){
    // Find its 2 children, 1 is guaranteed to be 1
    for(i = 0; i < max_node + 1; i++)
      if(parent_map[i] == 0 && i != 1) break;

    for(j = 0; j < degree[1]; j++)
      if(adj_list[1][j].dest == -1) adj_list[1][j].dest = i - 1;

    for(j = 0; j < degree[i]; j++)
      if(adj_list[i][j].dest == -1) adj_list[i][j].dest = 0;
  }

  mock = malloc(max_node * sizeof(BT_edge));
  for(i = 0;i < max_node; i++)
    mock[i] = &(adj_list[1 + i][0]);

  // Make the binary tree and return
  return tree_constructor(max_node, mock, &degree[1], &master_idx_map[1]);
}


int init(){
  memset(parent_map, -1, sizeof(parent_map));
  memset(adj_list, -1, sizeof(adj_list));
  memset(degree, 0, sizeof(degree));
  memset(master_idx_map, -1, sizeof(degree));
  return 0;
}

/* Helper function to check whether a node is in range
 * Input:   node to check
 * Output:  1 if the node is in range, 0 otherwise
 * Effect:  none
 */ 
int check(int node_a, BT * tree){
  if(!tree) return node_a >= 0; // no upper limit
  if(0 <= node_a && node_a < tree->n_node) return 1;
  else return 0;
}

// Delete a-c and make a-b, b-c 
int swap_adajcency(int node_a, int node_b, int node_c, BT* tree){
  int i, j;
  ASSERT(
      GENERAL_ERROR,
      F_CHECK_NODE_IN_SWAP_ADJ,
      check(node_a, tree) && check(node_b, tree) 
  );

  if(tree){
    for(i = 0; i < get_degree(tree, node_a); i++)
      if(get_adj(tree, node_a, i) == node_c){
        set_adj(tree, node_a, i, node_b);
                              // this is the only thing we need to change 
                              // since the index will be refreshed in the next 
                              // iteration and the leaf sample remains the same
        break; // preserve i for later use
      }
    for(j = 0; j < get_degree(tree, node_c); j++)
      if(get_adj(tree, node_c, j) == node_a){
        set_adj(tree, node_c, j, node_b);
        break; // preserve j for later use
      }

    // Assuming b is new so its fields are empty
    set_adj(tree, node_b, get_degree(tree, node_b), node_a);
    set_edge_sample(
        tree, 
        node_b, 
        get_degree(tree, node_b), 
        get_edge_sample(tree, node_c, j)
    );
    incr_deg(tree, node_b);

    set_adj(tree, node_b, get_degree(tree, node_b), node_c);
    set_edge_sample(
        tree, 
        node_b,
        get_degree(tree, node_b),
        get_edge_sample(tree, node_a, i)
    );
    incr_deg(tree, node_b);
  }
  return 0;
}

/* Helper function to make a node adjacent to another in the adjacency matrix
 * Input:   2 nodes
 * Output:  0 on success, ERROR otherwise
 * Effect:  none
 */ 
int make_adjacent(int node_a, int node_b, BT * tree){
  ASSERT(
      GENERAL_ERROR,
      F_CHECK_NODE_IN_SWAP_ADJ,
      check(node_a, tree) && check(node_b, tree) 
  );

  if(tree){
    set_adj(tree, node_a, get_degree(tree, node_a), node_b);
    incr_deg(tree, node_a);

    set_adj(tree, node_b, get_degree(tree, node_b), node_a);
    incr_deg(tree, node_b);
  } else {
    adj_list[node_a][degree[node_a]++].dest = node_b - 1;
    adj_list[node_b][degree[node_b]++].dest = node_a - 1;
  }

  return 0;
}

/* Helper function to make a node a parent of another 
 * Input:   2 nodes (order matters, the parent node comes first)
 * Output:  0 on success, ERROR otherwise
 * Effect:  none
 */ 
int make_parent(int node_p, int node_c, BT * tree){
  FCAL(
      GENERAL_ERROR,
      F_MK_ADJ_IN_MK_PAR,
      make_adjacent(node_p, node_c, tree)
  ); 

  parent_map[node_c] = node_p;
  return 0;
}

/* Helper function to save name for a numeric node and map it back to the 
 *    master ordering
 * Input:   node number and its name 
 * Output:  0 on success
 * Effect:  none
 */ 
int save_name(int cur_node, char * name, MAP_GRP * map, int tree_idx){
  int i; // store the master index
  for(i = 0; i < map->n_taxa; i++){
    if(strcmp(name, map->master_to_name[i]) != 0){
      if(i == map->n_taxa - 1) return 0;  // no name matches, which is fine for 
                                          // internal nodes
    } else break;
  }
  if(map){
    map->master_to_ctree[i]     = tree_idx;
    map->master_to_cidx[i]      = cur_node - 1;
  }

  master_idx_map[cur_node]    = i;

  STR_CLR(name);

  return 0;
}

void petite_dfs(
    int node, 
    int parent, 
    char ** name_map, 
    char * builder, 
    BT * tree)
{
  int i;
  if(tree->degree[node] == 1) {
    strcat(builder, name_map[tree->master_idx_map[node]]);
    return;
  }

  strcat(builder, "(");
  for(i = 0; i < tree->degree[node]; i++){
    if(tree->adj_list[node][i].dest == parent) continue;

    petite_dfs(tree->adj_list[node][i].dest, node, name_map, builder, tree);
    if(i != tree->degree[node] - 1) strcat(builder, ",");
  }

  strcat(builder, ")");
  return;
}




