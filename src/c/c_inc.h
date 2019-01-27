// File in inc_ml, created by Thien Le in July 2018

#ifndef C_INC_H
#define C_INC_H

#include "utilities.h"

typedef enum CTREE_METHOD{
  C_NO,
  C_SUBTREE,
  C_RAXML,
  C_FASTTREE,
  C_NJ,
  C_FASTME
} CTREE_METHOD;

typedef enum QTREE_METHOD{
  Q_FPM,
  Q_SUBTREE,
  Q_RAXML,
  Q_ML
} QTREE_METHOD;

typedef enum DIST_MOD{
  D_JC,
  D_LOGDET,
  D_K2P,
  D_P
} DIST_MOD;

typedef enum ITREE_METHOD{
  I_FASTTREE,
  I_RAXML,
  I_NJ,
  I_FASTME
} ITREE_METHOD;

typedef struct msa {
  int   N;        // size of one sequence
  int   num_seq;  // number of sequences
  char** msa;
  char** name;
} msa_t;


// Trees
typedef struct m{ // this is an adjacency list edge, meaning that it is 
                  // directed and that the source is know when the 
                  // structure is accessed

  int   dest;      // destination 

  int   master_idx;   // indexing into the vote array (only for the 
                      // voting), it maybe possible to remove this field

  int   sample;     // leaf sample at the destination
  double weight;
} BT_edge;

typedef struct tree{
  int n_node;           // number of nodes in the tree
  BT_edge ** adj_list;      // n_node array of adjacent nodes
  int * degree;          // n_node array of degee

  int * master_idx_map;      // n_node indexing from the constrained 
                  // tree's scheme to the master scheme, 
                  // this is not a surjection
} BT;

typedef struct rmq{
  // Statistics
  int     n;
  int     num_blk;
  int     num_exhaust;
  int     blk_sz;

  // Arrays and mappings
  int     * a;
  int     * blk_idx_to_config;

  // Tables
  int     ** sparse_table;    // num_blk x log2(num_blk) array,  
                    // where A[i][j] is the index of the 
                    // block with the smallest value from 
                    // the jth block to the j + 2^i block

  int     *** exhaustive_table;  // config_num x blk_sz x blk_sz array 
                    // where A[i][j][k] is the relative 
                    // index of the smallest value from the 
                    // jth position to the kth position 
                    // for the ith configuration
} RMQ_T;

typedef struct lca_t{
  // LCA stuff
  int   * euler_tour;
  int   * level;
  int   * first_occurence;
  int   * visited;
  BT   * tree;
  RMQ_T  * RMQ;

  // Smallest distance stuff
  double * d_from_root;
} LCA_T;

typedef struct ml_options{
  char * input_alignment;   // currently only accepts path
  char * output_prefix;    // currently only accepts path
  char * init_tree_name;   // currently only accepts path 
  char * init_d_name;     // currently only accepts path
  char * guide_tree_name;   // currently only accepts path

  int num_trees;
  char ** tree_names ;

  // Use distance matrix 
  int use_distance_matrix;

  // Constraint trees settings
  int recompute_constraint_trees;

  // Methods
  ITREE_METHOD itree_method;
  CTREE_METHOD ctree_method; 
  QTREE_METHOD qtree_method; 
  DIST_MOD distance_model; 

  // Use the intial tree as the spanning tree
  int use_initial_tree_as_spanning_tree; 
  LCA_T *   LCA; 

  // PASTA decomposition setting
  int ss_threshold; 
} ml_options;
// Whole graphs 
typedef struct list{
  int     dest;
  struct list *  next;
} ADJ_LIST;

typedef struct graph{
  int n;
  ADJ_LIST ** last_node; // O(1) insertion
  ADJ_LIST ** adjacency_list;
} GRAPH;

// Groups
typedef struct inc_grp{
  int     n_taxa;
  int     n_ctree;

  BT *    gtree;
  BT **    ctree;   // n_ctree-length array of constraint tree

  int *    visited;  // n_taxa-lenght array of visted taxa in the 
              // building tree, this is in the master indexing 
              // scheme
            
  float **  d;     // n_taxa by n_taxa array of distance matrix; 
              // in the event that the matrix is not used, 
              // d[0][0] stores twice the maximum defined value

  float **  dm;     // n_taxa by n_taxa array of distance matrix on the 
              // FastTree

  double   correction;

  msa_t *   msa;    // only when distance matrix is not used 

  ml_options * master_ml_options;

} INC_GRP;

typedef struct mapping{
  int   n_taxa;

  int *  master_to_midx;   // n_taxa length array mapping the master 
                // indexing scheme to the index within the 
                // binary constraint tree

  int *  master_to_ctree;  // n_taxa length array mapping the master 
                // indexing scheme to the index of the binary 
                // constraint tree

  int *  master_to_cidx;   // n_taxa length array mapping the master 
                // indexing scheme to the index within the 
                // binary constraint tree

  int *  master_to_gidx;
  char ** master_to_name;   // n_taxa length array mapping the master 
                // indexing scheme to its name
} MAP_GRP;

typedef struct mst{
  int   n_taxa;
  int *  prim_ord;    // n_taxa length array denoting the Prim's 
              // insertion order in building the mst

  int *  prim_par;    // n_taxa length array denoting the parent of 
              // each node in the mst (this can be viewed as the
              // tree itself)

  float  max_w;     // maximum edge weight in the mst
} MST_GRP;

// Special edges 
typedef struct s_edge{
  int   c;     // child node, usually the starting node
  int   p;     // parent node, usually the avoiding node
} S_EDGE;

typedef struct vote{
  BT * tree;

  // WARNING: redo initialization if the struct changes
  S_EDGE valid_st;    // the valid 'component', indexed in growing tree 
              // scheme

  S_EDGE c_lca;     // lca of taxa in the induced constraint tree, 
              // this follows the constraint's indexing scheme

  S_EDGE st_lca;     // lca of taxa in the growing tree following the 
              // first bipartition, this follows the growing 
              // tree's indexing scheme

  S_EDGE nd_lca;     // lca of taxa in the growing tree following the 
              // second bipartition, this follow the growing 
              // tree's indexing scheme

  S_EDGE ins;      // insertion edge (edge with the most vote), 
              // this follws the growing tree's indexing scheme

  int   ctree_idx; 
} VOTE_GRP;

// typedef struct options{
//   int num_options;
//   int num_trees;

//   int input_index;
//   char * input_name;

//   int output_index;
//   char * output_name;

//   int tree_index;
//   char ** tree_names ;
// } option_t;

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

#define N_CTREE_METHOD  6

static const char CTREE_FLAG[N_CTREE_METHOD][MAX_FLAG_SZ] = {
  "no",
  "subtree",
  "raxml",
  "fasttree",
  "nj",
  "me",
};

#define N_QTREE_METHOD  4

static const char QTREE_FLAG[N_QTREE_METHOD][MAX_FLAG_SZ] = {
  "fpm",
  "subtree",
  "raxml",
  "ml",
};

#define N_ITREE_METHOD 4

static const char ITREE_FLAG[N_ITREE_METHOD][MAX_FLAG_SZ] = {
  "fasttree",
  "raxml",
  "nj",
  "fastme",
};

// extern const int N_DIST_MOD = 4;
#define N_DIST_MOD  4

static const char DIST_MOD_FLAG[N_DIST_MOD][MAX_FLAG_SZ] = {
  "JC",
  "logDet",
  "K2P",
  "P",
};

extern int constraint_inc_main(
    int argc, 
    char ** argv, 
    ml_options * master_ml_options
);

// String defaults 
static const char DEFAULT_FT_SUF[] 
  = "first_tree.tree";

static const char DEFAULT_DIST_SUF[] 
  = "c_inc_input";

// static const char STATE_PASTA_DEC[]
//   = "starting PASTA decompositions...\n";

static const char F_RD_CMD_ARG_IN_CINC[] 
  = "read_cmd_arg failed in constraint_inc\n";

// static const char STATE_PARSE_MAT_TREE[]
//   = "parsing distance matrix and trees...\n";

static const char F_PARSE_MAT_TREE_IN_CINC[]
  = "parsing matrix and tree failed in constraint inc\n";

static const char F_PARSE_INPUT_IN_CINC[]
  = "parsing input matrices failed in constraint inc\n";

static const char F_INIT_META_IN_CINC[]
  = "initialize meta failed in constraint inc\n";

// static const char STATE_MST[]
//   = "constructing the mst...\n";

static const char F_PRIM_IN_CINC[]
  = "prim failed in constraint inc\n";

static const char F_PARSE_INIT_AS_MST_IN_CINC[]
  = "parse initial matrices as mst in constraint inc\n";

static const char F_PARSE_TREE_IN_CINC[]
  = "parse initial tree failed in constraint inc\n";

// static const char STATE_INIT_GTREE[]
//   = "initialize growing tree...\n";

static const char F_INIT_GTREE_IN_CINC[]
  = "initializing growing tree failed in constraint_inc\n";

// static const char STATE_BUILD_TREE[]
//   = "building the tree...\n";

static const char F_SERIAL_MAIN_LOOP_IN_CINC[]
  = "serial main loop failed in constraint inc\n";

// static const char STATE_OUT_TREE[] 
//   = "outputting the tree...\n";

static const char F_WRITE_TREE_IN_CINC[]
  = "write tree failed in constraint inc\n";

// static const char STATE_CLEAN[]
//   = "cleaning up...\n";

// static const char ITER_COUNT[] = "current iter is";

static const char F_INIT_VOTE_IN_CINC[] 
  = "initializing vote failed in constraint inc \n";

static const char F_FIND_BIPART_IN_CINC[]
  = "finding bipartition failed in constraint inc\n";

static const char F_FIND_VALID_ST_IN_CINC[]
  = "finding valid subtree failed in constraint inc\n";

static const char F_BFS_VOTE_IN_CINC[]
  = "bfs voting failed in constraint inc\n";

static const char F_ATTACH_IN_CINC[]
  = "attaching leaves failed in constraint inc\n";

// static const char STATE_NO_DIST[] 
//   = "without distance matrix...\n";

#endif