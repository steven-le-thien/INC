// File in inc_ml, created by Thien Le in July 2018

#ifndef TRAVERSAL_H
#define TRAVERSAL_H

#include "c_inc.h"

#define INV_SQR  2
#define INV      1
#define ALL_QUARTET 1
#define VALID_QUARTET 0
#define IS_REVOTING 1
#define NO_REVOTING 0

extern int init_vote(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote, int i);
extern int find_bipartition(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote, int i);
extern int find_valid_subtree(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote);
extern int bfs_vote(INC_GRP * meta,  MAP_GRP * map,MST_GRP * mst, VOTE_GRP * vote, int i);

static const char F_FPM_IN_BFS[]
  = "fpm failed in bfs\n";

static const char F_FPM_MAT_IN_BFS[]
  = "subtree failed in bfs\n";

static const char F_RAXML_Q_IN_BFS[]
  = "quartet with raxml failed in bfs\n";

static const char F_ML_Q_IN_BFS[]
  = "ml quartet (VP6) failed in bfs\n";

static const char F_DFS_PRE_IN_FIND_BIPART[]
  = "dfs preorder failed in finding bipartition\n";

static const char F_DFS_LCA_IMP_IN_FIND_BIPART[]
  = "dfs_lca_implementation failed in finding bipartition\n";

static const char MISSING_BIP[]
  = "missing bipartition\n";

static const char F_DFS_LCA_IMP_IN_FIND_VALID[]
  = "dfs_lca_implementation failed in find_valid_subtree\n";

static const char F_DFS_BACKTRACK_IN_FIND_VALID[] 
  = "dfs_backtrack failed in find_valid_subtree\n";

static const char F_BFS_VOTE_IMPL_IN_BFS_VOTE[]
  = "bfs_vote_implementation failed in bfs_vote\n";

static const char F_REC_DFS_LCA[] 
  = "recursive dfs_lca failed in dfs_lca_implementation\n";

static const char F_REC_DFS_BACK[]
  = "recursive dfs_backtrack failed\n"; 

static const char F_REC_DFS_PRR[]
  = "recursive dfs_preorder failed\n";

static const char F_INIT_BFS_IN_BFS[] 
  = "initial bfs failed in bfs\n";

static const char F_INIT_IN_REVOTE_IN_BFS[]
  = "init failed in revote in bfs\n";

static const char F_COMPUTE_M_IN_BFS[]
  = "compute m failed in bfs\n";

static const char F_RUN_QMETHOD_IN_DO_Q[] 
  = "run_qmethod failed in do_q\n"; 

static const char F_UPDATE_EDGE_IN_BFS[]
  = "update_edge failed in bfs\n";

#endif 