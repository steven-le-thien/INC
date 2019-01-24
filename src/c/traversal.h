// File in inc_ml, created by Thien Le in July 2018

#ifndef TRAVERSAL_H
#define TRAVERSAL_H

#include "c_inc.h"

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

#endif 