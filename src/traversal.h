#ifndef TRAVERSAL_H
#define TRAVERSAL_H

#include "c_inc.h"

extern int init_vote(INC_GRP * meta, VOTE_GRP * vote);
extern int find_bipartition(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote, int i);
extern int find_valid_subtree(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote);
extern int bfs_vote(INC_GRP * meta,  MAP_GRP * map,MST_GRP * mst, VOTE_GRP * vote, int i);
extern int find_insertion_edge(INC_GRP * meta, VOTE_GRP * vote);
#endif 