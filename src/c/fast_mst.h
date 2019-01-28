#ifndef FAST_MST_H
#define FAST_MST_H

#include "c_inc.h"

extern int fast_mst(char ** data, int n, DIST_MOD distance_model, int seed, MST_GRP * mst, int *** disjoint_subset);

static const char F_RANDOM_CENTROID_IN_NN_FROM_SEQ[]  
  = "random centroid failed in nn_from_seq\n"; 

static const char STATE_DONE_CENTROID[] 
  = "done with centroid...\n";

static const char F_EXT_SUBSET_IN_NN_FROM_SEQ[]
  = "extending subset failed in nn from seq\n";

static const char STATE_DONE_EXTEND_SS[]
  = "done with extending subsets...\n";

static const char F_MAKE_CLUSTER_CLIQUES_IN_NN_FROM_SEQ[]
  = "make clusters cliques faield in nn from seq\n";

static const char STATE_DONE_MAKE_CLUSTER_CLIQUES[]
  = "done with making clusters cliques\n";

static const char F_MAKE_CENTROID_CLIQUES_IN_NN_FROM_SEQ[]
  = "make centroid cliques failed in nn_from_seq\n";

static const char F_MAKE_COMPLETE_BP_IN_NN_FROM_SEQ[]
  = "make complete bipartite failed in nn_from_seq\n";

static const char F_INIT_GRAPH_IN_FAST_MST[] 
  = "init graph failed in fast mst\n";

static const char F_NN_FROM_SEQ_IN_FAST_MST[]
  = "nn from seq failed in fast mst\n";

static const char STATE_DONE_NN[]
  = "done with nn...\n";

static const char F_PRIM_ON_SMALL_GRAPH[]
  = "prim on small graph failed\n";

static const char STATE_DONE_PRIM[]
  = "done with prim...\n";

#endif 