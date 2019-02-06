// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "inc_ml.h"
#include "msa.h"
#include "c_inc.h"
#include "tree.h"
#include "utilities.h"
#include "options.h"
#include "tools.h"
#include "dist.h"
#include "prim.h"
#include "traversal.h"
#include "fast_mst.h"

#define SEED 12345

int serial_main_loop(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst);
int init_meta_with_msa(msa_t * msa, INC_GRP * meta, MAP_GRP * map);


// Implementation of constrained version of the INC algorithm. 
// Command line argument is as followed: constraint_inc -i <alignment file> 
//      -t <tree1> <tree2> ... 
int constraint_inc_main(int argc, char ** argv, ml_options * master_ml_options){
  // Meta variables
  INC_GRP     meta;
  MAP_GRP     map;
  MST_GRP     mst;

  // No distance matrix tmp variables
  msa_t msa;
  // int ** disjoint_subset; 
  meta.master_ml_options = master_ml_options;

  // Parse options
  printf(STATE_RD_OPTIONS_IN_CINC);
  FCAL(
      GENERAL_ERROR,
      F_RD_CMD_ARG_IN_CINC,
      read_cmd_arg(argc, argv, master_ml_options)
  );

  if(!meta.master_ml_options->use_initial_tree_as_spanning_tree){
    if(meta.master_ml_options->use_distance_matrix){
      // Getting the distance matrix
      printf(STATE_PARSE_MAT_TREE);
      FCAL(
          GENERAL_ERROR,
          F_PARSE_MAT_TREE_IN_CINC,
          parse_distance_matrix(&meta, &map, master_ml_options)
      );
    } else {    
      printf(STATE_NO_DIST);
      FCAL(
          GENERAL_ERROR,
          F_PARSE_INPUT_IN_CINC,
          parse_input(&msa, meta.master_ml_options->input_alignment)
      );
      FCAL(
          GENERAL_ERROR,
          F_INIT_META_IN_CINC,
          init_meta_with_msa(&msa, &meta, &map)
      );
    }

    // Compute the MST 
    printf(STATE_MST);
    FCAL(
        GENERAL_ERROR, 
        F_PRIM_IN_CINC,
        prim(&meta, &mst)
    );
  } else
    FCAL(
        GENERAL_ERROR, 
        F_PARSE_INIT_AS_MST_IN_CINC,
        parse_initial_tree_as_mst(&meta, &mst)
    );

  // At this stage, all tree names should be in options->tree_name. 
  // It should be ok to parse them all at once since together they 
  //  have at most 4M nodes
  FCAL(
      GENERAL_ERROR, 
      F_PARSE_TREE_IN_CINC,
      parse_tree(&meta, &map, master_ml_options)
  );
  
  // Initialize growing tree using the first 3 taxa in the ordering
  printf(STATE_INIT_GTREE);
  FCAL(
      GENERAL_ERROR, 
      F_INIT_GTREE_IN_CINC,
      init_growing_tree(&meta, &map, &mst)
  );

  // Loop through Prim's ordering
  printf(STATE_BUILD_TREE);
  FCAL(
      GENERAL_ERROR, 
      F_SERIAL_MAIN_LOOP_IN_CINC,
      serial_main_loop(&meta, &map, &mst)
  );

  // Report the growing tree
  printf(STATE_OUT_TREE);
  FCAL(
      GENERAL_ERROR, 
      F_WRITE_TREE_IN_CINC,
      write_newick(
          meta.gtree, 
          master_ml_options->output_prefix, 
          map.master_to_name
      )
  );

  printf(STATE_CLEAN);
  // Clean up
  return 0;
}

int serial_main_loop(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst){
  int i, j; //loop counter
  // Init voting structure
  VOTE_GRP    vote;
  printf(ITER_COUNT);
  for(i = 3; i < meta->n_taxa; i++){
    print_inline_iteration(i, j, meta->n_taxa, 3);

    FCAL(
        GENERAL_ERROR,
        F_INIT_VOTE_IN_CINC, 
        init_vote(meta, map, mst, &vote, i)
    );

    FCAL(
        GENERAL_ERROR,
        F_FIND_BIPART_IN_CINC, 
        find_bipartition(meta, map, mst, &vote, i)
    );

    FCAL(
        GENERAL_ERROR,
        F_FIND_VALID_ST_IN_CINC,
        find_valid_subtree(meta, map, mst, &vote)
    );

    FCAL(
        GENERAL_ERROR,
        F_BFS_VOTE_IN_CINC,
        bfs_vote(meta, map, mst, &vote, i)
    );

    FCAL(
        GENERAL_ERROR,
        F_ATTACH_IN_CINC,
        attach_leaf_to_edge(meta, map, mst, &vote, i)
    );
  }
  printf("\n");
  return 0;
}

int init_meta_with_msa(msa_t * msa, INC_GRP * meta, MAP_GRP * map)
{
  int i;

  meta->n_taxa = map->n_taxa = msa->num_seq;
  meta->msa = msa;
  map->master_to_name         = malloc(msa->num_seq * sizeof(char*));
  for(i = 0; i < msa->num_seq; i++)
    map->master_to_name[i]  = msa->name[i];
  return 0;
}
