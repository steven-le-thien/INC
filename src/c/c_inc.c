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

int serial_main_loop(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote);
int init_meta_with_msa(msa_t * msa, INC_GRP * meta, MAP_GRP * map);


// Implementation of constrained version of the INC algorithm. 
// Command line argument is as followed: constraint_inc -i <alignment file> 
//      -t <tree1> <tree2> ... 
int constraint_inc_main(int argc, char ** argv, ml_options * master_ml_options){
  // Meta variables
  INC_GRP     meta;
  MAP_GRP     map;
  MST_GRP     mst;

  // Init voting structure
  VOTE_GRP    vote;

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

  FCAL(
      GENERAL_ERROR,
      "voting initialization failed in main",
      pass_vote_options(master_ml_options, &vote)
  );

  // Loop through Prim's ordering
  printf(STATE_BUILD_TREE);
  FCAL(
      GENERAL_ERROR, 
      F_SERIAL_MAIN_LOOP_IN_CINC,
      serial_main_loop(&meta, &map, &mst, &vote)
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

int serial_main_loop(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote){
  int i, j; //loop counter

  printf(ITER_COUNT);
  for(i = 3; i < meta->n_taxa; i++){
    print_inline_iteration(i, j, meta->n_taxa, 3);

    FCAL(
        GENERAL_ERROR,
        F_INIT_VOTE_IN_CINC, 
        init_vote(meta, map, mst, vote, i)
    );

    FCAL(
        GENERAL_ERROR,
        F_FIND_BIPART_IN_CINC, 
        find_bipartition(meta, map, mst, vote, i)
    );

    FCAL(
        GENERAL_ERROR,
        F_FIND_VALID_ST_IN_CINC,
        find_valid_subtree(meta, map, mst, vote)
    );

    FCAL(
        GENERAL_ERROR,
        F_BFS_VOTE_IN_CINC,
        bfs_vote(meta, map, mst, vote, i)
    );

    FCAL(
        GENERAL_ERROR,
        F_ATTACH_IN_CINC,
        attach_leaf_to_edge(meta, map, mst, vote, i)
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


///////////////////
//  Main function for INC (Zhang, Rao, Warnow)
//
//
/////////////////
#if INC_CMPL

char F_MAIN[] = "inc failed in main\n";

// Commands are of the form inc -i <distance_matrix> -o <output_path>
int main(int argc, char ** argv){
  int i;
  int argc_in = 5; 
  char * argv_in[argc_in];
  ml_options master_ml_options;

  // There should be an odd number of arguments
  ASSERT(
      GENERAL_ERROR, 
      F_RD_CMD_ARG_IN_CINC,
      argc % 2 == 1
  );

  FCAL(
      GENERAL_ERROR,
      F_INIT_ML_OPT,
      init_ml_options(&master_ml_options)
  ); 

  // Parse user input
  for(i = 1; i < argc; i += 2){ // ignore the first word
    ASSERT(
        GENERAL_ERROR,
        F_RD_CMD_ARG_IN_CINC,
        argv[i][0] == '-'
    );

    switch(argv[i][1]){
      case 'i':
        master_ml_options.init_d_name = argv[i + 1]; 
        argv_in[1] = argv[i];
        argv_in[2] = argv[i + 1];
        break;
      case 'o':
        master_ml_options.output_prefix = argv[i + 1];
        argv_in[3] = argv[i];
        argv_in[4] = argv[i + 1];
        break;
      default:
        break;
    }
  }

  argv[0] = NULL;

  master_ml_options.qtree_method = Q_FPM;

  FCAL(
      GENERAL_ERROR,
      F_MAIN,
      constraint_inc_main(
          argc_in,
          argv_in,
          &master_ml_options
      )
  );

  return 0;
}


#endif //INC_CMPL


///////////////////
//  Main function for cosntraint INC (Zhang, Rao, Warnow)
//  Commands are of the form inc -i <distance_matrix> -o <output_path> 
//      -t <trees>.. -q <quartet_type> -g <guide_tree>
//  where -q is either `fpm' or `subtree' and -g must be provided if it is 
//  quartet type is `subtree'
//
/////////////////
#if CINC_CMPL
char F_MAIN[] = "cinc failed in main\n";
char F_ARGC_TOO_MANY[] = "too many arguments in the command\n";
char F_SUBTREE_NO_GUIDE[] = "using subtree quartet trees without guide tree\n";

int main(int argc, char ** argv){
  int argc_in; // capped at 10000 arguments
  char * argv_in[MAX_NUM_FLAG];
  int i, guide_tree_idx;
  int just_seen_tree_flag = 0;
  ml_options master_ml_options;

  // Prepare master_ml_options
  FCAL(
      GENERAL_ERROR,
      F_INIT_ML_OPT,
      init_ml_options(&master_ml_options)
  ); 

  // Check if subtree or fpm
  for(i = 0; i < argc; i++){
    if(argv[i][0] == '-' && argv[i][1] == 'q'){
      if(STR_EQ(argv[i + 1], "subtree"))
        master_ml_options.qtree_method = Q_SUBTREE;
      else if(STR_EQ(argv[i + 1], "fpm"))
        master_ml_options.qtree_method = Q_FPM;
      else{
        printf("wrong flag inputted into -q\n");
        return 0;
      }
    }
  }

  // If it is subtree then make sure that there is a guide tree that follows
  if(master_ml_options.qtree_method == Q_SUBTREE){
    for(i = 0; i < argc; i++)
      if(argv[i][0] == '-' && argv[i][1] == 'g')
        master_ml_options.guide_tree_name = argv[i + 1];

    ASSERT(GENERAL_ERROR, F_SUBTREE_NO_GUIDE, master_ml_options.guide_tree_name);

    for(guide_tree_idx = 0; guide_tree_idx < argc; guide_tree_idx++){
      if(argv[guide_tree_idx][0] == '-' && argv[guide_tree_idx][1] == 'g'){
        guide_tree_idx++;
        break;
      }
    }
  }

  // Set up argc_in and argv_in
  ASSERT(
      GENERAL_ERROR,
      F_ARGC_TOO_MANY,
      argc < MAX_NUM_FLAG
  );
  argc_in = 1;
  argv_in[0] = NULL;
  for(i = 1; i < argc; i++){
    if(STR_EQ(argv[i - 1], "-t")) just_seen_tree_flag = 1;
    if(STR_EQ(argv[i], "-q") || STR_EQ(argv[i], "-g"))
      i++; //skip these flags
    else if(master_ml_options.qtree_method == Q_SUBTREE 
        && just_seen_tree_flag == 1)
      argv_in[argc_in++] = argv[guide_tree_idx];
    else
      argv_in[argc_in++] = argv[i];
    if(just_seen_tree_flag == 1) just_seen_tree_flag = 2;
  }

  //never seen it
  if(master_ml_options.qtree_method == Q_SUBTREE && !just_seen_tree_flag){
    argv_in[argc_in] = malloc(3);
    argv_in[argc_in][0] = '-'; argv_in[argc_in][1] = 't'; argv_in[argc_in][2] = 0;
    argc_in++;
    argv_in[argc_in++] = argv[guide_tree_idx]; 
  }

  // Prepare argv to pass into

  argv[0] = NULL;
  FCAL(
      GENERAL_ERROR,
      F_MAIN,
      constraint_inc_main(
          argc_in,
          argv_in,
          &master_ml_options
      )
  );

  return 0;
}

#endif // CINC_CMPL




