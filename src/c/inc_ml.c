// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utilities.h"
#include "options.h"
#include "tools.h"
#include "inc_ml.h"
#include "msa.h"

// Some extern initialization
char TMP_FILE1[]                = "./tool_tmp_file1";
char TMP_FILE2[]                = "./tool_tmp_file2";

// Internal functions
int set_up_trees_to_inc(ml_options * master_ml_options, int * num_ctree);
int set_up_itree(ml_options * master_ml_options);
int make_constraint_trees(int * num_ctree,  ml_options  * master_ml_options);
int set_up_dist(ml_options * master_ml_options);

int main(int argc, char ** argv){
  // Pipe in a bunch of programs together
  int         num_ctree;
  ml_options master_ml_options;

  FCAL(
      GENERAL_ERROR,
      F_RD_CMD_ARG,
      read_ml_cmd_arg(argc, argv, &master_ml_options)
  );

  FCAL(
      GENERAL_ERROR,
      F_SET_UP_TREES_IN_MAIN,
      set_up_trees_to_inc(&master_ml_options, &num_ctree)
  );

  // Check if the initial distance matrix is set up or need to set up
  FCAL(
      GENERAL_ERROR,
      F_SET_UP_DIST_IN_MAIN,
      set_up_dist(&master_ml_options)
  );
  
  // Piping into constrained_inc code
  FCAL( 
      GENERAL_ERROR,
      F_CINC_IN_MAIN,
      constraint_inc(num_ctree, &master_ml_options)
  );

  return 0; 
}   

int set_up_trees_to_inc(ml_options * master_ml_options, int * num_ctree){
  FCAL(
      GENERAL_ERROR,
      F_SET_UP_I_TREE_IN_SET_UP_TREE,
      set_up_itree(master_ml_options)
  );

  // Making constraint trees 
  if(master_ml_options->ctree_method != C_NO)
    FCAL(
        GENERAL_ERROR,
        F_MAKE_CTREE_IN_SET_UP_TREE, 
        make_constraint_trees(num_ctree, master_ml_options)
    );
  else (*num_ctree) = 0;
  return 0;
}

int set_up_itree(ml_options * master_ml_options){
  FILE        *f;
  int         (*itree_job) (char * , char * , ml_options *, int);     
  char        name[MAX_BUFFER_SIZE];

  // Piping into fasttree 
  printf(STATE_INIT_TREE);
  if(!master_ml_options->init_tree_name){

    sprintf(
        master_ml_options->init_tree_name 
            = SAFE_MALLOC(MAX_BUFFER_SIZE * sizeof(char)), 
        "%s%s", 
        master_ml_options->output_prefix, 
        DEFAULT_FT_SUF
    );
    
    f = fopen(name, "r");
    if(!f){ // Call the tree constructor
      switch(master_ml_options->itree_method){
        case I_FASTTREE:
          itree_job = make_fasttree_constraint;
          break;
        case I_RAXML:
          itree_job = make_raxml_constraint;
          break; 
        case I_NJ:
          itree_job = make_nj_constraint;
          break; 
        case I_FASTME:
          itree_job = make_fastme_constraint;
          break;
      }

      FCAL(
          GENERAL_ERROR, 
          F_ITREE_IN_MAIN,
          itree_job(
              master_ml_options->input_alignment,
              master_ml_options->init_tree_name,
              master_ml_options,
              0
          )
      );
    } else fclose(f);
  }
  return 0; 
}

// int skip_counter = 0;
int make_constraint_trees(int * num_ctree,  ml_options  * master_ml_options){
  FILE *      f;
  FILE *      p;
  int (*ctree_job) (char * , char * , ml_options *, int); 
 
  char        in_name[MAX_BUFFER_SIZE];
  char        out_name[MAX_BUFFER_SIZE];
  char        msa_name[MAX_BUFFER_SIZE];
  msa_t       msa;

  // Recomputing the constraint trees if necessary
  if(master_ml_options->recompute_constraint_trees){

    printf(STATE_PASTA_DEC);

    FCAL(
        GENERAL_ERROR,
        F_MAKE_SS_LBL_IN_MAKE_CTREE,
        make_subset_label(
            master_ml_options->init_tree_name,
            master_ml_options->output_prefix, 
            master_ml_options
        )
    ); 

    if(master_ml_options->ctree_method != C_SUBTREE)
      FCAL(
          GENERAL_ERROR,
          F_PARSE_IN_IN_MAKE_CTREE,
          parse_input(
            &msa, 
            master_ml_options->input_alignment
          )
      );   
  }

  *num_ctree = 0;
  while(1){
    sprintf(
        in_name, 
        "%s_ctree%d.lab", 
        master_ml_options->output_prefix, 
        *num_ctree
    );

    sprintf(
        out_name, 
        "%s_ctree%d.tree", 
        master_ml_options->output_prefix, 
        *num_ctree
    );

    sprintf(
        msa_name, 
        "%s_ctree%d.msa", 
        master_ml_options->output_prefix, 
        *num_ctree
    );

    f = fopen(in_name, "r");
    p = fopen(out_name, "r");

    // Checking if files already exists
    if(!f && !p)
      break;
    else{
      if(f) fclose(f);
      if(p) fclose(p);

      if(master_ml_options->recompute_constraint_trees){
        if(master_ml_options->ctree_method != C_SUBTREE)
          FCAL(
              GENERAL_ERROR,
              F_SS_MSA_IN_MAKE_CTREE,
              subset_msa(in_name, msa_name, &msa)
          );

        switch(master_ml_options->ctree_method){
          case C_SUBTREE:
            ctree_job = make_subtree;
            break;
          case C_RAXML:
            ctree_job = make_raxml_constraint;
            break;
          case C_FASTTREE:
            ctree_job = make_fasttree_constraint;
            break;
          case C_NJ:
            ctree_job = make_nj_constraint;
            break;
          case C_FASTME:
            ctree_job = make_fastme_constraint;
            break;
          default:
            break;
        }

        FCAL(
            GENERAL_ERROR, 
            F_CTREE_IN_MAIN,
            ctree_job(
                master_ml_options->ctree_method == C_SUBTREE ? in_name : msa_name,
                out_name,
                master_ml_options,
                1
            )
        );
      }
      (*num_ctree)++;
    }
  }

  return 0;
}

int set_up_dist(ml_options * master_ml_options){
  FILE * f;
  if(!master_ml_options->init_d_name && 
      !master_ml_options->use_initial_tree_as_spanning_tree && 
      master_ml_options->use_distance_matrix){

    sprintf(
        master_ml_options->init_d_name 
            = SAFE_MALLOC(MAX_BUFFER_SIZE * sizeof(char)), 
        "%s%s", 
        master_ml_options->output_prefix, 
        DEFAULT_DIST_SUF
    );

    f = fopen(master_ml_options->init_d_name, "r");
    if(!f){
      printf(
          "%s%s\n",  
          STATE_DIST, 
          master_ml_options->init_d_name 
      );
      FCAL(
          GENERAL_ERROR,
          F_DIST_MAT_IN_MAIN,
          distance_matrix_job(
              master_ml_options->tmp_folder,
              master_ml_options->distance_model, 
              master_ml_options->input_alignment,
              master_ml_options->init_d_name
          )
      ); 
    } else fclose(f);        
  }

  return 0;
}


