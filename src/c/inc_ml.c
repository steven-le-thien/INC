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
char TMP_FILE1[]                = "tool_tmp_file1";
char TMP_FILE2[]                = "tool_tmp_file2";

// int skip_counter = 0;
int make_constraint_trees(int * num_ctree,  ml_options  * master_ml_options){
  FILE *      f;
  FILE *      p;
  char        in_name[10000];
  char        out_name[10000];
  char        msa_name[10000];
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
            FCAL(
                GENERAL_ERROR,
                F_MAKE_ST_IN_MAKE_CTREE,
                make_subtree(
                    in_name,
                    out_name,
                    master_ml_options->init_tree_name
                )
            );
            break;

          case C_RAXML:
            FCAL(
                GENERAL_ERROR,
                F_MAKE_RAXML_IN_MAKE_CTREE,
                make_raxml_constraint(
                    msa_name,
                    out_name,
                    master_ml_options
                )
            );
            break;

          case C_FASTTREE:
            FCAL(
                GENERAL_ERROR,
                F_MAKE_FT_IN_MAKE_CTREE,
                make_fasttree_constraint(
                    msa_name,
                    out_name,
                    master_ml_options
                )
            );
            break;

          case C_NJ:
            FCAL(
                GENERAL_ERROR,
                F_MAKE_NJ_IN_MAKE_CTREE,
                make_nj_constraint(
                    msa_name,
                    out_name,
                    master_ml_options
                )
            );
            break;  

          case C_FASTME:
            FCAL(
                GENERAL_ERROR,
                F_MAKE_FASTME_IN_MAKE_CTREE,
                make_fastme_constraint(
                    msa_name,
                    out_name,
                    master_ml_options
                )
            );
            break;

          default:
            break;
        }
      }
      (*num_ctree)++;

    }
  }
  return 0;
}


int main(int argc, char ** argv){
  // Pipe in a bunch of programs together
  FILE *      f;
  char        name[MAX_BUFFER_SIZE];
  // char        name2[MAX_BUFFER_SIZE];
  int         num_ctree;

  ml_options master_ml_options;
  // option_t    tmp_options;


  if(read_ml_cmd_arg(argc, argv, &master_ml_options)    != SUCCESS)         PRINT_AND_EXIT("read_cmd_arg failed in main\n", GENERAL_ERROR);
    // Check if the initial tree is set up 
    if(master_ml_options.ctree_method != C_NO){
      // Piping into fasttree 
      printf("checking for initial tree...\n");
      if(!master_ml_options.init_tree_name){
        sprintf(name, "%sfirst_tree.tree", master_ml_options.output_prefix);
        master_ml_options.init_tree_name = name;
        
        f = fopen(name, "r");
        if(!f){
          // tmp_options.input_name = master_ml_options.input_alignment;
          // tmp_options.tree_names = malloc(sizeof(char*));
          // tmp_options.tree_names[0] = master_ml_options.init_tree_name;
    

          // if(fasttree_job(&tmp_options, &master_ml_options)           != SUCCESS)         PRINT_AND_EXIT("fasttree_job failed in main\n", GENERAL_ERROR);
          // if(nj_job(&tmp_options, &master_ml_options) != SUCCESS) PRINT_AND_EXIT("n_job failed in main\n", GENERAL_ERROR);
          // if(fastme_job(&tmp_options, &master_ml_options) != SUCCESS) PRINT_AND_EXIT("n_job failed in main\n", GENERAL_ERROR);
          // if(fasttree_initial_tree_job(&tmp_options, &master_ml_options)           != SUCCESS)         PRINT_AND_EXIT("fasttree_job failed in main\n", GENERAL_ERROR);
          // if(make_raxml_constraint(master_ml_options.input_alignment, master_ml_options.init_tree_name, &master_ml_options) != SUCCESS) PRINT_AND_EXIT("raxml job failedi n main\n", GENERAL_ERROR);
        } else fclose(f);
        
      }
      // Making constraint trees 
      if(make_constraint_trees(&num_ctree, &master_ml_options) != SUCCESS)           PRINT_AND_EXIT("make constraint trees failed in main\n", GENERAL_ERROR);
    } else num_ctree = 0;

    // return 0;

    // Check if the initial distance matrix is set up or need to set up
    if(!master_ml_options.init_d_name && !master_ml_options.use_initial_tree_as_spanning_tree && master_ml_options.use_distance_matrix){
      master_ml_options.init_d_name = malloc(1000);
      sprintf(master_ml_options.init_d_name, "%sc_inc_input", master_ml_options.output_prefix);
      f = fopen(master_ml_options.init_d_name, "r");
      if(!f){
        printf("writing distance matrix using PAUP* to %sc_inc_input...\n",  master_ml_options.output_prefix);
        if(distance_matrix_job(
          master_ml_options.distance_model, 
          master_ml_options.input_alignment,
          master_ml_options.init_d_name
        ) != SUCCESS)         PRINT_AND_EXIT("fasttree_job failed in main\n", GENERAL_ERROR);
      } else fclose(f);        
    }

    // Piping into constrained_inc code
    FCAL( 
        GENERAL_ERROR,
        "constrained_inc failed in main\n",
        constraint_inc(num_ctree, &master_ml_options)
    );
    // if(constraint_inc(num_ctree, &master_ml_options) != SUCCESS)        PRINT_AND_EXIT("constraint_inc failed in main", GENERAL_ERROR);
  
  return 0; 
}   


