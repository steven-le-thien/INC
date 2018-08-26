// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utilities.h"
#include "options.h"
#include "tools.h"
#include "msa.h"



ml_options master_ml_options;
// int skip_counter = 0;

int make_constraint_trees(int * num_ctree){
    FILE *      f;
    FILE *      p;
    char        in_name[10000];
    char        out_name[10000];
    char        msa_name[10000];
    msa_t       msa;
    
    // Recomputing the constraint trees if necessary
    if(master_ml_options.recompute_constraint_trees){
        printf("performing PASTA decomposition\n");
        if(make_subset_label(master_ml_options.init_tree_name, master_ml_options.output_prefix, &master_ml_options) != SUCCESS)        PRINT_AND_RETURN("make_subset_label failed in main\n", GENERAL_ERROR); 
        if(!(master_ml_options.use_subtree_for_constraint_trees)){
            parse_input(&msa, master_ml_options.input_alignment);
        }
    }

    *num_ctree = 0;
    while(1){
        sprintf(in_name, "%s_ctree%d.lab", master_ml_options.output_prefix, *num_ctree);
        sprintf(out_name, "%s_ctree%d.tree", master_ml_options.output_prefix, *num_ctree);

        f = fopen(in_name, "r");
        p = fopen(out_name, "r");
        if(!f && !p)
            break;
        else {
            fclose(f);
            fclose(p);
            if(master_ml_options.recompute_constraint_trees){
                if(master_ml_options.use_subtree_for_constraint_trees){
                    if(make_subtree(in_name, out_name, master_ml_options.init_tree_name) != SUCCESS) PRINT_AND_RETURN("make subtree faield in main\n", GENERAL_ERROR);
                } else {
                    sprintf(msa_name, "%s_ctree%d.msa", master_ml_options.output_prefix, *num_ctree);
                    if(subset_msa(in_name, msa_name, &msa) != SUCCESS) PRINT_AND_RETURN("make subset msa failed in main\n", GENERAL_ERROR);

                    if(master_ml_options.use_raxml_for_constraint_trees){
                        if(make_raxml_constraint(msa_name, out_name, &master_ml_options) != SUCCESS) PRINT_AND_RETURN("make raxml constraint failed in main\n", GENERAL_ERROR);
                    } else if (master_ml_options.use_fasttree_for_constraint_trees){
                        if(make_fasttree_constraint(msa_name, out_name,  &master_ml_options) != SUCCESS) PRINT_AND_RETURN("make fasttree constraint failed in main \n", GENERAL_ERROR);
                    }
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
    char        name2[MAX_BUFFER_SIZE];
    option_t    tmp_options;
    int         num_ctree;

    if(read_ml_cmd_arg(argc, argv, &master_ml_options)    != SUCCESS)         PRINT_AND_EXIT("read_cmd_arg failed in main\n", GENERAL_ERROR);

    // Check if the initial tree is set up 
    if(master_ml_options.use_constraint){
        // Piping into fasttree 
        printf("checking for initial tree...\n");
        if(!master_ml_options.init_tree_name){
            sprintf(name, "%sfirst_tree.tree", master_ml_options.output_prefix);
            master_ml_options.init_tree_name = name;
            
            f = fopen(name, "r");
            if(!f){
                tmp_options.input_name = master_ml_options.input_alignment;
                tmp_options.tree_names = malloc(sizeof(char*));
                tmp_options.tree_names[0] = master_ml_options.init_tree_name;
                if(fasttree_job(&tmp_options, &master_ml_options)           != SUCCESS)         PRINT_AND_EXIT("fasttree_job failed in main\n", GENERAL_ERROR);
                // if(fasttree_initial_tree_job(&tmp_options, &master_ml_options)           != SUCCESS)         PRINT_AND_EXIT("fasttree_job failed in main\n", GENERAL_ERROR);
                // if(make_raxml_constraint(master_ml_options.input_alignment, master_ml_options.init_tree_name, &master_ml_options) != SUCCESS) PRINT_AND_EXIT("raxml job failedi n main\n", GENERAL_ERROR);
            } else fclose(f);
            
        }
        // Making constraint trees 
        if(make_constraint_trees(&num_ctree) != SUCCESS)           PRINT_AND_EXIT("make constraint trees failed in main\n", GENERAL_ERROR);
    } else num_ctree = 0;

    // Check if the initial distance matrix is set up
    if(!master_ml_options.init_d_name){
        sprintf(name2, "%sc_inc_input", master_ml_options.output_prefix);
        master_ml_options.init_d_name = name2;
        
        f = fopen(name2, "r");
        if(!f){
            printf("writing distance matrix using PAUP*...\n");
            tmp_options.input_name = master_ml_options.input_alignment;
            tmp_options.output_name = master_ml_options.init_d_name;
            if(distance_matrix_job(&tmp_options, &master_ml_options)          != SUCCESS)         PRINT_AND_EXIT("fasttree_job failed in main\n", GENERAL_ERROR);
        } else fclose(f);        
    }

    // Piping into constrained_inc code
    if(constraint_inc(num_ctree, &master_ml_options) != SUCCESS)        PRINT_AND_EXIT("constraint_inc failed in main", GENERAL_ERROR);
    return 0; 
}   

