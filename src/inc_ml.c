// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utilities.h"
#include "options.h"
#include "tools.h"
#include "msa.h"

#define MODE 1

int         use_raxml = 0;
// char        stock_path[] = "/Users/lethien_96/inc_ml_result/1000L1threshold_test";
char        stock_init_tree_name[] = "first_tree.tree";

int make_constraint_trees(int * num_ctree, option_t * options, msa_t * msa){
    FILE *      f;
    char        in_name[CMD_BUFFER_SIZE];
    char        out_name[CMD_BUFFER_SIZE];
    char        msa_name[CMD_BUFFER_SIZE];
    option_t    tmp_options;


    *num_ctree = 0;
    while(1){
        sprintf(in_name, "%s_ctree%d.lab", options->output_name, *num_ctree);
        f = fopen(in_name, "r");
        if(!f) break;
        else {
            sprintf(out_name, "%s_ctree%d.tree", options->output_name, *num_ctree);
            if(MODE){// 2.1 Pruning the tree to get the subtree
                tmp_options.input_name = options->tree_index == -1 ? stock_init_tree_name : options->tree_names[0];
                tmp_options.output_name = options->tree_index == -1 ? stock_init_tree_name : options->tree_names[0];
                if(rm_label_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("remove label failed in main", GENERAL_ERROR);

                tmp_options.input_name = in_name;
                tmp_options.output_name = out_name;
                tmp_options.tree_names = malloc(sizeof(char *));
                tmp_options.tree_names[0] = options->tree_index == -1 ? stock_init_tree_name : options->tree_names[0];
                // printf("%s %s\n", tmp_options.tree_names[0], options->tree_names[0]);
                if(nw_utils_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("nw_utils failed in main", GENERAL_ERROR);
            } else {// 2.2 Getting the msa and uses FastTree from them
                sprintf(msa_name, "%s_ctree%d.msa", options->output_name, *num_ctree);
                subset_msa(in_name, msa_name, msa);

                if(use_raxml){
                    tmp_options.input_name = msa_name;

                    // HARDCODING 
                    tmp_options.output_name = malloc(1000 * sizeof(char));
                    sprintf(tmp_options.output_name, "200_ctree%d.raxmltree", *num_ctree);
                    tmp_options.tree_names = malloc(sizeof(char *));
                    tmp_options.tree_names[0] = malloc(1000 * sizeof(char));
                    strncpy(tmp_options.tree_names[0], options->output_name, strlen(options->output_name - 3));
                    tmp_options.tree_names[0][strlen(options->output_name - 3)] = 0;
                    // END HARDCODING

                    // if(raxml_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("raxml_job failed in main", GENERAL_ERROR);
                } else {
                    tmp_options.input_name = msa_name;
                    tmp_options.output_name = malloc(sizeof(char));
                    tmp_options.output_name[0] = 0;
                    tmp_options.tree_names = malloc(sizeof(char *));
                    tmp_options.tree_names[0] = out_name;

                    if(fasttree_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("fasttree failed in main", GENERAL_ERROR);
                }

                tmp_options.input_name = out_name;
                tmp_options.output_name = out_name;
                if(rm_label_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("remove label failed in main", GENERAL_ERROR);
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
    option_t    options;
    msa_t       msa;
    float **    d;
    int         num_ctree;

    if(init_options(&options)               != SUCCESS)         PRINT_AND_EXIT("init_options failed in main\n", GENERAL_ERROR);
    if(read_cmd_arg(argc, argv, &options)   != SUCCESS)         PRINT_AND_EXIT("read_cmd_arg failed in main\n", GENERAL_ERROR);

    // Piping into fasttree 
    printf("checking for initial tree\n");
    if(options.tree_index == -1){
        options.tree_names = malloc(sizeof(char *));
        options.tree_names[0] = malloc(MAX_BUFFER_SIZE * sizeof(char));
        sprintf(options.tree_names[0], "%sfirst_tree.tree", options.output_name);
        if(fasttree_job(&options)           != SUCCESS)         PRINT_AND_EXIT("fasttree_job failed in main\n", GENERAL_ERROR);
    }

    printf("parsing input.. \n");
    parse_input(&msa, options.input_name);


    sprintf(name, "%sc_inc_input", options.output_name);
    f = fopen(name, "r");
    if(!f){
        printf("computing distance ...\n");
        compute_k2p_distance(&msa, &d);

        printf("writing matrix ...\n");
        write_distance_matrix(d, &options, &msa);
    }

    // Making constraint trees using PASTA code
    if(subset_job(&options)                 != SUCCESS)         PRINT_AND_EXIT("subset job failed in main\n", GENERAL_ERROR);
    make_constraint_trees(&num_ctree, &options, &msa);

    // Piping into constrained_inc code
    constraint_inc(num_ctree, &options);

    return 0; 
}   

