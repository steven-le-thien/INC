// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utilities.h"
#include "options.h"
#include "tools.h"
#include "msa.h"

int make_constraint_trees(int * num_ctree, option_t * options){
    FILE *      f;
    FILE *      p;
    char        in_name[10000];
    char        out_name[10000];

    
    // Recomputing the constraint trees if necessary
#if recompute_constraint_trees
    printf("performing PASTA decomposition\n");
    if(make_subset_label(options->tree_names[0], options->output_name) != SUCCESS)        PRINT_AND_RETURN("make_subset_label failed in main\n", GENERAL_ERROR); 
#if !(use_subtree_for_constraint_trees)
    char        msa_name[10000];
    msa_t       msa;

    parse_input(&msa, options->input_name);
#endif
#endif

    *num_ctree = 0;
    while(1){
        sprintf(in_name, "%s_ctree%d.lab", options->output_name, *num_ctree);
        sprintf(out_name, "%s_ctree%d.tree", options->output_name, *num_ctree);

        f = fopen(in_name, "r");
        p = fopen(out_name, "r");
        if(!f && !p)
            break;
        else {
            fclose(f);
            fclose(p);
#if recompute_constraint_trees
            // sprintf(out_name, "%s_ctree%d.tree", options->output_name, *num_ctree);
    #if use_subtree_for_constraint_trees
            // Just find the subtree
            if(make_subtree(in_name, out_name, options->tree_names[0]) != SUCCESS) PRINT_AND_RETURN("make subtree faield in main\n", GENERAL_ERROR);
    #else
            // Find the subalignment 
            sprintf(msa_name, "%s_ctree%d.msa", options->output_name, *num_ctree);
            if(subset_msa(in_name, msa_name, &msa) != SUCCESS) PRINT_AND_RETURN("make subset msa failed in main\n", GENERAL_ERROR);

        #if use_raxml_for_constraint_trees
            if(make_raxml_constraint(options->output_name, msa_name, out_name) != SUCCESS) PRINT_AND_RETURN("make raxml constraint failed in main\n", GENERAL_ERROR);

        #elif use_fasttree_for_constraint_trees
            if(make_fasttree_constraint(msa_name, out_name) != SUCCESS) PRINT_AND_RETURN("make fasttree constraint failed in main \n", GENERAL_ERROR);

        #else 
        #endif  // use_raxml_for_constraint_trees
    #endif  // use_subtree_for_constraint_trees         
#endif // recompute_constraint_trees
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
    int         num_ctree;

    if(init_options(&options)               != SUCCESS)         PRINT_AND_EXIT("init_options failed in main\n", GENERAL_ERROR);
    if(read_cmd_arg(argc, argv, &options)   != SUCCESS)         PRINT_AND_EXIT("read_cmd_arg failed in main\n", GENERAL_ERROR);

    if(use_constraint){
        // Piping into fasttree 
        printf("checking for initial tree...\n");
        if(options.tree_index == -1){
            sprintf(name, "%sfirst_tree.tree", options.output_name);
            options.tree_names = malloc(sizeof(char*));
            options.tree_names[0] = name;
            printf("%s\n", options.tree_names[0]);
            f = fopen(name, "r");
            if(!f){
                if(fasttree_job(&options)           != SUCCESS)         PRINT_AND_EXIT("fasttree_job failed in main\n", GENERAL_ERROR);
            }
            
        }
        // Making constraint trees 
        make_constraint_trees(&num_ctree, &options);
    } else num_ctree = 0;
    

    sprintf(name, "%sc_inc_input", options.output_name);
    f = fopen(name, "r");
    if(!f){
        printf("writing distance matrix using PAUP*...\n");
        distance_matrix_job(&options);
    }

    // Piping into constrained_inc code
    constraint_inc(num_ctree, &options);

    return 0; 
}   

