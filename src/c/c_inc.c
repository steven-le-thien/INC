// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

#define COMPILE_ALONE 0
#define SEED 12345

int serial_main_loop(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst);

// int skip_counter = 0;

int init_meta_with_msa(msa_t * msa, INC_GRP * meta, MAP_GRP * map){
    int i;

    meta->n_taxa = map->n_taxa = msa->num_seq;
    meta->msa = msa;
    map->master_to_name         = malloc(msa->num_seq * sizeof(char*));
    for(i = 0; i < msa->num_seq; i++)
        map->master_to_name[i]  = msa->name[i];
    return 0;
}


// Implementation of constrained version of the INC algorithm. 
// Command line argument is as followed: constraint_inc -i <alignment file> -t <tree1> <tree2> ... 
int constraint_inc_main(int argc, char ** argv, ml_options * master_ml_options){
    // Meta variables
    // option_t    options;

    INC_GRP     meta;
    MAP_GRP     map;
    MST_GRP     mst;

    // No distance matrix tmp variables
    msa_t msa;
    // int ** disjoint_subset; 

    meta.master_ml_options = master_ml_options;

    // Parse options
    printf("reading in options...\n");

    if(read_cmd_arg(argc, argv, master_ml_options)   != SUCCESS)         PRINT_AND_EXIT("read_cmd_arg failed in main\n", GENERAL_ERROR);
    if(!meta.master_ml_options->use_initial_tree_as_spanning_tree){
        if(meta.master_ml_options->use_distance_matrix){
            // Getting the distance matrix
            printf("parsing initial matrix and tree...\n");
            if(parse_distance_matrix(&meta, &map, master_ml_options) 
                                                    != SUCCESS)         PRINT_AND_EXIT("parse_distance_matrix failed in main\n", GENERAL_ERROR);

            
        } else {    
            printf("without distance matrix...\n");
            if(parse_input(&msa, meta.master_ml_options->input_alignment)
                                                    != SUCCESS)         PRINT_AND_EXIT("read data failed in main\n", GENERAL_ERROR);

            if(init_meta_with_msa(&msa, &meta, &map)!= SUCCESS)         PRINT_AND_EXIT("init meta with msa failed in main\n", GENERAL_ERROR);

            // Used to do FastMST her
            // printf("doing fast_mst...\n");
            // if(fast_mst(msa.msa, meta.n_taxa, meta.master_ml_options->distance_modelmeta.master_ml_options->distance_model, SEED, &mst, &disjoint_subset)
                                                    // != SUCCESS)         PRINT_AND_EXIT("fast_mst failed in main\n", GENERAL_ERROR); 

            // printf("doing constraint trees...\n");
            // if(make_constraint_trees_from_disjoint_subsets(meta.n_taxa, &msa, disjoint_subset, meta.master_ml_options) 
                                                    // != SUCCESS)         PRINT_AND_EXIT("make_constraint_trees_from_disjont_subsets failed\n", GENERAL_ERROR);
        }

        // Compute the MST 
        printf("computing the mst...\n");
        if(prim(&meta, &mst)                        != SUCCESS)         PRINT_AND_EXIT("prim algorithm failed in main\n", GENERAL_ERROR);
    } else {
        if(parse_initial_tree_as_mst(&meta, &mst)   != SUCCESS)         PRINT_AND_EXIT("parse_initial_tree_as_mst failed in main\n", GENERAL_ERROR);
    }
    

    
                                                                                                // printf("debug: the following print out prim's ctree\n");

                                                                                                // for(int i = 0; i < meta.n_taxa; i++)
                                                                                                //     printf("%d ", map.master_to_ctree[mst.prim_ord[i]]);
                                                                                                // printf("\n");

    // At this stage, all tree names should be in options->tree_name. It should be ok to parse them all at once since together they have at most 4M nodes
    if(parse_tree(&meta, &map, master_ml_options)    != SUCCESS)         PRINT_AND_EXIT("parse tree failed in main\n", GENERAL_ERROR);

    
    // Initialize growing tree using the first 3 taxa in the ordering
    printf("initializing the growing tree...\n");
    if(init_growing_tree(&meta, &map, &mst) != SUCCESS)         PRINT_AND_EXIT("init_growing_tree failed in main\n", GENERAL_ERROR);


                                                                                            #if 1 
    // int i, j;
                                                                                                // for(i = 0; i < 500; i++){
                                                                                                //     for(j = 0; j < 500; j++){
                                                                                                //         printf("%f ", meta.dm[i][j]);
                                                                                                //     }
                                                                                                //     printf("\n");
                                                                                                // }
                                                                                                // printf("\n"); 
                                                                                            #endif


    printf("building the tree...\n");
    // Loop through Prim's ordering
    if(serial_main_loop(&meta, &map, &mst)  != SUCCESS)          PRINT_AND_EXIT("serial_main_loop failed in main\n", GENERAL_ERROR);

    printf("outputing the tree...\n");
    // Report the growing tree
    write_newick(meta.gtree, master_ml_options->output_prefix, map.master_to_name);

    printf("done, cleaning up\n");
    // Clean up
    return 0;
}

int serial_main_loop(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst){
    int i, j; //loop counter

    // Init voting structure
    VOTE_GRP    vote;
    printf("current iteration is 2...");
    for(i = 3; i < meta->n_taxa; i++){
        print_inline_iteration(i, j, meta->n_taxa, 3);
        // printf("i is %d\n", i);

        if(init_vote(meta, map, mst, &vote, i)   != SUCCESS)             PRINT_AND_EXIT("init_vote failed in main loop\n", GENERAL_ERROR);

        // Determine the bipartition in the constraint tree
        if(find_bipartition(meta, map, mst, &vote, i) 
                                                != SUCCESS)             PRINT_AND_EXIT("find_bipartiiton failed in main loop\n", GENERAL_ERROR);

        // Determine the starting and ending of the valid subtree in the growing tree
        if(find_valid_subtree(meta, map, mst, &vote)
                                                != SUCCESS)             PRINT_AND_EXIT("find_valid_subtree failed in main loop\n", GENERAL_ERROR);   
                                                                                            // #if 1
                                                                                             //     printf("the vald subtree is %d %d %d %d\n", vote.st_lca.p, vote.st_lca.c, vote.nd_lca.p, vote.nd_lca.c); 
                                                                                             // if(i == 4 ) while(1);
                                                                                            // #endif
        // Vote!
                                                                                             // printf("i is %d\n", i);
        if(bfs_vote(meta, map, mst, &vote, i)  
                                                != SUCCESS)             PRINT_AND_EXIT("bfs_vote failed in main loop\n", GENERAL_ERROR);

            // printf("thje attachng edge is %d %d\n", vote.ins.p, vote.ins.c);

        // Attach to the growing subtree
        if(attach_leaf_to_edge(meta, map, mst, &vote, i)
                                                != SUCCESS)             PRINT_AND_EXIT("attach_leaf_to_edge failed in main loop\n", GENERAL_ERROR);
        
    }
    printf("\n");
    // printf("sip is %d\n", skip_counter);
    return 0;
}

// Driver for when the constraint inc is called without any constraint trees, using the default settings
#if COMPILE_ALONE
int main(int argc, char ** argv){
    int i; //loop counter
    ml_options options;

    options.init_d_name        = malloc(sizeof(char) * GENERAL_BUFFER_SIZE);
    options.init_d_name[0] = '\0';
    for(i = 0; i < argc; i++)
        if(strcmp(argv[i], "-i") == 0)
            strcpy(options.init_d_name, argv[i + 1]); 
        
    options.output_prefix      = malloc(sizeof(char) * GENERAL_BUFFER_SIZE);
    options.output_prefix[0] = '\0';
    for(i = 0; i < argc; i++)
        if(strcmp(argv[i], "-o") == 0) 
            strcpy( options.output_prefix, argv[i + 1]);    

    options.init_tree_name     = NULL;
    options.input_alignment    = NULL;
    options.guide_tree_name    = NULL;

    options.use_constraint                     = -1; // these shouldn't be used at all
    options.recompute_constraint_trees         = -1;    
    options.use_subtree_for_constraint_trees   = -1;
    options.use_raxml_for_constraint_trees     = -1;
    options.use_fasttree_for_constraint_trees  = -1;

    options.use_four_point_method_with_distance    = 1;
    options.use_four_point_method_with_tree        = 0;
    options.use_new_quartet_raxml                  = 0;
    options.use_ml_method                          = 0;

    options.distance_model                         = "logDet";
    options.ss_threshold                           = -1;

    if(constraint_inc_main(argc, argv, &options) != SUCCESS) PRINT_AND_EXIT("constraint_inc_main failed in main", GENERAL_ERROR);

    return 0;
} 
#endif


