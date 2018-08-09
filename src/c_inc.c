// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tree.h"
#include "utilities.h"
#include "options.h"
#include "tools.h"
#include "dist.h"
#include "prim.h"
#include "traversal.h"

int serial_main_loop(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst);

int skip_counter = 0;

// Implementation of constrained version of the INC algorithm. 
// Command line argument is as followed: constraint_inc -i <alignment file> -t <tree1> <tree2> ... 
int main(int argc, char ** argv){
    // Meta variables
    option_t    options;

    INC_GRP     meta;
    MAP_GRP     map;
    MST_GRP     mst;

    // Parse options
    printf("reading in options...\n");

    if(init_options(&options)               != SUCCESS)         PRINT_AND_EXIT("init_options failed in main\n", GENERAL_ERROR);
    if(read_cmd_arg(argc, argv, &options)   != SUCCESS)         PRINT_AND_EXIT("read_cmd_arg failed in main\n", GENERAL_ERROR);

    // Getting the distance matrix
    printf("parsing initial matrix and tree...\n");
    if(parse_distance_matrix(&meta, &map, &options) 
                                            != SUCCESS)         PRINT_AND_EXIT("parse_distance_matrix failed in main\n", GENERAL_ERROR);

    // At this stage, all tree names should be in options->tree_name. It should be ok to parse them all at once since together they have at most 4M nodes
    if(parse_tree(&meta, &map, &options)    != SUCCESS)         PRINT_AND_EXIT("parse tree failed in main\n", GENERAL_ERROR);

    // Compute the MST 
    printf("computing the mst...\n");
    if(prim(&meta, &mst)                    != SUCCESS)         PRINT_AND_EXIT("prim algorithm failed in main\n", GENERAL_ERROR);
                                                                                                // printf("debug: the following print out prim's ctree\n");

                                                                                                // for(int i = 0; i < meta.n_taxa; i++)
                                                                                                //     printf("%d ", map.master_to_ctree[mst.prim_ord[i]]);
                                                                                                // printf("\n");
    // Initialize growing tree using the first 3 taxa in the ordering
    printf("initializing the growing tree...\n");
    if(init_growing_tree(&meta, &map, &mst) != SUCCESS)         PRINT_AND_EXIT("init_growing_tree failed in main\n", GENERAL_ERROR);

    printf("building the tree...\n");
    // Loop through Prim's ordering
    if(serial_main_loop(&meta, &map, &mst)  != SUCCESS)          PRINT_AND_EXIT("serial_main_loop failed in main\n", GENERAL_ERROR);

    printf("outputing the tree...\n");
    // Report the growing tree
    write_newick(meta.gtree, options.output_name, map.master_to_name);

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
                                                                                                 // printf("the vald subtree is %d %d %d %d\n", vote.st_lca.p, vote.st_lca.c, vote.nd_lca.p, vote.nd_lca.c); 
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


