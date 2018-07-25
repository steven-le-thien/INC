// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "msa.h"
#include "tree.h"
#include "utilities.h"
#include "options.h"
#include "tools.h"
#include "stat.h"
#include "dist.h"
#include "prim.h"
#include "c_inc.h"


// Implementation of constrained version of the INC algorithm. 
// Command line argument is as followed: constrained_inc -i <alignment file> -t <tree1> <tree2> ... 
int main(int argc, char ** argv){
    // Meta variables
    option_t    options;

    INC_GRP     meta;
    MAP_GRP     map;
    MST_GRP     mst;

    // Loop counter variables
    int         i; //loop counter
    int         j;
    int         bipartition_counter;
    int         constraint_tree_idx;

    // Other variables
    int         placeholder; 

    // Parse options
    printf("reading in options...\n");

    if(init_options(&options)               != SUCCESS)         PRINT_AND_EXIT("init_options failed in main\n", GENERAL_ERROR);
    if(read_cmd_arg(argc, argv, &options)   != SUCCESS)         PRINT_AND_EXIT("read_cmd_arg failed in main\n", GENERAL_ERROR);

    // Getting the distance matrix
    printf("parsing initial matrix and tree...\n");
    if(parse_distance_matrix(&meta, &map, &options) 
                                            != SUCCESS)         PRINT_AND_EXIT("parse_distance_matrix failed in main\n", GENERAL_ERROR);

    // At this stage, all tree names should be in options->tree_name. It should be ok to parse them all at once since together they have at most 4M nodes
    if(parse_trees(&meta, &map, &options)   != SUCCESS)         PRINT_AND_EXIT("parse tree failed in main\n", GENERAL_ERROR);

    // Compute the MST 
    printf("computing the mst...\n");
    if(prim(&meta, &mst)                    != SUCCESS)         PRINT_AND_EXIT("prim algorithm failed in main\n", GENERAL_ERROR);

    // Initialize growing tree using the first 3 taxa in the ordering
    printf("initializing the growing tree\n");
    if(init_growing_tree(&meta, &mst)       != SUCCESS)         PRINT_AND_EXIT("init_growing_tree failed in main\n", GENERAL_ERROR);

    printf("building the tree, current iteration is 2");
    // Loop through Prim's ordering
    if(serial_main_loop(&meta, &map, &mst)  != SUCCESS)          PRINT_AND_EXIT("serial_main_loop failed in main\n", GENERAL_ERROR);

    printf("\noutputing the tree...\n");
    // Report the growing tree
    write_newick(&meta, &options);

    printf("done, cleaning up\n");
    // Clean up
    destroy_options(&options);
}

int serial_main_loop(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst){
    int i; //loop counter
    int printer_counter;

    // Init voting structure
    VOTE_GRP    vote;
    if(init_vote(meta, vote)                != SUCCESS)         PRINT_AND_EXIT("init_vote failed in main loop\n", GENERAL_ERROR);

    for(i = 3; i < meta->n_taxa; i++){
        // Print current iteration to terminal, this is slow so turn off for large dataset
        for(printer_counter = 0; printer_counter < (int)(log10(i - 1)); printer_counter++)
            printf("\b");
        printf("%d", i);

        // Clear voter array
        if(clear_vote(&vote)                != SUCCESS)         PRINT_AND_EXIT("clear_vote failed in main loop\n", GENERAL_ERROR);

        // Identify the constrained tree
        vote.ctree_idx = map->master_to_ctree[mst->prim_ord[i]];

        if(vote.ctree_idx != -1 && meta->ctree[vote->ctree_idx]->n_node == 1){ // if the constraint tree is non-trivial
            // Determine the bipartition in the constraint tree
            if(find_bipartition(meta, map, mst, &vote, i) 
                                            != SUCCESS)         PRINT_AND_EXIT("find_bipartiiton failed in main loop\n", GENERAL_ERROR);

            // Determine the starting and ending of the valid subtree in the growing tree
            if(find_valid_subtree(meta, map, mst, &vote)
                                            != SUCCESS)         PRINT_AND_EXIT("find_valid_subtree failed in main loop\n", GENERAL_ERROR);
            
        } else // set to start at some arbitray leaf 
            vote.st_lca.c = 0; 

        // BFS to find the legal subtree and update the vote. We will index into an edge array using its order in a bfs traversal starting from sequence 0
        if(bfs_vote(meta, mst, &vote, i)  
                                            != SUCCESS)         PRINT_AND_EXIT("bfs_vote failed in main loop\n", GENERAL_ERROR);
        
        // Find the edge to add the new taxon into 
        if(find_insertion_edge(meta, &vote) 
                                            != SUCCESS)         PRINT_AND_EXIT("find_insertion_edge failed in main loop\n", GENERAL_ERROR);

        // Attach to the growing subtree
        if(attach_leaf_to_edge(meta, mst, vote, i)
                                            != SUCCESS)         PRINT_AND_EXIT("attach_leaf_to_edge failed in main loop\n", GENERAL_ERROR);
    }
    return 0;
}
