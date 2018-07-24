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
    VOTE_GRP    vote;

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
    if(prim(adjacency_mat, num_sequence, &max_mst_weight, prim_ordering, adj_in_mst)
                                            != SUCCESS)         PRINT_AND_EXIT("prim algorithm failed in main\n", GENERAL_ERROR);

    // Initialize growing tree using the first 3 taxa in the ordering
    // TODO: check if there are less than 3 taxa 
    printf("initializing the growing tree\n");
    if(init_in_building(in_building, num_sequence) 
                                            != SUCCESS)         PRINT_AND_EXIT("init_in_building failed in main\n", GENERAL_ERROR); 
    if(init_growing_tree(growing_tree, prim_ordering, in_building, adj_in_mst[2])
                                            != SUCCESS)         PRINT_AND_EXIT("init_growing_tree failed in main\n", GENERAL_ERROR);

    printf("building the tree, current iteration is 2");
    // Loop through Prim's ordering

    if(serial_main_loop()                  != SUCCESS)          PRINT_AND_EXIT("serial_main_loop failed in main\n", GENERAL_ERROR);
    

    printf("\noutputing the tree...\n");
    // Report the growing tree
    write_newick(growing_tree, options->output_name);

    printf("done, cleaning up\n");
    // Clean up
    destroy_options(&options);
}

int serial_main_loop(int num_sequence, int * in_building, BT * growing_tree, int * prim_ordering){
    int i; //loop counter

    // Voting variable
    int * vote;

    // Other misc. variables
    int printer_counter;

    for(i = 3; i < num_sequence; i++){
        // Print current iteration to terminal, this is slow so turn off for large dataset
        for(printer_counter = 0; printer_counter < (int)(log10(i - 1)); printer_counter++)
            printf("\b");
        
        printf("%d", i);

        // Clear latent variables
        memset(in_building, 0, num_sequence * sizeof(int));
        memset(vote, 0, (num_sequence - 1) * sizeof(int));

        // Identify the constrained tree
        constraint_tree_idx = master_to_ctree_map[prim_ordering[i]];

        if(constraint_tree_idx != -1){
            // LCA for existing nodes
            lca_in_constraint = dfs_lca_change(master_to_ctreeindex[prim_ordering[i]], -1, constraint_trees[constraint_tree_idx], &rubbish, in_building, &lca_parent_in_constraint, 0);

            // DFS to identify bipartition (this should flip the corresponding leaves in the growing tree)
            bipartition_counter = 0;
            for(j = 0; j < 2; j++){
                if(constraint_trees[constraint_tree_idx]->adj_list[j] != lca_parent_in_constraint){
                    bipartition_counter++;
                    dfs_preorder(constraint_trees[constraint_tree_idx]->adj_list[j], lca_in_constraint, constraint_trees[constraint_tree_idx], in_building, bipartition_counter);
                }
            }

            // LCA in the growing tree, part 1. If the first taxon is in bipartition K, we will do bipartition 3 - K first
            first_lca_in_growing = dfs_lca_change(0, -1, growing_tree, in_building, &first_lca_parent_in_growing, in_building[0]);

            // LCA in the growing tree, part 2
            second_lca_in_growing = dfs_lca_change(first_lca_in_growing, first_lca_parent_in_growing, in_building, &second_lca_parent_in_growing, 3 - in_building[0]);
        } else {
            // The whole growing tree is valid
            first_lca_in_growing = 0;
            second_lca_in_growing = -1;
            first_lca_parent_in_growing = -1;
            second_lca_parent_in_growing = -1;
        }
        
        // BFS to find the legal subtree and update the vote. We will index into an edge array using its order in a bfs traversal starting from sequence 0
        bfs_vote(growing_tree, first_lca_in_growing, second_lca_in_growing, first_lca_parent_in_growing, vote, edge_child, edge_parent, num_sequence, adjacency_mat, prim_ordering[i], max_mst_weight);
        
        // Find the edge to add the new taxon into 
        find_addition_edge(vote, edge_child, edge_parent, &additional_edge_child, &addition_edge_parent, growing_tree->n - 1);

        // Attach to the growing subtree
        attach_leaf_to_edge(growing_tree, prim_ordering[i], addition_edge_parent, additional_edge_child, adj_in_mst[i]);
    }
    return 0;
}
