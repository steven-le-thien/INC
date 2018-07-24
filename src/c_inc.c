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


// Implementation of constrained version of the INC algorithm. 
// Command line argument is as followed: constrained_inc -i <alignment file> -t <tree1> <tree2> ... 
int main(int argc, char ** argv){
    // Meta variables
    option_t    options;
    BT **       constraint_trees;
    BT *        growing_tree;
    float **    adjacency_mat; // distance matrix
    int num_sequence;


    // Mapping variables
    int *       master_to_ctree_map;
    int *       master_to_ctreeindex;
    char **     master_name_map;

    // Prim mst variable
    int *       prim_ordering;
    int *       adj_in_mst;
    int max_mst_weight;

    // Tree growing variables
    int *       in_building; // for O(1) look up: 0 if the node is not in the growing tree, 
                             // 1 if the node is in the first bp, 2 if the node is in the second bipartition
    
    // Voting variables
    int *       vote; 
    int *       vote_first_endpoint;
    int *       vote_second_endpoint;

    int         lca_in_constraint;
    int         lca_parent_in_constraint;

    int         first_lca_in_growing;
    int         first_lca_parent_in_growing;

    int         second_lca_in_growing;
    int         second_lca_parent_in_growing;

    // Taxon insertion variables
    int         addition_child;
    int         addition_parent;

    // Loop counter variables
    int         i; //loop counter
    int         j;
    int         bipartition_counter;
    int         constraint_tree_idx;

    // Other variables
    int         placeholder; 
    int         printer_counter;

    // Parse options
    printf("reading in options...\n");
    init_options(&options);
    read_cmd_arg(argc, argv, &options);

    // Getting the distance matrix
    printf("parsing initial matrix and tree...\n");
    parse_distance_matrix(master_name_map, adjacency_mat, &options, &num_sequence);

    // At this stage, all tree names should be in options->tree_name. It should be ok to parse them all at once since together they have at most 4M nodes
    parse_trees(constraint_trees, &options, master_to_ctree_map, master_to_ctreeindex, num_sequence, master_name_map);

    // Compute the MST 
    printf("computing the mst...\n");
    prim(adjacency_mat, num_sequence, &max_mst_weight, prim_ordering, adj_in_mst);

    // Initialize growing tree using the first 3 taxa in the ordering
    // TODO: check if there are less than 3 taxa 
    printf("initializing the growing tree\n");
    init_in_building(in_building, num_sequence);
    init_growing_tree(growing_tree, prim_ordering, in_building, adj_in_mst[2]);

    printf("building the tree, current iteration is 2");
    // Loop through Prim's ordering
    for(i = 3; i < num_sequence; i++){
        // Print current iteration to terminal, this is slow so turn off for large dataset
        for(printer_counter = 0; printer_counter < (int)(log10(i - 1)); printer_counter++){
            printf("\b");
        }
        printf("%d", i);
        // Clear latent
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

    printf("\noutputing the tree...\n");
    // Report the growing tree
    write_newick(growing_tree, options->output_name);

    printf("done, cleaning up\n");
    // Clean up
    destroy_options(&options);
}
