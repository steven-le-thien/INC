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
    option_t    options;
    BT **       constraint_trees;
    BT *        growing_tree;
    char **     master_name_map;
    int *       prim_ordering;
    float **    adjacency_mat; // for Prim's algorithm
    int *       master_to_ctree_map;
    int *       master_to_ctreeindex;
    int *       in_building; // for O(1) look up: 0 if the node is not in the growing tree, 1 if the node is in the first bp, 2 if the node is in the second bipartition
    int *       vote; 
    int         lca_in_constraint;
    int         lca_parent_in_constraint;

    int         first_lca_in_growing;
    int         first_lca_parent_in_growing;
    int         second_lca_in_growing;
    int         second_lca_parent_in_growing;
    adj_list ** mst;

    int max_mst_weight;
    int num_sequence;
    int i; //loop counter
    int j;
    int bipartition_counter;
    int constraint_tree_idx;

    int placeholder; // rubbish variable

    // Parse options
    init_options(&options);
    read_cmd_arg(argc, argv, &options);

    // Getting the distance matrix
    parse_distance_matrix(master_name_map, adjacency_mat, &options, &num_sequence);

    // At this stage, all tree names should be in options->tree_name. It should be ok to parse them all at once since together they have at most 4M nodes
    parse_trees(constraint_trees, &options, master_to_ctree_map, master_to_ctreeindex, num_sequence, master_name_map);

    // Compute the MST 
    prim(adjacency_mat, num_sequence, &max_mst_weight, mst, prim_ordering);

    // Initialize growing tree using the first 3 taxa in the ordering
    // TODO: check if there are less than 3 taxa 
    init_in_building(in_building, num_sequence);
    init_growing_tree(growing_tree, prim_ordering, in_building);

    // Loop through Prim's ordering
    for(i = 3; i < num_sequence; i++){
        // Identify the constrained tree
        constraint_tree_idx = master_to_ctree_map[prim_ordering[i]];

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

        // DFS to find the legal subtree and update the vote

        // Attach to the subtree
    }

    // Report the growing tree
    write_newick(growing_tree, options->output_name);

    // Clean up
    destroy_options(&options);
}
