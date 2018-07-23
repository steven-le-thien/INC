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
    adj_list ** mst;

    int max_mst_weight;
    int num_sequence;
    int i; //loop counter

    // Parse options
    init_options(&options);
    read_cmd_arg(argc, argv, &options);

    // At this stage, all tree names should be in options->tree_name. It should be ok to parse them all at once since together they have at most 4M nodes
    parse_trees(constraint_trees, &options);

    // Getting the distance matrix
    parse_distance_matrix(master_name_map, adjacency_mat, &options, &num_sequence);

    // Compute the MST 
    prim(adjacency_mat, num_sequence, &max_mst_weight, mst, prim_ordering);

    // Initialize growing tree using the first 3 taxa in the ordering
    // TODO: check if there are less than 3 taxa 
    init_growing_tree(growing_tree, prim_ordering);

    // Loop through Prim's ordering
    for(i = 3; i < num_sequence; i++){

    }
        // Base case

        // Identify the constrained tree

        // LCA for existing nodes

        // DFS to identify bipartition (this should flip the corresponding leaves in the growing tree)

        // LCA in the growing tree, part 1

        // LCA in the growing tree, part 2

        // DFS to find the legal subtree and update the vote

        // Check for ending condition

    // Report the growing tree

    // Clean up
    destroy_options(&options);
}
