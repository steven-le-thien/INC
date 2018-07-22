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


// Implementation of constrained version of the INC algorithm. 
// Command line argument is as followed: constrained_inc -i <alignment file> -t <tree1> <tree2> ... 
int main(int argc, char ** argv){
    // Parse options
    option_t options;
    BT ** constrained_trees;

    init_options(&options);

    read_cmd_arg(argc, argv, &options);

    // At this stage, all tree names should be in options->tree_name. It should be ok to parse them all at once since together they have at most 4M nodes
    parse_trees(constrained_trees, &options);

    // Getting the distance matrix

    // Compute the MST 

    // Initialize the growing tree

    // Loop through Prim's ordering
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
