// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "test.h"
#include "tree.h"
#include "utilities.h"
#include "options.h"
#include "tools.h"
#include "dist.h"
#include "prim.h"
#include "traversal.h"

int serial_main_loop(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst);

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
    if(parse_tree(&meta, &map, &options)   != SUCCESS)         PRINT_AND_EXIT("parse tree failed in main\n", GENERAL_ERROR);

    // Compute the MST 
    printf("computing the mst...\n");
    if(prim(&meta, &mst)                    != SUCCESS)         PRINT_AND_EXIT("prim algorithm failed in main\n", GENERAL_ERROR);

    // Initialize growing tree using the first 3 taxa in the ordering
    printf("initializing the growing tree...\n");
    if(init_growing_tree(&meta, &mst)       != SUCCESS)         PRINT_AND_EXIT("init_growing_tree failed in main\n", GENERAL_ERROR);
    map.master_to_gidx[mst.prim_ord[0]] = 0;
    map.master_to_gidx[mst.prim_ord[1]] = 1;
    map.master_to_gidx[mst.prim_ord[2]] = 2;

    printf("building the tree...\n");
    // Loop through Prim's ordering
    if(serial_main_loop(&meta, &map, &mst)  != SUCCESS)          PRINT_AND_EXIT("serial_main_loop failed in main\n", GENERAL_ERROR);

    printf("\noutputing the tree...\n");
    // Report the growing tree
    write_newick(meta.gtree, "outo", map.master_to_name);

    printf("done, cleaning up\n");
    // Clean up
    return 0;
}

int serial_main_loop(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst){
    int i, j; //loop counter
    int printer_counter;
    int cap_counter;

    // Init voting structure
    VOTE_GRP    vote;
    VOTE_GRP    vote1;
    VOTE_GRP    vote2;
    VOTE_GRP    vote3;
    if(init_vote(meta, &vote)                != SUCCESS)         PRINT_AND_EXIT("init_vote failed in main loop\n", GENERAL_ERROR);
    for(i = 3; i < meta->n_taxa; i++){
                if(i % 1 == 0) printf("current iteration is %d...\n", i);

                                                                                            #if DEBUG 
                                                                                               
                                                                                                printf("debug: high level, the following print out the number of nodes with degree 1 and 3 and their sum compared to total nodes\n");
                                                                                                int j,k1, k3, k4, kt;
                                                                                                kt = 0;
                                                                                                // for(i = 0; i < meta->n_ctree; i++){
                                                                                                    k1 = 0;
                                                                                                    k3 = 0;
                                                                                                    k4 = 0;
                                                                                                    for(j = 0; j < meta->gtree->n_node; j++){
                                                                                                        k1 += meta->gtree->degree[j] == 1;
                                                                                                        k3 += meta->gtree->degree[j] == 3;
                                                                                                        k4 += (meta->gtree->degree[j] != 1) && ( meta->gtree->degree[j] != 3);

                                                                                                    }
                                                                                                    printf("for growing tree node with deg 1 is %d, 3 is %d, total is %d, n_node is %d, non13 is %d\n", k1, k3, k1 + k3, meta->gtree->n_node, k4);
                                                                                                    kt += k1;
                                                                                                // }
                                                                                                printf("debug: total number of leaf is %d\n", kt);

                                                                                                printf("debug: the following tests actual adjacency\n");
                                                                                                int flag = 0; //this sholuld be 0 all the way until the end
                                                                                                // for(i = 0; i < meta->n_ctree; i++){
                                                                                                    BT * tree = meta->gtree;

                                                                                                    for(int k = 0; k < tree->n_node; k++){
                                                                                                        for(j = 0; j < tree->degree[k]; j++){
                                                                                                            int iflag = 1;  
                                                                                                            int ad = tree->adj_list[k][j].dest;
                                                                                                            for(k1 = 0; k1 < tree->degree[ad]; k1++){
                                                                                                                if(tree->adj_list[ad][k1].dest == k){
                                                                                                                    if(iflag) iflag = 0;
                                                                                                                    else
                                                                                                                        flag++;
                                                                                                                }
                                                                                                            }
                                                                                                            if(iflag) flag++;
                                                                                                        }
                                                                                                    }
                                                                                                // }
                                                                                                if(!flag) printf("passed test\n");
                                                                                                else {
                                                                                                    printf("failed test with flag %d\n", flag);
                                                                                                    while(1);
                                                                                                }

                                                                                                printf("debug: this test tests mapping from ctree to master index\n");
                                                                                                flag = 0;
                                                                                                // for(i = 0; i< meta->n_ctree; i++){
                                                                                                    for(j = 0; j < tree->n_node; j++){
                                                                                                        if((tree->master_idx_map[j] < 0 && tree->degree[j] == 1) || 
                                                                                                            (tree->master_idx_map[j] >= 0 && tree->degree[j] != 1))
                                                                                                            flag++; 
                                                                                                    }

                                                                                                    for(j = 0; j < tree->n_node; j++)
                                                                                                        for(k1 = 0; k1 < tree->n_node; k1++)
                                                                                                            if(tree->master_idx_map[j] != -1 && j != k1 && tree->master_idx_map[j] == tree->master_idx_map[k1]){
                                                                                                                printf("failure at %d %d %d %d\n", j, k1, tree->master_idx_map[j], tree->master_idx_map[k1]);
                                                                                                                flag++;
                                                                                                            }
                                                                                                // }


                                                                                                if(!flag) printf("passed test\n");
                                                                                                else printf("failed test\n");
                                                                                            #endif

        cap_counter = 0;

        // Clear voter array
        memset(vote.vote, 0, 4 * vote.n_taxa * sizeof(int));
        memset(vote.edge_c, -1, 4 * vote.n_taxa * sizeof(int));
        memset(vote.edge_p, -1, 4 * vote.n_taxa * sizeof(int));
        memset(&vote.valid_st, -1, 10 * sizeof(int));

        // Clear visited array from bipartition
        for(j = 0; j < meta->n_taxa; j++)
            if(meta->visited[j] > 0) 
                meta->visited[j] = -1;

        // Identify the constraint tree
        vote.ctree_idx = map->master_to_ctree[mst->prim_ord[i]];
                                                                                            #if DEBUG 
                                                                                                printf("debug: node of interest is %d constraint tree is %d, name is %s\n", mst->prim_ord[i], vote.ctree_idx, map->master_to_name[mst->prim_ord[i]]); 
                                                                                            #endif
        if(vote.ctree_idx != -1 && meta->ctree[vote.ctree_idx]->n_node != 1){ // if the constraint tree is non-trivial
                                                                                            #if DEBUG 
                                                                                                printf("debug: constraint tree is non-trivial\n"); 
                                                                                            #endif
            // Determine the bipartition in the constraint tree
            if(find_bipartition(meta, map, mst, &vote, i) 
                                            != SUCCESS)         PRINT_AND_EXIT("find_bipartiiton failed in main loop\n", GENERAL_ERROR);

            // Here we may have the case where the intersection of constraint tree and growing tree taxa is less than 3
            for(j = 0; j < meta->n_taxa; j++)
                cap_counter += meta->visited[j] > 0;

            if(cap_counter < 3){
                vote.st_lca.p = meta->gtree->adj_list[0][0].dest; 
                vote.st_lca.c = 0;
            } else {
                                                                                            #if DEBUG 
                                                                                                 printf("DEBUG: i'm findng valid subtree\n"); 
                                                                                            #endif
                // Determine the starting and ending of the valid subtree in the growing tree
                if(find_valid_subtree(meta, map, mst, &vote)
                                                != SUCCESS)         PRINT_AND_EXIT("find_valid_subtree failed in main loop\n", GENERAL_ERROR);   
            }
        } else { // set to start at some arbitray leaf 
            vote.st_lca.p = meta->gtree->adj_list[0][0].dest; 
            vote.st_lca.c = 0;
        }

        // If there are only one edge then we don't need to do any voting'
        if(vote.st_lca.p == vote.nd_lca.c && vote.st_lca.c == vote.nd_lca.p){
            vote.ins.c = vote.st_lca.c;
            vote.ins.p = vote.st_lca.p;
        } else {
                                                                                            #if DEBUG 
                                                                                                 printf("DEBUG: i'm doing the votin\n"); 
                                                                                            #endif
            if(bfs_vote(meta, map, mst, &vote, i)  
                                                != SUCCESS)         PRINT_AND_EXIT("bfs_vote failed in main loop\n", GENERAL_ERROR);

            // Find the edge to add the new taxon into 
            // if(find_insertion_edge(meta, &vote) 
            //                                     != SUCCESS)         PRINT_AND_EXIT("find_insertion_edge failed in main loop\n", GENERAL_ERROR);
        }
                                                                                            #if DEBUG 
                                                                                                 printf("i is %d\n", i); 
                                                                                            #endif
        // BFS to find the legal subtree and update the vote. We will index into an edge array using its order in a bfs traversal starting from sequence 0
                                                                                            #if DEBUG 
                                                                                                printf("debug: the edge i'm attaching to is (%d %d) with i %d\n", vote.ins.c, vote.ins.p, i);
                                                                                                // if(i == 5) while(1);
                                                                                            #endif
        // Attach to the growing subtree
        if(attach_leaf_to_edge(meta, map, mst, &vote, i)
                                            != SUCCESS)         PRINT_AND_EXIT("attach_leaf_to_edge failed in main loop\n", GENERAL_ERROR);
        map->master_to_gidx[mst->prim_ord[i]] = meta->gtree->n_node - 1;


    }
    return 0;
}
