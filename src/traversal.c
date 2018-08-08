// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "traversal.h"
#include "utilities.h"
#include "options.h"
#include "quartet.h"

extern int skip_counter;

// All trees should be BTs 
int dfs_lca_implementation(int node, int parent, BT * tree, int * dp, int * in_building, int * lca_parent, int * lca_child, int mode);
int dfs_preorder(int node, int parent, BT * tree, int * in_building, int flag);
int bfs_vote_implementation(BT * tree, int valid_start, int valid_end, int valid_start_parent, int * mapping, int * edge_child, int * edge_parent, int n, float ** d, int x, float q0, int revote_power, int ** revote_map, int all_quartets, int revoting, char ** name_map, char * master_tree);
int find_insertion_edge_implementation(int * vote, int * edge_child, int *edge_parent, int * additional_edge_child, int * addition_edge_parent, int num_edges);
void dfs_backtrack(int node, int parent, int findme, BT * tree, int * next_to_start, int * found);
int less_than_3_shared_taxa(INC_GRP * meta);
void all_valid(INC_GRP * meta, VOTE_GRP * vote);

/* Function to initialize variables relevant to voting to a blank state so that a new round 
 *      of voting can be done 
 * Input:   meta,   including the array of visited taxa in the growing tree, the growing tree as well as the constraint tree
 *          map     including bookkeeping arrays for fast lookups
 *          mst     including the prim ordering to get the next taxon
 *          i       indexing into the prim ordering
 * Output:  0 on success, ERROR otherwise
 * Effet:   run memset to clear fields, reset meta's visited array, identify the correct constraint tree
 */
int init_vote(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote, int i){
    int j;
    // Set all other special edges to 0 at one go. WARNING: redo this when the struct changes   
    memset(&(vote->valid_st), -1, 10 * sizeof(int));

    // Clear visited array from bipartition (but still keep the visted flags)
    for(j = 0; j < meta->n_taxa; j++)
        if(meta->visited[j] > 0) 
            meta->visited[j] = -1;

    // Identify the constraint tree
    vote->ctree_idx = map->master_to_ctree[mst->prim_ord[i]];
                                                                                        #if DEBUG 
                                                                                            printf("debug: node of interest is %d constraint tree is %d, name is %s\n", mst->prim_ord[i], vote->ctree_idx, map->master_to_name[mst->prim_ord[i]]); 
                                                                                        #endif

                                                                                        #if DEBUG 
                                                                                            printf("debug: high level, the following print out the number of nodes with degree 1 and 3 and their sum compared to total nodes\n");
                                                                                            int k1, k3, k4, kt;
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
                                                                                            printf("debug: done init vote %d\n", vote->ctree_idx);
                                                                                        #endif    


    return 0;
}

/* Assuming that the correct constraint tree is found and meta's visited array has the correct 
 *      record for leaves in the growing tree, this function 
 *      1. identify the intersection A between the leaf set of the growing tree and the correct constraint tree
 *      2. determine the bipartition of the intersection set within the constraint tree by removing the
 *              component containing a certain taxon in t_c^A
 * Exact implementation detail is as followed:
 *      1. Given L, the set of leaves in the growing tree, root the tree at the query taxon and
 *          run a tree dp (/a post-order dfs) recursion that, at each internal node n, count the 
 *          number of leaves in the subtree rooted at n contained in L, this is the dp value
 *      2. As we move up the tree, we also keep track of the last time dp changes (by storing the parent and child
 *          of the edge on which the value change). This gives the LCA of taxa in L within the constraint tree 
 *      3. The LCA, when removed, divides the constraint tree into 3 subtrees, 1 containing the query taxon and 
 *          no leaves in L; another containing a subset of L_1 leaves in L and the last contain the rest of L.
 *          Note that we want to mark L_1 and L/L_1 
 *      4. 2 simple preorder dfs, each starting at 2 children of the LCA found previously, suffice to identify and 
 *          mark out this bipartition of L
 * Input:   all meta groups
 * Output:  0 on success, ERROR otherwise
 * Effect:  set some fields in vote
 */ 
int find_bipartition(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote, int i){
    int placeholder = 0;    // extra variable for dp
    int bp_counter = 0;     // bipartition counter
    int cur_node;           // loop counters 
    int j; 

    
    if(vote->ctree_idx < 0 || meta->ctree[vote->ctree_idx]->n_node < 3) return 0; // trivial tree, don't have any bipartition (from lca)
                                                                                                #if DEBUG 
                                                                                                printf("debug: after bipartition\n"); 
                                                                                                printf("debug: value of visited that is non zeor\n");
                                                                                                for(j = 0; j < meta->n_taxa; j++)
                                                                                                    if(meta->visited[j])
                                                                                                        printf("(%d) %d in %d (%s), ", meta->visited[j], j, map->master_to_ctree[j], map->master_to_name[j]);
                                                                                                printf("\n");

                                                                                                printf("debug: this shoi;d ne all -1\n");
                                                                                                for(j = 0; j < meta->ctree[vote->ctree_idx]->n_node; j++)
                                                                                                    if(meta->ctree[vote->ctree_idx]->degree[j] == 1)
                                                                                                        printf("%d(%s) ", meta->visited[meta->ctree[vote->ctree_idx]->master_idx_map[j]], map->master_to_name[meta->ctree[vote->ctree_idx]->master_idx_map[j]]);
                                                                                                printf("\n");
                                                                                            #endif

                                                                                            #if DEBUG 
                                                                                                printf("debug: input to dfs_lca_implementation is: starting node %d with name %s, master index %d, tree %d\n", map->master_to_cidx[mst->prim_ord[i]], map->master_to_name[meta->ctree[vote->ctree_idx]->master_idx_map[map->master_to_cidx[mst->prim_ord[i]]]], meta->ctree[vote->ctree_idx]->master_idx_map[map->master_to_cidx[mst->prim_ord[i]]], vote->ctree_idx); 
                                                                                            #endif
    if(dfs_lca_implementation(map->master_to_cidx[mst->prim_ord[i]],    /* strarting node   root the tree at the query taxon */
                                -1,                                     /* avoiding node    no parent since it's a leaf */
                                meta->ctree[vote->ctree_idx],           /* tree             get the correct constraint tree */
                                &placeholder,                           /* dp value         we won't use this (except for checking whether we are getting back the correct value) */
                                meta->visited,                          /* array mapping bipartion */
                                &(vote->c_lca.p),                       /* the address to write the parent and child in the resulting process to*/
                                &(vote->c_lca.c), 
                                0) != SUCCESS)              PRINT_AND_RETURN("dfs_lca_implementation failed in wrapper\n", GENERAL_ERROR);
                                                                                            #if DEBUG 
                                                                                                    printf("placehiolder is %d\n", placeholder);

                                                                                                printf("debug: output of dfs_loca is %d %d (%d %d) in master index with dp%d\n", vote->c_lca.p, vote->c_lca.c, meta->ctree[vote->ctree_idx]->master_idx_map[vote->c_lca.p], meta->ctree[vote->ctree_idx]->master_idx_map[vote->c_lca.c], placeholder); 
                                                                                            #endif
    // If there are less than 3 shared common taxa then don't bother
    if(placeholder < 3) return 0;

    // Identify the 2 children of the LCA found above and mark the correct bipartition
    for(j = 0; j < meta->ctree[vote->ctree_idx]->degree[vote->c_lca.c]; j++){
        cur_node = meta->ctree[vote->ctree_idx]->adj_list[vote->c_lca.c][j].dest;
        if(cur_node != vote->c_lca.p)
            if(dfs_preorder(cur_node,                                   /* starting node: child of the lca */
                            vote->c_lca.c,                              /* starting parent: the lca (this is important so that there's no backedge in the dfs) */
                            meta->ctree[vote->ctree_idx],               /* the constraint tree */
                            meta->visited,                              /* the array to record the bipartition */
                            ++bp_counter) != SUCCESS)       PRINT_AND_RETURN("dfs_preorder failed in wrapper\n", GENERAL_ERROR);

    }

    return 0;
}

/* Assuming the correct bipartition of A (LeafSet(t_c) \cap LeafSet(t)) where t_c is the constraint tree and t is the growing tree
 *      has been identified, this function find the valid subtree (the 'component in t_c') in which the query taxon may be added. 
 * The exact implementation detail is as followed:
 *      1. Let the biparition be {A1, A2}. There guarantees to be a node in A1 (i.e, A1 is nonempty). Find that node, root the tree there and find the LCA of A2, call this l2
 *      2. Likewise, find some node in A2, root the tree there and find the LCA of A1, call this l1
 *      3. Run any graph traversal (here, a dfs equiped with backtracking) from l1 and stops when reaching l2. This attempts to find the correct orientation of l1 and l2 in the valid subtree
 *          3.1 Here we denote the parent of l1 is the neighbor of l1 that leads to l2. This is unique since there's a unqiue path from any 2 nodes in a tree 
 * Input: all meta variables
 * Output: 0 on success, ERROR otherwise
 * Effect: set fields in vote
 */
int find_valid_subtree(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote){
    int i;
    int placeholder;
    
    if(less_than_3_shared_taxa(meta)){ // if less than 3 shared taxa then everything is valid
        all_valid(meta, vote);
        return 0;
    }

    for(i = 0; i < meta->n_taxa; i++)
        if(meta->visited[i] == 2)
            break; // there guaranteed to be 1
                                                                                            #if DEBUG 
                                                                                                printf("debug: input to find_valid_subtree st dfs is: starting node %d with flag %d\n", map->master_to_gidx[i], meta->visited[mst->prim_ord[0]] ? 3 - meta->visited[mst->prim_ord[0]] : 1); 
                                                                                            #endif

    // Run from the first leaf. We root the tree at the first taxon. If the taxon is in the constraint tree, we find the other bipartition; else we find the 1st biparition
    if(dfs_lca_implementation(map->master_to_gidx[i],       /* i is the taxon that is in the secon partition. here we obtain its index in the growing tree */
                        -1,                                 /* it is a leaf in the growing tree so we don't need to avoid anything */
                        meta->gtree,                        /* the growing tree */
                        &placeholder,                       /* some field used in the dfs (to hold the dp value) */
                        meta->visited,                      /* the array containing the biparition */
                        &(vote->st_lca.p),                  /* the parent isn't really important but the general implementation needs it */
                        &(vote->st_lca.c),                  /* the lca of leaves in the first partition */
                        1)  
                                != SUCCESS)     PRINT_AND_RETURN("dfs_lca_implementation failed in find_valid_subtree\n", GENERAL_ERROR);
                                                                                            #if DEBUG 
                                                                                                printf("debug: output of st dfs_lca is %d %d (%d %d) (p, c) in master index \n", vote->st_lca.p, vote->st_lca.c, meta->gtree->master_idx_map[vote->st_lca.p], meta->gtree->master_idx_map[vote->st_lca.c]); 
                                                                                            #endif
    for(i = 0; i < meta->n_taxa; i++)
        if(meta->visited[i] == 1)
            break; // there guaranteed to be 1
    // Find lca of the other bipartition that was not found in the first phase, by starting from the first phase lca, avoiding its own parent ()
    if(dfs_lca_implementation(map->master_to_gidx[i],       /* th index of some leaf in the first partition in the growing tree */
                        -1,                                 /* it is a leaf so we dont need to avoid anything */
                        meta->gtree,                        /* the growing tree */
                        &placeholder,                       /* some field used in the dfs (to hold the dp value) */
                        meta->visited,                      /* the array containing the biparition */
                        &(vote->nd_lca.p),                  /* the parent isn't really important but the general implementation needs it */
                        &(vote->nd_lca.c),                  /* the lca of leaves in the second partition */
                        2) 
                                != SUCCESS)     PRINT_AND_RETURN("dfs_lca_implementation failed in find_valid_subtree\n", GENERAL_ERROR);
                                                                                            #if DEBUG 
                                                                                                printf("debug: output of nd dfs_lca is %d %d (%d %d) (p, c) in master index \n", vote->nd_lca.p, vote->nd_lca.c, meta->gtree->master_idx_map[vote->nd_lca.p], meta->gtree->master_idx_map[vote->nd_lca.c]); 
                                                                                            #endif
    // Find correct orientation
    dfs_backtrack(vote->st_lca.c,
                    -1,                                     /* find in all direction from the starting point */
                    vote->nd_lca.c,                         /* search ends when we find the second lca */
                    meta->gtree,                            /* the growing tree */
                    &(vote->st_lca.p),                      /* storing the correct orientation */
                    &placeholder);                          /* not used here */
                                                                                            #if DEBUG 
                                                                                                printf("debug: output of backtrack is %d %d (%d %d); %d %d (%d %d) (p, c) in master index \n", vote->st_lca.p, vote->st_lca.c, meta->gtree->master_idx_map[vote->st_lca.p], meta->gtree->master_idx_map[vote->st_lca.c], vote->nd_lca.p, vote->nd_lca.c, meta->gtree->master_idx_map[vote->nd_lca.p], meta->gtree->master_idx_map[vote->nd_lca.c]); 
                                                                                            #endif
    return 0;
}

/* Assuming that the correct valid subtree is identified and the orientation of its endpoint is correct, this runs a bfs from one of the endpoints 
 *      to update vote counting relative to the starting edge. This is described in Theorem 7 of the original paper. Implementation detail is as followed
 *      1. Recall that we determine the orientation of the 2 endpoints l1, l2 by the neighbor of l1 that leads to l2. Call this vertex lp, the algo in Theorem 7 starts with (l1, lp)
 *      2. At each vertex visited in the bfs order, we are only interested in the internal nodes (since we want to update values on the edge). Thus seeing a leaf or seeing l2 ends the run.
 *      3. At each internal node, we identify the `parent' in the bfs order. Theorem 7 gives a way to update the 2 children edge based on the parents' value; the validity of the voting quartet
 *          and results from 4 point method. 
 *      4. Node with higher vote than its parent is immediate checked whether it is the node with the most vote encountered so far.
 * Input: all meta variables
 * Output: 0 on success, ERROR otherwise
 * Effect: set fields in vote
 */
int bfs_vote(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote, int i){
    int * revote_map;
    if(vote->st_lca.p == vote->nd_lca.c && vote->st_lca.c == vote->nd_lca.p){ // only 1 edge is valid
        vote->ins.c = vote->st_lca.c;
        vote->ins.p = vote->st_lca.p;
                                                                                            #if DEBUG 
                                                                                                printf("debug: the edge i'm attaching to is (%d %d) with i %d\n", vote->ins.c, vote->ins.p, i);
                                                                                            #endif
                                                                                                skip_counter++;
        return 0;
    }

    if(bfs_vote_implementation(meta->gtree,             // growing tree
                            vote->st_lca.p,             // starting nodes and parents, notice the reverse of ordering since the subtree is within the growing tree
                            vote->nd_lca.c,             // ending node
                            vote->st_lca.c,
                            meta->gtree->master_idx_map,// mapping to get the correct indexing from the distance matrix
                            &(vote->ins.c),             // store the best edge 
                            &(vote->ins.p),
                            meta->gtree->n_node,                
                            meta->d,                    // distance matrix 
                            mst->prim_ord[i],           // query taxon
                            mst->max_w,                 // q0
                            0,                          // power, 0 for 0-1, 1 or 2 for weighted voting
                            &revote_map,                // map
                            0,                          // are all quartets voting ?
                            0,                          // is revoting ? 
                            map->master_to_name,
                            meta->mtree)
                                != SUCCESS)     PRINT_AND_RETURN("bfs_vote_implementation failed in bfs_voite\n", GENERAL_ERROR);

    // if(bfs_vote_implementation(meta->gtree,             // growing tree
    //                         vote->st_lca.p,             // starting nodes and parents, notice the reverse of ordering since the subtree is within the growing tree
    //                         vote->nd_lca.c,             // ending node
    //                         vote->st_lca.c,
    //                         meta->gtree->master_idx_map,// mapping to get the correct indexing from the distance matrix
    //                         &(vote->ins.c),             // store the best edge 
    //                         &(vote->ins.p),
    //                         meta->gtree->n_node,                
    //                         meta->d,                    // distance matrix 
    //                         mst->prim_ord[i],           // query taxon
    //                         mst->max_w,                 // q0
    //                         3,                          // power, 0 for 0-1, 1 or 2 for weighted voting
    //                         &revote_map,                // map
    //                         1,                          // are all quartets voting ?
    //                         1)                          // is revoting ?
    //                             != SUCCESS)     PRINT_AND_RETURN("bfs_vote_implementation failed in bfs_voite\n", GENERAL_ERROR);

    // if(bfs_vote_implementation(meta->gtree,             // growing tree
    //                         vote->st_lca.p,             // starting nodes and parents, notice the reverse of ordering since the subtree is within the growing tree
    //                         vote->nd_lca.c,             // ending node
    //                         vote->st_lca.c,
    //                         meta->gtree->master_idx_map,// mapping to get the correct indexing from the distance matrix
    //                         &(vote->ins.c),             // store the best edge 
    //                         &(vote->ins.p),
    //                         meta->gtree->n_node,                
    //                         meta->d,                    // distance matrix 
    //                         mst->prim_ord[i],           // query taxon
    //                         mst->max_w,                 // q0
    //                         2,                          // power, 0 for 0-1, 1 or 2 for weighted voting
    //                         &revote_map,                // map
    //                         1,                          // are all quartets voting ?
    //                         1)                          // is revoting ?
    //                             != SUCCESS)     PRINT_AND_RETURN("bfs_vote_implementation failed in bfs_voite\n", GENERAL_ERROR);

                                                                                            #if DEBUG
                                                                                                printf("debug: the edge i'm attaching to is (%d %d) with i %d\n", vote->ins.c, vote->ins.p, i);
                                                                                            #endif
    free(revote_map);
    return 0;
}



// All trees should be BTs 
int dfs_lca_implementation(int node, int parent, BT * tree, int * dp, int * in_building, int * lca_parent, int * lca_child, int mode){
    int i;
    int child_dp; 
    int child_lca_parent;
    int child_lca_child;
    int child_dfs_return;
    int flag;
    int ret;
                                                                                            #if DEBUG && DEBUG_REC
                                                                                                printf("debug: in recursion, node is %d(%d), parent is %d\n", node, tree->master_idx_map[node], parent); 
                                                                                            #endif

    if(tree->degree[node] == 1 && parent != -1){
        if((!mode && in_building[tree->master_idx_map[node]]) ||
            (mode && in_building[tree->master_idx_map[node]] == mode)){
                                                                                            #if DEBUG && DEBUG_REC 
                                                                                                printf("debug: base case 1 (pos) reached, node is %d(%d), parent is %d\n", node, tree->master_idx_map[node], parent); 
                                                                                            #endif
            // printf("heee%d\n", node);
            *dp = 1;
            *lca_parent = parent;
            *lca_child = node;
            return 0;
        } else {
                                                                                            #if DEBUG && DEBUG_REC 
                                                                                                printf("debug: base case 2 (neg) reached, node is %d(%d), parent is %d\n", node, tree->master_idx_map[node], parent); 
                                                                                            #endif
            *dp = 0;
            *lca_child = -1;
            *lca_parent = -1;
            return 0;
        }
    }

    // Initialization
    *dp = 0;
    ret = -1;
    flag = 0;
    child_dp = -1;
    child_lca_parent    = -1;
    child_lca_child     = -1;
    child_dfs_return    = -1;

    for(i = 0; i < tree->degree[node]; i++){
        if(tree->adj_list[node][i].dest == parent) continue;
        child_dfs_return = dfs_lca_implementation(tree->adj_list[node][i].dest, node, tree, &child_dp, in_building, &child_lca_parent, &child_lca_child, mode);
                                                                                            
        if(child_dp){
            *dp += child_dp;
            if(!flag)           // haven't seen any positive child so far, so current is at least as good as this child
                {*lca_child = child_lca_child; *lca_parent = child_lca_parent; flag = 1;}
            else                                // have seen one positive child, so current is definitely better 
                {*lca_child = node;  *lca_parent = parent;}
        }
                                                                                            #if DEBUG && DEBUG_REC 
                                                                                                printf("in node %d, after child %d, child dp is %d, my dp is %d, my lca is (%d %d)\n", node, tree->adj_list[node][i].dest, child_dp, *dp, *lca_child, *lca_parent); 
                                                                                            #endif      
    }
    // If no child has any dp then the fields remain as -1 upon returning
                                                                                            #if DEBUG && DEBUG_REC 
                                                                                                printf("debug: at the end %d\n", node); 
                                                                                            #endif
    return 0;
}

void dfs_backtrack(int node, int parent, int findme, BT * tree, int * next_to_start, int * found){
    int i;
    int child_found = 0;
    if(found) return; // someone found it first (but nothing should be able to reach here since it's a dfs)
    if(node == findme){ // gottem 
        *found = 1;
        *next_to_start = node;
    } else { // search next
        for(i = 0; i < tree->degree[node]; i++){
            if(tree->adj_list[node][i].dest == parent) continue;

            dfs_backtrack(tree->adj_list[node][i].dest, node, findme, tree, next_to_start, &child_found);
            if(child_found){
                *next_to_start = node; // propagate the neighbor node all the way to the start
                break;
            }
        }
    }
}

int dfs_preorder(int node, int parent, BT * tree, int * in_building, int flag){
    int i;
    if(tree->degree[node] == 1)
        if(in_building[tree->master_idx_map[node]])
            in_building[tree->master_idx_map[node]] = flag;

    for(i = 0; i < tree->degree[node]; i++){
        if(tree->adj_list[node][i].dest == parent) continue;
        dfs_preorder(tree->adj_list[node][i].dest, node, tree, in_building, flag);
    }

    return 0;
}

int bfs_vote_implementation(BT * tree, int valid_start, int valid_end, int valid_start_parent, int * mapping, int * edge_child, int * edge_parent, int n, float ** d, int x, float q0, int revote_power, int ** revote_map, int all_quartets, int revoting, char ** name_map, char * master_tree){
    // Initialize the first edge to some arbitrary array
    int edge_count;
    int qs, qe; // mock queue start/end pointers
    int queue[n]; // mock queue
    float vote[n];
    int parent_map[n];
    int parent_idx = -1;
    int child_idx;

    // State
    int cur_vertex;
    int i, j;

    // Quartet method
    int parent_sample;
    int child_sample[2];
    int child_idx_a[2];
    int child_counter;
    int up;
    int u[2];
    int quartet_result;
    float M;

    float max_vote;
    int best_c, best_p;

    // Stat tracker
    max_vote = revoting ? (revote_map[0] ? 0.0 : -10000.0) : 0.0;
    best_c = valid_start;
    best_p = valid_start_parent;

    // Initialization
    int edge_counter = 0;

    memset(vote, 0, n * sizeof(float));
    memset(parent_map, -1, n * sizeof(int));
    memset(queue, -1, n * sizeof(int));
    if(!revoting) {
        *revote_map = malloc(n * sizeof(int));
        for(i = 0; i < n; i++)
            (*revote_map)[i] = 1; // initially, all internal vertices are 'ties'
    }

    parent_map[valid_start] = valid_start_parent;   
    for(i = 0; i < tree->degree[valid_start]; i++)
        if(tree->adj_list[valid_start][i].dest == valid_start_parent)
            tree->adj_list[valid_start][i].master_idx = 0;
    
    for(i = 0; i < tree->degree[valid_start_parent]; i++)
        if(tree->adj_list[valid_start_parent][i].dest == valid_start)
            tree->adj_list[valid_start_parent][i].master_idx = 0;

    edge_count = 1;
    qs = qe = 0;

    // Check valid pointer
                                                                                            #if 1
                                                                                                // printf("debug: assert that all valid start are inner nodes\n"); 
                                                                                                if(valid_start != -1 && tree->degree[valid_start] != 3){
                                                                                                    printf("failed, valid start is %d with deg %d\n", valid_start, tree->degree[valid_start]);
                                                                                                    while(1);
                                                                                                }
                                                                                            #endif

    // First vertex in the BFS
    queue[qe++] = valid_start; 
                                                                                            #if DEBUG && DEBUG_BFS 
                                                                                                 printf("debug: voting array is with max cote %f\n", max_vote);
                                                                                                 for(i = 0; i < n - 1; i++){
                                                                                                    printf("%f ", vote[i]);
                                                                                                 } 
                                                                                                 printf("\n");
                                                                                            #endif
    while(qs != qe){
        cur_vertex = queue[qs++];

        // We only deal with inner edges
        if(tree->degree[cur_vertex] == 1 || cur_vertex == valid_end || cur_vertex == valid_start_parent) continue;
                                                                                            #if DEBUG && DEBUG_BFS 
                                                                                                 printf("debug: voting array is with max cote %f\n", max_vote);
                                                                                                 for(i = 0; i < n - 1; i++){
                                                                                                    printf("%f ", vote[i]);
                                                                                                 } 
                                                                                                 printf("\n");
                                                                                            #endif
        // Work on better logic for this part
        for(i = 0; i < 3; i++){
            if(tree->adj_list[cur_vertex][i].dest == parent_map[cur_vertex]){
                parent_idx = i;
                parent_sample = tree->adj_list[cur_vertex][i].sample;
            }
        }

        child_counter = 0;
        for(i = 0; i < 3; i++)
            if(i != parent_idx){
                child_idx_a[child_counter] = i; 
                child_sample[child_counter] = tree->adj_list[cur_vertex][i].sample;
                child_counter++;


                for(j = 0; j < tree->degree[tree->adj_list[cur_vertex][i].dest]; j++)
                    if(tree->adj_list[tree->adj_list[cur_vertex][i].dest][j].dest == cur_vertex)
                        break;
                
                tree->adj_list[cur_vertex][i].master_idx = edge_count;
                tree->adj_list[tree->adj_list[cur_vertex][i].dest][j].master_idx = edge_count;
                edge_count++;
            }

        // 4 point
        up = tree->adj_list[cur_vertex][parent_idx].sample;
        u[0] = tree->adj_list[cur_vertex][child_idx_a[0]].sample;
        u[1] = tree->adj_list[cur_vertex][child_idx_a[1]].sample;


                                                                                                            #if DEBUG && DEBUG_BFS 
                                                                                                                printf("debug: current-node is %d\n", cur_vertex);
                                                                                                                printf("debug: in bfs, indexing is %d %d %d for %d %d %d\n", up, u[0], u[1], tree->adj_list[cur_vertex][child_idx_a[0]].master_idx, tree->adj_list[cur_vertex][child_idx_a[1]].master_idx,tree->adj_list[cur_vertex][parent_idx].master_idx); 
                                                                                                                printf("debug: in bfs, mapped is %d %d %d\n", mapping[up], mapping[u[0]], mapping[u[1]]);
                                                                                                                printf("test all mem %f %f %f %f %f %f\n", d[mapping[up]][mapping[u[0]]], d[mapping[up]][mapping[u[1]]], d[mapping[up]][x], d[mapping[u[0]]][mapping[u[1]]], d[mapping[u[0]]][x],d[mapping[u[1]]][x]);
                                                                                                            #endif


        M = MAX(MAX(MAX(MAX(MAX(d[mapping[up]][mapping[u[0]]], d[mapping[up]][mapping[u[1]]]), d[mapping[up]][x]), d[mapping[1]][mapping[u[0]]]), d[mapping[1]][x]), d[x][mapping[u[0]]]); 

        for(i = 0; i < 2; i++) // base case, invalid vote
            vote[tree->adj_list[cur_vertex][child_idx_a[i]].master_idx] = vote[tree->adj_list[cur_vertex][parent_idx].master_idx];
                                                                                                            #if  DEBUG && DEBUG_BFS 
                                                                                                                printf("index is %d %d, vote is %f %f\n", tree->adj_list[cur_vertex][child_idx_a[i]].master_idx, tree->adj_list[cur_vertex][parent_idx].master_idx, vote[tree->adj_list[cur_vertex][child_idx_a[i]].master_idx], vote[tree->adj_list[cur_vertex][parent_idx].master_idx]); 
                                                                                                            #endif
        if(all_quartets || M - 8.0 * q0 <= EPS){ // if all quartets are voting, we don't care about validity
            // Do quartet method here
#if use_four_point_method 
            four_point_method(d, mapping[up], mapping[u[0]], mapping[u[1]], x, &quartet_result);
#elif use_induced_quartet


            induced_quartets(name_map[mapping[up]], name_map[mapping[u[0]]], name_map[mapping[u[1]]], name_map[x], &quartet_result, master_tree);
#elif use_ml_method

#endif
            if(quartet_result == 0)  // parent wins
                for(i = 0; i < 2; i++)
                    vote[tree->adj_list[cur_vertex][child_idx_a[i]].master_idx] -= POWER(1.0 / M, revote_power);

            else if (quartet_result == 1)// u1 wins
                vote[tree->adj_list[cur_vertex][child_idx_a[0]].master_idx] += POWER(1.0 / M, revote_power);

            else if (quartet_result == 2)
                vote[tree->adj_list[cur_vertex][child_idx_a[1]].master_idx] += POWER(1.0 / M, revote_power);
        }
            
        
        for(i = 0; i < 2; i++){
            if((*revote_map)[tree->adj_list[cur_vertex][child_idx_a[i]].master_idx] // if we are revoting, we make sure that only latest best are eligible
                && vote[tree->adj_list[cur_vertex][child_idx_a[i]].master_idx] - max_vote > EPS){
                    max_vote = vote[tree->adj_list[cur_vertex][child_idx_a[i]].master_idx];
                    best_c = tree->adj_list[cur_vertex][child_idx_a[i]].dest;
                    best_p = cur_vertex;
            }
        }
        for(i = 0; i < tree->degree[cur_vertex]; i++){
            child_idx = tree->adj_list[cur_vertex][i].dest;

            // Stopping condition
            if(child_idx == parent_map[cur_vertex]) continue;

            parent_map[child_idx] = cur_vertex;
            edge_counter++;
            queue[qe++] = tree->adj_list[cur_vertex][i].dest;
        }
    }
                                                                                                            #if DEBUG && DEBUG_BFS 
                                                                                                                printf("debug: fin bf %d\n", n); 
                                                                                                            #endif
    *edge_child = best_c;
    *edge_parent = best_p;

    if(!revoting) // first round, write down revote
        for(i = 0; i < n - 1; i++)
            (*revote_map)[i] = vote[i] - max_vote < EPS && vote[i] - max_vote > -EPS;
    else  // others
        for(i = 0; i < n - 1; i++)
            (*revote_map)[i] = ((*revote_map)[i]) && (vote[i] - max_vote < EPS && vote[i] - max_vote > -EPS);
                                                                                            #if 0
                                                                                                 printf("debug: voting array is with max cote %f\n", max_vote);
                                                                                                 for(i = 0; i < n - 1; i++){
                                                                                                    printf("%f ", vote[i]);
                                                                                                 } 
                                                                                                 printf("\n");
                                                                                            #endif
    // printf("edge counter is %d\n", edge_counter);
    return 0;
}

int less_than_3_shared_taxa(INC_GRP * meta){
    int cap_counter = 0;
    int j;

    cap_counter = 0;
    for(j = 0; j < meta->n_taxa; j++)
        cap_counter += meta->visited[j] > 0;

    return cap_counter < 3; 
}

void all_valid(INC_GRP * meta, VOTE_GRP * vote){
    vote->st_lca.p = meta->gtree->adj_list[0][0].dest; 
    vote->st_lca.c = 0;
}


