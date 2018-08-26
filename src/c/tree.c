// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tree.h"
#include "utilities.h"
#include "options.h"
#include "tools.h"

// Global variable for quickly intializing a tree
BT_edge adj_list[MAXN][3];
int degree[MAXN];
int parent_map[MAXN];
int master_idx_map[MAXN];

// Private functions
int attach_leaf_to_edge_impl(BT * growing_tree, int x, int addition_edge_parent, int additional_edge_child, int adjacent_in_mst);
BT * tree_constructor(int n, BT_edge ** adj_list, int * degree, int * master_idx_map);
BT * read_newick(MAP_GRP * map, char * filename, int tree_idx);
int init();
int check(int node_a, BT * tree);
int swap_adajcency(int nod_a, int node_b, int node_c, BT* tree);
int make_adjacent(int node_a, int node_b, BT * tree);
int make_parent(int node_p, int node_c, BT * tree);
int save_name(int cur_node, char * name, MAP_GRP * map, int tree_idx);
void petite_dfs(int node, int parent, char ** name_map, char * builder, BT * tree);

/* Read all constraint trees at once into memory
 * Input:       meta        meta variables, including the pointers to constraint trees in memory
 *              map         map variables, including the indexing map and name map
 *              options     user input options, including the name of constraint trees
 * Output: 0 on success, ERROR otherwise
 * Effect: allocate some memories, build some trees, open some files, init 2 arrays in map and 1 in meta 
 */
int parse_tree(INC_GRP * meta, MAP_GRP * map, option_t * options){
    int i;  // loop variable

    // Init meta's field
    if(meta->master_ml_options->use_four_point_method_with_tree){
        meta->n_ctree = options->num_trees - 1;
        map->master_to_midx    = malloc(meta->n_taxa  * sizeof(int));
        meta->dm = malloc(meta->n_taxa * sizeof(float*));
        if(!meta->dm)           PRINT_AND_RETURN("malloc failed to construct dm in parse_tree\n", MALLOC_ERROR);
        for(i = 0; i < meta->n_taxa; i++){
            meta->dm[i] = malloc(meta->n_taxa * sizeof(float));
        }

        construct_unweighted_matrix_job(options->tree_names[0], options->output_name, meta->dm, map->master_to_name, map->master_to_midx);
    } else 
        meta->n_ctree = options->num_trees;
 

    // Malloc sequence
    meta->ctree             = malloc(meta->n_ctree * sizeof(BT *));
    map->master_to_ctree    = malloc(meta->n_taxa  * sizeof(int));
    map->master_to_cidx     = malloc(meta->n_taxa  * sizeof(int));
    map->master_to_gidx     = malloc(meta->n_taxa  * sizeof(int));

    if(!meta->ctree)            PRINT_AND_RETURN("malloc failed for constraint trees in parse_tree\n", MALLOC_ERROR);
    if(!map->master_to_ctree)   PRINT_AND_RETURN("malloc failed for master_to_ctree in parse_tree\n", MALLOC_ERROR);
    if(!map->master_to_cidx)    PRINT_AND_RETURN("malloc failed for master_to_cidx in parse_tree\n", MALLOC_ERROR);
    if(!map->master_to_gidx)    PRINT_AND_RETURN("malloc failed for master_to_gidx in parse_tree\n", MALLOC_ERROR);

    // Init
    for(i = 0; i < meta->n_taxa; i++){
        map->master_to_ctree[i]     = -1;
        map->master_to_cidx[i]      = -1;
        map->master_to_gidx[i]      = -1;        
    }

    // Call the newick reader
    for(i = 0; i < meta->n_ctree; i++){
        printf("i is %d\n", i);
        if(meta->master_ml_options->use_four_point_method_with_tree)
            meta->ctree[i]      = read_newick(map, options->tree_names[i + 1], i);
        else 
            meta->ctree[i]      = read_newick(map, options->tree_names[i], i);
    }

                                                                                            #if DEBUG 
                                                                                                printf("debug: number of taxa is %d, number of constraint trees is %d\n", meta->n_taxa, meta->n_ctree); 
                                                                                                printf("debug: high level read newick tree debug, the following print out number of nodes of each constraint trees and their sum\n");
                                                                                                int sum = 0;
                                                                                                for(i = 0; i < meta->n_ctree; i++){
                                                                                                    // printf("%d ", meta->ctree[i]->n_node);
                                                                                                    sum += meta->ctree[i]->n_node;
                                                                                                }
                                                                                                printf("\n");
                                                                                                printf("sum is %d\n", sum);


                                                                                                printf("debug: high level, the following print out the number of nodes with degree 1 and 3 and their sum compared to total nodes\n");
                                                                                                int j,k1, k3, k4, kt;
                                                                                                kt = 0;
                                                                                                for(i = 0; i < meta->n_ctree; i++){
                                                                                                    k1 = 0;
                                                                                                    k3 = 0;
                                                                                                    k4 = 0;
                                                                                                    for(j = 0; j < meta->ctree[i]->n_node; j++){
                                                                                                        k1 += meta->ctree[i]->degree[j] == 1;
                                                                                                        k3 += meta->ctree[i]->degree[j] == 3;
                                                                                                        k4 += (meta->ctree[i]->degree[j] != 1) && ( meta->ctree[i]->degree[j] != 3);
                                                                                                    }
                                                                                                    // printf("for tree %d, node with deg 1 is %d, 3 is %d, total is %d, n_node is %d, non13 is %d\n", i, k1, k3, k1 + k3, meta->ctree[i]->n_node, k4);
                                                                                                    kt += k1;
                                                                                                }
                                                                                                printf("debug: total number of leaf is %d\n", kt);

                                                                                                printf("debug: the following tests actual adjacency\n");
                                                                                                int flag = 0; //this sholuld be 0 all the way until the end
                                                                                                for(i = 0; i < meta->n_ctree; i++){
                                                                                                    BT * tree = meta->ctree[i];

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
                                                                                                }
                                                                                                if(!flag) printf("passed test\n");
                                                                                                else printf("failed test with flag %d\n", flag);

                                                                                                printf("debug: this test tests mapping from ctree to master index\n");
                                                                                                flag = 0;
                                                                                                for(i = 0; i< meta->n_ctree; i++){
                                                                                                    BT * tree = meta->ctree[i];
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
                                                                                                }


                                                                                                if(!flag) printf("passed test\n");
                                                                                                else printf("failed test\n");
                                                                                             

                                                                                                // printf("debug: checking 2 master maps, all should be non empty, the first oen should be unique\n");
                                                                                                // for(i = 0; i < meta->n_taxa; i++){
                                                                                                //     printf("%d ", map->master_to_ctree[i]);
                                                                                                // }
                                                                                                // printf("\n");
                                                                                                // for(i = 0; i < meta->n_taxa; i++){
                                                                                                //     printf("%d ", map->master_to_cidx[i]);
                                                                                                // }
                                                                                                // printf("\n");
                                                                                                // while(1);
                                                                                            #endif
    
                                                                                                    
    return 0;
}


/* Creating a 3-leaf star from the firt 3 nodes visited in the Prim's ordering of the mst
 * Input:       meta        meta variables, including the pointers to growing tree in memory
 *              mst         mst variables, including the prim ordering and the mst itself
 * Output: 0 on success, ERROR otherwise
 * Effect: allocate some memories, build some trees, init the growing tree and visited array in meta
 */
int init_growing_tree(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst){ 
    int i;

    // Update master->idx mapping
    map->master_to_gidx[mst->prim_ord[0]] = 0;
    map->master_to_gidx[mst->prim_ord[1]] = 1;
    map->master_to_gidx[mst->prim_ord[2]] = 2;

    // Check that there must be at least 3 leaves
    if(meta->n_taxa < 3)                            PRINT_AND_RETURN("there are less than 3 taxa, the tree is trivial\n", GENERAL_ERROR);

    // Meta initialization
    meta->visited = malloc(meta->n_taxa * sizeof(int));
    if(!meta->visited)                              PRINT_AND_RETURN("malloc failed for meta->visited in init_growing_tree\n", GENERAL_ERROR);
    memset(meta->visited, 0, meta->n_taxa * sizeof(int));

    // Growing tree initialization
    meta->gtree = tree_constructor(4 * meta->n_taxa,  NULL, NULL, NULL); // allocate the maximum number of nodes for mallocation to reserve enough space for tree addition
    if(!meta->gtree)                                PRINT_AND_RETURN("tree construction failed in init_growing_tree\n", GENERAL_ERROR);

    // Override tree logistics since there are only 4 real nodes in the tree, the rest can be viewed as singletons
    meta->gtree->n_node = 4;

    if(make_adjacent(0, 3, meta->gtree) != SUCCESS) PRINT_AND_RETURN("make_adjacent failed in init_growing tree\n", GENERAL_ERROR);
    if(make_adjacent(1, 3, meta->gtree) != SUCCESS) PRINT_AND_RETURN("make_adjacent failed in init_growing tree\n", GENERAL_ERROR);
    if(make_adjacent(2, 3, meta->gtree) != SUCCESS) PRINT_AND_RETURN("make_adjacent failed in init_growing tree\n", GENERAL_ERROR);

    // Init tree variables  
    for(i = 0; i < 3; i++){
        meta->gtree->master_idx_map[i] = mst->prim_ord[i];
        meta->visited[mst->prim_ord[i]] = -1;
        meta->gtree->adj_list[3][i].sample = i;
    }

    // This is a bit tricky, we need to find adjacent nodes in the mst
    meta->gtree->adj_list[0][0].sample = 1; // fist vertex in prim is obviously connected to second
    meta->gtree->adj_list[1][0].sample = 0; // vice versa
    meta->gtree->adj_list[2][0].sample = map->master_to_gidx[mst->prim_par[mst->prim_ord[2]]]; // check if the parent of 3rd vertex in mst is adjacent to 0 or 1
                                                                                            #if DEBUG 
                                                                                                printf("debug: testing modules for initializing the tree\n"); 

                                                                                                printf("debug: n_node is %d (4)\n", meta->gtree->n_node);
                                                                                                printf("debug: degree of 4 node is %d %d %d %d (1, 1, 1, 3)\n", meta->gtree->degree[0], meta->gtree->degree[1], meta->gtree->degree[2], meta->gtree->degree[3]);
                                                                                                printf("debug: master index of 4 node is %d %d %d %d (%d, %d, %d, -1)\n", meta->gtree->master_idx_map[0], meta->gtree->master_idx_map[1], meta->gtree->master_idx_map[2], meta->gtree->master_idx_map[3], mst->prim_ord[0], mst->prim_ord[1], mst->prim_ord[2]);

                                                                                                printf("debug: sample of 03 13 23 and reverse is %d %d %d, %d %d %d\n", meta->gtree->adj_list[0][0].sample,meta->gtree->adj_list[1][0].sample, meta->gtree->adj_list[2][0].sample, meta->gtree->adj_list[3][0].sample, meta->gtree->adj_list[3][1].sample, meta->gtree->adj_list[3][2].sample);
                                                                                                printf("debug: dest for upper sampel is %d %d %d, %d %d %d\n", meta->gtree->adj_list[0][0].dest,meta->gtree->adj_list[1][0].dest, meta->gtree->adj_list[2][0].dest, meta->gtree->adj_list[3][0].dest, meta->gtree->adj_list[3][1].dest, meta->gtree->adj_list[3][2].dest);
                                                                                                printf("debug: the following print out non 0 places in inbuilding\n");
                                                                                                for(i = 0; i < meta->n_taxa; i++){
                                                                                                    if(meta->visited[i] != 0) printf("%d(%d) ", i, meta->visited[i]);
                                                                                                }
                                                                                                printf("\n");
                                                                                                // while(1);
                                                                                            #endif  
                                                                                            #if DEBUG 
                                                                                                BT * growing_tree = meta->gtree;
                                                                                                 // printf("debug: adjacent_in_mst is %d\n", adjacent_in_mst); 
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[0][0].sample, growing_tree->adj_list[0][1].sample, growing_tree->adj_list[0][2].sample);
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[1][0].sample, growing_tree->adj_list[1][1].sample, growing_tree->adj_list[1][2].sample);
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[2][0].sample, growing_tree->adj_list[2][1].sample, growing_tree->adj_list[2][2].sample);
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[3][0].sample, growing_tree->adj_list[3][1].sample, growing_tree->adj_list[3][2].sample);
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[4][0].sample, growing_tree->adj_list[4][1].sample, growing_tree->adj_list[4][2].sample);
                                                                                            #endif
    return 0;
}

int attach_leaf_to_edge(INC_GRP * meta,  MAP_GRP * map,MST_GRP * mst, VOTE_GRP * vote, int i){
    // Set the new edge to be present in the growing tree
    meta->visited[mst->prim_ord[i]] = -1;
                                                                                            #if DEBUG 
                                                                                                 // printf("debug: the following print out prim's tree (parenting)\n");
                                                                                                 // for(i = 0; i < meta->n_taxa; i++)
                                                                                                 //     printf("%d ", mst->prim_par[i]);
                                                                                                printf("debug: parent in mater is %d, in growing is %d\n", mst->prim_par[mst->prim_ord[i]], map->master_to_gidx[mst->prim_par[mst->prim_ord[i]]]);
                                                                                            #endif
    // Actual attachment
    if(attach_leaf_to_edge_impl(meta->gtree,
                                    mst->prim_ord[i], 
                                    vote->ins.p,
                                    vote->ins.c,
                                    map->master_to_gidx[mst->prim_par[mst->prim_ord[i]]])
                                        != SUCCESS)                 PRINT_AND_RETURN("attach leaf to edge failed in wrapper\n", GENERAL_ERROR);

    // Set the correct mapping of the index of the new edge (the internal tree mapping is set in the function above)
    map->master_to_gidx[mst->prim_ord[i]] = meta->gtree->n_node - 1;
    return 0;

}


int write_newick(BT * tree, char * filename, char ** name_map){
    char buffer[MAX_BUFFER_SIZE * 100];
    FILE * f;

    // Subdivide an edge and make it the root
    tree->n_node++;

    swap_adajcency(0, tree->n_node - 1, tree->adj_list[0][0].dest, tree);
    strclr(buffer);
    petite_dfs(tree->n_node - 1, -1, name_map, buffer, tree);
    strcat(buffer, ";");

    f = fopen(filename, "w");
    fprintf(f, "%s", buffer);
    fclose(f);

    return 0;
}

// INTERNAL FUNCTION IMPLEMENTATIONS
int attach_leaf_to_edge_impl(BT * growing_tree, int x, int additional_edge_parent, int additional_edge_child, int adjacent_in_mst){
    int i;
    // Changing tree logistics (if this fails or terminate in the middle of some process then we are screwed)
    // Because of the way we mallocate the growing tree, there are already memory for the extra 2 nodes 
    growing_tree->n_node += 2;

                                                                                            #if DEBUG 
                                                                                                 printf("debug: adjacent_in_mst is %d\n", adjacent_in_mst); 
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[0][0].sample, growing_tree->adj_list[0][1].sample, growing_tree->adj_list[0][2].sample);
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[1][0].sample, growing_tree->adj_list[1][1].sample, growing_tree->adj_list[1][2].sample);
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[2][0].sample, growing_tree->adj_list[2][1].sample, growing_tree->adj_list[2][2].sample);
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[3][0].sample, growing_tree->adj_list[3][1].sample, growing_tree->adj_list[3][2].sample);
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[4][0].sample, growing_tree->adj_list[4][1].sample, growing_tree->adj_list[4][2].sample);
                                                                                            #endif
    // Subdivide the edge
    if(swap_adajcency(additional_edge_parent, growing_tree->n_node - 2, additional_edge_child, growing_tree)
                    != SUCCESS) PRINT_AND_RETURN("first swap_adajcency failed in attach_leaf_to_edge_impl\n", GENERAL_ERROR);

    if(make_adjacent(growing_tree->n_node - 1, growing_tree->n_node - 2, growing_tree) 
                    != SUCCESS) PRINT_AND_RETURN("make_adjacent failed in attach_leaf_to_edge_impl\n", GENERAL_ERROR);
    // Set up the leaf sample properly
    for(i = 0; i < 3; i++)
        if(growing_tree->adj_list[growing_tree->n_node - 2][i].dest == growing_tree->n_node - 1)
            growing_tree->adj_list[growing_tree->n_node - 2][i].sample = growing_tree->n_node - 1;
        
    growing_tree->adj_list[growing_tree->n_node - 1][0].sample = adjacent_in_mst;
                                                                                            #if DEBUG 
                                                                                                 printf("debug: adjacent_in_mst is %d\n", adjacent_in_mst); 
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[0][0].sample, growing_tree->adj_list[0][1].sample, growing_tree->adj_list[0][2].sample);
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[1][0].sample, growing_tree->adj_list[1][1].sample, growing_tree->adj_list[1][2].sample);
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[2][0].sample, growing_tree->adj_list[2][1].sample, growing_tree->adj_list[2][2].sample);
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[3][0].sample, growing_tree->adj_list[3][1].sample, growing_tree->adj_list[3][2].sample);
                                                                                                 printf("%d %d %d\n", growing_tree->adj_list[4][0].sample, growing_tree->adj_list[4][1].sample, growing_tree->adj_list[4][2].sample);
                                                                                                // while(1);
                                                                                            #endif
    growing_tree->master_idx_map[growing_tree->n_node - 1] = x;

    return 0;
}

BT * tree_constructor(int n, BT_edge ** adj_list, int * degree, int * master_idx_map){
    // The arrays mayeb null (in case of null construction)
    int i, j;

    BT * tree =     malloc(sizeof(BT)); 

    if(!tree)           PRINT_AND_RETURN("malloc error in tree_constructor\n", NULL);

    tree->n_node                     = n;

    tree->adj_list              = malloc(n * sizeof(BT_edge *));
    tree->degree                = malloc(n * sizeof(int));
    tree->master_idx_map        = malloc(n * sizeof(int));

    for(i = 0; i < n; i++){
        tree->adj_list[i] = malloc(3 * sizeof(BT_edge));
    }
    // Initialization 
    for(i = 0; i < n; i++){
        if(adj_list && adj_list[i]) memcpy(tree->adj_list[i], adj_list[i], 3 * sizeof(BT_edge));
        else 
            for(j = 0; j < 3; j++){
                tree->adj_list[i][j].dest = -1;
                tree->adj_list[i][j].sample = -1;
                tree->adj_list[i][j].master_idx = -1;
            }

        tree->master_idx_map[i]
                                    = -1;
        tree->degree[i]             = 0; 
    }

    if(degree)  
        memcpy(tree->degree, degree, n * sizeof(int));

    if(master_idx_map)
        memcpy(tree->master_idx_map, master_idx_map, n * sizeof(int));

    return tree;
}

#define READ_NAME_STATE 0
#define OTHER_STATE     1

#define incr_node(node, max_node)               node = ++max_node
#define incr_level(node_p, node_c, max_node)    do{node_p = node_c; incr_node(node_c, max_node);} while(0)
#define decr_level(node_p, node_c)              do{node_c = node_p; node_p = parent_map[node_c];} while(0)

BT * read_newick(MAP_GRP * map, char * filename, int tree_idx){
    FILE *  f;
    int i, j;

    BT_edge ** mock;

    // FSM variables
    int     cur_state;

    int     cur_node;
    int     max_node;
    int     parent_node;

    char    cur_char;
    char    cur_name[MAX_NAME_SIZE];

    // File opening
    f = fopen(filename, "r");
    if(!f)                      PRINT_AND_RETURN("fail to open in read_newick file\n",    NULL);

    // Intialization
    cur_state       = READ_NAME_STATE;

    cur_node        = 0;
    max_node        = 0;
    parent_node     = -1;

    cur_char        = 0;
    strclr(cur_name);

    // Initialize global variables
    init();

    while(fscanf(f, "%c", &cur_char) == 1){
        switch(cur_char){
            case '(': //start of a new level and start of a new node
                incr_level(parent_node, cur_node, max_node);
                if(make_parent(parent_node, cur_node, NULL)     != SUCCESS) PRINT_AND_RETURN("make parent failed in read newick\n", NULL);
                cur_state = READ_NAME_STATE;
                break;
            case ',': // start of new node end of old node
                if(save_name(cur_node, cur_name, map, tree_idx) != SUCCESS) PRINT_AND_RETURN("save name failed in read_newick\n", NULL);
                incr_node(cur_node, max_node);
                if(make_parent(parent_node, cur_node, NULL)     != SUCCESS) PRINT_AND_RETURN("make parent failed in read newick\n", NULL);
                cur_state = READ_NAME_STATE;
                break;
            case ')': // end of level
                if(save_name(cur_node, cur_name, map, tree_idx) != SUCCESS) PRINT_AND_RETURN("save name failed in read_newick\n", NULL);
                decr_level(parent_node, cur_node); 
                cur_state = OTHER_STATE;
                break;
            case ':': 
                cur_state = OTHER_STATE;
                break;
            default:
                if(cur_state == READ_NAME_STATE){
                    cat_c_to_a(cur_char, cur_name);
                }
        }
    }

    fclose(f);
    // Deal with the case where the root has degree 2
    if(degree[0] == 2){
        // Find its 2 children, 1 is guaranteed to be 1
        for(i = 0; i < max_node + 1; i++)
            if(parent_map[i] == 0 && i != 1) break;

        for(j = 0; j < degree[1]; j++)
            if(adj_list[1][j].dest == -1) adj_list[1][j].dest = i - 1;

        for(j = 0; j < degree[i]; j++)
            if(adj_list[i][j].dest == -1) adj_list[i][j].dest = 0;
    }

    mock = malloc(max_node * sizeof(BT_edge));
    for(i = 0;i < max_node; i++)
        mock[i] = &(adj_list[1 + i][0]);
    // Make the binary tree and return
    return tree_constructor(max_node, mock, &degree[1], &master_idx_map[1]);
}


int init(){
    memset(parent_map, -1, sizeof(parent_map));
    memset(adj_list, -1, sizeof(adj_list));
    memset(degree, 0, sizeof(degree));
    memset(master_idx_map, -1, sizeof(degree));
    return 0;
}

/* Helper function to check whether a node is in range
 * Input:   node to check
 * Output:  1 if the node is in range, 0 otherwise
 * Effect:  none
 */ 
int check(int node_a, BT * tree){
    if(!tree) return node_a >= 0; // no upper limit
    if(0 <= node_a && node_a < tree->n_node) return 1;
    else return 0;
}

// Delete a-c and make a-b, b-c 
int swap_adajcency(int node_a, int node_b, int node_c, BT* tree){
    int i, j;
    if(!check(node_a, tree) || !check(node_b, tree))  PRINT_AND_RETURN("node out of range in swap_adajcency", GENERAL_ERROR);

    if(tree){
        for(i = 0; i < tree->degree[node_a]; i++)
            if(tree->adj_list[node_a][i].dest == node_c){
                tree->adj_list[node_a][i].dest = node_b; // this is the only thing we need to change since the index will be refreshed in the next iteration
                                                         // and the leaf sample remains the same
                break; // preserve i for later use
            }
        for(j = 0; j < tree->degree[node_c]; j++)
            if(tree->adj_list[node_c][j].dest == node_a){
                tree->adj_list[node_c][j].dest = node_b;
                break; // preserve j for later use
            }

        // Assuming b is new so its fields are empty
        tree->adj_list[node_b][tree->degree[node_b]].dest = node_a;
        tree->adj_list[node_b][tree->degree[node_b]].sample = tree->adj_list[node_c][j].sample; // we can do this since we know b is not a leaf
        tree->degree[node_b]++;

        tree->adj_list[node_b][tree->degree[node_b]].dest = node_c;
        tree->adj_list[node_b][tree->degree[node_b]].sample = tree->adj_list[node_a][i].sample; // we can do this since we know b is not a leaf
        tree->degree[node_b]++;

    }
    return 0;
}

/* Helper function to make a node adjacent to another in the adjacency matrix
 * Input:   2 nodes
 * Output:  0 on success, ERROR otherwise
 * Effect:  none
 */ 
int make_adjacent(int node_a, int node_b, BT * tree){
    if(!check(node_a, tree) || !check(node_b, tree))  PRINT_AND_RETURN("node out of range in make_adjacent", GENERAL_ERROR);

    if(tree){
        tree->adj_list[node_a][tree->degree[node_a]++].dest = node_b;
        tree->adj_list[node_b][tree->degree[node_b]++].dest = node_a;
    } else {
        adj_list[node_a][degree[node_a]++].dest = node_b - 1;
        adj_list[node_b][degree[node_b]++].dest = node_a - 1;
    }
    

    return 0;
}

/* Helper function to make a node a parent of another 
 * Input:   2 nodes (order matters, the parent node comes first)
 * Output:  0 on success, ERROR otherwise
 * Effect:  none
 */ 
int make_parent(int node_p, int node_c, BT * tree){
    if(make_adjacent(node_p, node_c, tree) != SUCCESS)    PRINT_AND_RETURN("node out of range in make_parent", GENERAL_ERROR);
    parent_map[node_c] = node_p;
    return 0;
}

/* Helper function to save name for a numeric node and map it back to the master ordering
 * Input:   node number and its name 
 * Output:  0 on success
 * Effect:  none
 */ 
int save_name(int cur_node, char * name, MAP_GRP * map, int tree_idx){
    int i; // store the master index
    for(i = 0; i < map->n_taxa; i++){
        if(strcmp(name, map->master_to_name[i]) != 0){
            // printf("%s, %s\n", name, map->master_to_name[i]);
            if(i == map->n_taxa - 1) return 0; // no name matches, which is fine for internal nodes
        } else break;
    }
    if(map){
        map->master_to_ctree[i]     = tree_idx;
        map->master_to_cidx[i]      = cur_node - 1;
    }

    master_idx_map[cur_node]    = i;

    strclr(name);

    return 0;
}

void petite_dfs(int node, int parent, char ** name_map, char * builder, BT * tree){
    int i;
    if(tree->degree[node] == 1) {
        strcat(builder, name_map[tree->master_idx_map[node]]);
        return;
    }

    strcat(builder, "(");
    for(i = 0; i < tree->degree[node]; i++){
        if(tree->adj_list[node][i].dest == parent) continue;

        petite_dfs(tree->adj_list[node][i].dest, node, name_map, builder, tree);
        if(i != tree->degree[node] - 1) strcat(builder, ",");
    }

    strcat(builder, ")");
}


