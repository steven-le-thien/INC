// File in HMMDecompositionDecision, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tree.h"
#include "utilities.h"
#include "options.h"


// Global variable for quickly intializing a tree
int adj_list[maxN][3];
int degree[maxN];
int parent_map[maxN];
int is_in_building_tree[maxN];
char * name_map[maxN][maxNameSize];

// Private functions
int size_dfs(int node, int parent);
int centroid_search(int node, int parent);
int check(int node_a);
int make_adjacent(int node_a, int node_b);
int make_parent(int node_p, int node_c);
int save_name(int cur_node, char * name);
int init();

/* Brute force algorithm to test whether a node is a leaf (has degree exactly 1 in a binary tree)
 * Input:   node to be tested
 * Output:  1 if the node is a leave, 0 otherwise
 * Effect:  none
 */ 
int is_leave(int cur_node, BT * tree){
    if(!check(cur_node, BT * tree))        PRINT_AND_RETURN("cur_node does not exists in the tree in is_leave", GENERAL_ERROR);
    return degree[cur_node] == 1;
}

BT * tree_constructor(int n, int root, int ** adj_list, int * degree, int * parent_map, int * is_in_building_tree; char ** name_map){
    // The arrays mayeb null (in case of null construction)

    int i, j;
    BT * tree = malloc(sizeof(BT)); 

    if(!tree)           PRINT_AND_RETURN("malloc error in tree_constructor", NULL);

    tree->n                     = n;
    tree->root                  = root;

    tree->adj_list              = malloc(n * sizeof(int *));
    tree->degree                = malloc(n * sizeof(int));
    tree->parent_map            = malloc(n * sizeof(int)); 
    tree->is_in_building_tree   = malloc(n * sizeof(int));
    tree->name_map              = malloc(n * sizeof(char *))

    for(i = 0; i < n; i++){
        tree->adj_list[i] = malloc(3 * sizeof(int));
        tree->name_map[i] = malloc(maxNameSize * sizeof(char)); 
    }

    // Initialization to default values; this may or may not be override in other constructors
    for(i = 0; i < n; i++){
        if(adj_list[i]) memcpy(tree->adj_list[i], adj_list[i], 3 * sizeof(int));
        else 
            for(j = 0; j < 3; j++)
                tree->adj_list[i][j] = -1;

        if(name_map[i]) memcpy(tree->name_map[i], name_map[i], maxNameSize * sizeof(char));
        else name_map[i][0] = 0;

        tree->degree[i]             = 0; 
        tree->parent_map[i]         = -1;
        tree->is_in_building_tree   = 0;
    }

    if(degree)  
        memcpy(tree->degree, degree, n * sizeof(int));

    if(parent_map)  
        memcpy(tree->parent_map, parent_map, n * sizeof(int));

    if(is_in_building_tree)  
        memcpy(tree->is_in_building_tree, is_in_building_tree, n * sizeof(int));

    return tree;
}

void tree_destructor(BT * tree){
    int i; // loop counter

    if(tree->adj_list){                  
        for(i = 0; i < tree->n; i++)
            if(tree->adj_list[i])
                free(tree->adj_list[i]);

        free(tree->adj_list);
    }

    if(tree->degree)                    free(tree->degree[i]);
    if(tree->parent_map)                free(tree->parent_map[i]);
    if(tree->is_in_building_tree)       free(tree->is_in_building_tree[i]);

    free(tree);
}

#define READ_NAME_STATE 0
#define OTHER_STATE     1

#define incr_level(node_p, node_c, max_node) do{node_p = node_c; incr_node(node_c, max_node);} while(0)
#define decr_level(node_p, node_c, max_node) do{node_c = node_p; node_p = parent_map[node_c];} while(0)
#define incr_node(node, max_node)               node = ++max_node

BT * read_newick(char * filename){
    FILE *  f;
    BT * tree;

    // FSM variabless
    int     node_counter;
    int     cur_state;
    int     cur_node;
    int     max_node;
    int     parent_node;
    char    cur_char;
    char    cur_name[maxNameSize];

    strclr(cur_name);

    // File opening
    f = freopen(filename, "r", stdin);
    if(!f) PRINT_AND_RETURN("fail to open read_newick file",    OPEN_ERROR);

    // FSM for reading newick format
    cur_state = READ_NAME_STATE;
    node_counter = 0; 
    cur_node = 0;
    max_node = 0;
    parent_node = -1;
    cur_char = 0;
    cur_name[maxNameSize];

    init();

    while(scanf("%c", &cur_char) == 1){
        // printf("cur char is %c cur node is %d cur par is %d cur state is %d\n", cur_char, cur_node, parent_node, cur_state);
        switch(cur_char){
            case '(': //start of a new level and start of a new node
                incr_level(parent_node, cur_node, max_node);
                make_parent(parent_node, cur_node);
                cur_state = READ_NAME_STATE;
                node_counter++;
                break;
            case ',': // start of new node end of old node
                save_name(cur_node, cur_name);
                incr_node(cur_node, max_node);
                make_parent(parent_node, cur_node);
                cur_state = READ_NAME_STATE;
                node_counter++;
                break;
            case ')': // end of level
                save_name(cur_node, cur_name);
                decr_level(parent_node, cur_node, max_node, tree);
                break;
            case ':': 
                cur_state = OTHER_STATE;
                break;
            default:
                if(cur_state == READ_NAME_STATE){
                    char buf[2];
                    buf[0] = cur_char;
                    buf[1] = 0;
                    strcat(cur_name, buf);
                }
        }
    }

    fclose(f);

    // Make the binary tree and return
    tree = tree_constructor(node_counter, 0, adj_list, degree, parent_map, is_in_building_tree; name_map);
    return tree;
}

int parse_tree(BT ** constrained_trees, options_t * options){
    int i;

    constrained_trees = malloc(options->num_tree, sizeof(BT *));

    if(!constrained_trees)      PRINT_AND_RETURN("malloc error in parse_tree", MALLOC_ERROR);

    for(i = 0; i < options->num_tree; i++){
        constrained_trees[i] = read_newick(options->tree_name[i]);

        if(!constrained_trees[i]){
            i--;
            for(; i >= 0; i--) free(constrained_trees[i]);
            free(constrained_trees);
            PRINT_AND_RETURN("malloc error in parse_tree", MALLOC_ERROR);
        }
    }
    return 0;
}

int init_growing_tree(int)

// INTERNAL FUNCTION IMPLEMENTATIONS
int init(){
    int i, j;

    for(i = 0; i < maxN; i++){
        for(j = 0; j < 3; j++){
            adj_list[i][j] = - 1;
        }

        degree[i] = 0; 
        parent_map[i] = -1;
        is_in_building_tree[i] = -1;
        name_map[i][0] = 0;
    }
}

/* Helper function to check whether a node is in range
 * Input:   node to check
 * Output:  1 if the node is in range, 0 otherwise
 * Effect:  none
 */ 
int check(int node_a, BT * tree){
    if(0 <= node_a && node_a < tree->n) return 1;
    else return 0;
}

/* Helper function to make a node adjacent to another in the adjacency matrix
 * Input:   2 nodes
 * Output:  0 on success, ERROR otherwise
 * Effect:  none
 */ 
int make_adjacent(int node_a, int node_b, BT * tree){
    if(!check(node_a) || !check(node_b))  PRINT_AND_RETURN("node out of range in make_adjacent", GENERAL_ERROR);

    adj_list[node_a][degree[node_a]++] = node_b;
    adj_list[node_b][degree[node_b]++] = node_a;

    return 0;
}

/* Helper function to make a node a parent of another 
 * Input:   2 nodes (order matters, the parent node comes first)
 * Output:  0 on success, ERROR otherwise
 * Effect:  none
 */ 
int make_parent(int node_p, int node_c, BT * tree){
    // printf("make paret %d is paret f %d\n", node_p, node_c);
    if(make_adjacent(node_p, node_c, tree) != SUCCESS)    PRINT_AND_RETURN("node out of range in make_parent", GENERAL_ERROR);

    parent_map[node_c] = node_p;
    return 0;
}

/* Helper function to save name for a numeric node
 * Input:   node number and its name 
 * Output:  0 on success
 * Effect:  none
 */ 
int save_name(int cur_node, char * name){
    // printf("save_name %d %s\n", cur_node, name);
    if(strempty(name)){
        char buf[1000];
        snprintf(buf, 1000, "%d", cur_node);
        strcpy(name_map[cur_node], buf);
    } else {
        strcpy(name_map[cur_node], name);
    }
    strclr(name); 

    return 0;
}


