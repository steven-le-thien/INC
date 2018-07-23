#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tree.h"
#include "traversal.h"
#include "utilities.h"
#include "options.h"

// All trees should be BTs 
int dfs_lca_change(int node, int parent, BT * tree, int * dp, int * in_building, int * lca_parent, int mode){
	int i;
	int child_dp; 
	int child_lca_parent;
	int child_dfs_return;
	int ret;


	if(tree->degree[node] == 1){
		if((!mode && in_building[tree->index_in_master_name_map[node]]) ||
			(mode && in_building[tree->index_in_master_name_map[node]] == mode)){
			
			*dp = 1;
			*lca_parent = parent;
			return node;
		} else {
			*dp = 0;
			return -1;
		}
	}

	*dp = 0;
	ret = -1;
	for(i = 0; i < tree->degree[node]; i++){
		if(tree->adj_list[i] == parent) continue;

		child_dfs_return = dfs_counting(tree->adj_list[i], node, goal, tree, &child_dp, &child_lca_parent);

		if(child_dp){
			*dp += child_dp;
			if(ret == -1) 			// haven't seen any positive child so far, so current is at least as good as this child
				{ret = child_dfs_return; *lca_parent = child_lca_parent;}
			else if(ret != node) 	// have seen one positive child, so current is definitely better 
				{ret = node;  *lca_parent = parent;}
		}
	}

	return ret;
}

void dfs_preorder(int node, int parent, BT * tree, int * in_building, int flag){
	int i;

	if(tree->degree[node] == 1)
		if(in_building[tree->index_in_master_name_map[node]])
			in_building[tree->index_in_master_name_map[node]] = flag;

	for(i = 0; i < tree->degree[node]; i++){
		if(tree->adj_list[i] == parent) continue;
		dfs_preorder(tree->adj_list[i], node, tree);
	}
}
