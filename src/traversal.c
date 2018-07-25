#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tree.h"
#include "traversal.h"
#include "utilities.h"
#include "options.h"


int init_vote(INC_GRP * meta, VOTE_GRP * vote){
	// Init vote's field
	vote->n_taxa 		= meta->n_taxa;

	vote->vote 			= malloc(4 * vote->n_taxa * sizeof(int));
	if(!vote->vote)	 		PRINT_AND_RETURN("malloc failure in init_vote\n", GENERAL_ERROR);

	// Set all other special edges to 0 at one go. WARNING: redo this when the struct changes	
	memset(&(vote->valid_st), -1, sizeof(int) * 10);
	return 0;
}

int dfs_lca_wrapper(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote, int i){
	int placeholder;
	int bp_counter = 0;
	int cur_node;
	int j; 

	if(dfs_lca_implementation(map->master_to_cidx[mst->prim_ord[i]],	/* strarting node 	root the tree at the query taxon */
								-1, 									/* avoiding node 	no parent since it's a leaf :( */
								meta->ctree[vote->ctree_idx], 			/* tree 			get the correct constraint tree */
								&placeholder,							/* dp value 		we won't use this (except for checking whether we are getting back the correct value) */
								meta->visited,							/* array mapping bipartion */
								&(vote->c_lca.p), 						/* the address to write the parent of the resuling process to*/
								&(vote->c_lca.c), 
								0) != SUCCESS) 				PRINT_AND_RETURN("dfs_lca_implementation failed in wrapper\n", GENERAL_ERROR);

	for(j = 0; j < 3; j++){
		cur_node = meta->ctree[vote->ctree_idx]->adj_list[vote->c_lca.c][j];
	    if(cur_node != vote->c_lca.p)
	        if(dfs_preorder(cur_node, 
	        				vote->c_lca.c, 
	        				meta->ctree[vote->ctree_idx], 
	        				meta->visited, 
	        				++bp_counter) != SUCCESS) 		PRINT_AND_RETURN("dfs_preorder failed in wrapper\n", GENERAL_ERROR);

	}
	return 0;
}

int find_valid_subtree(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote){
	int placeholder;

	// Run from the first leaf. We root the tree at the first taxon. If the taxon is in the constraint tree, we find the other bipartition; else we find the 1st biparition
	if(dfs_lca_implementation(0, 
						-1, 
						meta->gtree, 
						&placeholder,
						meta->visited,
						&(vote->st_lca.p),
						&(vote->st_lca.c),
						meta->visited[mst->prim_ord[0]] ? 3 - meta->visited[mst->prim_ord[0]] : 1)  
								!= SUCCESS)		PRINT_AND_RETURN("dfs_lca_implementation failed in find_valid_subtree\n", GENERAL_ERROR);

	// Find lca of the other bipartition that was not found in the first phase, by starting from the first phase lca, avoiding its own parent ()
	if(dfs_lca_implementation(vote->st_lca.p, 		// notice the order switching here, since the parent will actualy be the node that leads to the other bipartition lca
						vote->st_lca.c, 
						meta->gtree, 
						&placeholder,
						meta->visited,
						&(vote->nd_lca.p),
						&(vote->nd_lca.c), 
						meta->visited[mst->prim_ord[0]] ? meta->visited[mst->prim_ord[0]] : 2) 
								!= SUCCESS)		PRINT_AND_RETURN("dfs_lca_implementation failed in find_valid_subtree\n", GENERAL_ERROR);
}

int bfs_vote(INC_GRP * meta, MST_GRP * mst, VOTE_GRP * vote, int i){
	if(bfs_vote_implementation(meta->gtree,  			// growing tree
							vote->st_lca.p, 			// starting nodes and parents, notice the reverse of ordering since the subtree is within the growing tree
							vote->st_lca.c,
							vote->nd_lca.p,				// ending nodes and parents, notice the reverse of ordering
							vote->nd_lca.c,
							vote->vote,					// voting arrays
							vote->edge_c, 			
							vote->edge_c,
							vote->n_taxa, 				
							meta->d,					// distance matrix and its limit
							mst->prim_ord[i], 			// query taxon
							mst->max_w); 				// q0
		 						!= SUCCESS) 	PRINT_AND_RETURN("bfs_vote_implementation failed in bfs_voite\n", GENERAL_ERROR);

	return 0;
}

int find_insertion_edge(INC_GRP * meta, VOTE_GRP * vote){
	if(find_insertion_edge_implementation(vote->vote,
											vote->edge_c,
											vote->edge_p,
											&(vote->ins.c),
											&(vote->ins.p),
											meta->gtree->n_node -1)
								!= SUCCESS) 	PRINT_AND_RETURN("find_insertion_edge_implementation faied in wrapper\n", GENERAL_ERROR);

	return 0;
}




// All trees should be BTs 
int dfs_lca_implementation(int node, int parent, BT * tree, int * dp, int * in_building, int * lca_parent, int * lca_child, int mode){
	int i;
	int child_dp; 
	int child_lca_parent;
	int child_lca_child;
	int child_dfs_return;
	int ret;

	// Base case
	if(tree->degree[node] == 1){
		if((!mode && in_building[tree->master_idx_map[node]]) ||
			(mode && in_building[tree->master_idx_map[node]] == mode)){

			*dp = 1;
			*lca_parent = parent;
			*lca_child = node;
			return 0;
		} else {
			*dp = 0;
			return 0;
		}
	}

	// Initialization
	*dp = 0;
	ret = -1;

	child_lca_parent 	= -1;
	child_lca_child 	= -1;
	child_dfs_return 	= -1;

	for(i = 0; i < tree->degree[node]; i++){
		if(tree->adj_list[node][i].dest == parent) continue;

		child_dfs_return = dfs_counting(tree->adj_list[node][i].dest, node, goal, tree, &child_dp, &child_lca_parent, &child_lca_child, mode);

		if(child_dp){
			*dp += child_dp;
			if(child_lca_child == -1) 			// haven't seen any positive child so far, so current is at least as good as this child
				{*lca_child = child_lca_child; *lca_parent = child_lca_parent;}
			else 							 	// have seen one positive child, so current is definitely better 
				{*lca_child = node;  *lca_parent = parent;}
		}
	}
	// If no child has any dp then the fields remain as -1 upon returning

	return ret;
}

void dfs_preorder(int node, int parent, BT * tree, int * in_building, int flag){
	int i;

	if(tree->degree[node] == 1)
		if(in_building[tree->index_in_master_name_map[node]])
			in_building[tree->index_in_master_name_map[node]] = flag;

	for(i = 0; i < tree->degree[node]; i++){
		if(tree->adj_list[node][i].dest == parent) continue;
		dfs_preorder(tree->adj_list[node][i].dest, node, tree);
	}
}

void bfs_vote_implementation(BT * tree, int valid_start, int valid_end, int valid_start_parent, int * vote, int * edge_child, int * edge_parent, int n, float ** d, int x, float q0){
	// Initialize the first edge to some arbitrary array
	int edge_count;
	int qs, qe; // mock queue start/end pointers
	int queue[n]; // mock queue
	int parent_map[n];
	int parent_idx;
	int child_idx;

	// State
	int cur_vertex;
	int i;

	// 4 point
	int parent_sample;
	int child_sample[2];
	int child_idx[2];
	int child_counter;
	int up, u1, u2;
	float up1, up2, upx, m;

	memset(vote, 0, (n - 1) * sizeof(int));
	memset(edge_child, -1, (n - 1) * sizeof(int));
	memset(edge_parent, -1, (n - 1) * sizeof(int));
	memset(parent_map, -1, n * sizeof(int));
	memset(queue, -1, n * sizeof(int));
	edge_count = 0;
	qs = qe = 0;

	// Check valid pointer

	// First vertex in the BFS
	visited[0] = 1;
	edge_count++;
	queue[qe++] = valid_start; 

	while(qs != qe){
		cur_vertex = queue[qs++];

		if(tree->degree[cur_vertex] == 1 && cur_vertex != 0) continue;

		// Work on better logic for this part
		for(i = 0; i < 3; i++)
			if(tree->adj_list[cur_vertex][i].dest == parent_map[cur_vertex]){
				parent_idx = i;
				parent_sample = tree->adj_list[cur_vertex][i].dest_sample;
			}

		child_counter = 0;
		for(i = 0; i < 3; i++)
			if(i != parent_idx){
				child_idx[child_counter] = i; 
				child_sample[child_counter] = tree->adj_list[cur_vertex][i].dest_sample;
				child_counter++;
			}

		// 4 point
		up = tree->adj_list[node][parent_idx].dest_sample;
		u1 = tree->adj_list[node][0].dest_sample;
		u2 = tree->adj_list[node][1].dest_sample;

		upu1 = d[up][u1] + d[u2][x];
		upu2 = d[up][u2] + d[u1][x];
		upx = d[up][x] + d[u1][u2];
		m = max(max(u1u2, u1u3), u1x);

		// TODO: work on better logic for this part
		if(m <= 8.0 * q0){ 	// vote is valid
			if(m == upx){ 	// parent wins
				for(i = 0; i < 2; i++)
					vote[tree->adj_list[cur_vertex][child_idx[i]].idx] = vote[tree->adj_list[cur_vertex][parent_idx].idx] - 1;
			} else if (m == upu1){ // u1 wins
				vote[tree->adj_list[cur_vertex][child_idx[0]].idx] = vote[tree->adj_list[cur_vertex][parent_idx].idx] + 1;
				vote[tree->adj_list[cur_vertex][child_idx[1]].idx] = vote[tree->adj_list[cur_vertex][parent_idx].idx];
			}
		} else // vote is not valid, the children has the same vote as parents 	
			for(i = 0; i < 2; i++)
				vote[tree->adj_list[cur_vertex][child_idx[i]].idx] = vote[tree->adj_list[cur_vertex][parent_idx].idx];

		
		for(i = 0; i < tree->degree[cur_vertex]; i++){
			child_idx = tree->adj_list[cur_vertex][i].dest;

			// Stopping condition
			if(child_idx == parent_map[cur_vertex] || 
				child_idx == valid_start_parent ||
				child_idx == valid_end ||
				tree->degree[child_idx] == 1) continue;

			edge_child[*edge_count] = child_idx;
			edge_parent[*edge_count] = cur_vetex;
			tree->adj_list[cur_vertex][i].idx = (*edge_count)++;
			parent_map[child_idx] = cur_vertex;

			queue[qe++] = tree->adj_list[cur_vertex][i].dest;
		}
	}
}

int find_addition_edge(int * vote, int * edge_child, int *edge_parent, int * additional_edge_child, int * addition_edge_parent, int num_edges){
	int i;
	int max_vote;
	int max_vote_index;

	max_vote = vote[0]; //check for null pointers first
	for(i = 0; i < num_edges; i++){
		if(vote[i] > max_vote){
			max_vote = vote[i];
			max_vote_index = i;
		}
	}

	*addition_edge_parent 	= edge_parent[max_vote_index];
	*additional_edge_child 	= edge_child[max_vote_index];

	return 0;
}


