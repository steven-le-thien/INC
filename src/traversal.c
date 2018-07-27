#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "traversal.h"
#include "utilities.h"
#include "options.h"

#define DEBUG_REC 0
#define DEBUG_BFS 1


// All trees should be BTs 
int dfs_lca_implementation(int node, int parent, BT * tree, int * dp, int * in_building, int * lca_parent, int * lca_child, int mode);
int dfs_preorder(int node, int parent, BT * tree, int * in_building, int flag);
int bfs_vote_implementation(BT * tree, int valid_start, int valid_end, int valid_start_parent, int *, int * edge_child, int * edge_parent, int n, float ** d, int x, float q0);
int find_insertion_edge_implementation(int * vote, int * edge_child, int *edge_parent, int * additional_edge_child, int * addition_edge_parent, int num_edges);
void dfs_backtrack(int node, int parent, int findme, BT * tree, int * next_to_find_me, int * next_to_start, int * found);

int init_vote(INC_GRP * meta, VOTE_GRP * vote){
	// Init vote's field
	vote->n_taxa 		= meta->n_taxa;

	vote->vote 			= malloc(4 * vote->n_taxa * sizeof(int));
	vote->edge_c 		= malloc(4 * vote->n_taxa * sizeof(int));
	vote->edge_p 		= malloc(4 * vote->n_taxa * sizeof(int));
	if(!vote->vote)	 		PRINT_AND_RETURN("malloc failure in init_vote\n", GENERAL_ERROR);
	// Set all other special edges to 0 at one go. WARNING: redo this when the struct changes	
	memset(&(vote->valid_st), -1, sizeof(int) * 10);

	memset(vote->vote, 0, 4 * vote->n_taxa * sizeof(int));
	memset(vote->edge_c, -1, 4 * vote->n_taxa * sizeof(int));
	memset(vote->edge_p, -1, 4 * vote->n_taxa * sizeof(int));
																							#if DEBUG 
 																								printf("debug: init vote, n_taxa is %d\n", vote->n_taxa); 
																							#endif
	return 0;
}

int find_bipartition(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote, int i){
	int placeholder = 0;
	int bp_counter = 0;
	int cur_node;
	int j; 
																							#if DEBUG 
 																								printf("debug: input to dfs_lca_implementation is: starting node %d with name %s, master index %d, tree %d\n", map->master_to_cidx[mst->prim_ord[i]], map->master_to_name[meta->ctree[vote->ctree_idx]->master_idx_map[map->master_to_cidx[mst->prim_ord[i]]]], meta->ctree[vote->ctree_idx]->master_idx_map[map->master_to_cidx[mst->prim_ord[i]]], vote->ctree_idx); 
																							#endif
	if(dfs_lca_implementation(map->master_to_cidx[mst->prim_ord[i]],	/* strarting node 	root the tree at the query taxon */
								-1, 									/* avoiding node 	no parent since it's a leaf :( */
								meta->ctree[vote->ctree_idx], 			/* tree 			get the correct constraint tree */
								&placeholder,							/* dp value 		we won't use this (except for checking whether we are getting back the correct value) */
								meta->visited,							/* array mapping bipartion */
								&(vote->c_lca.p), 						/* the address to write the parent of the resuling process to*/
								&(vote->c_lca.c), 
								0) != SUCCESS) 				PRINT_AND_RETURN("dfs_lca_implementation failed in wrapper\n", GENERAL_ERROR);
																							#if DEBUG 
 																								printf("debug: output of dfs_loca is %d %d (%d %d) in master index \n", vote->c_lca.p, vote->c_lca.c, meta->ctree[vote->ctree_idx]->master_idx_map[vote->c_lca.p], meta->ctree[vote->ctree_idx]->master_idx_map[vote->c_lca.c]); 
																							#endif
	// If there are less than 3 nodes then don't bother
 	if(placeholder < 3) return 0;
	for(j = 0; j < meta->ctree[vote->ctree_idx]->degree[vote->c_lca.c]; j++){
		cur_node = meta->ctree[vote->ctree_idx]->adj_list[vote->c_lca.c][j].dest;
	    if(cur_node != vote->c_lca.p)
	        if(dfs_preorder(cur_node, 
	        				vote->c_lca.c, 
	        				meta->ctree[vote->ctree_idx], 
	        				meta->visited, 
	        				++bp_counter) != SUCCESS) 		PRINT_AND_RETURN("dfs_preorder failed in wrapper\n", GENERAL_ERROR);

	}
																							#if DEBUG 
 																								printf("debug: after bipartition\n"); 
 																								printf("debug: value of visited that is non zeor\n");
 																								for(j = 0; j < meta->n_taxa; j++)
 																									if(meta->visited[j])
 																										printf("(%d) %d, ", meta->visited[j], j);
 																								printf("\n");
																							#endif
	// meta->visited[mst->prim_ord[i]] = 3;
	// printf("checking, lca is %d %d, visited is %d %d %d %d %d %d\n", vote->c_lca.p, vote->c_lca.c, meta->visited[0], meta->visited[1], meta->visited[2], meta->visited[3], meta->visited[4], meta->visited[5]);

	return 0;
}

int find_valid_subtree(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote){
	int i;
	int placeholder;

	for(i = 0; i < meta->n_taxa; i++)
		if(meta->visited[i] == 2)
			break; // there guaranteed to be 1
																							#if DEBUG 
 																								printf("debug: input to find_valid_subtree st dfs is: starting node %d with flag %d\n", map->master_to_gidx[i], meta->visited[mst->prim_ord[0]] ? 3 - meta->visited[mst->prim_ord[0]] : 1); 
																							#endif

	// Run from the first leaf. We root the tree at the first taxon. If the taxon is in the constraint tree, we find the other bipartition; else we find the 1st biparition
	if(dfs_lca_implementation(map->master_to_gidx[i], 
						-1, 
						meta->gtree, 
						&placeholder,
						meta->visited,
						&(vote->st_lca.p),
						&(vote->st_lca.c),
						1)  
								!= SUCCESS)		PRINT_AND_RETURN("dfs_lca_implementation failed in find_valid_subtree\n", GENERAL_ERROR);
																							#if DEBUG 
 																								printf("debug: output of st dfs_lca is %d %d (%d %d) (p, c) in master index \n", vote->st_lca.p, vote->st_lca.c, meta->gtree->master_idx_map[vote->st_lca.p], meta->gtree->master_idx_map[vote->st_lca.c]); 
																							#endif
	for(i = 0; i < meta->n_taxa; i++)
		if(meta->visited[i] == 1)
			break; // there guaranteed to be 1
	// Find lca of the other bipartition that was not found in the first phase, by starting from the first phase lca, avoiding its own parent ()
	if(dfs_lca_implementation(map->master_to_gidx[i], 		// notice the order switching here, since the parent will actualy be the node that leads to the other bipartition lca
						-1, 
						meta->gtree, 
						&placeholder,
						meta->visited,
						&(vote->nd_lca.p),
						&(vote->nd_lca.c), 
						2) 
								!= SUCCESS)		PRINT_AND_RETURN("dfs_lca_implementation failed in find_valid_subtree\n", GENERAL_ERROR);
	// printf("checking, lca is %d %d, lca2 is %d %d, visited is %d %d %d %d %d %d\n", vote->st_lca.p, vote->st_lca.c, vote->nd_lca.p, vote->nd_lca.c, meta->visited[0], meta->visited[1], meta->visited[2], meta->visited[3], meta->visited[4], meta->visited[5]);
																							#if DEBUG 
 																								printf("debug: output of nd dfs_lca is %d %d (%d %d) (p, c) in master index \n", vote->nd_lca.p, vote->nd_lca.c, meta->gtree->master_idx_map[vote->nd_lca.p], meta->gtree->master_idx_map[vote->nd_lca.c]); 
																							#endif
	// while(1);
 	// Find correct orientation
 	dfs_backtrack(vote->st_lca.c,
 					-1, 					// find in all direction from the starting point
 					vote->nd_lca.c, 
 					meta->gtree,
 					&(vote->st_lca.p),
 					&(vote->nd_lca.p),
 					&placeholder);
																							#if DEBUG 
 																								printf("debug: output of backtrack is %d %d (%d %d); %d %d (%d %d) (p, c) in master index \n", vote->st_lca.p, vote->st_lca.c, meta->gtree->master_idx_map[vote->st_lca.p], meta->gtree->master_idx_map[vote->st_lca.c], vote->nd_lca.p, vote->nd_lca.c, meta->gtree->master_idx_map[vote->nd_lca.p], meta->gtree->master_idx_map[vote->nd_lca.c]); 
																							#endif
	return 0;
}

int bfs_vote(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst, VOTE_GRP * vote, int i){
	int j;
	if(bfs_vote_implementation(meta->gtree,  			// growing tree
							vote->st_lca.p, 			// starting nodes and parents, notice the reverse of ordering since the subtree is within the growing tree
							vote->nd_lca.c,				// ending nodes notice the reverse of ordering
							vote->st_lca.c,
							// vote->vote,
							meta->gtree->master_idx_map,					// voting arrays
							&(vote->ins.c), 			
							&(vote->ins.p),
							meta->gtree->n_node, 				
							meta->d,					// distance matrix and its limit
							mst->prim_ord[i], 			// query taxon
							mst->max_w) 				// q0
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
	int flag;
	int ret;
																							#if DEBUG && DEBUG_REC
 																								printf("debug: in recursion, node is %d(%d), parent is %d\n", node, tree->master_idx_map[node], parent); 
																							#endif
	// printf("node is %d %d\n", node, parent);

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
	child_lca_parent 	= -1;
	child_lca_child 	= -1;
	child_dfs_return 	= -1;

	for(i = 0; i < tree->degree[node]; i++){
		if(tree->adj_list[node][i].dest == parent) continue;
		child_dfs_return = dfs_lca_implementation(tree->adj_list[node][i].dest, node, tree, &child_dp, in_building, &child_lca_parent, &child_lca_child, mode);
																							
		if(child_dp){
			*dp += child_dp;
			if(!flag) 			// haven't seen any positive child so far, so current is at least as good as this child
				{*lca_child = child_lca_child; *lca_parent = child_lca_parent; flag = 1;}
			else 							 	// have seen one positive child, so current is definitely better 
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

void dfs_backtrack(int node, int parent, int findme, BT * tree, int * next_to_find_me, int * next_to_start, int * found){
	int i;
	int child_found = 0;
	if(found) return; // someone found it first (but nothing should be able to reach here since it's a dfs)
	if(node == findme){ // gottem 
		*found = 1;
		*next_to_find_me = parent;
		*next_to_start = node;
	} else { // search next
		for(i = 0; i < tree->degree[node]; i++){
			if(tree->adj_list[node][i].dest == parent) continue;

			dfs_backtrack(tree->adj_list[node][i].dest, node, findme, tree, next_to_find_me, next_to_start, &child_found);
			if(child_found){
				*next_to_start = node;
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

int bfs_vote_implementation(BT * tree, int valid_start, int valid_end, int valid_start_parent, int * mapping, int * edge_child, int * edge_parent, int n, float ** d, int x, float q0){
	// Initialize the first edge to some arbitrary array
	int edge_count;
	int qs, qe; // mock queue start/end pointers
	int queue[n * 4]; // mock queue
	int vote[n * 4];
	int parent_map[n * 4];
	int parent_idx = -1;
	int child_idx;

	// State
	int cur_vertex;
	int i, j;
	// 4 point
	int parent_sample;
	int child_sample[2];
	int child_idx_a[2];
	int child_counter;
	int up, u1, u2;
	float up1, up2, upx, m;

	int max_vote, best_c, best_p;

	max_vote = 0;
	best_c = valid_start;
	best_p = valid_start_parent;

	memset(vote, 0, n * 4 * sizeof(int));
	// memset(edge_child, -1, (n - 1) * sizeof(int));

	// memset(edge_parent, -1, (n - 1) * sizeof(int));
				// printf("vaid start is %d valid end is %d vald start parent is %d n is %d x is %d float q0 %f\n", valid_start, valid_end, valid_start_parent, n, x, q0);

	memset(parent_map, -1, n * 4 * sizeof(int));

	memset(queue, -1, n * 4 * sizeof(int));

	parent_map[valid_start] = valid_start_parent;

	edge_count = 0;
	qs = qe = 0;

	// Check valid pointer
																							#if DEBUG && DEBUG_BFS 
 																								printf("debug: assert that all valid start are inner nodes\n"); 
 																								if(valid_start != -1 && tree->degree[valid_start] != 3){
 																									printf("failed, valid start is %d with deg %d\n", valid_start, tree->degree[valid_start]);
 																									while(1);
 																								}
 																								// if(valid_end != -1 &&tree->degree[valid_end] != 3){
 																								// 	printf("failed, valid valid_end is %d with deg %d\n", valid_end, tree->degree[valid_end]);
 																								// 	while(1);
 																								// }
																							#endif

	// First vertex in the BFS
	queue[qe++] = valid_start; 

	while(qs != qe){
		cur_vertex = queue[qs++];

		// We only deal with inner edges
		if(tree->degree[cur_vertex] == 1 || cur_vertex == valid_end || cur_vertex == valid_start_parent) continue;

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

				// edge_child[edge_count] = tree->adj_list[cur_vertex][i].dest;
				// edge_parent[edge_count] = cur_vertex;

				for(j = 0; j < tree->degree[tree->adj_list[cur_vertex][i].dest]; j++)
					if(tree->adj_list[tree->adj_list[cur_vertex][i].dest][j].dest == cur_vertex)
						break;
				
				tree->adj_list[cur_vertex][i].master_idx = edge_count;
				tree->adj_list[tree->adj_list[cur_vertex][i].dest][j].master_idx = edge_count;
				edge_count++;
			}

		// 4 point
		up = tree->adj_list[cur_vertex][parent_idx].sample;
		u1 = tree->adj_list[cur_vertex][child_idx_a[0]].sample;
		u2 = tree->adj_list[cur_vertex][child_idx_a[1]].sample;

		// printf("cur vertex is %d parent idx is %d parent is %d\n", cur_vertex, parent_idx, tree->adj_list[cur_vertex][parent_idx].dest);
		// printf("%d %d %d %d\n", up, u1, u2, x);

																											#if DEBUG && DEBUG_BFS 
 																								printf("debug: in bfs, indexing is %d %d %d for %d %d %d\n", up, u1, u2, tree->adj_list[cur_vertex][child_idx_a[0]].master_idx, tree->adj_list[cur_vertex][child_idx_a[1]].master_idx,tree->adj_list[cur_vertex][parent_idx].master_idx); 
 																								printf("debug: in bfs, mapped is %d %d %d\n", mapping[up], mapping[u1], mapping[u2]);
 																								printf("test all mem %f %f %f %f %f %f\n", d[mapping[up]][mapping[u1]], d[mapping[up]][mapping[u2]], d[mapping[up]][x], d[mapping[u1]][mapping[u2]], d[mapping[u1]][x],d[mapping[u2]][x]);
																							#endif
		up1 = d[mapping[up]][mapping[u1]] + d[mapping[u2]][x];
		up2 = d[mapping[up]][mapping[u2]] + d[mapping[u1]][x];
		upx = d[mapping[up]][x] + d[mapping[u1]][mapping[u2]];
						// printf("dw\n");

 		m = max(max(up1, up2), upx);

		// TODO: work on better logic for this part
		if(m <= 8.0 * q0){ 	// vote is valid
			if(m == upx){ 	// parent wins
				for(i = 0; i < 2; i++)
					vote[tree->adj_list[cur_vertex][child_idx_a[i]].master_idx] = vote[tree->adj_list[cur_vertex][parent_idx].master_idx] - 1;
			} else if (m == up1){ // u1 wins
				vote[tree->adj_list[cur_vertex][child_idx_a[0]].master_idx] = vote[tree->adj_list[cur_vertex][parent_idx].master_idx] + 1;
				vote[tree->adj_list[cur_vertex][child_idx_a[1]].master_idx] = vote[tree->adj_list[cur_vertex][parent_idx].master_idx];
			} else if (m == up2) {
				vote[tree->adj_list[cur_vertex][child_idx_a[1]].master_idx] = vote[tree->adj_list[cur_vertex][parent_idx].master_idx] + 1;
				vote[tree->adj_list[cur_vertex][child_idx_a[0]].master_idx] = vote[tree->adj_list[cur_vertex][parent_idx].master_idx];
			}
		} else // vote is not valid, the children has the same vote as parents 	
			for(i = 0; i < 2; i++)
				vote[tree->adj_list[cur_vertex][child_idx_a[i]].master_idx] = vote[tree->adj_list[cur_vertex][parent_idx].master_idx];

		for(i = 0; i < 2; i++){
			if(vote[tree->adj_list[cur_vertex][child_idx_a[i]].master_idx] > max_vote){
				max_vote = vote[tree->adj_list[cur_vertex][child_idx_a[i]].master_idx];
				best_c = tree->adj_list[cur_vertex][child_idx_a[i]].dest;
				best_p = cur_vertex;
			}
		}
		for(i = 0; i < tree->degree[cur_vertex]; i++){
			child_idx = tree->adj_list[cur_vertex][i].dest;

			// Stopping condition
			if(child_idx == parent_map[cur_vertex]) continue;

			// printf("%d %d %d www\n", child_idx, cur_vertex, edge_count);
			// if(child_idx == valid_start_parent ||
			// 	child_idx == valid_end ||
			// 	tree->degree[child_idx] == 1) continue;

			parent_map[child_idx] = cur_vertex;

			queue[qe++] = tree->adj_list[cur_vertex][i].dest;
		}
	}
																							#if DEBUG && DEBUG_BFS 
 																								printf("debug: fin bf \n"); 
																							#endif
	// for(int i = 0; i < 3; i++)
	// 	printf("vote is %d\n", vote[i]);
	*edge_child = best_c;
	*edge_parent = best_p;
	return 0;
}

int find_insertion_edge_implementation(int * vote, int * edge_child, int *edge_parent, int * additional_edge_child, int * addition_edge_parent, int num_edges){
	int i;
	int max_vote;
	int max_vote_index;

	max_vote = vote[0]; //check for null pointers first
	max_vote_index = 0;
																							#if DEBUG 
																								printf("debug: num edge is %d\n", num_edges);
 																								printf("debug: this prints out all vtes\n"); 
 																								for(i = 0; i < num_edges; i++)
 																									printf("%d ", vote[i]);
 																								printf("\n");
																							#endif
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


