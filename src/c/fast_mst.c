#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "fast_mst.h"
#include "dist.h"
#include "prim.h"

int max_subset_size_mult = 10;

int int_swap(int* a, int x, int y){
	int z;

	z = a[x];
	a[x] = a[y];
	a[y] = z;

	return 0;
}

// Internal functions 
int init_graph(GRAPH * graph, int n){
	int i;

	graph->n = n;
	graph->adjacency_list = malloc(n * sizeof(ADJ_LIST*));

	for(i = 0; i < n; i++)
		graph->adjacency_list[i] = NULL; // no edges at the beginning

	graph->last_node = malloc(n * sizeof(ADJ_LIST*));
	for(i = 0; i < n; i++)
		graph->last_node[i] = NULL;

	return 0;
}

int add_edge(GRAPH * graph, int u, int v){
	if(!graph->adjacency_list[u]){
		graph->adjacency_list[u] = malloc(sizeof(ADJ_LIST));
		graph->adjacency_list[u]->dest = v;
		graph->adjacency_list[u]->next = NULL;
		graph->last_node[u] = graph->adjacency_list[u]; 

		graph->adjacency_list[v] = malloc(sizeof(ADJ_LIST));
		graph->adjacency_list[v]->dest = u;
		graph->adjacency_list[v]->next = NULL;
		graph->last_node[v] = graph->adjacency_list[v]; 
	} else {
		graph->last_node[u]->next = malloc(sizeof(ADJ_LIST));
		graph->last_node[u]->next->dest = v;
		graph->last_node[u]->next->next = NULL;
		graph->last_node[u] = graph->last_node[u]->next;

		graph->last_node[v]->next = malloc(sizeof(ADJ_LIST));
		graph->last_node[v]->next->dest = u;
		graph->last_node[v]->next->next = NULL;
		graph->last_node[v] = graph->last_node[u]->next;
	}
	return 0;
}

int random_centroids(char ** data, char * distance_model, int n, int seed, int ** disjoint_subset, int * centroid_arr, double ** small_dist_mat, int ** overlapping_subset){
	int * permu; 
	int * centroid_size_count;
	int * fast_disjoint_subset;
	int condition_flag;
	int sqrt_n; // useful precalculation
	int i, j;
	int num_site;
	char * tmp_dist_data[2];

	sqrt_n = (int) sqrt(1.0 * n);
	num_site = strlen(data[0]);

	centroid_size_count = malloc(sqrt_n * sizeof(int));
	for(i = 0; i < sqrt_n; i++)
		centroid_size_count[i] = 0;

	disjoint_subset = malloc(n * sizeof(int*)); // we don't actually construct the subset but maintain an array of indices
	for(i = 0; i < n; i++){
		disjoint_subset[i] = malloc(sqrt_n * sizeof(int));
		for(j = 0; j < sqrt_n; j++)
			disjoint_subset[i][j] = 0; 
	}

	fast_disjoint_subset = malloc(n * sizeof(int));
	for(i = 0; i < n; i++)
		fast_disjoint_subset[i] = -1;

	permu = malloc(n * sizeof(int));
	for(i = 0; i < n; i++)
		permu[i] = i;

	small_dist_mat = malloc(n * sizeof(double*));
	for(i = 0; i < n; i++)
		small_dist_mat[i] = malloc(sqrt_n * sizeof(double));	

	srand(seed); // use the random seed
	while(1){
		// Select the centroid at random by shuffling the permutation array. 
		for(i = n - 1; i > 0; i--)
			int_swap(permu, i, rand() % (i + 1)); // Fisher-Yates 

		for(i = 0; i < sqrt_n; i++)
			for(j = 0; j < n; j++)
				if(j == permu[i]) 
					small_dist_mat[j][i] = 0.0;
				else{
					if(fast_disjoint_subset[j] == -1){
						fast_disjoint_subset[j] = 0;
						tmp_dist_data[0] = data[permu[0]];
						tmp_dist_data[1] = data[j];
						small_dist_mat[j][0] = (strcmp(distance_model, "JC") == 0) ? compute_jc_distance(tmp_dist_data, num_site) : compute_logdet_distance(tmp_dist_data, num_site);
					}
					else {
						tmp_dist_data[0] = data[permu[j]];
						tmp_dist_data[1] = data[permu[i]];
						small_dist_mat[j][i] = (strcmp(distance_model, "JC") == 0) ? compute_jc_distance(tmp_dist_data, num_site) : compute_logdet_distance(tmp_dist_data, num_site);
						if(small_dist_mat[j][i] < small_dist_mat[j][fast_disjoint_subset[j]]){
							fast_disjoint_subset[j] = i;
							centroid_size_count[i]++;
						}
					}
				}

		// Check condition
		condition_flag = 1;

		for(i = 0; i < sqrt_n; i++)
			if(centroid_size_count[i] > max_subset_size_mult * sqrt_n){
				condition_flag = 0;
				break;
			}
		if(condition_flag) break;
	}

	for(i = 0; i < n; i++)
		disjoint_subset[i][fast_disjoint_subset[i]] = 1;

	centroid_arr = malloc(sqrt_n * sizeof(int));
	for(i = 0; i < sqrt_n; i++)
		centroid_arr[i] = permu[i];

	// Centroids are in their own clusters
	for(i = 0; i < sqrt_n; i++)
		for(j = 0; j < sqrt_n; j++)
			if(i == j)
				disjoint_subset[centroid_arr[i]][centroid_arr[j]] = 1;
			else 
				disjoint_subset[centroid_arr[i]][centroid_arr[j]] = 0; 

	overlapping_subset = malloc(n * sizeof(int*)); // we don't actually construct the subset but maintain an array of indices
	for(i = 0; i < n; i++){
		overlapping_subset[i] = malloc(sqrt_n * sizeof(int));
		for(j = 0; j < sqrt_n; j++)
			overlapping_subset[i][j] = disjoint_subset[i][j]; 
	}

	return 0;
}

int extend_subset(int n, int ** disjoint_subset, int p, double ** small_dist_mat){ 
	int i, j, k; 
	int sqrt_n;
	int min_idx;

	sqrt_n = (int) sqrt(1.0 * n);

	for(i = 0; i < n; i++){ // the first closest point is already chosen, we only need the remaining p - 1 point
		for(j = 0; j < p - 1; j++){
			min_idx = -1;
			for(k = 0; k < sqrt_n; k++)
				if(disjoint_subset[i][k]) continue;
				else if(min_idx == -1 || small_dist_mat[i][k] < small_dist_mat[i][min_idx])
					min_idx = k;
			
			disjoint_subset[i][min_idx] = 1; 
		}
	}

	return 0; 
}

int make_cluster_cliques(GRAPH * graph, int n, int ** disjoint_subset, int * centroid_arr){
	int * clique_v;
	int i, j, k;
	int sqrt_n;
	int clique_count;

	sqrt_n = (int) sqrt(1.0 * n);
	clique_v = malloc(max_subset_size_mult * sqrt_n * sizeof(int) + sizeof(int)); // one more for the centroid

	for(i = 0; i < sqrt_n; i++){
		clique_count = 0;
		for(j = 0; j < n; j++){
			if(disjoint_subset[j][i]) 
				clique_v[clique_count++] = j;
		}
		clique_v[clique_count++] = centroid_arr[i];
		for(j = 0; j < clique_count; j++)
			for(k = j + 1; k < clique_count; k++)
				add_edge(graph, clique_v[j], clique_v[k]); // need to check when there are overlaps
	}

	return 0;
}

int make_centroid_clique(GRAPH * graph, int n, int * centroid_arr){
	int i, j;
	int sqrt_n = (int) sqrt(1.0 * n);

	for(i = 0; i < sqrt_n; i++)
		for(j = i+1; j < sqrt_n; j++)
			add_edge(graph, centroid_arr[i], centroid_arr[j]);

	return 0; 
}

int make_complete_bipartite(GRAPH * graph, int n, int * centroid_arr){
	int * centroid_map;
	int sqrt_n;
	int i, j;

	sqrt_n = (int) sqrt(1.0 * n);
	centroid_map = malloc(n * sizeof(int)); 

	for(i = 0; i < n; i++)
		centroid_map[i] = 0; 

	for(i = 0; i < sqrt_n; i++) 
		centroid_map[centroid_arr[i]] = 1;

	for(i = 0; i < sqrt_n; i++)
		for(j = 0; j < n; j++)
			if(j != centroid_arr[i] && !centroid_map[j])
				add_edge(graph, j, centroid_arr[i]);

	return 0;
}

int nn_from_sequences(GRAPH * graph, int ** disjoint_subset, char ** data, char * distance_model, int n, int seed){
	int ** overlapping_subset = NULL;
	int * centroid_arr = NULL;
	double ** small_dist_mat = NULL;

	int p = 10;

	// Randomly select the centroid and create the X_i's, this also copies disjoint subset to overlapping one
	random_centroids(data, distance_model, n, seed, disjoint_subset, centroid_arr, small_dist_mat, overlapping_subset);

	// Create the Y_i's by adding additional close centroids
	extend_subset(n, overlapping_subset, p, small_dist_mat);

	// Make cliques from the Y_i's
	make_cluster_cliques(graph, n, overlapping_subset, centroid_arr);

	// Make centroid clique
	make_centroid_clique(graph, n, centroid_arr);

	// Make complete bipartite and centroid clique
	make_complete_bipartite(graph, n, centroid_arr); 

	return 0;
}

// This is an approximate mst that does not depend on the distance matrix and runs in time O(n^1.5)
// Input is a set of sequences in a char array 
// Output is an mst as defined in prim.c (since the last step of this function is always to call the prim algorithm)
int fast_mst(char ** data, int n, char * distance_model, int seed, MST_GRP * mst, int ** disjoint_subset){
	GRAPH * small_graph = NULL; 

	// Compute NN graph
	init_graph(small_graph, n);
	nn_from_sequences(small_graph, disjoint_subset, data, distance_model, n, seed);

	// Call Prim's algorthm on the small graph
	prim_on_small_graph(n, small_graph, mst, distance_model, data);

	return 0; 
}