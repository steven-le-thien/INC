#ifndef TRAVERSAL_H
#define TRAVERSAL_H

extern int dfs_lca_change(int node, int parent, BT * tree, int * dp, int * in_building, int mode);
extern void dfs_preorder(int node, int parent, BT * tree, int * in_building, int flag);
extern void bfs_vote(BT * tree, int valid_start, int valid_end, int valid_start_parent, int * vote, int * edge_child, int * edge_parent, int n, float ** d, int x, float q0);
extern int find_addition_edge(int * vote, int * edge_child, int *edge_parent, int * additional_edge_child, int * addition_edge_parent, int num_edges);
#endif 