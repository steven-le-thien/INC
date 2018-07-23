#ifndef TRAVERSAL_H
#define TRAVERSAL_H

extern int dfs_lca_change(int node, int parent, BT * tree, int * dp, int * in_building, int mode);
extern void dfs_preorder(int node, int parent, BT * tree, int * in_building, int flag);

#endif 