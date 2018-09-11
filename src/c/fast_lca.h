#ifndef FAST_LCA_H
#define FAST_LCA_H

#include "c_inc.h"

typedef struct lca_t{
	// LCA stuff
	int 	* euler_tour;
	int 	* level;
	int 	* first_occurence;
	int 	* visited;
	BT 		* tree;
	RMQ_T 	* RMQ;

	// Smallest distance stuff
	double 	* d_from_root;
} LCA_T;

extern int fast_lca_init(LCA_T * LCA);
extern int fpm_on_tree(LCA_T * LCA, int p, int u1, int u2, int x, int * ret, double * weight);

#endif 