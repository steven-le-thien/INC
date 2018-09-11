#ifndef FAST_RMQ_H
#define FAST_RMQ_H

typedef struct rmq{
    // Statistics
    int         n;
    int         num_blk;
    int         num_exhaust;
    int         blk_sz;

    // Arrays and mappings
    int         * a;
    int         * blk_idx_to_config;

    // Tables
    int         ** sparse_table;        // num_blk x log2(num_blk) array, where A[i][j] is the index of the block with the smallest value from the jth block to the j + 2^i block
    int         *** exhaustive_table;   // config_num x blk_sz x blk_sz array where A[i][j][k] is the relative index of the smallest value from the jth position to the kth position for the ith configuration
} RMQ_T;

extern int fast_rmq_init(int n, int * a,  RMQ_T * RMQ);
extern int fast_rmq(int qs, int qe, int * min_idx, RMQ_T * RMQ);

#endif