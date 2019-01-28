#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utilities.h"
#include "fast_rmq.h"

#define TEST_RMQ 0

// Internal functions
int num_to_config(int blk_sz, int num, int * config);
int config_to_num(int blk_sz, int * config, int * num);
int config_to_full(int blk_sz, int * config, int * full);
int full_to_config(int blk_sz, int * full, int * config);
int compare_with_sparse_table(RMQ_T * RMQ, int level, int * idx);
int rmq_with_sparse_table(RMQ_T * RMQ, int bs, int be);
int rmq_in_blk_with_exhaustive_table(RMQ_T * RMQ, int blk_idx, int qs, int qe);

// Public functions
int fast_rmq_init(int n, int * a,  RMQ_T * RMQ){
  int i, j, k, l;
  int * config, * full; // a config is a bit array (0, 1)
  int tmp_a[2];
  int st_levels;

  // 1. Partition the array into sets of size logn/2 
  RMQ->n              = n;
  RMQ->a              = a;
  RMQ->blk_sz         = ((int) log(1.0 * n) / 2) > 0 ? 
                          (int) log(1.0 * n) / 2 : 1;
  // this is the real number of block, including the last one
  RMQ->num_blk        = n / RMQ->blk_sz + (n % RMQ->blk_sz != 0); 
  RMQ->num_exhaust    = 1 << (RMQ->blk_sz - 1);

  // 2. Set up exhaustive table
  config              = SAFE_MALLOC((RMQ->blk_sz - 1) * sizeof(int));
  full                = SAFE_MALLOC(RMQ->blk_sz * sizeof(int));
  full[0] = 0;

  RMQ->exhaustive_table   = SAFE_MALLOC(RMQ->num_exhaust * sizeof(int **));
  for(i = 0; i < RMQ->num_exhaust; i++){
    FCAL(
        GENERAL_ERROR,
        F_NUM_TO_CFG_IN_FAST_RMQ_INIT,
        num_to_config(RMQ->blk_sz, i, config)
    );
    FCAL(
        GENERAL_ERROR,
        F_NUM_TO_CFG_IN_FAST_RMQ_INIT,
        config_to_full(RMQ->blk_sz, config, full)
    );

    RMQ->exhaustive_table[i] = SAFE_MALLOC(RMQ->blk_sz * sizeof(int *));
    for(j = 0; j < RMQ->blk_sz; j++){
      RMQ->exhaustive_table[i][j] = SAFE_MALLOC(RMQ->blk_sz * sizeof(int));
      for(k = 0; k < RMQ->blk_sz; k++){
        RMQ->exhaustive_table[i][j][k] = j;
        for(l = j + 1; l <= k; l++)
          if(full[l] < full[RMQ->exhaustive_table[i][j][k]]) 
            RMQ->exhaustive_table[i][j][k] = l;
      }
    }
  }

  // 3. Set up block index to config
  RMQ->blk_idx_to_config = SAFE_MALLOC(RMQ->num_blk * sizeof(int));
  for(i = 0; i < RMQ->num_blk; i++){
    // Copy it to the full array
    for(j = 0; j < RMQ->blk_sz; j++)  // for the last block, we may go out of 
                                      // range, in which case we ony need to 
                                      // keep increasing the level
      full[j] = i * RMQ->blk_sz + j < n ? 
          RMQ->a[i * RMQ->blk_sz + j] : 
          full[j - 1] + 1; 

    // Make configuration array from full array and then convert to the 
    // enumeration
    FCAL(
        GENERAL_ERROR,
        F_NUM_TO_CFG_IN_FAST_RMQ_INIT,
        full_to_config(RMQ->blk_sz, full, config)
    );
    FCAL(
        GENERAL_ERROR,
        F_NUM_TO_CFG_IN_FAST_RMQ_INIT,
        config_to_num(RMQ->blk_sz, config, &RMQ->blk_idx_to_config[i])
    );
  }

  // 4. Set up a sparse table for the minimum of each block
  st_levels = (int) LOG2(1.0 * RMQ->num_blk) + 1;

  RMQ->sparse_table = malloc(st_levels * sizeof(int *));
  for(i = 0; i < st_levels; i++){
    RMQ->sparse_table[i] = malloc(RMQ->num_blk * sizeof(int));
    for(j = 0; j < RMQ->num_blk; j++){
      RMQ->sparse_table[i][j] = -1;
    }
  }
  // First row / baes case, i = 0
  for(j = 0; j < RMQ->num_blk; j++) // the rows are exactly the array 
                                    // (range 1 segments)
    RMQ->sparse_table[0][j] = j;    // RMQ->exhaustive_table
                                    //    [RMQ->blk_idx_to_config[j]]

  tmp_a[0] = j;
  tmp_a[1] = j + (1 << (i - 1));
  // Run the rest of DP
  for(i = 1; i < st_levels; i++)
    for(j = 0; j < RMQ->num_blk; j++)
      RMQ->sparse_table[i][j] = j + (1 << i) >= RMQ->num_blk ?
          -1 : 
          compare_with_sparse_table(RMQ, i - 1, tmp_a);

  return 0;
}

// Input is an array of integer point values; 2 indices. Output is the min 
//    value and index within the 2 indices.
// This implementation follows Bender and Farach-Colton
// This function must only be called after initialization (preproc) is done
int fast_rmq(int qs, int qe, int * min_idx, RMQ_T * RMQ){
  // Block index variables
  int         i, m;
  int         bs, be;
  int         ms_idx, mm_idx, me_idx;
  int         abs_idx[3];
  int         **tmp_blk;

  // 2. Get the index of the block after qs and that before qe
  bs = qs / RMQ->blk_sz + 1;
  be = qe / RMQ->blk_sz - 1;
  mm_idx = 
      bs > be ? 
      -1 : 
      rmq_with_sparse_table(RMQ, bs, be);

  ms_idx = 
      bs > be + 1 ? 
      -1 : 
      rmq_in_blk_with_exhaustive_table(
          RMQ, 
          bs - 1, 
          qs % RMQ->blk_sz,
          RMQ->blk_sz - 1
      );

  me_idx = 
      bs > be + 1 ? 
      -1 : 
      rmq_in_blk_with_exhaustive_table(
          RMQ, 
          be + 1, 
          0, 
          qe % RMQ->blk_sz
      );

  // if this happens the ms_idx and me_idx must also be -1
  if(ms_idx == -1 && mm_idx == -1) 
    *min_idx =  
        (bs - 1) * RMQ->blk_sz + 
        rmq_in_blk_with_exhaustive_table(
            RMQ, 
            bs - 1, 
            qs % RMQ->blk_sz, 
            qe % RMQ->blk_sz
        );

  else if(mm_idx == -1 && ms_idx != -1){
    abs_idx[0] = (bs - 1) * RMQ->blk_sz + ms_idx;
    abs_idx[2] = (be + 1) * RMQ->blk_sz + me_idx;
    *min_idx = 
        RMQ->a[abs_idx[0]] > RMQ->a[abs_idx[2]] ? 
        abs_idx[2] : 
        abs_idx[0];
  }
  else {
    tmp_blk = RMQ->exhaustive_table[RMQ->blk_idx_to_config[mm_idx]];

    abs_idx[0] = (bs - 1) * RMQ->blk_sz + ms_idx;
    abs_idx[2] = (be + 1) * RMQ->blk_sz + me_idx;
    abs_idx[1] = 
        mm_idx * RMQ->blk_sz + tmp_blk[0][RMQ->blk_sz - 1];

    m = (int) 1e9; 
    for(i = 0; i < 3; i++)
      if(RMQ->a[abs_idx[i]] < m){
        m = RMQ->a[abs_idx[i]];
        *min_idx = abs_idx[i];
      }
  }
  return 0;
}

// Internal implementation
int num_to_config(int blk_sz, int num, int * config){
  int i;
  for(i = 0; i < blk_sz - 1; i++){
    config[i] = num % 2;
    num /= 2;
  }
  return 0;
}

int config_to_num(int blk_sz, int * config, int * num){
  int i;
  *num = 0; 
  for(i = 0; i < blk_sz - 1; i++)
    *num += config[i] << i;
  return 0;
}

int config_to_full(int blk_sz, int * config, int * full){
  int i;
  for(i = 1; i < blk_sz; i++)
    full[i] = full[i - 1] + (config[i - 1] ? 1 : -1);
  
  return 0;
}

int full_to_config(int blk_sz, int * full, int * config){
  int i;

  for(i = 0; i < blk_sz - 1; i++)
    config[i] = full[i + 1] - full[i] == 1 ? 1 : 0;

  return 0;
}

// This takes in an RMQ struct, a level, 2 blocks index and return the index of 
//    the block with the smaller minimum
int compare_with_sparse_table(RMQ_T * RMQ, int level, int idx[2]){
  int i;

  // DP variables
  int blk_cand_idx[2];
  int rl_min_idx[2];
  int abs_min_idx[2];
  int** exhaustive_blk[2];

  for(i = 0; i < 2; ++i){
    blk_cand_idx[i] = RMQ->sparse_table[level][idx[i]];
    exhaustive_blk[i] = 
        RMQ->exhaustive_table[RMQ->blk_idx_to_config[blk_cand_idx[i]]];
    rl_min_idx[i] = exhaustive_blk[i][0][RMQ->blk_sz - 1];
    abs_min_idx[i] = blk_cand_idx[i] * RMQ->blk_sz + rl_min_idx[i];
  }

  return RMQ->a[abs_min_idx[0]] < RMQ->a[abs_min_idx[1]] ? 
      blk_cand_idx[0] : 
      blk_cand_idx[1]; 
}

int rmq_with_sparse_table(RMQ_T * RMQ, int bs, int be){
  int level = (int) LOG2(1.0 * (be - bs + 1));
  int tmp_a[] = {bs, be - (1 << level) + 1};
  return compare_with_sparse_table(RMQ, level, tmp_a);
}

// Give a relative index within a block and the block index, do rmq within the 
// range and return the relative index of the minimum
int rmq_in_blk_with_exhaustive_table(RMQ_T * RMQ, int blk_idx, int qs, int qe){
  return RMQ->exhaustive_table[RMQ->blk_idx_to_config[blk_idx]][qs][qe];
}

#if TEST_RMQ

// Unit test
int main(){
  int n, q;
  int qs, qe;
  int a[1000000];
  int res1, res2;
  int t;

  int i, j;
  RMQ_T RMQ;

  scanf("%d %d", &n, &q);
  for(i = 0; i < n; i++){
    scanf("%d", &a[i]);
  }
  fast_rmq_init(n, &a[0], &RMQ);

  for(i = 0; i < q; i++){
    scanf("%d %d", &qs, &qe);

    // Our method
    fast_rmq(qs, qe, &res1, &RMQ);
    res1 = a[res1];

    // Naive method
    res2 = a[qs];
    for(j = qs; j <= qe; j++){
      if(a[j] < res2) res2 = a[j];
    }

    // Compare
    if(res1 != res2){
      printf("wrong!!!!!\n");
      while(1);
    }
  }
    
  return 0;
}

#endif