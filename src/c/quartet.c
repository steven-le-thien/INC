// File in inc_ml, created by Thien Le in July 2018
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utilities.h"
#include "tools.h"
#include "quartet.h"

char STOCK_MSA_NAME[] = "tmp.msa";
char STOCK_QUARTET_NAME[] = "tmp.quartet";

int match_and_print_msa(char * aln, char ** name_map, int * u);
int make_quartet(DIST_MOD distance_model);
int process_quartet(char ** name_map, int * u, int * res);
int do_quartet_ll(
    char ** name_map, 
    int * u, 
    DIST_MOD distance_model, 
    int * res, 
    double * M);

// Implementation of FPM. Input is a distance matrix and a quartet, u1 - u3 are 
//    known leaves, x is query taxon, results shows which leave is the sibling 
//    of the query taxon
int four_point_method(float ** d, int * u, int * res){
  int i;
  float tmp, m = 1e9;

  for(i = 0; i < 3; i++){
    tmp = d[u[i]][u[3]] + ((i == 1) ? d[u[0]][u[2]] : d[u[1]][u[2 - i]]);
    if(m > tmp){
      m = tmp;
      *res = i;
    }
  }

  return 0;
}

int new_quartets_raxml(
    char ** name_map, 
    int * u, 
    int * res, 
    ml_options * master_ml_options)
{
  FCAL(
      GENERAL_ERROR,
      F_MATCH_PRINT_IN_NEW_Q_RXML,
      match_and_print_msa(master_ml_options->input_alignment, name_map, u) 
  ); 

  FCAL(
      GENERAL_ERROR,
      F_MK_Q_IN_NEW_Q_RXML,
      make_quartet(master_ml_options->distance_model) 
  );

  FCAL(
      GENERAL_ERROR,
      F_PROCESS_Q_IN_NEW_Q_RXML,
      process_quartet(name_map, res, u)
  );
  
  SYSCAL(
      GENERAL_ERROR, 
      ERR_RM, 
      "rm %s %s %s", 
      STOCK_QUARTET_NAME, 
      STOCK_MSA_NAME, 
      "RAxML_*"
  ); 

  return 0;
} 

int ml_quartet(
    char ** name_map, 
    int * u, 
    int * res, 
    ml_options * master_ml_options, 
    double * M)
{
  FCAL(
      GENERAL_ERROR,
      F_MATCH_PRINT_IN_ML_Q,
      match_and_print_msa(master_ml_options->input_alignment, name_map, u) 
  ); 

  FCAL(
      GENERAL_ERROR,
      F_DO_Q_LL_IN_ML_Q,
      do_quartet_ll(name_map, u, master_ml_options->distance_model, res, M)
  );

  SYSCAL(
      GENERAL_ERROR, 
      ERR_RM, 
      "rm %s %s %s", 
      STOCK_QUARTET_NAME, 
      STOCK_MSA_NAME, 
      "RAxML_*"
  );
  return 0;
}

int match_and_print_msa(char * aln, char ** name_map, int * u){
  FILE * f, * p;
  int i;

  char buf[GENERAL_BUFFER_SIZE];

  f = fopen(STOCK_MSA_NAME, "w");
  p = SAFE_FOPEN_RD(aln);

  while(fscanf(p, "%s", buf) >= 0){
    for(i = 0; i < QUAD + 1; i++)
      if(i == QUAD) // no match
        fscanf(p, "%s", buf);
      else if(STR_EQ(&buf[1], name_map[u[i]])){
        fprintf(f, "%s\n", buf);
        fscanf(p, "%s", buf);
        fprintf(f, "%s\n", buf);
      }
  }
  fclose(f);
  fclose(p);

  return 0;
}

int make_quartet(DIST_MOD distance_model){
  char buf[GENERAL_BUFFER_SIZE];
  FILE * f = SAFE_FOPEN_RD(STOCK_QUARTET_NAME);

  if(f){
    fclose(f);
    ASSERT(GENERAL_ERROR, W_TMP_Q_EXIST, 0);
  } 

  FCAL(
      GENERAL_ERROR,
      F_RXML_JOB_IN_MAKE_Q,
      raxml_job(
          distance_model, 
          STOCK_QUARTET_NAME, 
          STOCK_MSA_NAME, 
          NULL
      )
  );

  sprintf(buf, "RAxML_bestTree.%s", STOCK_QUARTET_NAME);

  FCAL(
      GENERAL_ERROR,
      F_RM_LBL_JOB_IN_MAKE_Q,
      rm_label_job(buf, STOCK_QUARTET_NAME)
  );
  return 0;
}




int process_quartet(char ** name_map, int * u, int * res){
  char sib[2][GENERAL_BUFFER_SIZE];
  FILE * f = SAFE_FOPEN_RD(STOCK_QUARTET_NAME);
  int counter = 0, yes = 0;
  int idx[2], i;
  char cur_char;

  while(fscanf(f, "%c", &cur_char) >= 0){
    if(cur_char == '(') {
      yes = 1;
      counter = 0;
      continue;
    }
    if(yes == 1){
      if(cur_char == ',') {
        sib[0][counter] = 0;
        yes = 2;
        counter = 0;
      } else {
        sib[0][counter++] = cur_char;
      }
    } else if(yes == 2){
      if(cur_char == ')'){
        sib[1][counter] = 0;
        break;
      } else {
        sib[1][counter++] = cur_char;
      }
    }
  }
  fclose(f);

  for(i = 0; i < 2; ++i)
    for(idx[i] = 0; idx[i] < QUAD; ++idx[i])
      if(STR_EQ(sib[i], name_map[u[idx[i]]]))
        break;

  if(idx[0] == 3)
    *res = idx[1];
  else if(idx[1] == 3)
    *res = idx[0];
  else
    *res = (0 + 1 + 2 + 3) - idx[0] - idx[1] - 3;

  return 0;
}


int do_quartet_ll(
    char ** name_map, 
    int * u, 
    DIST_MOD distance_model, 
    int * res, 
    double * M)
{
  int i;
  double ll[3], m;
  char cq[GENERAL_BUFFER_SIZE], buf[GENERAL_BUFFER_SIZE];
  FILE * f = fopen(STOCK_QUARTET_NAME, "r"); 


  if(f){
    fclose(f);
    ASSERT(GENERAL_ERROR, W_TMP_Q_EXIST, 0);
  }

  STR_CLR(cq);
  for(i = 0; i < 3; ++i){
    sprintf(
        buf,
        "((%s,%s),(%s,%s))\n", 
        name_map[u[i]], 
        name_map[u[3]], 
        name_map[u[(i + 1) % 3]], 
        name_map[u[(i + 2) % 3]]
    );
    strcat(cq, buf);
  }

  FCAL(
      GENERAL_ERROR,
      F_GET_LL_IN_DO_Q_LL,
      get_ll_from_raxml(
          distance_model,
          STOCK_QUARTET_NAME,
          STOCK_MSA_NAME,
          NULL, 
          cq, 
          ll
      )
  );

  m = -1e9;
  for(i = 0; i < 3; ++i)
    if(ll[i] - m > -EPS){
      m = ll[i]; 
      *res = i;
      *M = 1.0 / 
          (1.0 + exp(1ll * (-ABS(ll[i] - ll[(i + 1) % 3]))) +
          exp(1ll * (-ABS(ll[i] - ll[(i + 2) % 3]))));
    }

  if(*M < 0.001) *M = 0.001;
  *M = 1.0/(*M);

  return 0;
}
