// File in inc_ml, created by Thien Le in July 2018
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utilities.h"
#include "tools.h"


// Implementation of FPM. Input is a distance matrix and a quartet, u1 - u3 are known leaves, x is query taxon, results shows which leave is the sibling of the query taxon
int four_point_method(float ** d, int up, int u1, int u2, int x, int * res){
	float px, p1, p2, m;

	p2 = d[up][u1] + d[u2][x];
	p1 = d[up][u2] + d[u1][x];
	px = d[up][x] + d[u1][u2];

	m = MIN(MIN(p1, p2), px);

	if(ABS(m - px) < EPS) 
		*res = 0; // parent is sibling
	else if(ABS(m - p1) < EPS) 
		*res = 1; // first child is sibling
	else 
		*res = 2; // second child is sibling
	// printf("px is %f p1 is %f p2 is %f, res is %d\n", px, p1, p2, *res);
	return 0;
}

int new_quartets_raxml(char * up, char * u1, char * u2, char * x, int * res, ml_options * master_ml_options){
	char sibling1[MAX_NAME_SIZE];
	char sibling2[MAX_NAME_SIZE];
	char cur_char;
	int counter;
	int yes;

	FILE * f, * p;

	char stock_msa_name[] = "tmp.msa";
	char stock_quartet_name[] = "tmp.quartet";
	// char stock_msa_name = "tmp.msa";

	f = fopen(stock_msa_name, "w");
	p = fopen(master_ml_options->input_alignment, "r");
	// printf("hre-1, %s %s %s %s %s\n", up, u1, u2, x, master_ml_options->input_alignment);
	while(fscanf(p, "%s", sibling1) >= 0){
		// printf("curent is %s\n", sibling1);
		if(!strcmp(&sibling1[1], up) || !strcmp(&sibling1[1], u1) || !strcmp(&sibling1[1], u2) || !strcmp(&sibling1[1], x)){
			// printf("insderd %s\n", sibling1);
			fprintf(f, "%s\n", sibling1);
			fscanf(p, "%s", sibling1);
			fprintf(f, "%s\n", sibling1);
		} else 
			fscanf(p, "%s", sibling1);
	}
	fclose(f);
	fclose(p);

	// Make subtree
	f = fopen(stock_quartet_name, "r"); 
	if(f){
		fclose(f);
		PRINT_AND_RETURN("tmp.quartet already exits in working directory\n", GENERAL_ERROR);
	} 


// int raxml_job(char * dist_model,  int l_distance_model,
//               char * out_pfx,     int l_out_pfx,
//               char * in_seq,      int l_in_seq,
//               char * wd_sfx,      int l_wd_sfx)

        if(raxml_job(master_ml_options->distance_model, stock_quartet_name, stock_msa_name, NULL) != SUCCESS) PRINT_AND_RETURN("raxml_job failed in main", GENERAL_ERROR);
    // if(fasttree_job(&tmp_options, master_ml_options) != SUCCESS) PRINT_AND_RETURN("raxml_job failed in main", GENERAL_ERROR);

        sprintf(sibling1, "RAxML_bestTree.%s", stock_quartet_name);

       if(rm_label_job(sibling1, stock_quartet_name) != SUCCESS) PRINT_AND_RETURN("remove label failed in main", GENERAL_ERROR);

    // tmp_options.out

	// f = fopen("tmp.quartet", "r");  // the tree looks like ((A, B), (C, D));
	// fscanf(f, "%s", sibling1);
	// fclose(f);
	// printf("tree is %s\n", sibling1);

	// tmp_options.output_name = NULL;
	// tmp_options.input_name = stock_msa_name;
 //    tmp_options.tree_names = malloc(sizeof(char *));
 //    tmp_options.tree_names[0] = stock_quartet_name;
 //    // if(raxml_job(&tmp_options, master_ml_options) != SUCCESS) PRINT_AND_RETURN("raxml_job failed in main", GENERAL_ERROR);
 //    if(fasttree_job(&tmp_options, master_ml_options) != SUCCESS) PRINT_AND_RETURN("raxml_job failed in main", GENERAL_ERROR);

 //    // sprintf(sibling1, "RAxML_bestTree.%s", stock_quartet_name);
 //    // tmp_options.input_name = sibling1;
 //    tmp_options.input_name = stock_quartet_name;
 //    tmp_options.output_name = stock_quartet_name;
 //    if(rm_label_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("remove label failed in main", GENERAL_ERROR);

        f = fopen("tmp.quartet", "r");  // the tree looks like ((A, B), (C, D));
	fscanf(f, "%s", sibling1);
	fclose(f);
	// printf("tree is %s\n", sibling1);

	f = fopen("tmp.quartet", "r");  // the tree looks like ((A, B), (C, D));
	counter = 0;
	yes = 0;
	while(fscanf(f, "%c", &cur_char) >= 0){
		if(cur_char == '(') {
			yes = 1;
			counter = 0;
			continue;
		}
		if(yes == 1){
			if(cur_char == ',') {
				sibling1[counter] = 0;
				yes = 2;
				counter = 0;
			} else {
				sibling1[counter++] = cur_char;
 			}
		} else if(yes == 2){
			if(cur_char == ')'){
				sibling2[counter] = 0;
				break;
			} else {
				sibling2[counter++] = cur_char;
			}
		}
	}
	fclose(f);

	// Set result
	if((strcmp(sibling1, up) == 0 && strcmp(sibling2, u1) == 0) ||
		(strcmp(sibling2, up) == 0 && strcmp(sibling1, u1) == 0) ||
		(strcmp(sibling1, u2) == 0 && strcmp(sibling2, x) == 0) ||
		(strcmp(sibling2, u2) == 0 && strcmp(sibling1, x) == 0)) *res = 2;

	else if((strcmp(sibling1, up) == 0 && strcmp(sibling2, u2) == 0) ||
		(strcmp(sibling2, up) == 0 && strcmp(sibling1, u2) == 0) ||
		(strcmp(sibling1, u1) == 0 && strcmp(sibling2, x) == 0) ||
		(strcmp(sibling2, u1) == 0 && strcmp(sibling1, x) == 0)) *res = 1;

	else *res = 0;
	// printf("siblin %s %s %d\n", sibling1, sibling2, *res);
	// Remmove tmp files (careful, this is destructive)
	system("rm tmp.msa tmp.quartet RAxML_*");

	return 0;
} 

int ml_quartet(char * up, char * u1, char * u2, char * x, int * res, ml_options * master_ml_options, double * M,  int * revote_power){
	// option_t tmp_options;
	char cq[GENERAL_BUFFER_SIZE];
	double ll_p, ll_1, ll_2;
	FILE * f, * p;

	char stock_msa_name[] = "tmp.msa";
	char stock_quartet_name[] = "tmp.quartet";
	// char stock_msa_name = "tmp.msa";

	f = fopen(stock_msa_name, "w");
	p = fopen(master_ml_options->input_alignment, "r");
	// printf("hre-1, %s %s %s %s %s\n", up, u1, u2, x, master_ml_options->input_alignment);
	while(fscanf(p, "%s", cq) >= 0){
		// printf("curent is %s\n", sibling1);
		if(!strcmp(&cq[1], up) || !strcmp(&cq[1], u1) || !strcmp(&cq[1], u2) || !strcmp(&cq[1], x)){
			// printf("insderd %s\n", sibling1);
			fprintf(f, "%s\n", cq);
			fscanf(p, "%s", cq);
			fprintf(f, "%s\n", cq);
		} else 
			fscanf(p, "%s", cq);
	}
	fclose(f);
	fclose(p);

	// Make subtree
	f = fopen(stock_quartet_name, "r"); 
	if(f){
		fclose(f);
		PRINT_AND_RETURN("tmp.quartet already exits in working directory\n", GENERAL_ERROR);
	} 

	// tmp_options.output_name = NULL;
	// tmp_options.input_name = stock_msa_name;
 //    tmp_options.tree_names = malloc(sizeof(char *));
 //    tmp_options.tree_names[0] = stock_quartet_name;

	// parent 
	sprintf(cq, "((%s,%s),(%s,%s));\n((%s,%s),(%s,%s));\n((%s,%s),(%s,%s));\n", up, x, u1, u2,  u1, x, up, u2,  u2, x, u1, up);
	get_ll_from_raxml(
		master_ml_options->distance_model,
		stock_quartet_name,
		stock_msa_name,
		NULL, 
		cq, 
		&ll_p, 
		&ll_1, 
		&ll_2
	);	
	// sprintf(cq, "((%s,%s),(%s,%s));\n", u1, x, up, u2);
	// get_ll_from_raxml(&tmp_options, master_ml_options, cq, &ll_1);	

	// sprintf(cq, "((%s,%s),(%s,%s));\n", u2, x, u1, up);
	// get_ll_from_raxml(&tmp_options, master_ml_options, cq, &ll_2);
	// printf("%lf\n", ll_p);
	// while(1);


	if(ll_p - ll_1 > -EPS && ll_p - ll_2 > -EPS){
		*res = 0;
		*M = 1.0 / (1.0 + exp(1ll * (-ABS(ll_p - ll_1))) + exp(1ll * (-ABS(ll_p - ll_2))));
	} else if (ll_1 - ll_p > -EPS && ll_1 - ll_2 > -EPS){
		*res = 1;
		*M = 1.0 / (1.0 + exp(1ll * (-ABS(ll_1 - ll_p))) + exp(1ll * (-ABS(ll_1 - ll_2))));
	} else {
		*res = 2;
		*M = 1.0 / (1.0 + exp(1ll * (-ABS(ll_2 - ll_1))) + exp(1ll * (-ABS(ll_2 - ll_p))));
	}
	if(*M < 0.001) *M = 0.001;
	*M = 1.0/(*M);
	*revote_power = 1; 
	// if(*M < 0.001) *M = 0.001;
	// printf("checking , p %lf, 1 %lf, 2 %lf, res %d M %lf\n", ll_p, ll_1, ll_2, *res, *M);

	system("rm tmp.msa");

	// while(1);s
	return 0;
}
