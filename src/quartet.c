// File in inc_ml, created by Thien Le in July 2018
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
		*res = 1; // first choild is sibling
	else 
		*res = 2; // second child is sibling

	return 0;
}

int induced_quartets(char * up, char * u1, char * u2, char * x, int * res, char * tree_name){
	char sibling1[MAX_NAME_SIZE];
	char sibling2[MAX_NAME_SIZE];
	char cur_char;
	int counter;

	FILE * f;

	// Make label
	f = fopen("tmp.label", "r"); 
	if(!f){
		// we good
		f = fopen("tmp.label", "w");
	} else {
		fclose(f);
		PRINT_AND_RETURN("tmp.label already exits in working directory\n", GENERAL_ERROR);
	}
																							// #if 1 
 																						// 		printf("%s\n%s\n%s\n%s\n", up, u1, u2, x);
 																						// 		while(1);
																							// #endif
	fprintf(f, "%s\n%s\n%s\n%s\n", up, u1, u2, x);
	fclose(f);

	// Make subtree
	f = fopen("tmp.quartet", "r"); 
	if(f){
		fclose(f);
		PRINT_AND_RETURN("tmp.quarte already exits in working directory\n", GENERAL_ERROR);
	}
	make_subtree("tmp.label", "tmp.quartet", tree_name);
	// Read subtree
	f = fopen("tmp.quartet", "r");  // the tree looks like ((A, B), (C, D));
	counter = 0;
	fscanf(f, "%c", &cur_char); // skip first 2 characters 
	cur_char = ' ';
	while(cur_char != '(') // look for second open bracket 
		fscanf(f, "%c", &cur_char);

	fscanf(f, "%c", &cur_char);
	while(cur_char != ','){
		sibling1[counter++] = cur_char;
		fscanf(f, "%c", &cur_char);
	}
	sibling1[counter] = 0; // null terminal

	counter = 0;
	fscanf(f, "%c", &cur_char);
	while(cur_char != ')'){
		sibling2[counter++] = cur_char;
		fscanf(f, "%c", &cur_char);
	}
	sibling2[counter] = 0; // null terminal

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

	// Remmove tmp files
	system("rm tmp.label tmp.quartet");

	return 0;
} 