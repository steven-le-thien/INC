// File in HMMDecompositionDecision, created by Thien Le in July 2018

#ifndef UTILITIES_H
#define UTILITIES_H

#define maxN 			4000000 // these are number of nodes, not just number of leaves
#define maxNameSize 	10
#define maxBufferSize 	1000000

#define GENERAL_BUFFER_SIZE                     10000

// Error values
#define GENERAL_ERROR                           -1
#define MALLOC_ERROR                            -2
#define OPEN_ERROR                              -3

#define IS_ERROR(val)                           val < 0

#define SUCCESS                                 0 

// String utils
#define strclr(string)                          string[0] = 0
#define strempty(string)                        string[0] == 0
#define str_start_with(string, c)               string[0] == c
#define str_add_space(string)                   strcat(string, " ")

// Print utils
#define PRINT_AND_RETURN(p, r)                  do{printf("%s\n", p); return r;}while(0) 
#define PRINT_AND_EXIT(p, r, a)                 do{printf("%s\n", p); clean_up(a); return r;}while(0) 

// Math utils
#define max(a, b) (a > b ? a : b)


#endif