// File in inc_ml, created by Thien Le in July 2018

#ifndef UTILITIES_H
#define UTILITIES_H

#define DEBUG       0

#define MAXN                                    4000000 // these are number of nodes, not just number of leaves
#define MAX_NAME_SIZE                           10
#define MAX_BUFFER_SIZE                         1000000

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
#define cat_c_to_a(c, a)                        do{\
                                                    char buf[2];\
                                                    buf[0] = c;\
                                                    buf[1] = 0;\
                                                    strcat(a, buf);\
                                                } while(0)

// Print utils
#define PRINT_AND_RETURN(p, r)                  do{printf("%s\n", p); return r;}while(0) 
#define PRINT_AND_EXIT(p, r)                    do{printf("%s\n", p); return r;}while(0) 
#define print_inline_iteration(i, j, n, s)      do{\
                                                    if(i % ((n) / 10) == 0 || i == s){ \
                                                        for(j = 0; j < (int)log10(i - 1) + 4; j++) \
                                                           printf("\b");\
                                                        printf("%d...", i);\
                                                        fflush(stdout);\
                                                    }\
                                                } while(0)


// Math utils
#define max(a, b) (a > b ? a : b)
#define INT_MAX                                 10000000
#define LN2                                     1.4426950408
#define EPS                                     1e-7                                           


#endif