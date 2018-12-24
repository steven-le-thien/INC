// File in inc_ml, created by Thien Le in July 2018

#ifndef UTILITIES_H
#define UTILITIES_H


#define DEBUG                                   0
#define LARGE_DEBUG                             0
#define DEBUG_REC                               0
#define DEBUG_BFS                               0

#define REVOTE_POWER                            2

// Limits 
#define MAXN                                    4000000 // these are number of nodes, not just number of leaves
#define MAX_NAME_SIZE                           1000
#define MAX_BUFFER_SIZE                         10000
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
                                                    if(i % (n / 10) == 0){\
                                                        for(j = 0; j < (i <= n / 10 ? 0 : (int)log10(i - n / 10)) + 4; j++) \
                                                           printf("\b");\
                                                        printf("%d...", i);\
                                                        fflush(stdout);\
                                                    }\
                                                } while(0)


// Math utils
#define MAX(a, b)                               ((a) > (b) ? (a) : (b))
#define MIN(a, b)                               ((a) < (b) ? (a) : (b))

#ifndef INT_MAX
#define INT_MAX                                 10000000
#endif 

#define LN2                                     1.4426950408
#define EPS                                     1e-7              
#define POWER(a, b)                             ((b) == 0 ? 1 : ((b) == 1 ? (a) : ((b) == 2 ? (a) * (a) : ((b) == 3 ? (a) * (a) * (a) : ((b) == 4 ? (a) * (a) * (a) * (a) : -1)))))
#define ABS(a)                                  ((a) > 0 ? (a) : -(a))
#define LOG2(a)                                 (log(a) / log(2))


// Options util
#define RECIP_WEIGHT                            1 
#define SQUARE_RECIP_WEIGHT                     2


#endif
