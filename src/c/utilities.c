#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utilities.h"

// Jim Baiter's 
void* safe_malloc(size_t n)
{
  void* p = malloc(n);
  if (!p)
  {
    fprintf(stderr, "Out of memory(%lu bytes)\n", (unsigned long) n);
    exit(EXIT_FAILURE);
  }
  return p;
}

FILE* safe_open_read(char * name)
{
  if(!name) return NULL;
  FILE * f = fopen(name, "r");
  if(!f) {
    fprintf(stderr, "Open error: %s\n", name);
    exit(EXIT_FAILURE);
  }
  else return f;
}

FILE* safe_reopen_read(char * name, FILE * stream)
{
  if(!name) return NULL;
  FILE * f = freopen(name, "r", stream);
  if(!f) {
    fprintf(stderr, "ReOpen error: %s \n", name);
    exit(EXIT_FAILURE);
  }
  else return f;
}


