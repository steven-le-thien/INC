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
