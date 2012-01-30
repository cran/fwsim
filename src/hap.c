#include "hap.h"

int* get_haplotype_copy(int* haplotype, int r) {
  int i;
  int* new_haplotype;
  
  new_haplotype = malloc(r*sizeof(int));
  
  if (!new_haplotype) error("Could not allocate memory for new haplotype");
 
  for (i = 0; i < r; i++) {
    new_haplotype[i] = haplotype[i];
  }

  return new_haplotype;
}

