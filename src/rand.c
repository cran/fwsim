#include "rand.h"

double random_poisson(double mean) {
  return rpois(mean);
}

int random_locus(int r) {
  return (int)((double)r*runif(0.0, 1.0));
}

int random_mutation() {
  return (runif(0.0, 1.0) < 0.5) ? -1 : 1;
}

/*
  Returns n of r loci marked, each with probability 1/r
*/
int* random_loci(int r, int n) {
  int* loci = malloc(r*sizeof(int));
  int i, chosen = 0;
  double comp = 1.0 / (double)r;
  
  for (i = 0; i < r; i++) {
    loci[i] = 0;
  }
  
  if (n <= 0) {
    return loci;
  }
  
  if (n >= r) {
    for (i = 0; i < r; i++) {
      loci[i] = 1;
    }
    
    return loci;
  }
  
  while (chosen < n) {    
    i = random_locus(r);
    
    while (loci[i] == 1) {
      i = random_locus(r);
    }
    
    loci[i] = 1;
    chosen++;
  };
  
  return loci;  
}

