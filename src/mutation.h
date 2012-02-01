#ifndef _MUTATION_H_
#define _MUTATION_H_

#include "common.h"
#include "hap.h"

struct mutation_category {
  int r;
  int d;
  
  int simple_nrow;
  int** simple_table;
  double* simple_not_others_probs;
    
  int extended_nrow;
  int** extended_table;
  double* extended_probs;
  double* extended_normalised_probs;
  
  double prob_sum;
};

#include "print.h"

int* get_mutated_haplotype(struct mutation_category* mut_cat, int* haplotype, int mutation_type);
struct mutation_category* create_mutation_category(int r, int d, double** mutation_prob_table);
void free_mutation_category(struct mutation_category* mut_cat);

#endif /* _MUTATION_H_ */

