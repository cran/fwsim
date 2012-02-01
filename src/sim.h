#ifndef _SIM_H_
#define _SIM_H_

#include "common.h"
#include "rand.h"
#include "print.h"
#include "hap.h"
#include "kdtree.h"
#include "tree-to-list.h"
#include "clean.h"
#include "mutation.h"

int population_to_matrix(struct kdnode* head, int r, int*** pop, int* nrowptr);
struct kdtree** simulate_generations(int g, int k, int r, double* alpha, double** mutation_prob_table, int* pop_sizes, int trace);

#endif /* _SIM_H_ */

