#ifndef _SIM_H_
#define _SIM_H_

#include "common.h"
#include "sim.h"
#include "rand.h"
#include "print.h"
#include "hap.h"
#include "kdtree.h"
#include "tree-to-list.h"
#include "clean.h"

int population_to_matrix(struct kdnode* head, int r, int*** pop, int* nrowptr);
struct kdtree** simulate_generations(int g, int k, int r, double alpha, double mu, int* pop_sizes, int trace, int trace_loc_mut);

#endif /* _SIM_H_ */

