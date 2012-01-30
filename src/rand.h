#ifndef _RAND_H_
#define _RAND_H_

#include "common.h"

double random_poisson(double mean);
int random_locus(int r);
int random_mutation();
int* random_loci(int r, int n);

#endif /* _RAND_H_ */

