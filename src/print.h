#ifndef _PRINT_H_
#define _PRINT_H_

#include "common.h"
#include "kdtree.h"

void print_h(const int* h, int r);
void print_node(const struct kdnode *node, const int r);
void print_tree(const struct kdtree *tree);
int print_list(struct kdnode* head, int r);
void print_generations(struct kdtree** generations, int g, int k, int r);

#endif /* _PRINT_H_ */

