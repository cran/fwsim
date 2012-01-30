#include "print.h"

void print_h(const int* h, int r) {
  int i;

  Rprintf("(");
      
  for (i = 0; i < r; i++) {
    if (i > 0) Rprintf(", ");
    Rprintf("%02d", h[i]);
  }
  
  Rprintf(")");
}

void print_node(const struct kdnode *node, const int r) {  
  if (!node) {
    return;
  }
  
  print_h(node->pos, r);
  Rprintf(" = %d\n", node->count);
  
  print_node(node->left, r);
  print_node(node->right, r);
}

void print_tree(const struct kdtree *tree) {
  if (tree->root) {
    print_node(tree->root, tree->dim);
  } else {
    Rprintf("(EMPTY TREE)\n");
  }
}

int print_list(struct kdnode* head, int r) {
  struct kdnode* current;
  int count;
  current = head;
  
  count = 0;
  
  while (current != NULL) {
    count += current->count;
    Rprintf("\t");
    print_h(current->pos, r);
    Rprintf(" = %d\n", current->count);
    current = current->right;
    if (current == head) break;
  }
  
  return count;
}

void print_generations(struct kdtree** generations, int g, int k, int r) {
  int i;
  int count;
  
  for (i = 0; i <= g; i++) {  
    Rprintf("Generation %4d:\n", i);
    count = print_list(generations[i]->root, r);
    Rprintf("\tPopulation size = %d\n\n", count);
  }
}

