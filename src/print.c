#include "print.h"

void print_mu(double** mutation_prob_table, int r) {
  int i;

  Rprintf("        ");
  for (i = 0; i < r; i++) {
    Rprintf("L%-5d ", (i+1));
  }
  
  Rprintf("\n");
    
  Rprintf("  Up:   ");
  for (i = 0; i < r; i++) Rprintf("%6.4f ", mutation_prob_table[i][0]); 
  Rprintf("\n");

  Rprintf("  Down: ");
  for (i = 0; i < r; i++) Rprintf("%6.4f ", mutation_prob_table[i][1]); 
  Rprintf("\n");
  
  Rprintf("  Sum:  ");
  for (i = 0; i < r; i++) Rprintf("%6.4f ", mutation_prob_table[i][0]+mutation_prob_table[i][1]); 
  Rprintf("\n");
}

void print_alpha(const double* alpha, int g) {
  int i;
  int printed = 0;
  int skipped = 0;
  double last_alpha = -1;

  Rprintf("(");

  for (i = 0; i < g; i++) {
    if (alpha[i] == last_alpha) {
      skipped++;
      continue;
    }

    if (skipped == 0 && i > 0) Rprintf(", ");    
    
    last_alpha = alpha[i];

    printed++;
    
    if (skipped > 0) {
      Rprintf(" x %d, ", skipped + 1);
      printed++;
      skipped = 0;
    }
    
    Rprintf("%.3f", alpha[i]);
  }

  if (skipped > 0) {
    Rprintf(" x %d", skipped + 1);
  }
    
  Rprintf(")");
}

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

void print_mutation_category(struct mutation_category* mut_cat) {
  int i, j, k;
  
  Rprintf("  %-2d mutations has probability %e\n", mut_cat->d, mut_cat->prob_sum);
  
  if (!mut_cat->simple_nrow) {
    return;
  }
  
  Rprintf("    Number of possible locus combinations    = %d\n", mut_cat->simple_nrow);
  Rprintf("    Number of possible mutation combinations = %d\n", mut_cat->extended_nrow);
  
  #ifdef VERBOSE
  for (i = 0; i < mut_cat->simple_nrow; i++) {
    for (j = 0; j < mut_cat->d; j++) {
      Rprintf("%d ", mut_cat->simple_table[i][j]);
    }
    Rprintf(" P(not mut other loci) = %f\n", mut_cat->simple_not_others_probs[i]);
  }  
  #endif

  #ifdef VERBOSE
  Rprintf("extended table:\n");

  for (i = 0; i < mut_cat->extended_nrow; i++) {
    for (j = 0; j < mut_cat->d; j++) {
      Rprintf("%d [%c] ", mut_cat->extended_table[i][j], mut_cat->extended_table[i][mut_cat->d + j] == 0 ? '+' : '-');
    }
    Rprintf(" P(row) = %f, norm prob = %f\n", mut_cat->extended_probs[i], mut_cat->extended_normalised_probs[i]);
  }
  #endif
}

