#include "sim.h"

int population_to_matrix(struct kdnode* head, int r, int*** pop, int* nrowptr) {
  int** matrix;
  int nrow = 0;
  int i, j;
  struct kdnode* current;
  
  current = head;
  while (current != NULL) {
    nrow++;
    current = current->right;
    if (current == head) break;
  }
  
  matrix = malloc(nrow * sizeof(int*));
  if (!matrix) error("Could not allocate memory for population matrix");

  i = 0;
  current = head;
  while (current != NULL) {
    if (i > nrow) error("Some counting went wrong in filling population matrix");
    
    matrix[i] = malloc((r+1)*sizeof(int));
    
    for (j = 0; j < r; j++) {
      matrix[i][j] = current->pos[j];
    }
    
    matrix[i][r] = current->count;
    
    current = current->right;
    i++;
    
    if (current == head) break;
  }
    
  *nrowptr = nrow;
  *pop = matrix;
  
  return 0;
}

static struct kdtree* generate_initial_population(int k, int r) {
  int i;
  int res;
  int* h;
  struct kdtree *tree;
  tree = kd_create(r);
  
  h = malloc(r*sizeof(int));
  if (!h) error("Could not allocate memory for new haplotype");
  
  for (i = 0; i < r; i++) {
    h[i] = 0;
  }
  
  res = kd_insert(tree, h, k);
  if (res != 0) error("Could not insert new node in kd-tree");
  
  /*
  Rprintf("Initial pop generated\n");
  */
  
  return tree;
}

static struct kdtree* simulate_generation(struct kdtree* fathers, int k, int r, 
  struct mutation_category** mut_cats, double alpha, int* popsize, int trace) 
{
  int size = 0;
  int res, l, j;
  int* new_haplotype;
  struct kdnode* fathers_list;
  struct kdnode* current;
  struct kdtree* children;
  struct mutation_category* mut_cat;
  double mj;
  int* nj;
  int* rN;
  
  fathers_list = tree_to_list(fathers, r);
  children = kd_create(r);
 
  current = fathers_list;
  
  while (current != NULL) {
    mj = (double)current->count;
    nj = malloc((r+1)*sizeof(int));
    if (!nj) error("Could not allocate memory for nj");

    /*
    if (trace == 1) {
      Rprintf("  Haplotype ");
      print_h(current->pos, r);
      Rprintf(":\n");
    }
    */
    
    for (l = 0; l <= r; l++) {
      nj[l] = (int)random_poisson(mut_cats[l]->prob_sum*alpha*mj);
      
      /*
      if (trace == 1 && nj[l] > 0) {
        Rprintf("    %d muts: %5d haplotypes\n", l, nj[l]);
      }
      */
    }

    /*
    No mutations
    */
    if (nj[0] > 0) {
      size += nj[0];
      res = kd_insert_or_update_count(children, current->pos, nj[0]);
      if (res != 0) error("Could not insert new node in kd-tree");
    }
    
    /*
    Multinomial mutation:
      void rmultinom(int n, double* prob, int K, int* rN)
      Return vector rN[1:K] {K := length(prob)}
    Binomial:
      double rbinom(double nin, double pp)
    */
    for (l = 1; l <= r; l++) {
      if (nj[l] <= 0) {
        continue;
      }
      
      size += nj[l];
      mut_cat = mut_cats[l];
      
      rN = malloc(mut_cat->extended_nrow*sizeof(int));
      if (!rN) error("Could not allocate memory for rN");
      
      rmultinom(nj[l], mut_cat->extended_normalised_probs, mut_cat->extended_nrow, rN);
      
      /*
      Rprintf("  %2d mutations:\n", mut_cat->d);
      */
      
      for (j = 0; j < mut_cat->extended_nrow; j++) {
        if (rN[j] == 0) {
          continue;
        }
        
        /*
        Rprintf("    %d of index %d\n", rN[j], j);
        */
        int* new_hap = get_mutated_haplotype(mut_cat, current->pos, j);
        res = kd_insert_or_update_count(children, new_hap, rN[j]);
        if (res != 0) error("Could not insert new node in kd-tree");
      }

      free(rN);
    }
    
    free(nj);
    
    current = current->right;

    if (current == fathers_list) {
      break;
    }
  }
  
  *popsize = size;
  
  return children;
}

struct kdtree** simulate_generations(int g, int k, int r, double* alpha, double** mutation_prob_table, int* pop_sizes, int trace) {
  int size = 0;
  int i;
  struct kdtree** generations;
  struct kdtree* fathers;
  struct kdtree* children;
  struct mutation_category** mut_cats;
  int mut_cats_threshold = (MUTATION_CATEGORY_THRESHOLD > 0 && MUTATION_CATEGORY_THRESHOLD < r) ? MUTATION_CATEGORY_THRESHOLD : r;
  double total_mut_prob_sum = 0;

  #ifdef VERBOSE
  Rprintf("\nMutation categories threshold = %d\n", mut_cats_threshold);
  #endif
  
  /*
  if (trace) {
    Rprintf("Mutation parameters in Poisson(f(d; r, mu)*alpha*N):\n");
  }  
  */
  generations = malloc((g+1)*sizeof(struct kdtree*));
  if (!generations) error("Could not allocate memory for generations");
    
  mut_cats = malloc((mut_cats_threshold + 1)*sizeof(struct mutation_category*));
  if (!mut_cats) error("Could not allocate memory for mut_cats");
  
  for (i = 0; i <= mut_cats_threshold; i++) {
    mut_cats[i] = malloc(sizeof(struct mutation_category));
    if (!mut_cats[i]) error("Could not allocate memory for mut_cats[i]");
    mut_cats[i] = create_mutation_category(r, i, mutation_prob_table);
    total_mut_prob_sum += mut_cats[i]->prob_sum;
  }
  
  if (trace == 1) {
    Rprintf("Mutation categories:\n");
    
    for (i = 0; i <= mut_cats_threshold; i++) { 
      print_mutation_category(mut_cats[i]);
    }
  
    Rprintf("  Sum of mutation categories probabilities: %f\n", total_mut_prob_sum);    
  }

  if (fabs(total_mut_prob_sum - 1.0) > 1.0e-6) {
    error("Expected sum of mutation probabilities to be 1");
  }
  
  fathers = generate_initial_population(k, r);  
  
  if (trace) {
    Rprintf("\nInitial population of %d created.\nContinuing evolution:\n", k);
  }
    
  generations[0] = fathers;
  pop_sizes[0] = k;
  
  for (i = 1; i <= g; i++) {
    /*
    Rprintf("Evolving generation %5d\n", i);
    */
    fathers = generations[i-1];

    children = simulate_generation(fathers, k, r, mut_cats, alpha[i-1], &size, trace);

    generations[i] = children;
    pop_sizes[i] = size;

    free_tree(fathers);
    
    if (trace) {
      Rprintf("%5d/%d: done (growth rate = %.3f, population size = %9d, haplotypes = %7d)\n", i, g, alpha[i-1], size, children->size);
    }
  }
  
  /* Up to g-1 is done in loop, but not this */
  generations[g]->root = tree_to_list(generations[g], r);
  
  for (i = 0; i <= mut_cats_threshold; i++) {
    free_mutation_category(mut_cats[i]);
  }
  
  free(mut_cats);  

  return generations;
}

