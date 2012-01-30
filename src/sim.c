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
  
  matrix = malloc(nrow*sizeof(int*));
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

/*
  r:       number of loci
  k:       number of fathers
*/
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

static struct kdtree* simulate_generation(struct kdtree* fathers, int k, int r, double alpha, double mu, int* popsize, int trace_loc_mut) {
  int size = 0;
  int res, l, j;
  int* new_haplotype;
  struct kdnode* fathers_list;
  struct kdnode* current;
  struct kdtree* children;
  double mj;
  int* nj;
  
  fathers_list = tree_to_list(fathers, r);
  children = kd_create(r);

  current = fathers_list;
  
  while (current != NULL) {
    mj = (double)current->count;
    nj = malloc((r+1)*sizeof(int));

    for (l = 0; l <= r; l++) {
      nj[l] = (int)random_poisson(dbinom(l, r, mu, 0)*alpha*mj);
      //printf("\t%d muts: %5d haplotypes\n", l, nj[l]);
    }

    // No mutations
    if (nj[0] > 0) {
      size += nj[0];
      res = kd_insert_or_update_count(children, current->pos, nj[0]);
      if (res != 0) error("Could not insert new node in kd-tree");
    }
    
    for (l = 1; l <= r; l++) {
      if (nj[l] <= 0) {
        continue;
      }
      
      size += nj[l];
      
      int* loci = random_loci(r, l);      
      int* new_haplotype = get_haplotype_copy(current->pos, r);
      
      for (j = 0; j < r; j++) {
        if (loci[j] == 1) {
          new_haplotype[j] += random_mutation();
        }
      }
      
      res = kd_insert_or_update_count(children, new_haplotype, nj[l]);
      if (res != 0) error("Could not insert new node in kd-tree");

      if (l >= trace_loc_mut) {
        Rprintf("\t%d mutations occured from ", l);
        print_h(current->pos, r);
        Rprintf(" to ");
        print_h(new_haplotype, r);
        Rprintf("\n", l);
      }
                
      free(loci);
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

/*
  g:        number of generations (results in g+1 because initial is included)
  r:        number of loci
  k:        number of fathers
  alpha:    growth rate
  mu:       mutation rate
*/
struct kdtree** simulate_generations(int g, int k, int r, double alpha, double mu, int* pop_sizes, int trace, int trace_loc_mut) {
  int size = 0;
  int i;
  struct kdtree** generations;
  struct kdtree* fathers;
  struct kdtree* children;
  
  generations = malloc((g+1)*sizeof(struct kdtree*));
  if (!generations) error("Could not allocate memory for generations");
  
  fathers = generate_initial_population(k, r);  
  
  if (trace) {
    Rprintf("\nInitial population of %d created.\nContinuing evolution:\n", k);
  }
    
  generations[0] = fathers;
  pop_sizes[0] = k;
  
  for (i = 1; i <= g; i++) {
    fathers = generations[i-1];
    children = simulate_generation(fathers, k, r, alpha, mu, &size, trace_loc_mut);
    generations[i] = children;
    pop_sizes[i] = size;

    /* FIXME */
    free_tree(fathers);
    
    if (trace) {
      Rprintf("%5d/%d done (population size = %10ld)\n", i, g, size);
    }
  }
  
  /* Up to g-1 is done in loop, but not this */
  generations[g]->root = tree_to_list(generations[g], r);
  
  return generations;
}

