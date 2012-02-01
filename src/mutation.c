#include "mutation.h"

int* get_mutated_haplotype(struct mutation_category* mut_cat, int* haplotype, int mutation_type) {
  int k, locus, dir;
  int* h = get_haplotype_copy(haplotype, mut_cat->r);
  
  if (mutation_type >= mut_cat->extended_nrow) {
    error("Illegal mutation type");
  }

  for (k = 0; k < mut_cat->d; k++) {
    locus = mut_cat->extended_table[mutation_type][k];
    dir = mut_cat->extended_table[mutation_type][mut_cat->d + k];
    h[locus] += (dir == 0) ? 1 : -1;
  }

  /*
  print_h(haplotype, mut_cat->r);
  Rprintf(" -> ");
  print_h(h, mut_cat->r);
  Rprintf("\n");
  */
  
  return(h);
}

static void gen_all_subset(int* s, int p, int k, int* t, int q, int r, int** res, int* res_idx) {
  int i;
  
  if (q == k) {
    for (i = 0; i < k; i++) {
      res[*res_idx][i] = t[i];
    }

    *res_idx += 1;
  } else {
    for (i = r; i < p; i++) {
      t[q] = s[i];
      gen_all_subset(s, p, k, t, q+1, i+1, res, res_idx);
    }
  }
}

static void calculate_simple_table(struct mutation_category* mut_cat) {
  int i, j;
  int* t;
  int* s;
  int** res;
  int res_idx = 0;
  int m = choose(mut_cat->r, mut_cat->d);
  mut_cat->simple_nrow = m;
  
  if (mut_cat->d == 0) {
    return;
  }
  
  if (mut_cat->d == mut_cat->r) {
    /*
    Special treatment for d == r
    */
    res = malloc(sizeof(int*));
    if (!res) error("generate_subsets, d == r: Could not allocate memory for res");
    res[0] = malloc(mut_cat->r*sizeof(int));
    if (!res[0]) error("generate_subsets, d == r: Could not allocate memory for res[0]");
  
    for (j = 0; j < mut_cat->r; j++) {
      res[0][j] = j;
    }
    
    mut_cat->simple_table = res;
    return;
  }
  
  s = malloc(m * sizeof(int));
  if (!s) error("generate_subsets: Could not allocate memory for s");
  for (i = 0; i < mut_cat->r; i++) s[i] = i;
  
  res = malloc(m * sizeof(int*));
  if (!res) error("generate_subsets: Could not allocate memory for res");
  for (i = 0; i < m; i++) {
    res[i] = malloc(mut_cat->d * sizeof(int));
    if (!res[i]) error("generate_subsets: Could not allocate memory for res[i]");
  }
  
  t = malloc(m * sizeof(int));
  if (!t) error("generate_subsets: Could not allocate memory for t");
    
  gen_all_subset(s, mut_cat->r, mut_cat->d, t, 0, 0, res, &res_idx);
  
  free(s);
  free(t);
  
  mut_cat->simple_table = res;
  
//  return res;  
}

static void calculate_simple_not_others_probs(struct mutation_category* mut_cat, double** mutation_prob_table) {
  int i, j, k;
  double* probs_not_others;
  double prob_not;
  int* loci_mut;
  int include_locus;
  
  probs_not_others = malloc(mut_cat->simple_nrow * sizeof(double));
  if (!probs_not_others) error("calculate_simple_probs: Could not allocate memory for probs_not_others");
  
  for (i = 0; i < mut_cat->simple_nrow; i++) {
    loci_mut = mut_cat->simple_table[i];
    prob_not = 1;
    
    for (j = 0; j < mut_cat->r; j++) {
      include_locus = 1;
      
      for (k = 0; k < mut_cat->d; k++) {
        if (j == loci_mut[k]) {
          include_locus = 0;
          break;
        }
      }
      
      if (include_locus == 1) {
        prob_not *= 1 - mutation_prob_table[j][2];
      }
    }
    
    probs_not_others[i] = prob_not;
  }
  
  mut_cat->simple_not_others_probs = probs_not_others;
}

static void calculate_extended_table(struct mutation_category* mut_cat, double** mutation_prob_table) {
  int i, j, k, n = 0;
  int** res;
  double* probs;
  double* probs_normalised;
  double prob_sum = 0;  
  int combs = pow(2, mut_cat->d);
  double prob;
  int dir, locus;
  
  mut_cat->extended_nrow = mut_cat->simple_nrow * combs;
  
  res = malloc(mut_cat->extended_nrow * sizeof(int*));
  if (!res) error("calculate_extended_table: Could not allocate memory for res");
  for (i = 0; i < mut_cat->extended_nrow; i++) {
    res[i] = malloc(2 * mut_cat->d * sizeof(int));
    if (!res[i]) error("calculate_extended_table: Could not allocate memory for res[i]");
  }

  probs = malloc(mut_cat->extended_nrow * sizeof(double));
  if (!probs) error("calculate_extended_table: Could not allocate memory for probs");

  probs_normalised = malloc(mut_cat->extended_nrow * sizeof(double));
  if (!probs_normalised) error("calculate_extended_table: Could not allocate memory for probs_normalised");  
  
  for (i = 0; i < mut_cat->simple_nrow; i++) { 
    for (j = 0; j < combs; j++) {
      prob = mut_cat->simple_not_others_probs[i];
    
      for (k = 0; k < mut_cat->d; k++) {
        dir = IS_BIT_SET(j, k); /* dir == 0 => +, dir == 1 => -*/
        locus = mut_cat->simple_table[i][k];
        
        res[n][k] = locus;
        res[n][mut_cat->d + k] = dir;
        prob *= mutation_prob_table[locus][dir];
        
        #ifdef EXTRA_VERBOSE
        Rprintf("%d [%c] ", res[n][k], res[n][mut_cat->d + k] == 0 ? '+' : '-');
        #endif
      }
      
      #ifdef EXTRA_VERBOSE
      Rprintf(" => Prob = P(not mut others)");
      
      for (k = 0; k < mut_cat->d; k++) {
        dir = IS_BIT_SET(j, k); /* dir == 0 => +, dir == 1 => -*/
        locus = mut_cat->simple_table[i][k];
        Rprintf(" * mu[%d][%d]", locus, dir);
      }
      
      Rprintf(" = %f ", mut_cat->simple_not_others_probs[i]);  
      
      for (k = 0; k < mut_cat->d; k++) {
        dir = IS_BIT_SET(j, k); /* dir == 0 => +, dir == 1 => -*/
        locus = mut_cat->simple_table[i][k];
        Rprintf(" * %f", mutation_prob_table[locus][dir]);
      }
      
      Rprintf("%f", prob);      
      Rprintf("\n");
      #endif
    
      probs[n] = prob;
      prob_sum += prob;
      n++;
    }
  }
  
  for (i = 0; i < mut_cat->extended_nrow; i++) {
    probs_normalised[i] = probs[i] / prob_sum;
  }
  
  mut_cat->extended_table = res;
  mut_cat->extended_probs = probs;
  mut_cat->extended_normalised_probs = probs_normalised;
  mut_cat->prob_sum = prob_sum;
    
  /*  
  Rprintf("extended table:\n");

  for (i = 0; i < mut_cat->simple_nrow; i++) {
    for (j = 0; j < mut_cat->d; j++) {
      for (dir = 0; dir <= 1; dir++) {
        Rprintf("%d%c ", mut_cat->simple_table[i][j], (dir == 0) ? '-' : '+');
        c++;
      }
    }
    
    Rprintf("\n");
  }

  Rprintf("extended nrow = %d = %d\n", mut_cat->extended_nrow, c);
  */
}

static void fill_no_mutation_category(struct mutation_category* mut_cat, double** mutation_prob_table) {
  int i;
  double prob = 1;
  
  for (i = 0; i < mut_cat->r; i++) {
    prob *= 1 - mutation_prob_table[i][2];
  }
  
  mut_cat->prob_sum = prob;
}

struct mutation_category* create_mutation_category(int r, int d, double** mutation_prob_table) {
  struct mutation_category* mut_cat;
  int n = 0;
  
  mut_cat = malloc(sizeof(struct mutation_category));
  if (!mut_cat) error("Could not allocate memory for mutation category");

  mut_cat->r = r;
  mut_cat->d = d;
  mut_cat->simple_nrow = 0;  
  mut_cat->simple_table = 0;
  mut_cat->simple_not_others_probs = 0;
  mut_cat->extended_nrow = 0;
  mut_cat->extended_table = 0;
  mut_cat->extended_probs = 0;
  mut_cat->extended_normalised_probs = 0;
  mut_cat->prob_sum = 0;
  
  if (d == 0) {
    fill_no_mutation_category(mut_cat, mutation_prob_table);
  } else {     
    calculate_simple_table(mut_cat);
    calculate_simple_not_others_probs(mut_cat, mutation_prob_table);
    calculate_extended_table(mut_cat, mutation_prob_table);
  }
  
  return mut_cat;
}

void free_mutation_category(struct mutation_category* mut_cat) {
  if (!mut_cat) {
    return;
  }
  
  int i, j;
  int d = mut_cat->d;
  
  for (i = 0; i < mut_cat->simple_nrow; i++) {
    free(mut_cat->simple_table[i]);
  }
  
  free(mut_cat->simple_table);
  
  free(mut_cat);  
}


















static int* build_1mut_table(int r, int* n, double** prob) {
  int i, idx;
  int* tab;
  double* prob_tab;
  
  *n = r;
  
  tab = malloc(r*sizeof(int));
  if (!tab) error("build_1mut_table: Could not allocate memory for tab");

  prob_tab = malloc(r*sizeof(double));
  if (!prob_tab) error("build_1mut_table: Could not allocate memory for prob_tab");
  
  for (i = 0; i < r; i++) {
    prob_tab[i] = 1.0 / (double)r;
    tab[i] = i;
  }
  
  *prob = prob_tab;  

  return tab;
}

static int** build_2mut_table(int r, int* n, double** prob) {
  int i, j, m, idx;
  int** tab;
  double* prob_tab;
  
  if (r < 2) {
    r = 2;
  }
  
  m = choose(r, 2);
  *n = m;
  
  tab = malloc(m*sizeof(int*));
  if (!tab) error("build_2mut_table: Could not allocate memory for tab");

  prob_tab = malloc(m*sizeof(double));
  if (!prob_tab) error("build_2mut_table: Could not allocate memory for prob_tab");
  
  for (i = 0; i < m; i++) {
    prob_tab[i] = 1.0 / (double)m;
  }
  
  *prob = prob_tab;
  
  idx = 0;
  
  for (i = 0; i < r - 1; i++) {
    for (j = i+1; j < r; j++) {
      tab[idx] = malloc(2*sizeof(int));
      if (!tab[idx]) error("build_2mut_table: Could not allocate memory for tab[k]");
      
      tab[idx][0] = i;
      tab[idx][1] = j;

      idx++;
    }    
  }

  return tab;
}

static int** build_3mut_table(int r, int* n, double** prob) {
  int i, j, k, m, idx;
  int** tab;  
  double* prob_tab;
  
  if (r < 3) {
    r = 3;
  }
  
  m = choose(r, 3);
  *n = m;
  
  tab = malloc(m*sizeof(int*));
  if (!tab) error("build_3mut_table: Could not allocate memory for tab");
  
  prob_tab = malloc(m*sizeof(double));
  if (!prob_tab) error("build_3mut_table: Could not allocate memory for prob_tab");
  
  for (i = 0; i < m; i++) {
    prob_tab[i] = 1.0 / (double)m;
  }
  
  *prob = prob_tab;
  
  idx = 0;
  
  for (i = 0; i < r - 2; i++) {
    for (j = i+1; j < r - 1; j++) {
      for (k = j+1; k < r; k++) {
        tab[idx] = malloc(3*sizeof(int));
        if (!tab[idx]) error("build_3mut_table: Could not allocate memory for tab[k]");
        
        tab[idx][0] = i;
        tab[idx][1] = j;
        tab[idx][2] = k;
        
        idx++;
      }
    }    
  }

  return tab;
}


