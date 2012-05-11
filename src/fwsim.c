#include "fwsim.h"

#define IMATRIX(m, i, j) (INTEGER(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])

static double** generate_mutation_prob_table(int r, double* mu_par) {
  int i, j;
  double** mutation_prob_table;
    
  /*
    mutation_prob_table[i][0]: upwards mutation rate
    mutation_prob_table[i][1]: downwards mutation rate
    mutation_prob_table[i][2]: summed mutation rate
  */

  mutation_prob_table = malloc(r*sizeof(double*));
  if (!mutation_prob_table) error("Could not allocate memory for mutation_prop_table");
  for (i = 0; i < r; i++) {
    mutation_prob_table[i] = malloc(3*sizeof(double));
  }

  for (i = 0; i < 2; i++) {
    for (j = 0; j < r; j++) {
      double mu = mu_par[j*2 + i];
      if (!(mu >= 0 && mu <= 1)) error("Mutation rate: 0 <= mu <= 1 required");
      mutation_prob_table[j][i] = mu;
    }    
  }
  
  for (i = 0; i < r; i++) {
    mutation_prob_table[i][2] = mutation_prob_table[i][0] + mutation_prob_table[i][1];
  }
  
  return mutation_prob_table;  
}

static double expected_target_population_size(const double* alpha, const int k, const int g) {
  int i;
  double size = (double)k;
  
  for (i = 0; i < g; i++) {
    size = size*alpha[i];
  }
  
  return size;
}

static double* expected_population_sizes(const double* alpha, const int k, const int g) {
  int i;
  double* sizes = malloc((g+1)*sizeof(double));
  
  if (!sizes) error("Could not allocate memory for expected populuation sizes");
    
  sizes[0] = (double)k;
  
  for (i = 1; i <= g; i++) {
    sizes[i] = sizes[i-1]*alpha[i-1];
  }
  
  return sizes;
}

SEXP fwsim(SEXP param_g, SEXP param_k, SEXP param_r, 
  SEXP param_alpha, SEXP param_mu, 
  SEXP param_save_gs, 
  SEXP param_trace,
  SEXP param_alpha_length)
{
  int i, j, i2;
  
  struct kdtree* population;
  struct kdtree** generations;
  int* pop_sizes;
  double* expected_pop_sizes;
  int** pop;
  int nrow = 0;
  double** mutation_prob_table;

  SEXP R_pop;
  SEXP R_pop_intermediate;
  SEXP R_pop_intermediate_tmp;
  SEXP R_pop_sizes;
  SEXP R_expected_pop_sizes;
  SEXP R_list_names;
  SEXP R_ans;
  
  int g = INTEGER(param_g)[0];
  int k = INTEGER(param_k)[0];
  int r = INTEGER(param_r)[0];
  double* mu_par = REAL(param_mu);
  int* save_gs = INTEGER(param_save_gs);
  int trace = INTEGER(param_trace)[0];
  int alpha_length = INTEGER(param_alpha_length)[0];
    
  if (!(g >= 1)) error("Number of generations: g >= 1 required");
  if (!(k >= 1)) error("Size of initial population: k >= 1 required");
  if (!(r >= 1)) error("Number of loci: r >= 1 required");
  if (!(alpha_length == 1 || alpha_length == g)) error("Growth rate must be of length 1 or g");

  //double alpha = REAL(param_alpha)[0];
  double* alpha = malloc(g*sizeof(double));
  
  if (alpha_length == 1) {
    if (!(REAL(param_alpha)[0] > 0)) error("Growth rate: alpha > 0 required");
    
    for (i = 0; i < g; i++) {
      alpha[i] = REAL(param_alpha)[0];
    }
  } else {
    for (i = 0; i < g; i++) {
      if (!(REAL(param_alpha)[i] > 0)) error("Growth rate: alpha > 0 required");
      
      alpha[i] = REAL(param_alpha)[i];
    }
  }
  
  mutation_prob_table = generate_mutation_prob_table(r, mu_par);
  
  if (trace) {
    Rprintf("Number of generations:                        %d\n", g);
    Rprintf("Size of initial population:                   %d\n", k);
    Rprintf("Number of loci:                               %d\n", r);
    Rprintf("Mutation rates:\n");
    print_mu(mutation_prob_table, r);
    
    if (alpha_length == 1) {
      Rprintf("Growth rate:                                  %f\n", alpha[0]);
    } else {
      Rprintf("Growth rate:                                  ");
      print_alpha(alpha, g);
      Rprintf("\n");
    }

    Rprintf("Generations to record besides end population: ");
    
    if (do_save_gs(save_gs)) {
      for (i = 1; i < g; i++) {
        if (save_gs[i] == 1) Rprintf("%d ", i);
      }
      Rprintf("\n");
    } else {
      Rprintf("NONE\n");
    }
    
    Rprintf("Expected end population size:                 %.2f\n", 
      expected_target_population_size(alpha, k, g));
  }
  
  pop_sizes = malloc((g+1)*sizeof(int));
  if (!pop_sizes) error("Could not allocate memory for populuation sizes");
  
  GetRNGstate();
  generations = simulate_generations(g, k, r, alpha, mutation_prob_table, pop_sizes, save_gs, trace);
  if (!generations) error("Generations could not be simulated");
  PutRNGstate();
  
  for (i = 0; i < r; i++) free(mutation_prob_table[i]);
  free(mutation_prob_table);
  
  population = generations[g];  
 
  if (population == NULL) {
    error("Did not expect the population to be NULL");
  }
    
  if (population_to_matrix(population->root, r, &pop, &nrow) != 0) {
    error("Could not convert population matrix to R object");
  }

  if (nrow == 0) {
    R_pop = R_NilValue;
  } else {
    PROTECT(R_pop = allocMatrix(INTSXP, nrow, r + 1));
    for (i = 0; i < nrow; i++) {
     for (j = 0; j < r + 1; j++) {
       IMATRIX(R_pop, i, j) = pop[i][j];
     }
     free(pop[i]);
    }
    UNPROTECT(1);

    //free_tree(population);
    //generations[g] = NULL;
  }

  free(pop);


  if (do_save_gs(save_gs)) {
    PROTECT(R_pop_intermediate = allocVector(VECSXP, g-1));
    for (i2 = 1; i2 < g; i2++) {
      if (save_gs[i2] != 1) {
        SET_VECTOR_ELT(R_pop_intermediate, i2, R_NilValue);
        continue;
      }
    
      population = generations[i2];
      
      if (population == NULL) {
        /* Possibly because the parent population died out */
        /*error("Did not expect the save_gs population to be NULL");
        */
        SET_VECTOR_ELT(R_pop_intermediate, i2-1, R_NilValue);
        continue;
      }
     
      if (population_to_matrix(population->root, r, &pop, &nrow) != 0) {
        error("Could not convert population matrix to R object");
      }

      if (nrow == 0) {
        SET_VECTOR_ELT(R_pop_intermediate, i2-1, R_NilValue);
      } else {
        PROTECT(R_pop_intermediate_tmp = allocMatrix(INTSXP, nrow, r + 1));
        for (i = 0; i < nrow; i++) {
         for (j = 0; j < r + 1; j++) {
           IMATRIX(R_pop_intermediate_tmp, i, j) = pop[i][j];
         }
         free(pop[i]);
        }
        UNPROTECT(1);
        
        SET_VECTOR_ELT(R_pop_intermediate, i2-1, R_pop_intermediate_tmp);
      }
      
      free(pop);
      
      free_tree(population);
      generations[i2] = NULL; 
    }  
    UNPROTECT(1);
  } else {
    R_pop_intermediate = R_NilValue;
  }


  

  PROTECT(R_pop_sizes = allocVector(INTSXP, g + 1));
  for (i = 0; i < g+1; i++) {
   INTEGER(R_pop_sizes)[i] = pop_sizes[i];
  }
  UNPROTECT(1);


  expected_pop_sizes = expected_population_sizes(alpha, k, g);
  PROTECT(R_expected_pop_sizes = allocVector(REALSXP, g + 1));
  for (i = 0; i <= g; i++) {
   REAL(R_expected_pop_sizes)[i] = expected_pop_sizes[i];
  }
  UNPROTECT(1);
  free(expected_pop_sizes);


  PROTECT(R_list_names = allocVector(STRSXP, 4));
  SET_STRING_ELT(R_list_names, 0, mkChar("intermediate.haplotypes"));
  SET_STRING_ELT(R_list_names, 1, mkChar("haplotypes"));
  SET_STRING_ELT(R_list_names, 2, mkChar("sizes"));
  SET_STRING_ELT(R_list_names, 3, mkChar("expected.sizes"));
  UNPROTECT(1);


  PROTECT(R_ans = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(R_ans, 0, R_pop_intermediate);
  SET_VECTOR_ELT(R_ans, 1, R_pop);
  SET_VECTOR_ELT(R_ans, 2, R_pop_sizes);
  SET_VECTOR_ELT(R_ans, 3, R_expected_pop_sizes);
  setAttrib(R_ans, R_NamesSymbol, R_list_names);
  UNPROTECT(1);  
  
  for (i = 0; i <= g; i++) {
    if (!generations[i]) continue;
    free_tree(generations[i]);
  }

  free(alpha);  
  free(generations);
  free(pop_sizes);  
  
  return(R_ans);
}

