#include "fwsim.h"

#define IMATRIX(m, i, j) (INTEGER(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])

SEXP fwsim(SEXP param_g, SEXP param_k, SEXP param_r, SEXP param_alpha, SEXP param_mu, SEXP param_trace, SEXP param_trace_loc_mut)
{
  int i, j;
  
  struct kdtree* population;
  struct kdtree** generations;
	int* pop_sizes;
	int** pop;
	int nrow = 0;
	
	SEXP R_pop;
	SEXP R_pop_sizes;
	SEXP R_list_names;
	SEXP R_ans;
	
	int g = INTEGER(param_g)[0];
	int k = INTEGER(param_k)[0];
	int r = INTEGER(param_r)[0];
	double alpha = REAL(param_alpha)[0];
	double mu = REAL(param_mu)[0];
	int trace = INTEGER(param_trace)[0];
	int trace_loc_mut = INTEGER(param_trace_loc_mut)[0];
		
	if (!(g >= 1)) error("Number of generations: g >= 1 required");
	if (!(k >= 1)) error("Size of initial population: k >= 1 required");
	if (!(r >= 1)) error("Number of loci: r >= 1 required");
	if (!(alpha > 0)) error("Growth rate: alpha > 0 required");
	if (!(mu > 0 && mu < 1)) error("Mutation rate: 0 < mu < 1 required");
	if (!(trace_loc_mut > 0)) error("Number of loci trace: trace.loc.mut > 0 required");

  if (trace) {
    Rprintf("Number of generations:      %d\n", g);
    Rprintf("Size of initial population: %d\n", k);
    Rprintf("Number of loci:             %d\n", r);
    Rprintf("Growth rate:                %f\n", alpha);
    Rprintf("Mutation rate:              %f\n", mu);
  }
  
  pop_sizes = malloc((g+1)*sizeof(int));
  if (!pop_sizes) error("Could not allocate memory for populuation sizes");
	
	GetRNGstate();
  generations = simulate_generations(g, k, r, alpha, mu, pop_sizes, trace, trace_loc_mut);
  if (!generations) error("Generations could not be simulated");
  PutRNGstate();

  population = generations[g];  
 
  if (population_to_matrix(population->root, r, &pop, &nrow) != 0) {
    error("Could not convert population matrix to R object");
  }
 	      
 	PROTECT(R_pop = allocMatrix(INTSXP, nrow, r + 1));
 	for (i = 0; i < nrow; i++) {
 	  for (j = 0; j < r + 1; j++) {
 	    IMATRIX(R_pop, i, j) = pop[i][j];
 	  }
 	}
 	UNPROTECT(1);
 	
 	PROTECT(R_pop_sizes = allocVector(INTSXP, g + 1));
 	for (i = 0; i < g + 1; i++) {
 	  INTEGER(R_pop_sizes)[i] = pop_sizes[i];
 	}
 	UNPROTECT(1);
 	
	PROTECT(R_list_names = allocVector(STRSXP, 2));
	SET_STRING_ELT(R_list_names, 0, mkChar("haplotypes"));
	SET_STRING_ELT(R_list_names, 1, mkChar("sizes"));
	UNPROTECT(1);
	
	PROTECT(R_ans = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(R_ans, 0, R_pop);
	SET_VECTOR_ELT(R_ans, 1, R_pop_sizes);
	setAttrib(R_ans, R_NamesSymbol, R_list_names);
	UNPROTECT(1);  
	
	free_tree(population);
  free(generations);
  free(pop_sizes);
  
  return(R_ans);
}

