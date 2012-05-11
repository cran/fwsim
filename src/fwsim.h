#ifndef _FWSIM_H_
#define _FWSIM_H_

#include "common.h"
#include "math.h"
#include "sim.h"
#include "kdtree.h"
#include "clean.h"
#include "hap.h"
#include "print.h"

SEXP fwsim(SEXP param_g, SEXP param_k, SEXP param_r, 
  SEXP param_alpha, SEXP param_mu, 
  SEXP param_save_gs,
  SEXP param_trace,
  SEXP param_alpha_length);

#endif /* _FWSIM_H_ */

