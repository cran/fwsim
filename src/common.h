#ifndef _COMMON_H_
#define _COMMON_H_

#define VERBOSE_
#define EXTRA_VERBOSE_

#define MUTATION_CATEGORY_THRESHOLD -1
/*
#define IS_BIT_SET(var, pos) ((var) & (1 << (pos)))
*/
/*
  Double negation to only get 0 and 1
*/
#define IS_BIT_SET(var, pos) !!((var) & (1 << (pos)))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#if defined(WIN32) || defined(__WIN32__)
#include <malloc.h>
#endif

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

/*
#define MATHLIB_STANDALONE
#include <Rmath.h>
*/

#endif /* _COMMON_H_ */

