#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

SEXP test_binom_r(SEXP r_n, SEXP r_p) {
  size_t n_sample = (size_t)length(r_n);
  int *n = INTEGER(r_n);
  double *p = REAL(r_p);
  GetRNGstate();
  for (size_t i = 0; i < n_sample; ++i) {
    Rf_rbinom(n[i], p[i]);
  }
  PutRNGstate();
  return R_NilValue;
}
