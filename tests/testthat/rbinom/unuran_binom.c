#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <unuran.h>

double unuran_rbinom(double n, double p);

SEXP test_unuran_binom_r(SEXP r_n, SEXP r_p) {
  size_t n_sample = (size_t)length(r_n);
  int *n = INTEGER(r_n);
  double *p = REAL(r_p);
  // GetRNGstate();
  for (size_t i = 0; i < n_sample; ++i) {
    unuran_rbinom(n[i], p[i]);
  }
  // PutRNGstate();
  return R_NilValue;
}

double unuran_rbinom(double n, double p) {
   double x; /* will hold the random number */
   double parameter[2] = {n, p};

   /* Declare the three UNURAN objects. */
   UNUR_DISTR *distr; /* distribution object */
   UNUR_PAR *par; /* parameter object */
   UNUR_GEN *gen; /* generator object */

   distr = unur_distr_binomial(parameter, 2);
   // par = unur_dstd_new(distr);
   //unur_dstd_set_variant(par, UNUR_STDGEN_DEFAULT);
   par = unur_dsrou_new(distr);

   gen = unur_init(par);
   if (gen == NULL) {
      fprintf(stderr, "ERROR: cannot create generator object\n");
      exit (EXIT_FAILURE);
   }
   /* It is possible to reuse the distribution object to create */
   /* another generator object. If you do not need it any more, */
   /* it should be destroyed to free memory. */
   unur_distr_free(distr);
   /* Now you can use the generator object ‘gen’ to sample from */
   /* the standard Gaussian distribution. */
   /* Eg.: */
   x = unur_sample_discr(gen);
   /* When you do not need the generator object any more, you */
   /* can destroy it. */
   unur_free(gen);

   return(x);
}

