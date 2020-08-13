// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"

// basic.cpp
SEXP dust_basic_alloc(cpp11::list r_data, size_t step, size_t n_particles, size_t n_threads, size_t seed);
extern "C" SEXP _sircovid2_dust_basic_alloc(SEXP r_data, SEXP step, SEXP n_particles, SEXP n_threads, SEXP seed) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_alloc(cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_data), cpp11::as_cpp<cpp11::decay_t<size_t>>(step), cpp11::as_cpp<cpp11::decay_t<size_t>>(n_particles), cpp11::as_cpp<cpp11::decay_t<size_t>>(n_threads), cpp11::as_cpp<cpp11::decay_t<size_t>>(seed)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_run(SEXP ptr, size_t step_end);
extern "C" SEXP _sircovid2_dust_basic_run(SEXP ptr, SEXP step_end) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_run(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<size_t>>(step_end)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_set_index(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _sircovid2_dust_basic_set_index(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_set_index(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_set_state(SEXP ptr, SEXP r_state, SEXP r_step);
extern "C" SEXP _sircovid2_dust_basic_set_state(SEXP ptr, SEXP r_state, SEXP r_step) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_set_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_state), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_step)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_reset(SEXP ptr, cpp11::list r_data, size_t step);
extern "C" SEXP _sircovid2_dust_basic_reset(SEXP ptr, SEXP r_data, SEXP step) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_reset(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_data), cpp11::as_cpp<cpp11::decay_t<size_t>>(step)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_state(SEXP ptr, SEXP r_index);
extern "C" SEXP _sircovid2_dust_basic_state(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_index)));
  END_CPP11
}
// basic.cpp
size_t dust_basic_step(SEXP ptr);
extern "C" SEXP _sircovid2_dust_basic_step(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_step(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// basic.cpp
void dust_basic_reorder(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _sircovid2_dust_basic_reorder(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    dust_basic_reorder(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index));
    return R_NilValue;
  END_CPP11
}
// basic.cpp
SEXP dust_basic_rng_state(SEXP ptr);
extern "C" SEXP _sircovid2_dust_basic_rng_state(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_simulate(cpp11::sexp r_steps, cpp11::list r_data, cpp11::doubles_matrix r_state, cpp11::sexp r_index, const size_t n_threads, const size_t seed);
extern "C" SEXP _sircovid2_dust_basic_simulate(SEXP r_steps, SEXP r_data, SEXP r_state, SEXP r_index, SEXP n_threads, SEXP seed) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_simulate(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_steps), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_data), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles_matrix>>(r_state), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index), cpp11::as_cpp<cpp11::decay_t<const size_t>>(n_threads), cpp11::as_cpp<cpp11::decay_t<const size_t>>(seed)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_alloc(cpp11::list r_data, size_t step, size_t n_particles, size_t n_threads, size_t seed);
extern "C" SEXP _sircovid2_dust_carehomes_alloc(SEXP r_data, SEXP step, SEXP n_particles, SEXP n_threads, SEXP seed) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_alloc(cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_data), cpp11::as_cpp<cpp11::decay_t<size_t>>(step), cpp11::as_cpp<cpp11::decay_t<size_t>>(n_particles), cpp11::as_cpp<cpp11::decay_t<size_t>>(n_threads), cpp11::as_cpp<cpp11::decay_t<size_t>>(seed)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_run(SEXP ptr, size_t step_end);
extern "C" SEXP _sircovid2_dust_carehomes_run(SEXP ptr, SEXP step_end) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_run(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<size_t>>(step_end)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_set_index(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _sircovid2_dust_carehomes_set_index(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_set_index(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_set_state(SEXP ptr, SEXP r_state, SEXP r_step);
extern "C" SEXP _sircovid2_dust_carehomes_set_state(SEXP ptr, SEXP r_state, SEXP r_step) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_set_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_state), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_step)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_reset(SEXP ptr, cpp11::list r_data, size_t step);
extern "C" SEXP _sircovid2_dust_carehomes_reset(SEXP ptr, SEXP r_data, SEXP step) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_reset(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_data), cpp11::as_cpp<cpp11::decay_t<size_t>>(step)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_state(SEXP ptr, SEXP r_index);
extern "C" SEXP _sircovid2_dust_carehomes_state(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_index)));
  END_CPP11
}
// carehomes.cpp
size_t dust_carehomes_step(SEXP ptr);
extern "C" SEXP _sircovid2_dust_carehomes_step(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_step(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// carehomes.cpp
void dust_carehomes_reorder(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _sircovid2_dust_carehomes_reorder(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    dust_carehomes_reorder(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index));
    return R_NilValue;
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_rng_state(SEXP ptr);
extern "C" SEXP _sircovid2_dust_carehomes_rng_state(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_simulate(cpp11::sexp r_steps, cpp11::list r_data, cpp11::doubles_matrix r_state, cpp11::sexp r_index, const size_t n_threads, const size_t seed);
extern "C" SEXP _sircovid2_dust_carehomes_simulate(SEXP r_steps, SEXP r_data, SEXP r_state, SEXP r_index, SEXP n_threads, SEXP seed) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_simulate(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_steps), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_data), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles_matrix>>(r_state), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index), cpp11::as_cpp<cpp11::decay_t<const size_t>>(n_threads), cpp11::as_cpp<cpp11::decay_t<const size_t>>(seed)));
  END_CPP11
}

extern "C" {
/* .Call calls */
extern SEXP _sircovid2_dust_basic_alloc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid2_dust_basic_reorder(SEXP, SEXP);
extern SEXP _sircovid2_dust_basic_reset(SEXP, SEXP, SEXP);
extern SEXP _sircovid2_dust_basic_rng_state(SEXP);
extern SEXP _sircovid2_dust_basic_run(SEXP, SEXP);
extern SEXP _sircovid2_dust_basic_set_index(SEXP, SEXP);
extern SEXP _sircovid2_dust_basic_set_state(SEXP, SEXP, SEXP);
extern SEXP _sircovid2_dust_basic_simulate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid2_dust_basic_state(SEXP, SEXP);
extern SEXP _sircovid2_dust_basic_step(SEXP);
extern SEXP _sircovid2_dust_carehomes_alloc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid2_dust_carehomes_reorder(SEXP, SEXP);
extern SEXP _sircovid2_dust_carehomes_reset(SEXP, SEXP, SEXP);
extern SEXP _sircovid2_dust_carehomes_rng_state(SEXP);
extern SEXP _sircovid2_dust_carehomes_run(SEXP, SEXP);
extern SEXP _sircovid2_dust_carehomes_set_index(SEXP, SEXP);
extern SEXP _sircovid2_dust_carehomes_set_state(SEXP, SEXP, SEXP);
extern SEXP _sircovid2_dust_carehomes_simulate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid2_dust_carehomes_state(SEXP, SEXP);
extern SEXP _sircovid2_dust_carehomes_step(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_sircovid2_dust_basic_alloc",         (DL_FUNC) &_sircovid2_dust_basic_alloc,         5},
    {"_sircovid2_dust_basic_reorder",       (DL_FUNC) &_sircovid2_dust_basic_reorder,       2},
    {"_sircovid2_dust_basic_reset",         (DL_FUNC) &_sircovid2_dust_basic_reset,         3},
    {"_sircovid2_dust_basic_rng_state",     (DL_FUNC) &_sircovid2_dust_basic_rng_state,     1},
    {"_sircovid2_dust_basic_run",           (DL_FUNC) &_sircovid2_dust_basic_run,           2},
    {"_sircovid2_dust_basic_set_index",     (DL_FUNC) &_sircovid2_dust_basic_set_index,     2},
    {"_sircovid2_dust_basic_set_state",     (DL_FUNC) &_sircovid2_dust_basic_set_state,     3},
    {"_sircovid2_dust_basic_simulate",      (DL_FUNC) &_sircovid2_dust_basic_simulate,      6},
    {"_sircovid2_dust_basic_state",         (DL_FUNC) &_sircovid2_dust_basic_state,         2},
    {"_sircovid2_dust_basic_step",          (DL_FUNC) &_sircovid2_dust_basic_step,          1},
    {"_sircovid2_dust_carehomes_alloc",     (DL_FUNC) &_sircovid2_dust_carehomes_alloc,     5},
    {"_sircovid2_dust_carehomes_reorder",   (DL_FUNC) &_sircovid2_dust_carehomes_reorder,   2},
    {"_sircovid2_dust_carehomes_reset",     (DL_FUNC) &_sircovid2_dust_carehomes_reset,     3},
    {"_sircovid2_dust_carehomes_rng_state", (DL_FUNC) &_sircovid2_dust_carehomes_rng_state, 1},
    {"_sircovid2_dust_carehomes_run",       (DL_FUNC) &_sircovid2_dust_carehomes_run,       2},
    {"_sircovid2_dust_carehomes_set_index", (DL_FUNC) &_sircovid2_dust_carehomes_set_index, 2},
    {"_sircovid2_dust_carehomes_set_state", (DL_FUNC) &_sircovid2_dust_carehomes_set_state, 3},
    {"_sircovid2_dust_carehomes_simulate",  (DL_FUNC) &_sircovid2_dust_carehomes_simulate,  6},
    {"_sircovid2_dust_carehomes_state",     (DL_FUNC) &_sircovid2_dust_carehomes_state,     2},
    {"_sircovid2_dust_carehomes_step",      (DL_FUNC) &_sircovid2_dust_carehomes_step,      1},
    {NULL, NULL, 0}
};
}

extern "C" void R_init_sircovid2(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
