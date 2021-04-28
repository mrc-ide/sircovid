// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"

// basic.cpp
SEXP dust_basic_alloc(cpp11::list r_pars, bool pars_multi, size_t step, cpp11::sexp r_n_particles, size_t n_threads, cpp11::sexp r_seed, cpp11::sexp device_id);
extern "C" SEXP _sircovid_dust_basic_alloc(SEXP r_pars, SEXP pars_multi, SEXP step, SEXP r_n_particles, SEXP n_threads, SEXP r_seed, SEXP device_id) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_alloc(cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars), cpp11::as_cpp<cpp11::decay_t<bool>>(pars_multi), cpp11::as_cpp<cpp11::decay_t<size_t>>(step), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_n_particles), cpp11::as_cpp<cpp11::decay_t<size_t>>(n_threads), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_seed), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(device_id)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_run(SEXP ptr, size_t step_end, bool device);
extern "C" SEXP _sircovid_dust_basic_run(SEXP ptr, SEXP step_end, SEXP device) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_run(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<size_t>>(step_end), cpp11::as_cpp<cpp11::decay_t<bool>>(device)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_simulate(SEXP ptr, cpp11::sexp step_end);
extern "C" SEXP _sircovid_dust_basic_simulate(SEXP ptr, SEXP step_end) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_simulate(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(step_end)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_set_index(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _sircovid_dust_basic_set_index(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_set_index(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_set_state(SEXP ptr, SEXP r_state, SEXP r_step);
extern "C" SEXP _sircovid_dust_basic_set_state(SEXP ptr, SEXP r_state, SEXP r_step) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_set_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_state), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_step)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_reset(SEXP ptr, cpp11::list r_pars, size_t step);
extern "C" SEXP _sircovid_dust_basic_reset(SEXP ptr, SEXP r_pars, SEXP step) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_reset(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars), cpp11::as_cpp<cpp11::decay_t<size_t>>(step)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_state(SEXP ptr, SEXP r_index);
extern "C" SEXP _sircovid_dust_basic_state(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_index)));
  END_CPP11
}
// basic.cpp
size_t dust_basic_step(SEXP ptr);
extern "C" SEXP _sircovid_dust_basic_step(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_step(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// basic.cpp
void dust_basic_reorder(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _sircovid_dust_basic_reorder(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    dust_basic_reorder(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index));
    return R_NilValue;
  END_CPP11
}
// basic.cpp
SEXP dust_basic_resample(SEXP ptr, cpp11::doubles r_weights);
extern "C" SEXP _sircovid_dust_basic_resample(SEXP ptr, SEXP r_weights) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_resample(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_weights)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_set_pars(SEXP ptr, cpp11::list r_pars);
extern "C" SEXP _sircovid_dust_basic_set_pars(SEXP ptr, SEXP r_pars) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_set_pars(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_rng_state(SEXP ptr, bool last_only);
extern "C" SEXP _sircovid_dust_basic_rng_state(SEXP ptr, SEXP last_only) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<bool>>(last_only)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_set_rng_state(SEXP ptr, cpp11::raws rng_state);
extern "C" SEXP _sircovid_dust_basic_set_rng_state(SEXP ptr, SEXP rng_state) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_set_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::raws>>(rng_state)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_set_data(SEXP ptr, cpp11::list data);
extern "C" SEXP _sircovid_dust_basic_set_data(SEXP ptr, SEXP data) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_set_data(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(data)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_compare_data(SEXP ptr, bool device);
extern "C" SEXP _sircovid_dust_basic_compare_data(SEXP ptr, SEXP device) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_compare_data(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<bool>>(device)));
  END_CPP11
}
// basic.cpp
SEXP dust_basic_filter(SEXP ptr, bool save_trajectories, cpp11::sexp step_snapshot, bool device);
extern "C" SEXP _sircovid_dust_basic_filter(SEXP ptr, SEXP save_trajectories, SEXP step_snapshot, SEXP device) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_filter(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<bool>>(save_trajectories), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(step_snapshot), cpp11::as_cpp<cpp11::decay_t<bool>>(device)));
  END_CPP11
}
// basic.cpp
cpp11::sexp dust_basic_capabilities();
extern "C" SEXP _sircovid_dust_basic_capabilities() {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_capabilities());
  END_CPP11
}
// basic.cpp
void dust_basic_set_n_threads(SEXP ptr, int n_threads);
extern "C" SEXP _sircovid_dust_basic_set_n_threads(SEXP ptr, SEXP n_threads) {
  BEGIN_CPP11
    dust_basic_set_n_threads(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<int>>(n_threads));
    return R_NilValue;
  END_CPP11
}
// basic.cpp
int dust_basic_n_state(SEXP ptr);
extern "C" SEXP _sircovid_dust_basic_n_state(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_n_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// basic.cpp
cpp11::sexp dust_basic_device_info();
extern "C" SEXP _sircovid_dust_basic_device_info() {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_device_info());
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_alloc(cpp11::list r_pars, bool pars_multi, size_t step, cpp11::sexp r_n_particles, size_t n_threads, cpp11::sexp r_seed, cpp11::sexp device_id);
extern "C" SEXP _sircovid_dust_carehomes_alloc(SEXP r_pars, SEXP pars_multi, SEXP step, SEXP r_n_particles, SEXP n_threads, SEXP r_seed, SEXP device_id) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_alloc(cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars), cpp11::as_cpp<cpp11::decay_t<bool>>(pars_multi), cpp11::as_cpp<cpp11::decay_t<size_t>>(step), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_n_particles), cpp11::as_cpp<cpp11::decay_t<size_t>>(n_threads), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_seed), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(device_id)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_run(SEXP ptr, size_t step_end, bool device);
extern "C" SEXP _sircovid_dust_carehomes_run(SEXP ptr, SEXP step_end, SEXP device) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_run(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<size_t>>(step_end), cpp11::as_cpp<cpp11::decay_t<bool>>(device)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_simulate(SEXP ptr, cpp11::sexp step_end);
extern "C" SEXP _sircovid_dust_carehomes_simulate(SEXP ptr, SEXP step_end) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_simulate(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(step_end)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_set_index(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _sircovid_dust_carehomes_set_index(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_set_index(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_set_state(SEXP ptr, SEXP r_state, SEXP r_step);
extern "C" SEXP _sircovid_dust_carehomes_set_state(SEXP ptr, SEXP r_state, SEXP r_step) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_set_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_state), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_step)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_reset(SEXP ptr, cpp11::list r_pars, size_t step);
extern "C" SEXP _sircovid_dust_carehomes_reset(SEXP ptr, SEXP r_pars, SEXP step) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_reset(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars), cpp11::as_cpp<cpp11::decay_t<size_t>>(step)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_state(SEXP ptr, SEXP r_index);
extern "C" SEXP _sircovid_dust_carehomes_state(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_index)));
  END_CPP11
}
// carehomes.cpp
size_t dust_carehomes_step(SEXP ptr);
extern "C" SEXP _sircovid_dust_carehomes_step(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_step(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// carehomes.cpp
void dust_carehomes_reorder(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _sircovid_dust_carehomes_reorder(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    dust_carehomes_reorder(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index));
    return R_NilValue;
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_resample(SEXP ptr, cpp11::doubles r_weights);
extern "C" SEXP _sircovid_dust_carehomes_resample(SEXP ptr, SEXP r_weights) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_resample(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_weights)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_set_pars(SEXP ptr, cpp11::list r_pars);
extern "C" SEXP _sircovid_dust_carehomes_set_pars(SEXP ptr, SEXP r_pars) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_set_pars(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_rng_state(SEXP ptr, bool last_only);
extern "C" SEXP _sircovid_dust_carehomes_rng_state(SEXP ptr, SEXP last_only) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<bool>>(last_only)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_set_rng_state(SEXP ptr, cpp11::raws rng_state);
extern "C" SEXP _sircovid_dust_carehomes_set_rng_state(SEXP ptr, SEXP rng_state) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_set_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::raws>>(rng_state)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_set_data(SEXP ptr, cpp11::list data);
extern "C" SEXP _sircovid_dust_carehomes_set_data(SEXP ptr, SEXP data) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_set_data(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(data)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_compare_data(SEXP ptr, bool device);
extern "C" SEXP _sircovid_dust_carehomes_compare_data(SEXP ptr, SEXP device) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_compare_data(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<bool>>(device)));
  END_CPP11
}
// carehomes.cpp
SEXP dust_carehomes_filter(SEXP ptr, bool save_trajectories, cpp11::sexp step_snapshot, bool device);
extern "C" SEXP _sircovid_dust_carehomes_filter(SEXP ptr, SEXP save_trajectories, SEXP step_snapshot, SEXP device) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_filter(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<bool>>(save_trajectories), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(step_snapshot), cpp11::as_cpp<cpp11::decay_t<bool>>(device)));
  END_CPP11
}
// carehomes.cpp
cpp11::sexp dust_carehomes_capabilities();
extern "C" SEXP _sircovid_dust_carehomes_capabilities() {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_capabilities());
  END_CPP11
}
// carehomes.cpp
void dust_carehomes_set_n_threads(SEXP ptr, int n_threads);
extern "C" SEXP _sircovid_dust_carehomes_set_n_threads(SEXP ptr, SEXP n_threads) {
  BEGIN_CPP11
    dust_carehomes_set_n_threads(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<int>>(n_threads));
    return R_NilValue;
  END_CPP11
}
// carehomes.cpp
int dust_carehomes_n_state(SEXP ptr);
extern "C" SEXP _sircovid_dust_carehomes_n_state(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_n_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// carehomes.cpp
cpp11::sexp dust_carehomes_device_info();
extern "C" SEXP _sircovid_dust_carehomes_device_info() {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_carehomes_device_info());
  END_CPP11
}

extern "C" {
/* .Call calls */
extern SEXP _sircovid_dust_basic_alloc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_basic_capabilities();
extern SEXP _sircovid_dust_basic_compare_data(SEXP, SEXP);
extern SEXP _sircovid_dust_basic_device_info();
extern SEXP _sircovid_dust_basic_filter(SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_basic_n_state(SEXP);
extern SEXP _sircovid_dust_basic_reorder(SEXP, SEXP);
extern SEXP _sircovid_dust_basic_resample(SEXP, SEXP);
extern SEXP _sircovid_dust_basic_reset(SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_basic_rng_state(SEXP, SEXP);
extern SEXP _sircovid_dust_basic_run(SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_basic_set_data(SEXP, SEXP);
extern SEXP _sircovid_dust_basic_set_index(SEXP, SEXP);
extern SEXP _sircovid_dust_basic_set_n_threads(SEXP, SEXP);
extern SEXP _sircovid_dust_basic_set_pars(SEXP, SEXP);
extern SEXP _sircovid_dust_basic_set_rng_state(SEXP, SEXP);
extern SEXP _sircovid_dust_basic_set_state(SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_basic_simulate(SEXP, SEXP);
extern SEXP _sircovid_dust_basic_state(SEXP, SEXP);
extern SEXP _sircovid_dust_basic_step(SEXP);
extern SEXP _sircovid_dust_carehomes_alloc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_capabilities();
extern SEXP _sircovid_dust_carehomes_compare_data(SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_device_info();
extern SEXP _sircovid_dust_carehomes_filter(SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_n_state(SEXP);
extern SEXP _sircovid_dust_carehomes_reorder(SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_resample(SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_reset(SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_rng_state(SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_run(SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_set_data(SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_set_index(SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_set_n_threads(SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_set_pars(SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_set_rng_state(SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_set_state(SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_simulate(SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_state(SEXP, SEXP);
extern SEXP _sircovid_dust_carehomes_step(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_sircovid_dust_basic_alloc",             (DL_FUNC) &_sircovid_dust_basic_alloc,             7},
    {"_sircovid_dust_basic_capabilities",      (DL_FUNC) &_sircovid_dust_basic_capabilities,      0},
    {"_sircovid_dust_basic_compare_data",      (DL_FUNC) &_sircovid_dust_basic_compare_data,      2},
    {"_sircovid_dust_basic_device_info",       (DL_FUNC) &_sircovid_dust_basic_device_info,       0},
    {"_sircovid_dust_basic_filter",            (DL_FUNC) &_sircovid_dust_basic_filter,            4},
    {"_sircovid_dust_basic_n_state",           (DL_FUNC) &_sircovid_dust_basic_n_state,           1},
    {"_sircovid_dust_basic_reorder",           (DL_FUNC) &_sircovid_dust_basic_reorder,           2},
    {"_sircovid_dust_basic_resample",          (DL_FUNC) &_sircovid_dust_basic_resample,          2},
    {"_sircovid_dust_basic_reset",             (DL_FUNC) &_sircovid_dust_basic_reset,             3},
    {"_sircovid_dust_basic_rng_state",         (DL_FUNC) &_sircovid_dust_basic_rng_state,         2},
    {"_sircovid_dust_basic_run",               (DL_FUNC) &_sircovid_dust_basic_run,               3},
    {"_sircovid_dust_basic_set_data",          (DL_FUNC) &_sircovid_dust_basic_set_data,          2},
    {"_sircovid_dust_basic_set_index",         (DL_FUNC) &_sircovid_dust_basic_set_index,         2},
    {"_sircovid_dust_basic_set_n_threads",     (DL_FUNC) &_sircovid_dust_basic_set_n_threads,     2},
    {"_sircovid_dust_basic_set_pars",          (DL_FUNC) &_sircovid_dust_basic_set_pars,          2},
    {"_sircovid_dust_basic_set_rng_state",     (DL_FUNC) &_sircovid_dust_basic_set_rng_state,     2},
    {"_sircovid_dust_basic_set_state",         (DL_FUNC) &_sircovid_dust_basic_set_state,         3},
    {"_sircovid_dust_basic_simulate",          (DL_FUNC) &_sircovid_dust_basic_simulate,          2},
    {"_sircovid_dust_basic_state",             (DL_FUNC) &_sircovid_dust_basic_state,             2},
    {"_sircovid_dust_basic_step",              (DL_FUNC) &_sircovid_dust_basic_step,              1},
    {"_sircovid_dust_carehomes_alloc",         (DL_FUNC) &_sircovid_dust_carehomes_alloc,         7},
    {"_sircovid_dust_carehomes_capabilities",  (DL_FUNC) &_sircovid_dust_carehomes_capabilities,  0},
    {"_sircovid_dust_carehomes_compare_data",  (DL_FUNC) &_sircovid_dust_carehomes_compare_data,  2},
    {"_sircovid_dust_carehomes_device_info",   (DL_FUNC) &_sircovid_dust_carehomes_device_info,   0},
    {"_sircovid_dust_carehomes_filter",        (DL_FUNC) &_sircovid_dust_carehomes_filter,        4},
    {"_sircovid_dust_carehomes_n_state",       (DL_FUNC) &_sircovid_dust_carehomes_n_state,       1},
    {"_sircovid_dust_carehomes_reorder",       (DL_FUNC) &_sircovid_dust_carehomes_reorder,       2},
    {"_sircovid_dust_carehomes_resample",      (DL_FUNC) &_sircovid_dust_carehomes_resample,      2},
    {"_sircovid_dust_carehomes_reset",         (DL_FUNC) &_sircovid_dust_carehomes_reset,         3},
    {"_sircovid_dust_carehomes_rng_state",     (DL_FUNC) &_sircovid_dust_carehomes_rng_state,     2},
    {"_sircovid_dust_carehomes_run",           (DL_FUNC) &_sircovid_dust_carehomes_run,           3},
    {"_sircovid_dust_carehomes_set_data",      (DL_FUNC) &_sircovid_dust_carehomes_set_data,      2},
    {"_sircovid_dust_carehomes_set_index",     (DL_FUNC) &_sircovid_dust_carehomes_set_index,     2},
    {"_sircovid_dust_carehomes_set_n_threads", (DL_FUNC) &_sircovid_dust_carehomes_set_n_threads, 2},
    {"_sircovid_dust_carehomes_set_pars",      (DL_FUNC) &_sircovid_dust_carehomes_set_pars,      2},
    {"_sircovid_dust_carehomes_set_rng_state", (DL_FUNC) &_sircovid_dust_carehomes_set_rng_state, 2},
    {"_sircovid_dust_carehomes_set_state",     (DL_FUNC) &_sircovid_dust_carehomes_set_state,     3},
    {"_sircovid_dust_carehomes_simulate",      (DL_FUNC) &_sircovid_dust_carehomes_simulate,      2},
    {"_sircovid_dust_carehomes_state",         (DL_FUNC) &_sircovid_dust_carehomes_state,         2},
    {"_sircovid_dust_carehomes_step",          (DL_FUNC) &_sircovid_dust_carehomes_step,          1},
    {NULL, NULL, 0}
};
}

extern "C" void R_init_sircovid(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
