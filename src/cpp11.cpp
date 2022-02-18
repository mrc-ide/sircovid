// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"

// basic.cpp
cpp11::sexp dust_basic_capabilities();
extern "C" SEXP _sircovid_dust_basic_capabilities() {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_capabilities());
  END_CPP11
}
// basic.cpp
cpp11::sexp dust_basic_gpu_info();
extern "C" SEXP _sircovid_dust_basic_gpu_info() {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_basic_gpu_info());
  END_CPP11
}
// basic.cpp
SEXP dust_cpu_basic_alloc(cpp11::list r_pars, bool pars_multi, size_t step, cpp11::sexp r_n_particles, size_t n_threads, cpp11::sexp r_seed, bool deterministic, cpp11::sexp gpu_config);
extern "C" SEXP _sircovid_dust_cpu_basic_alloc(SEXP r_pars, SEXP pars_multi, SEXP step, SEXP r_n_particles, SEXP n_threads, SEXP r_seed, SEXP deterministic, SEXP gpu_config) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_alloc(cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars), cpp11::as_cpp<cpp11::decay_t<bool>>(pars_multi), cpp11::as_cpp<cpp11::decay_t<size_t>>(step), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_n_particles), cpp11::as_cpp<cpp11::decay_t<size_t>>(n_threads), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_seed), cpp11::as_cpp<cpp11::decay_t<bool>>(deterministic), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(gpu_config)));
  END_CPP11
}
// basic.cpp
SEXP dust_cpu_basic_run(SEXP ptr, size_t step_end);
extern "C" SEXP _sircovid_dust_cpu_basic_run(SEXP ptr, SEXP step_end) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_run(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<size_t>>(step_end)));
  END_CPP11
}
// basic.cpp
SEXP dust_cpu_basic_simulate(SEXP ptr, cpp11::sexp step_end);
extern "C" SEXP _sircovid_dust_cpu_basic_simulate(SEXP ptr, SEXP step_end) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_simulate(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(step_end)));
  END_CPP11
}
// basic.cpp
SEXP dust_cpu_basic_set_index(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _sircovid_dust_cpu_basic_set_index(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_set_index(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index)));
  END_CPP11
}
// basic.cpp
SEXP dust_cpu_basic_update_state(SEXP ptr, SEXP r_pars, SEXP r_state, SEXP r_step, SEXP r_set_initial_state);
extern "C" SEXP _sircovid_dust_cpu_basic_update_state(SEXP ptr, SEXP r_pars, SEXP r_state, SEXP r_step, SEXP r_set_initial_state) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_update_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_pars), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_state), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_step), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_set_initial_state)));
  END_CPP11
}
// basic.cpp
SEXP dust_cpu_basic_state(SEXP ptr, SEXP r_index);
extern "C" SEXP _sircovid_dust_cpu_basic_state(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_index)));
  END_CPP11
}
// basic.cpp
size_t dust_cpu_basic_step(SEXP ptr);
extern "C" SEXP _sircovid_dust_cpu_basic_step(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_step(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// basic.cpp
void dust_cpu_basic_reorder(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _sircovid_dust_cpu_basic_reorder(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    dust_cpu_basic_reorder(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index));
    return R_NilValue;
  END_CPP11
}
// basic.cpp
SEXP dust_cpu_basic_resample(SEXP ptr, cpp11::doubles r_weights);
extern "C" SEXP _sircovid_dust_cpu_basic_resample(SEXP ptr, SEXP r_weights) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_resample(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_weights)));
  END_CPP11
}
// basic.cpp
SEXP dust_cpu_basic_rng_state(SEXP ptr, bool first_only, bool last_only);
extern "C" SEXP _sircovid_dust_cpu_basic_rng_state(SEXP ptr, SEXP first_only, SEXP last_only) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<bool>>(first_only), cpp11::as_cpp<cpp11::decay_t<bool>>(last_only)));
  END_CPP11
}
// basic.cpp
SEXP dust_cpu_basic_set_rng_state(SEXP ptr, cpp11::raws rng_state);
extern "C" SEXP _sircovid_dust_cpu_basic_set_rng_state(SEXP ptr, SEXP rng_state) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_set_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::raws>>(rng_state)));
  END_CPP11
}
// basic.cpp
SEXP dust_cpu_basic_set_data(SEXP ptr, cpp11::list data);
extern "C" SEXP _sircovid_dust_cpu_basic_set_data(SEXP ptr, SEXP data) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_set_data(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(data)));
  END_CPP11
}
// basic.cpp
SEXP dust_cpu_basic_compare_data(SEXP ptr);
extern "C" SEXP _sircovid_dust_cpu_basic_compare_data(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_compare_data(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// basic.cpp
SEXP dust_cpu_basic_filter(SEXP ptr, SEXP step_end, bool save_trajectories, cpp11::sexp step_snapshot, cpp11::sexp min_log_likelihood);
extern "C" SEXP _sircovid_dust_cpu_basic_filter(SEXP ptr, SEXP step_end, SEXP save_trajectories, SEXP step_snapshot, SEXP min_log_likelihood) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_filter(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(step_end), cpp11::as_cpp<cpp11::decay_t<bool>>(save_trajectories), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(step_snapshot), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(min_log_likelihood)));
  END_CPP11
}
// basic.cpp
void dust_cpu_basic_set_n_threads(SEXP ptr, int n_threads);
extern "C" SEXP _sircovid_dust_cpu_basic_set_n_threads(SEXP ptr, SEXP n_threads) {
  BEGIN_CPP11
    dust_cpu_basic_set_n_threads(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<int>>(n_threads));
    return R_NilValue;
  END_CPP11
}
// basic.cpp
int dust_cpu_basic_n_state(SEXP ptr);
extern "C" SEXP _sircovid_dust_cpu_basic_n_state(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_basic_n_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// lancelot.cpp
cpp11::sexp dust_lancelot_capabilities();
extern "C" SEXP _sircovid_dust_lancelot_capabilities() {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_lancelot_capabilities());
  END_CPP11
}
// lancelot.cpp
cpp11::sexp dust_lancelot_gpu_info();
extern "C" SEXP _sircovid_dust_lancelot_gpu_info() {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_lancelot_gpu_info());
  END_CPP11
}
// lancelot.cpp
SEXP dust_cpu_lancelot_alloc(cpp11::list r_pars, bool pars_multi, size_t step, cpp11::sexp r_n_particles, size_t n_threads, cpp11::sexp r_seed, bool deterministic, cpp11::sexp gpu_config);
extern "C" SEXP _sircovid_dust_cpu_lancelot_alloc(SEXP r_pars, SEXP pars_multi, SEXP step, SEXP r_n_particles, SEXP n_threads, SEXP r_seed, SEXP deterministic, SEXP gpu_config) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_alloc(cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars), cpp11::as_cpp<cpp11::decay_t<bool>>(pars_multi), cpp11::as_cpp<cpp11::decay_t<size_t>>(step), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_n_particles), cpp11::as_cpp<cpp11::decay_t<size_t>>(n_threads), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_seed), cpp11::as_cpp<cpp11::decay_t<bool>>(deterministic), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(gpu_config)));
  END_CPP11
}
// lancelot.cpp
SEXP dust_cpu_lancelot_run(SEXP ptr, size_t step_end);
extern "C" SEXP _sircovid_dust_cpu_lancelot_run(SEXP ptr, SEXP step_end) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_run(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<size_t>>(step_end)));
  END_CPP11
}
// lancelot.cpp
SEXP dust_cpu_lancelot_simulate(SEXP ptr, cpp11::sexp step_end);
extern "C" SEXP _sircovid_dust_cpu_lancelot_simulate(SEXP ptr, SEXP step_end) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_simulate(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(step_end)));
  END_CPP11
}
// lancelot.cpp
SEXP dust_cpu_lancelot_set_index(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _sircovid_dust_cpu_lancelot_set_index(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_set_index(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index)));
  END_CPP11
}
// lancelot.cpp
SEXP dust_cpu_lancelot_update_state(SEXP ptr, SEXP r_pars, SEXP r_state, SEXP r_step, SEXP r_set_initial_state);
extern "C" SEXP _sircovid_dust_cpu_lancelot_update_state(SEXP ptr, SEXP r_pars, SEXP r_state, SEXP r_step, SEXP r_set_initial_state) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_update_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_pars), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_state), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_step), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_set_initial_state)));
  END_CPP11
}
// lancelot.cpp
SEXP dust_cpu_lancelot_state(SEXP ptr, SEXP r_index);
extern "C" SEXP _sircovid_dust_cpu_lancelot_state(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_index)));
  END_CPP11
}
// lancelot.cpp
size_t dust_cpu_lancelot_step(SEXP ptr);
extern "C" SEXP _sircovid_dust_cpu_lancelot_step(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_step(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// lancelot.cpp
void dust_cpu_lancelot_reorder(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _sircovid_dust_cpu_lancelot_reorder(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    dust_cpu_lancelot_reorder(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index));
    return R_NilValue;
  END_CPP11
}
// lancelot.cpp
SEXP dust_cpu_lancelot_resample(SEXP ptr, cpp11::doubles r_weights);
extern "C" SEXP _sircovid_dust_cpu_lancelot_resample(SEXP ptr, SEXP r_weights) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_resample(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_weights)));
  END_CPP11
}
// lancelot.cpp
SEXP dust_cpu_lancelot_rng_state(SEXP ptr, bool first_only, bool last_only);
extern "C" SEXP _sircovid_dust_cpu_lancelot_rng_state(SEXP ptr, SEXP first_only, SEXP last_only) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<bool>>(first_only), cpp11::as_cpp<cpp11::decay_t<bool>>(last_only)));
  END_CPP11
}
// lancelot.cpp
SEXP dust_cpu_lancelot_set_rng_state(SEXP ptr, cpp11::raws rng_state);
extern "C" SEXP _sircovid_dust_cpu_lancelot_set_rng_state(SEXP ptr, SEXP rng_state) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_set_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::raws>>(rng_state)));
  END_CPP11
}
// lancelot.cpp
SEXP dust_cpu_lancelot_set_data(SEXP ptr, cpp11::list data);
extern "C" SEXP _sircovid_dust_cpu_lancelot_set_data(SEXP ptr, SEXP data) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_set_data(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(data)));
  END_CPP11
}
// lancelot.cpp
SEXP dust_cpu_lancelot_compare_data(SEXP ptr);
extern "C" SEXP _sircovid_dust_cpu_lancelot_compare_data(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_compare_data(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// lancelot.cpp
SEXP dust_cpu_lancelot_filter(SEXP ptr, SEXP step_end, bool save_trajectories, cpp11::sexp step_snapshot, cpp11::sexp min_log_likelihood);
extern "C" SEXP _sircovid_dust_cpu_lancelot_filter(SEXP ptr, SEXP step_end, SEXP save_trajectories, SEXP step_snapshot, SEXP min_log_likelihood) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_filter(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(step_end), cpp11::as_cpp<cpp11::decay_t<bool>>(save_trajectories), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(step_snapshot), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(min_log_likelihood)));
  END_CPP11
}
// lancelot.cpp
void dust_cpu_lancelot_set_n_threads(SEXP ptr, int n_threads);
extern "C" SEXP _sircovid_dust_cpu_lancelot_set_n_threads(SEXP ptr, SEXP n_threads) {
  BEGIN_CPP11
    dust_cpu_lancelot_set_n_threads(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<int>>(n_threads));
    return R_NilValue;
  END_CPP11
}
// lancelot.cpp
int dust_cpu_lancelot_n_state(SEXP ptr);
extern "C" SEXP _sircovid_dust_cpu_lancelot_n_state(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_lancelot_n_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}

extern "C" {
/* .Call calls */
extern SEXP _sircovid_dust_basic_capabilities();
extern SEXP _sircovid_dust_basic_gpu_info();
extern SEXP _sircovid_dust_cpu_basic_alloc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_basic_compare_data(SEXP);
extern SEXP _sircovid_dust_cpu_basic_filter(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_basic_n_state(SEXP);
extern SEXP _sircovid_dust_cpu_basic_reorder(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_basic_resample(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_basic_rng_state(SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_basic_run(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_basic_set_data(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_basic_set_index(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_basic_set_n_threads(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_basic_set_rng_state(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_basic_simulate(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_basic_state(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_basic_step(SEXP);
extern SEXP _sircovid_dust_cpu_basic_update_state(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_alloc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_compare_data(SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_filter(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_n_state(SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_reorder(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_resample(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_rng_state(SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_run(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_set_data(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_set_index(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_set_n_threads(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_set_rng_state(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_simulate(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_state(SEXP, SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_step(SEXP);
extern SEXP _sircovid_dust_cpu_lancelot_update_state(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sircovid_dust_lancelot_capabilities();
extern SEXP _sircovid_dust_lancelot_gpu_info();

static const R_CallMethodDef CallEntries[] = {
    {"_sircovid_dust_basic_capabilities",         (DL_FUNC) &_sircovid_dust_basic_capabilities,         0},
    {"_sircovid_dust_basic_gpu_info",             (DL_FUNC) &_sircovid_dust_basic_gpu_info,             0},
    {"_sircovid_dust_cpu_basic_alloc",            (DL_FUNC) &_sircovid_dust_cpu_basic_alloc,            8},
    {"_sircovid_dust_cpu_basic_compare_data",     (DL_FUNC) &_sircovid_dust_cpu_basic_compare_data,     1},
    {"_sircovid_dust_cpu_basic_filter",           (DL_FUNC) &_sircovid_dust_cpu_basic_filter,           5},
    {"_sircovid_dust_cpu_basic_n_state",          (DL_FUNC) &_sircovid_dust_cpu_basic_n_state,          1},
    {"_sircovid_dust_cpu_basic_reorder",          (DL_FUNC) &_sircovid_dust_cpu_basic_reorder,          2},
    {"_sircovid_dust_cpu_basic_resample",         (DL_FUNC) &_sircovid_dust_cpu_basic_resample,         2},
    {"_sircovid_dust_cpu_basic_rng_state",        (DL_FUNC) &_sircovid_dust_cpu_basic_rng_state,        3},
    {"_sircovid_dust_cpu_basic_run",              (DL_FUNC) &_sircovid_dust_cpu_basic_run,              2},
    {"_sircovid_dust_cpu_basic_set_data",         (DL_FUNC) &_sircovid_dust_cpu_basic_set_data,         2},
    {"_sircovid_dust_cpu_basic_set_index",        (DL_FUNC) &_sircovid_dust_cpu_basic_set_index,        2},
    {"_sircovid_dust_cpu_basic_set_n_threads",    (DL_FUNC) &_sircovid_dust_cpu_basic_set_n_threads,    2},
    {"_sircovid_dust_cpu_basic_set_rng_state",    (DL_FUNC) &_sircovid_dust_cpu_basic_set_rng_state,    2},
    {"_sircovid_dust_cpu_basic_simulate",         (DL_FUNC) &_sircovid_dust_cpu_basic_simulate,         2},
    {"_sircovid_dust_cpu_basic_state",            (DL_FUNC) &_sircovid_dust_cpu_basic_state,            2},
    {"_sircovid_dust_cpu_basic_step",             (DL_FUNC) &_sircovid_dust_cpu_basic_step,             1},
    {"_sircovid_dust_cpu_basic_update_state",     (DL_FUNC) &_sircovid_dust_cpu_basic_update_state,     5},
    {"_sircovid_dust_cpu_lancelot_alloc",         (DL_FUNC) &_sircovid_dust_cpu_lancelot_alloc,         8},
    {"_sircovid_dust_cpu_lancelot_compare_data",  (DL_FUNC) &_sircovid_dust_cpu_lancelot_compare_data,  1},
    {"_sircovid_dust_cpu_lancelot_filter",        (DL_FUNC) &_sircovid_dust_cpu_lancelot_filter,        5},
    {"_sircovid_dust_cpu_lancelot_n_state",       (DL_FUNC) &_sircovid_dust_cpu_lancelot_n_state,       1},
    {"_sircovid_dust_cpu_lancelot_reorder",       (DL_FUNC) &_sircovid_dust_cpu_lancelot_reorder,       2},
    {"_sircovid_dust_cpu_lancelot_resample",      (DL_FUNC) &_sircovid_dust_cpu_lancelot_resample,      2},
    {"_sircovid_dust_cpu_lancelot_rng_state",     (DL_FUNC) &_sircovid_dust_cpu_lancelot_rng_state,     3},
    {"_sircovid_dust_cpu_lancelot_run",           (DL_FUNC) &_sircovid_dust_cpu_lancelot_run,           2},
    {"_sircovid_dust_cpu_lancelot_set_data",      (DL_FUNC) &_sircovid_dust_cpu_lancelot_set_data,      2},
    {"_sircovid_dust_cpu_lancelot_set_index",     (DL_FUNC) &_sircovid_dust_cpu_lancelot_set_index,     2},
    {"_sircovid_dust_cpu_lancelot_set_n_threads", (DL_FUNC) &_sircovid_dust_cpu_lancelot_set_n_threads, 2},
    {"_sircovid_dust_cpu_lancelot_set_rng_state", (DL_FUNC) &_sircovid_dust_cpu_lancelot_set_rng_state, 2},
    {"_sircovid_dust_cpu_lancelot_simulate",      (DL_FUNC) &_sircovid_dust_cpu_lancelot_simulate,      2},
    {"_sircovid_dust_cpu_lancelot_state",         (DL_FUNC) &_sircovid_dust_cpu_lancelot_state,         2},
    {"_sircovid_dust_cpu_lancelot_step",          (DL_FUNC) &_sircovid_dust_cpu_lancelot_step,          1},
    {"_sircovid_dust_cpu_lancelot_update_state",  (DL_FUNC) &_sircovid_dust_cpu_lancelot_update_state,  5},
    {"_sircovid_dust_lancelot_capabilities",      (DL_FUNC) &_sircovid_dust_lancelot_capabilities,      0},
    {"_sircovid_dust_lancelot_gpu_info",          (DL_FUNC) &_sircovid_dust_lancelot_gpu_info,          0},
    {NULL, NULL, 0}
};
}

extern "C" void R_init_sircovid(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
