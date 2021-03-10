#include <cpp11.hpp>
[[cpp11::register]]
SEXP dust_carehomes_alloc(cpp11::list r_pars, bool pars_multi, size_t step,
                         cpp11::sexp r_n_particles, size_t n_threads,
                         cpp11::sexp r_seed, cpp11::sexp device_id);

[[cpp11::register]]
SEXP dust_carehomes_run(SEXP ptr, size_t step_end, bool device);

[[cpp11::register]]
SEXP dust_carehomes_simulate(SEXP ptr, cpp11::sexp step_end);

[[cpp11::register]]
SEXP dust_carehomes_set_index(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_carehomes_set_state(SEXP ptr, SEXP r_state, SEXP r_step);

[[cpp11::register]]
SEXP dust_carehomes_reset(SEXP ptr, cpp11::list r_pars, size_t step);

[[cpp11::register]]
SEXP dust_carehomes_state(SEXP ptr, SEXP r_index);

[[cpp11::register]]
size_t dust_carehomes_step(SEXP ptr);

[[cpp11::register]]
void dust_carehomes_reorder(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_carehomes_resample(SEXP ptr, cpp11::doubles r_weights);

[[cpp11::register]]
SEXP dust_carehomes_set_pars(SEXP ptr, cpp11::list r_pars);

[[cpp11::register]]
SEXP dust_carehomes_rng_state(SEXP ptr, bool first_only);

[[cpp11::register]]
SEXP dust_carehomes_set_rng_state(SEXP ptr, cpp11::raws rng_state);

[[cpp11::register]]
SEXP dust_carehomes_set_data(SEXP ptr, cpp11::list data);

[[cpp11::register]]
SEXP dust_carehomes_compare_data(SEXP ptr);

[[cpp11::register]]
SEXP dust_carehomes_filter(SEXP ptr, bool save_history);

[[cpp11::register]]
cpp11::sexp dust_carehomes_capabilities();

[[cpp11::register]]
void dust_carehomes_set_n_threads(SEXP ptr, int n_threads);

[[cpp11::register]]
int dust_carehomes_n_state(SEXP ptr);

[[cpp11::register]]
cpp11::sexp dust_carehomes_device_info();
