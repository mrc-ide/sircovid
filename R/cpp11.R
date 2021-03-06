# Generated by cpp11: do not edit by hand

dust_basic_alloc <- function(r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, device_id) {
  .Call(`_sircovid_dust_basic_alloc`, r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, device_id)
}

dust_basic_run <- function(ptr, step_end, device) {
  .Call(`_sircovid_dust_basic_run`, ptr, step_end, device)
}

dust_basic_simulate <- function(ptr, step_end, device) {
  .Call(`_sircovid_dust_basic_simulate`, ptr, step_end, device)
}

dust_basic_set_index <- function(ptr, r_index) {
  .Call(`_sircovid_dust_basic_set_index`, ptr, r_index)
}

dust_basic_set_state <- function(ptr, r_state, r_step) {
  .Call(`_sircovid_dust_basic_set_state`, ptr, r_state, r_step)
}

dust_basic_reset <- function(ptr, r_pars, step) {
  .Call(`_sircovid_dust_basic_reset`, ptr, r_pars, step)
}

dust_basic_state <- function(ptr, r_index) {
  .Call(`_sircovid_dust_basic_state`, ptr, r_index)
}

dust_basic_step <- function(ptr) {
  .Call(`_sircovid_dust_basic_step`, ptr)
}

dust_basic_reorder <- function(ptr, r_index) {
  invisible(.Call(`_sircovid_dust_basic_reorder`, ptr, r_index))
}

dust_basic_resample <- function(ptr, r_weights) {
  .Call(`_sircovid_dust_basic_resample`, ptr, r_weights)
}

dust_basic_set_pars <- function(ptr, r_pars) {
  .Call(`_sircovid_dust_basic_set_pars`, ptr, r_pars)
}

dust_basic_rng_state <- function(ptr, first_only, last_only) {
  .Call(`_sircovid_dust_basic_rng_state`, ptr, first_only, last_only)
}

dust_basic_set_rng_state <- function(ptr, rng_state) {
  .Call(`_sircovid_dust_basic_set_rng_state`, ptr, rng_state)
}

dust_basic_set_data <- function(ptr, data) {
  .Call(`_sircovid_dust_basic_set_data`, ptr, data)
}

dust_basic_compare_data <- function(ptr, device) {
  .Call(`_sircovid_dust_basic_compare_data`, ptr, device)
}

dust_basic_filter <- function(ptr, save_trajectories, step_snapshot, device) {
  .Call(`_sircovid_dust_basic_filter`, ptr, save_trajectories, step_snapshot, device)
}

dust_basic_capabilities <- function() {
  .Call(`_sircovid_dust_basic_capabilities`)
}

dust_basic_set_n_threads <- function(ptr, n_threads) {
  invisible(.Call(`_sircovid_dust_basic_set_n_threads`, ptr, n_threads))
}

dust_basic_n_state <- function(ptr) {
  .Call(`_sircovid_dust_basic_n_state`, ptr)
}

dust_basic_device_info <- function() {
  .Call(`_sircovid_dust_basic_device_info`)
}

dust_carehomes_alloc <- function(r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, device_id) {
  .Call(`_sircovid_dust_carehomes_alloc`, r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, device_id)
}

dust_carehomes_run <- function(ptr, step_end, device) {
  .Call(`_sircovid_dust_carehomes_run`, ptr, step_end, device)
}

dust_carehomes_simulate <- function(ptr, step_end, device) {
  .Call(`_sircovid_dust_carehomes_simulate`, ptr, step_end, device)
}

dust_carehomes_set_index <- function(ptr, r_index) {
  .Call(`_sircovid_dust_carehomes_set_index`, ptr, r_index)
}

dust_carehomes_set_state <- function(ptr, r_state, r_step) {
  .Call(`_sircovid_dust_carehomes_set_state`, ptr, r_state, r_step)
}

dust_carehomes_reset <- function(ptr, r_pars, step) {
  .Call(`_sircovid_dust_carehomes_reset`, ptr, r_pars, step)
}

dust_carehomes_state <- function(ptr, r_index) {
  .Call(`_sircovid_dust_carehomes_state`, ptr, r_index)
}

dust_carehomes_step <- function(ptr) {
  .Call(`_sircovid_dust_carehomes_step`, ptr)
}

dust_carehomes_reorder <- function(ptr, r_index) {
  invisible(.Call(`_sircovid_dust_carehomes_reorder`, ptr, r_index))
}

dust_carehomes_resample <- function(ptr, r_weights) {
  .Call(`_sircovid_dust_carehomes_resample`, ptr, r_weights)
}

dust_carehomes_set_pars <- function(ptr, r_pars) {
  .Call(`_sircovid_dust_carehomes_set_pars`, ptr, r_pars)
}

dust_carehomes_rng_state <- function(ptr, first_only, last_only) {
  .Call(`_sircovid_dust_carehomes_rng_state`, ptr, first_only, last_only)
}

dust_carehomes_set_rng_state <- function(ptr, rng_state) {
  .Call(`_sircovid_dust_carehomes_set_rng_state`, ptr, rng_state)
}

dust_carehomes_set_data <- function(ptr, data) {
  .Call(`_sircovid_dust_carehomes_set_data`, ptr, data)
}

dust_carehomes_compare_data <- function(ptr, device) {
  .Call(`_sircovid_dust_carehomes_compare_data`, ptr, device)
}

dust_carehomes_filter <- function(ptr, save_trajectories, step_snapshot, device) {
  .Call(`_sircovid_dust_carehomes_filter`, ptr, save_trajectories, step_snapshot, device)
}

dust_carehomes_capabilities <- function() {
  .Call(`_sircovid_dust_carehomes_capabilities`)
}

dust_carehomes_set_n_threads <- function(ptr, n_threads) {
  invisible(.Call(`_sircovid_dust_carehomes_set_n_threads`, ptr, n_threads))
}

dust_carehomes_n_state <- function(ptr) {
  .Call(`_sircovid_dust_carehomes_n_state`, ptr)
}

dust_carehomes_device_info <- function() {
  .Call(`_sircovid_dust_carehomes_device_info`)
}
