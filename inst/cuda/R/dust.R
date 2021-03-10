carehomes <- R6::R6Class(
  "dust",
  cloneable = FALSE,

  private = list(
    pars_ = NULL,
    pars_multi_ = NULL,
    index_ = NULL,
    info_ = NULL,
    n_threads_ = NULL,
    n_particles_ = NULL,
    n_particles_each_ = NULL,
    shape_ = NULL,
    ptr_ = NULL,
    device_id_ = NULL,
    param_ = list(G_D_transmission = list(has_default = FALSE, default_value = NULL,
    rank = 0, min = -Inf, max = Inf, integer = FALSE), ICU_transmission = list(
    has_default = FALSE, default_value = NULL, rank = 0, min = -Inf,
    max = Inf, integer = FALSE), I_A_transmission = list(has_default = FALSE,
    default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE),
    I_C_1_transmission = list(has_default = FALSE, default_value = NULL,
        rank = 0, min = -Inf, max = Inf, integer = FALSE), I_C_2_transmission = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), I_P_transmission = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), N_tot_15_64 = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), N_tot_all = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), N_tot_over25 = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), N_tot_react = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), beta_step = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), exp_noise = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), hosp_transmission = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), index_dose = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), k_A = list(has_default = FALSE,
        default_value = NULL, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), k_C_1 = list(has_default = FALSE, default_value = NULL,
        rank = 0, min = -Inf, max = Inf, integer = FALSE), k_C_2 = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), k_E = list(has_default = FALSE,
        default_value = NULL, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), k_G_D = list(has_default = FALSE, default_value = NULL,
        rank = 0, min = -Inf, max = Inf, integer = FALSE), k_H_D = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), k_H_R = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), k_ICU_D = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), k_ICU_W_D = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), k_ICU_W_R = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), k_ICU_pre = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), k_P = list(has_default = FALSE,
        default_value = NULL, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), k_PCR_pos = list(has_default = FALSE,
        default_value = NULL, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), k_PCR_pre = list(has_default = FALSE,
        default_value = NULL, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), k_W_D = list(has_default = FALSE, default_value = NULL,
        rank = 0, min = -Inf, max = Inf, integer = FALSE), k_W_R = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), k_sero_pos = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), kappa_ICU = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), kappa_admitted = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), kappa_all_admission = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), kappa_death = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), kappa_death_carehomes = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), kappa_death_comm = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), kappa_death_hosp = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), kappa_death_non_hosp = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), kappa_diagnoses = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), kappa_general = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), kappa_hosp = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), kappa_pillar2_cases = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), m = list(has_default = FALSE,
        default_value = NULL, rank = 2, min = -Inf, max = Inf,
        integer = FALSE), p_C = list(has_default = FALSE, default_value = NULL,
        rank = 1, min = -Inf, max = Inf, integer = FALSE), p_G_D_step = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), p_H_D_step = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), p_H_step = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), p_ICU_D_step = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), p_ICU_step = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), p_NC = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), p_W_D_step = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), p_sero_pos = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), p_star_step = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), phi_ICU = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), phi_admitted = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), phi_all_admission = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), phi_death_carehomes = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), phi_death_comm = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), phi_death_hosp = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), phi_diagnoses = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), phi_general = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), phi_hosp = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), phi_pillar2_cases = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), pillar2_sensitivity = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), pillar2_specificity = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), psi_G_D = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), psi_H = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), psi_H_D = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), psi_ICU = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), psi_ICU_D = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), psi_W_D = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), psi_star = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), react_sensitivity = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), react_specificity = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), rel_infectivity = list(
        has_default = FALSE, default_value = NULL, rank = 2,
        min = -Inf, max = Inf, integer = FALSE), rel_p_hosp_if_sympt = list(
        has_default = FALSE, default_value = NULL, rank = 2,
        min = -Inf, max = Inf, integer = FALSE), rel_p_sympt = list(
        has_default = FALSE, default_value = NULL, rank = 2,
        min = -Inf, max = Inf, integer = FALSE), rel_susceptibility = list(
        has_default = FALSE, default_value = NULL, rank = 2,
        min = -Inf, max = Inf, integer = FALSE), rho_pillar2_tests = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), sero_sensitivity = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), sero_specificity = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), steps_per_day = list(
        has_default = FALSE, default_value = NULL, rank = 0,
        min = -Inf, max = Inf, integer = FALSE), strain_seed_step = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), strain_transmission = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), vaccine_dose_step = list(
        has_default = FALSE, default_value = NULL, rank = 3,
        min = -Inf, max = Inf, integer = FALSE), vaccine_progression_rate_base = list(
        has_default = FALSE, default_value = NULL, rank = 2,
        min = -Inf, max = Inf, integer = FALSE), waning_rate = list(
        has_default = FALSE, default_value = NULL, rank = 1,
        min = -Inf, max = Inf, integer = FALSE), gamma_A = list(
        has_default = TRUE, default_value = 0.1, rank = 0, min = -Inf,
        max = Inf, integer = FALSE), gamma_C_1 = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_C_2 = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_E = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_G_D = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_H_D = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_H_R = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_ICU_D = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_ICU_W_D = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_ICU_W_R = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_ICU_pre = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_P = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_PCR_pos = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_PCR_pre = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_U = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_W_D = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_W_R = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_sero_pos = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_sero_pre_1 = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), gamma_sero_pre_2 = list(has_default = TRUE,
        default_value = 0.1, rank = 0, min = -Inf, max = Inf,
        integer = FALSE), model_pcr_and_serology_user = list(
        has_default = TRUE, default_value = 1L, rank = 0, min = -Inf,
        max = Inf, integer = FALSE), p_sero_pre_1 = list(has_default = TRUE,
        default_value = 0.5, rank = 0, min = -Inf, max = Inf,
        integer = FALSE))
  ),

  public = list(
    ##' @description
    ##' Create a new model. Note that the behaviour of this object
    ##' created by this function will change considerably based on
    ##' whether the `pars_multi` argument is `TRUE`. If not (the
    ##' default) then we create `n_particles` which all share the same
    ##' parameters as specified by the `pars` argument. If `pars_multi`
    ##' is `TRUE` then `pars` must be an unnamed list, and each element
    ##' of it represents a different set of parameters. We will
    ##' create `length(pars)` *sets* of `n_particles` particles which
    ##' will be simulated together. These particles must have the same
    ##' dimension - that is, they must correspond to model state that
    ##' is the same size.
    ##'
    ##' @param pars Data to initialise your model with; a `list`
    ##' object, but the required elements will depend on the details of
    ##' your model. If `pars_multi` is `TRUE`, then this must be an
    ##' *unnamed* list of `pars` objects (see Details).
    ##'
    ##' @param step Initial step - must be nonnegative
    ##'
    ##' @param n_particles Number of particles to create - must be at
    ##' least 1
    ##'
    ##' @param n_threads Number of OMP threads to use, if `dust` and
    ##' your model were compiled with OMP support (details to come).
    ##' `n_particles` should be a multiple of `n_threads` (e.g., if you use 8
    ##' threads, then you should have 8, 16, 24, etc particles). However, this
    ##' is not compulsory.
    ##'
    ##' @param seed The seed to use for the random number generator. Can
    ##' be a positive integer, `NULL` (initialise with R's random number
    ##' generator) or a `raw` vector of a length that is a multiple of
    ##' 32 to directly initialise the generator (e..g., from the
    ##' [`dust`] object's `$rng_state()` method).
    ##'
    ##' @param pars_multi Logical, indicating if `pars` should be
    ##' interpreted as a set of different initialisations, and that we
    ##' should prepare `n_particles * length(pars)` particles for
    ##' simulation. This has an effect on many of the other methods of
    ##' the object.
    ##'
    ##' @param device_id Integer, indicating the device to use, where the
    ##' model has GPU support. If not given, then the default value of
    ##' `NULL` will fall back on the first found device if any are
    ##' available. An error is thrown if the device id given is larger
    ##' than those reported to be available (note that CUDA numbers devices
    ##' from 0, so that '0' is the first device, and so on). Negative values
    ##' disable the use of a device.' See the method `$device_info()` for
    ##' available device ids; this can be called before object creation as
    ##' `carehomes$public_methods$device_info()`
    initialize = function(pars, step, n_particles, n_threads = 1L,
                          seed = NULL, pars_multi = FALSE, device_id = NULL) {
      res <- dust_carehomes_alloc(pars, pars_multi, step, n_particles,
                        n_threads, seed, device_id)
      private$pars_ <- pars
      private$pars_multi_ <- pars_multi
      private$n_threads_ <- n_threads
      private$ptr_ <- res[[1L]]
      private$info_ <- res[[2L]]
      private$shape_ <- res[[3L]]
      private$device_id_ <- res[[4L]]
      private$n_particles_ <- prod(private$shape_)
      if (pars_multi) {
        private$n_particles_each_ <- private$n_particles_ / length(pars)
      } else {
        private$n_particles_each_ <- private$n_particles_
      }
    },

    ##' @description
    ##' Returns friendly model name
    name = function() {
      "carehomes"
    },

    ##' @description
    ##' Returns parameter information, if provided by the model. This
    ##' describes the contents of pars passed to the constructor or to
    ##' `reset` as the `pars` argument, and the details depend on the model.
    param = function() {
      private$param_
    },

    ##' @description
    ##' Run the model up to a point in time, returning the filtered state
    ##' at that point.
    ##'
    ##' @param step_end Step to run to (if less than or equal to the current
    ##' step(), silently nothing will happen)
    ##'
    ##' @param device **Experimental!**: This argument may allow running on
    ##' a GPU once support is finished, if the model supports it, and if
    ##' the model is compiled appropriately (and assuming you have a
    ##' suitable GPU). At present it exists for testing and will run
    ##' slower than running with `device = TRUE`. The interpretation of
    ##' this argument will likely change to allow selecting the GPU on
    ##' systems with more than one. In short, please leave this argument
    ##' alone unless you're developing dust.
    run = function(step_end, device = FALSE) {
      m <- dust_carehomes_run(private$ptr_, step_end, device)
      rownames(m) <- names(private$index_)
      m
    },

    ##' @description
    ##' Iterate all particles forward in time over a series of steps,
    ##' collecting output as they go. This is a helper around `$run()`
    ##' where you want to run to a series of points in time and save
    ##' output. The returned object will be filtered by your active index,
    ##' so that it has shape (`n_state` x `n_particles` x `length(step_end)`)
    ##' for single-parameter objects, and (`n_state` x `n_particles` x
    ##' `n_pars` x `length(step_end)`) for multiparameter objects. Note that
    ##' this method is very similar to `$run()` except that the rank of
    ##' the returned array is one less. For a scalar `step_end` you would
    ##' ordinarily want to use `$run()` but the resulting numbers would
    ##' be identical.
    ##'
    ##' @param step_end A vector of time points that the simulation should
    ##'   report output at. This the first time must be at least the same
    ##'   as the current time, and every subsequent time must be equal or
    ##'   greater than those before it (ties are allowed though probably
    ##'   not wanted).
    simulate = function(step_end) {
      m <- dust_carehomes_simulate(private$ptr_, step_end)
      rownames(m) <- names(private$index_)
      m
    },

    ##' @description
    ##' Set the "index" vector that is used to return a subset of pars
    ##' after using `run()`. If this is not used then `run()` returns
    ##' all elements in your state vector, which may be excessive and slower
    ##' than necessary. This method must be called after any
    ##' call to `reset()` as `reset()` may change the size of the state
    ##' and that will invalidate the index.
    ##'
    ##' @param index The index vector - must be an integer vector with
    ##' elements between 1 and the length of the state (this will be
    ##' validated, and an error thrown if an invalid index is given).
    set_index = function(index) {
      dust_carehomes_set_index(private$ptr_, index)
      private$index_ <- index
      invisible()
    },

    ##' @description
    ##' Returns the `index` as set by `$set_index`
    index = function() {
      private$index_
    },

    ##' @description
    ##' Returns the number of threads that the model was constructed with
    n_threads = function() {
      private$n_threads_
    },

    ##' @description
    ##' Returns the length of the per-particle state
    n_state = function() {
      dust_carehomes_n_state(private$ptr_)
    },

    ##' @description
    ##' Returns the number of particles
    n_particles = function() {
      private$n_particles_
    },

    ##' @description
    ##' Returns the number of particles per parameter set
    n_particles_each = function() {
      private$n_particles_each_
    },

    ##' @description
    ##' Returns the shape of the particles
    shape = function() {
      private$shape_
    },

    ##' @description
    ##' Set the "state" vector for all particles, overriding whatever your
    ##' models `initial()` method provides.
    ##'
    ##' @param state The state vector - can be either a numeric vector with the
    ##' same length as the model's current state (in which case the same
    ##' state is applied to all particles), or a numeric matrix with as
    ##' many rows as your model's state and as many columns as you have
    ##' particles (in which case you can set a number of different starting
    ##' states at once).
    ##'
    ##' @param step If not `NULL`, then this sets the initial step. If this
    ##' is a vector (with the same length as the number of particles), then
    ##' particles are started from different initial steps and run up to the
    ##' largest step given (i.e., `max(step)`)
    set_state = function(state, step = NULL) {
      dust_carehomes_set_state(private$ptr_, state, step)
      invisible()
    },

    ##' @description
    ##' Reset the model while preserving the random number stream state
    ##'
    ##' @param pars New pars for the model (see constructor)
    ##' @param step New initial step for the model (see constructor)
    reset = function(pars, step) {
      private$info_ <- dust_carehomes_reset(private$ptr_, pars, step)
      private$pars_ <- pars
      invisible()
    },

    ##' @description
    ##' Set the 'pars' element in a dust object while holding model state,
    ##' index, etc constant. In contrast to `$reset`, the old state must
    ##' be compatible with the new one (e.g., don't change model size), and
    ##' the index will remain valid.
    ##'
    ##' @param pars New pars for the model (see constructor)
    set_pars = function(pars) {
      private$info_ <- dust_carehomes_set_pars(private$ptr_, pars)
      private$pars_ <- pars
    },

    ##' @description
    ##' Return full model state
    ##' @param index Optional index to select state using
    state = function(index = NULL) {
      m <- dust_carehomes_state(private$ptr_, index)
      rownames(m) <- names(index)
      m
    },

    ##' @description
    ##' Return current model step
    step = function() {
      dust_carehomes_step(private$ptr_)
    },

    ##' @description
    ##' Reorder particles.
    ##' @param index An integer vector, with values between 1 and n_particles,
    ##' indicating the index of the current particles that new particles should
    ##' take.
    reorder = function(index) {
      storage.mode(index) <- "integer"
      dust_carehomes_reorder(private$ptr_, index)
      invisible()
    },

    ##' @description
    ##' Resample particles according to some weight.
    ##'
    ##' @param weights A numeric vector representing particle weights.
    ##' For a "multi-parameter" dust object this should be be a matrix
    ##' with the number of rows being the number of particles per
    ##' parameter set and the number of columns being the number of
    ##' parameter sets.
    ##' long as all particles or be a matrix.
    resample = function(weights) {
      invisible(dust_carehomes_resample(private$ptr_, weights))
    },

    ##' @description
    ##' Returns information about the pars that your model was created with.
    ##' Only returns non-NULL if the model provides a `dust_info` template
    ##' specialisation.
    info = function() {
      private$info_
    },

    ##' @description
    ##' Returns the `pars` object that your model was constructed with.
    pars = function() {
      private$pars_
    },

    ##' @description
    ##' Returns the state of the random number generator. This returns a
    ##' raw vector of length 32 * n_particles. This can be useful for
    ##' debugging or for initialising other dust objects.
    ##'
    ##' @param first_only Logical, indicating if we should return only the
    ##' *first* particle's random number state. If `FALSE` (the default)
    ##' all particles states are returned, being 32 bytes per particle.
    ##' If `TRUE` then we take just the first particle's state, which
    ##' will be a total of 32 bytes. Both forms are suitable for seeding
    ##' a new [`dust`] object as the shorter version will be used for
    ##' the first particle, followed by jumps for each subsequent particle.
    rng_state = function(first_only = FALSE) {
      dust_carehomes_rng_state(private$ptr_, first_only)
    },

    ##' @description Set the random number state for this model. This
    ##' replaces the RNG state that the model is using with a state of
    ##' your choosing, saved out from a different model object. This method
    ##' is designed to support advanced use cases where it is easier to
    ##' manipulate the state of the random number generator than the
    ##' internal state of the dust object.
    ##'
    ##' @param rng_state A random number state, as saved out by the
    ##' `$rng_state()` method. Note that unlike `seed` as passed to the
    ##' constructor, this *must* be a raw vector of the expected length.
    set_rng_state = function(rng_state) {
      dust_carehomes_set_rng_state(private$ptr_, rng_state)
      invisible()
    },

    ##' @description
    ##' Returns a logical, indicating if this model was compiled with
    ##' "OpenMP" support, in which case it will react to the `n_threads`
    ##' argument passed to the constructor. This method can also be used
    ##' as a static method by running it directly
    ##' as `carehomes$public_methods$has_openmp()`
    has_openmp = function() {
      dust_carehomes_capabilities()[["openmp"]]
    },

    ##' @description
    ##' Returns a logical, indicating if this model was compiled with
    ##' "CUDA" support, in which case it will react to the `device`
    ##' argument passed to the run method. This method can also be used
    ##' as a static method by running it directly
    ##' as `carehomes$public_methods$has_cuda()`
    has_cuda = function() {
      dust_carehomes_capabilities()[["cuda"]]
    },

    ##' @description
    ##' Returns the number of distinct pars elements required. This is `0`
    ##' where the object was initialised with `pars_multi = FALSE` and
    ##' an integer otherwise.  For multi-pars dust objects, Where `pars`
    ##' is accepted, you must provide an unnamed list of length `$n_pars()`.
    n_pars = function() {
      if (private$pars_multi_) length(private$pars_) else 0L
    },

    ##' @description
    ##' Change the number of threads that the dust object will use. Your
    ##' model must be compiled with "OpenMP" support for this to have an
    ##' effect. Returns (invisibly) the previous value.
    ##'
    ##' @param n_threads The new number of threads to use. You may want to
    ##'   wrap this argument in [dust::dust_openmp_threads()] in order to
    ##'   verify that you can actually use the number of threads
    ##'   requested (based on environment variables and OpenMP support).
    set_n_threads = function(n_threads) {
      prev <- private$n_threads_
      dust_carehomes_set_n_threads(private$ptr_, n_threads)
      private$n_threads_ <- n_threads
      invisible(prev)
    },

    ##' @description
    ##' Returns a logical, indicating if this model was compiled with
    ##' "compare" support, in which case the `set_data` and `compare_data`
    ##' methods are available (otherwise these methods will error). This
    ##' method can also be used as a static method by running it directly
    ##' as `carehomes$public_methods$has_compare()`
    has_compare = function() {
      dust_carehomes_capabilities()[["compare"]]
    },

    ##' @description
    ##' Set "data" into the model for use with the `$compare_data()` method.
    ##' This is not supported by all models, depending on if they define a
    ##' `data_t` type.  See [dust::dust_data()] for a helper function to
    ##' construct suitable data and a description of the required format. You
    ##' will probably want to use that here, and definitely if using multiple
    ##' parameter sets.
    ##'
    ##' @param data A list of data to set.
    set_data = function(data) {
      dust_carehomes_set_data(private$ptr_, data)
    },

    ##' @description
    ##' Compare the current model state against the data as set by
    ##' `set_data`. If there is no data set, or no data corresponding to
    ##' the current time then `NULL` is returned. Otherwise a numeric vector
    ##' the same length as the number of particles is returned. If model's
    ##' underlying `compare_data` function is stochastic, then each call to
    ##' this function may be result in a different answer.
    compare_data = function() {
      dust_carehomes_compare_data(private$ptr_)
    },

    ##' @description
    ##' Run a particle filter. The interface here will change a lot over the
    ##' next few versions. You *must* `$reset()` the filter before using
    ##' this method to get sensible values. We will tinker with this in
    ##' future versions to allow things like partial runs.
    ##'
    ##' @param save_history Logical, indicating if the filtered particle
    ##' trajectories should be saved. If `TRUE` then the `history` element
    ##' will be a 3d array (state x particles x time) containing the state
    ##' values, selected according to the index set with `$set_index()`.
    ##' If you have a multi-parameter dust object this will be a 4d array
    ##' (state x particles x parameter x time).
    filter = function(save_history = FALSE) {
      dust_carehomes_filter(private$ptr_, save_history)
    },

    ##' @description
    ##' **Experimental!** Return information about GPU devices, if the model
    ##' has been compiled with CUDA/GPU support. This can be called as a
    ##' static method by running `carehomes$public_methods$device_info()`
    device_info = function() {
      ret <- dust_carehomes_device_info()
      if (ret$has_cuda && exists("private", inherits = FALSE)) {
        ret$device_id <- private$device_id_
      }
      ret
    }
  ))
class(carehomes) <- c("dust_generator", class(carehomes))
