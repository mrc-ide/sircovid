##' Run a grid search of the particle filter over beta and start date.
##' This is parallelised, first run \code{plan(multiprocess)} to set
##' this up.
##'
##' @title Grid search of beta and start date
##'
##' @param min_beta Minimum value of beta in the search
##'
##' @param max_beta Maximum value of beta in the search
##'
##' @param beta_step Step to increment beta between min and max
##'
##' @param first_start_date Earliest start date as 'yyyy-mm-dd'
##'
##' @param last_start_date Latest start date as 'yyyy-mm-dd'
##'
##' @param day_step Step to increment date in days
##'
##' @param data Hosptial data to fit to. See \code{example.csv}
##'   and \code{particle_filter_data()}
##'
##' @param sircovid_model An odin model generator and comparison function.
##'
##' @param model_params Model parameters, from a call to
##'   \code{generate_parameters()}
##'
##' @param pars_obs list of parameters to use for the comparison function.
##'
##' @param n_particles Number of particles. Positive Integer. Default = 100
##'
##' @param scale_prior Set to use a gamma prior on beta, and report the
##'   posterior rather than the likelihood. Sets scale of gamma prior.
##'
##' @param shape_prior Set to use a gamma prior on beta, and report the
##'   posterior rather than the likelihood. Sets shape of gamma prior.
##'
##' @param tolerance The smallest difference from 0 that is acceptable
##' in the probability matrix.
##' @return List of beta and start date grid values, and
##'   normalised probabilities at each point
##'
##' @export
##' @import furrr
##' @importFrom stats dgamma
scan_beta_date <- function(
                           min_beta,
                           max_beta,
                           beta_step,
                           first_start_date,
                           last_start_date,
                           day_step,
                           data,
                           sircovid_model = basic_model(),
                           model_params,
                           pars_obs,
                           n_particles = 100,
                           scale_prior = NULL,
                           shape_prior = NULL,
                           tolerance = 1e-10) {

  # Parameter checks
  if (!is.null(scale_prior) || !is.null(shape_prior)) {
    if (!is.numeric(scale_prior) || !is.numeric(shape_prior)) {
      stop("If provided, both scale_prior and shape_prior must both be numeric")
    }
  }

  #
  # Set up parameter space to scan
  #
  beta_1D <- seq(min_beta, max_beta, beta_step)
  date_list <- seq(first_start_date, last_start_date, day_step)
  param_grid <- expand.grid(beta = beta_1D, start_date = date_list)

  #
  # Set up calls to simulator runs
  #
    if (length(model_params$beta_y) > 1) {
      stop("Set beta variation through generate_beta_func in sircovid_model, not model_params")
    }


  #
  # Multi-core futures with furrr (parallel purrr)
  #
  ## Particle filter outputs, extracting log-likelihoods
  pf_run_ll <- furrr::future_pmap_dbl(
    .l = param_grid, .f = run_grid_particle_filter,
    sircovid_model = sircovid_model, model_params = model_params, data = data,
    pars_obs = pars_obs, n_particles = n_particles,
    forecast_days = 0, save_particles = FALSE, return = "ll"
  )

  ## Construct a matrix with start_date as columns, and beta as rows
  ## order of return is set by order passed to expand.grid, above
  ## Returned column-major (down columns of varying beta) - set byrow = FALSE
  mat_log_ll <- matrix(
    pf_run_ll,
    nrow = length(beta_1D),
    ncol = length(date_list),
    byrow = FALSE
  )

  # Exponentiate elements and normalise to 1 to get probabilities
  prob_matrix <- exp(mat_log_ll)

  renorm_mat_LL <- prob_matrix / sum(prob_matrix)


  # Apply the prior, if provided
  if (!is.null(shape_prior) && !is.null(scale_prior)) {
    log_prior <- matrix(rep(dgamma(beta_1D,
      shape = shape_prior,
      scale = scale_prior,
      log = TRUE
    ), length(date_list)),
    ncol = length(date_list)
    )
    mat_log_ll <- mat_log_ll + log_prior
    exp_mat <- exp(mat_log_ll - max(mat_log_ll))
    renorm_mat_LL <- exp_mat / sum(exp_mat)
  }

  ## Check if the edges of the matrix in each dimension are close
  ## enough to 0.
  ## The first and last rows and the first and last columns
  ## Or in case of multidimensional arrays, edges in each dimension.

  close_enough <- zero_boundary(renorm_mat_LL, tolerance = tolerance)

  if (!close_enough) {
    warning("Edges of the probability matrix are not close enough to 0.")
  }


  results <- list(
    x = beta_1D,
    y = date_list,
    mat_log_ll = mat_log_ll,
    renorm_mat_LL = renorm_mat_LL,
    inputs = list(
      model = sircovid_model,
      model_params = model_params,
      pars_obs = pars_obs,
      data = data
    )
  )

  class(results) <- "sircovid_scan"
  results
}

##' @export
plot.sircovid_scan <- function(x, ..., what = "likelihood", title = NULL) {
  if (what == "likelihood") {
    graphics::image(
      x = x$x, y = x$y, z = x$mat_log_ll,
      xlab = "beta", ylab = "Start date", main = title
    )
  } else if (what == "probability") {
    graphics::image(
      x = x$x, y = x$y, z = x$renorm_mat_LL,
      xlab = "beta", ylab = "Start date", main = title
    )
  }
}

##' Particle filter outputs
##'
##' Helper function to run the particle filter with a
##' new beta and start date
##'
##' @noRd
run_grid_particle_filter <- function(beta, start_date,
                                      sircovid_model,
                                      model_params, data,
                                      pars_obs, n_particles,
                                      forecast_days = 0,
                                      save_particles = FALSE,
                                      return = "full") {
  # Edit beta in parameters
  new_beta <- update_beta(sircovid_model, 
                          beta_start = beta, 
                          beta_end = NULL, 
                          beta_pl = NULL,
                          start_date,
                          model_params$dt)
  model_params[names(new_beta)] <- new_beta

  X <- run_particle_filter(
    data, sircovid_model, model_params, start_date, pars_obs,
    pars_seeding = NULL, n_particles, forecast_days, save_particles, return
  )

  X
}
