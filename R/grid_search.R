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
##'   \code{generate_parameters()}. If NULL, uses defaults as
##'   in unit tests.
##'
##' @param pars_obs list of parameters to use for the comparison function.
##'   
##' @param n_particles Number of particles. Positive Integer. Default = 100
##' 
##' @return List of beta and start date grid values, and
##'   normalised probabilities at each point
##' 
##' @export
##' @import furrr
scan_beta_date <- function(
  min_beta, 
  max_beta, 
  beta_step,
  first_start_date, 
  last_start_date, 
  day_step,
  data,
  sircovid_model = basic_model(),
  model_params = NULL,
  pars_obs = NULL,
  n_particles = 100) {
  #
  # Set up parameter space to scan
  #
  beta_1D <- seq(min_beta, max_beta, beta_step)
  date_list <- seq(as.Date(first_start_date), as.Date(last_start_date), day_step)
  param_grid <- expand.grid(beta = beta_1D, start_date = date_list)
  
  #
  # Set up calls to simulator runs
  #
  if (is.null(model_params)) {
    time_steps_per_day <- 4
    model_params <- generate_parameters(
      sircovid_model = sircovid_model,
      transmission_model = "POLYMOD",
      beta = 0.1,
      beta_times = "2020-01-01",
      trans_profile = 1,
      trans_increase = 1,
      dt = 1/time_steps_per_day
    )
  } else {
    if (length(model_params$beta_y) > 1) {
      stop("Set beta variation through generate_beta, not model_params")
    }
  }

  if (is.null(pars_obs)) {
    pars_obs <- list(phi_general = 0.95,
                     k_general = 2,
                     phi_ICU = 0.95,
                     k_ICU = 2,
                     phi_death = 926 / 1019,
                     k_death = 2,
                     exp_noise = 1e6)
  }
  
  #
  # Multi-core futures with furrr (parallel purrr)
  #
  ## Particle filter outputs, extracting log-likelihoods
  pf_run_ll <- furrr::future_pmap_dbl(
    .l = param_grid, .f = beta_date_particle_filter,
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
    byrow = FALSE)

  # Exponentiate elements and normalise to 1 to get probabilities
  prob_matrix <- exp(mat_log_ll)
  renorm_mat_LL <- prob_matrix/sum(prob_matrix)
  
  results <- list(x = beta_1D, 
                  y = date_list,
                  mat_log_ll = mat_log_ll,
                  renorm_mat_LL = renorm_mat_LL,
                  inputs = list(
                    model = sircovid_model,
                    model_params = model_params,
                    pars_obs = pars_obs,
                    data = data))
  
  class(results) <- "sircovid_scan"
  results
}

##' @export
plot.sircovid_scan <- function(x, ..., what = "likelihood") {
  if (what == "likelihood") {
    graphics::image(x=x$x, y=x$y, z=x$mat_log_ll,
                    xlab="beta", ylab="Start date", main = "Log-likelihood")
  } else if (what == "probability") {
    graphics::image(x=x$x, y=x$y, z=x$renorm_mat_LL,
                    xlab="beta", ylab="Start date", main = "Probability")
  }
}

##' @export
plot.sample_grid_search <- function(x, ..., what = "ICU") {

  idx <- odin_index(x$inputs$model$odin_model(user = x$inputs$model_params,
                                              unused_user_action = "message"))

  # what are we plotting
  if (what == "ICU") {
    index <- c(idx$I_ICU_D,idx$I_ICU_R) - 1L
    ylab <- "ICU"
    particles <- vapply(seq_len(dim(x$trajectories)[3]), function(y) {
      rowSums(x$trajectories[,index,y], na.rm = TRUE)}, 
      FUN.VALUE = numeric(dim(x$trajectories)[1]))
    plot_particles(particles, ylab = ylab)
    points(as.Date(x$inputs$data$date), x$inputs$data$itu / x$inputs$pars_obs$phi_ICU, pch = 19)
    
  } else if (what == "general") {
    
    index <- c(idx$I_triage,idx$I_hosp_R,idx$I_hosp_D,idx$R_stepdown) - 1L
    ylab <- "General beds"
    particles <- vapply(seq_len(dim(x$trajectories)[3]), function(y) {
      rowSums(x$trajectories[,index,y], na.rm = TRUE)}, 
      FUN.VALUE = numeric(dim(x$trajectories)[1]))
    plot_particles(particles, ylab = ylab)
    points(as.Date(x$inputs$data$date), x$inputs$data$general / x$inputs$pars_obs$phi_general, pch = 19)
    
  }
  
  else if(what == "deaths") {
    
    index <- c(idx$D) - 1L
    ylab <- "Deaths"
    particles <- vapply(seq_len(dim(x$trajectories)[3]), function(y) {
      out <- c(0,diff(rowSums(x$trajectories[,index,y], na.rm = TRUE)))
      names(out)[1] <- rownames(x$trajectories)[1]
      out}, 
      FUN.VALUE = numeric(dim(x$trajectories)[1]))
    plot_particles(particles, ylab = ylab)
    points(as.Date(x$inputs$data$date), 
           x$inputs$data$deaths/ x$inputs$pars_obs$phi_death, pch = 19)
    
  } else {
    
    stop("Requested what must be one of 'ICU' or 'Deaths'")
    
  } 
  
}

##' Particle filter outputs
##' 
##' Helper function to run the particle filter with a
##' new beta and start date
##' 
##' @noRd
beta_date_particle_filter <- function(beta, start_date,
                                      sircovid_model,
                                      model_params, data, 
                                      pars_obs, n_particles,
                                      forecast_days = 0,
                                      save_particles = FALSE,
                                      return = "full") {
  # Edit beta in parameters
  new_beta <- sircovid_model$generate_beta(beta, start_date)
  beta_t <- normalise_beta(new_beta$beta_times, model_params$dt)
  
  model_params$beta_y <- new_beta$beta
  model_params$beta_t <- beta_t

  X <- run_particle_filter(data, sircovid_model, model_params, start_date, pars_obs,
                           n_particles, forecast_days, save_particles, return)
  
  X
}
