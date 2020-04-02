##' Run a grid search of the particle filter over beta and start date
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
##' @param model_params Model parameters, from a call to
##'   \code{generate_parameters()}. If NULL, uses defaults as
##'   in unit tests.
##'   
##' @param pars_obs list of parameters to use in comparison
##'   with \code{compare_icu}.
##'   
##' @return List of beta and start date grid values, and
##'   normalised probabilities at each point
##' 
##' @export
##' @import furrr
##' @import future
scan_beta_date <- function(
  min_beta, 
  max_beta, 
  beta_step, 
  first_start_date, 
  last_start_date, 
  day_step,
  data, 
  model_params = NULL,
  pars_obs = list(phi_ICU = 0.95,
                  k_ICU = 2,
                  phi_death = 926 / 1019,
                  k_death = 2,
                  exp_noise = 1e6)) {
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
      transmission_model = "POLYMOD",
      progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 1, hosp = 2, ICU = 2, rec = 2),
      gammas = list(E = 1/2.5, asympt = 1/2.09, mild = 1/2.09, ILI = 1/4, hosp = 2/1, ICU = 2/5, rec = 2/5),
      beta = 0.1,
      beta_times = "2020-01-01",
      hosp_transmission = 0,
      ICU_transmission = 0,
      trans_profile = 1,
      trans_increase = 1,
      dt = 1/time_steps_per_day
    )
  } else {
    if (length(model_params$beta_y)) {
      stop("Time-varying beta not yet supported by grid search")
    }
  }
  
  #
  # Multi-core futures with furrr (parallel purrr)
  #
  ## Particle filter outputs
  plan(multiprocess)
  pf_run_outputs <- furrr::future_pmap(
    param_grid,
    function(beta, start_date) {
      # Edit beta in parameters
      model_params$beta_y <- beta
      model_params$beta_t <- 0
      
      # Objects used by particle filter run
      pf_data <- particle_filter_data(data, start_date, steps_per_day= 1/model_params$dt)
      model_run <- sircovid(params = model_params)
      compare_func <- compare_icu(model_run, pars_obs, pf_data)
      
      X <- particle_filter(pf_data, model_run, compare_func, n_particles=100)
    }
  )
  ## Extract Log-likelihoods
  log_ll <- furrr::future_map(pf_run_outputs, ~ .$log_likelihood) 
  
  ## Construct a matrix with start_date as columns, and beta as rows
  ## order of return is set by order passed to expand.grid, above
  ## Returned column-major (down columns of varying beta) - set byrow = FALSE
  mat_log_ll <- matrix(
    as.numeric(log_ll), 
    nrow = length(beta_1D),
    ncol = length(date_list),
    byrow = FALSE)

  # Exponentiate elements and normalise to 1 to get probabilities
  prob_matrix <- exp(mat_log_ll)
  renorm_mat_LL <- prob_matrix/sum(prob_matrix)
  
  # Display likelihood/probability surfaces across grid
  graphics::image(x=beta_1D, y=date_list, z=mat_log_ll,
        xlab="beta", ylab="Start date", main = "Log-likelihood")
  graphics::image(x=beta_1D, y=date_list, z=renorm_mat_LL,
        xlab="beta", ylab="Start date", main = "Probability")
  
  results <- list(beta_vals=param_grid$beta, 
                  start_dates=param_grid$start_date,
                  renorm_mat_LL=renorm_mat_LL)
  
  results
}
