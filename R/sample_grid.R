##' Take a grid search produced by \code{\link{scan_beta_date}} and 
##' sample \code{n_sample_pairs} from the parameter grid uses based
##' on their probaility. For each parameter pair chosen, run particle
##' filter with \code{num_particles} and sample 1 trajectory
##' 
##' @title Sample Grid Scan
##' 
##' @param scan_results Output of \code{\link{scan_beta_date}}.
##' 
##' @param n_sample_pairs Number of parameter pairs to be sampled. Integer.
##'   Default = 10.
##'   
##' @param n_particles Number of particles. Positive Integer. Default = 100
##'   
##' @return Array of trajectories. First dimension is time, second is the state
##'   and the third dimension is the index of the trajectory
##' 
##' @export
##' @import furrr
sample_grid_scan <- function(scan_results, 
                             n_sample_pairs = 10, 
                             n_particles = 100) {
  
  # checks on args
  assert_custom_class(scan_results, "sircovid_scan")
  assert_pos_int(n_sample_pairs)
  assert_pos_int(n_particles)
  
  # grab the pobability matrix
  prob <- scan_results$renorm_mat_LL
  nr <- nrow(prob)
  nc <- ncol(prob)
  
  # construct what the x and y dimensions look like
  x_grid <- matrix(scan_results$x, nrow = nr, ncol = nc)
  y_grid <- matrix(as.character(scan_results$y), nrow = nr, ncol = nc, byrow = TRUE)
  
  # draw which grid pairs are chosen
  pairs <- sample(x =  length(prob), size = n_sample_pairs, 
                  replace = TRUE, prob = prob)
  
  # what are related beta and dates
  beta <- x_grid[pairs]
  dates <- y_grid[pairs]
  
  # recreate parameters for re running
  param_grid <- data.frame("beta" = beta, "start_date" = dates)
  model_params <- scan_results$inputs$model_params
  pars_obs <- scan_results$inputs$pars_obs
  data <- scan_results$inputs$data
  
  
  # Multi-core futures with furrr (parallel purrr)
  
  ## Particle filter outputs
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
      
      X <- particle_filter(pf_data, model_run, compare_func, 
                           n_particles = n_particles, save_particles = TRUE)
    }
  )
  
  # sample 1 particle from each run
  traces <- lapply(pf_run_outputs, function(x) {
    x$states[,,sample(n_particles, 1)]
  })
  
  
  # collapse into an array of trajectories
  trajectories <- array(unlist(traces), 
                        dim = c(nrow(traces[[1]]), ncol(traces[[1]]), length(traces)),
                        dimnames = list(rownames(traces[[1]]), NULL, NULL))
  
  return(trajectories)
  
}