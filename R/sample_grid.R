##' Take a grid search produced by \code{\link{scan_beta_date}} and 
##' sample \code{n_sample_pairs} from the parameter grid uses based
##' on their probability. For each parameter pair chosen, run particle
##' filter with \code{num_particles} and sample 1 trajectory
##' 
##' @title Sample Grid Scan
##' 
##' @param scan_results Output of \code{\link{scan_beta_date}}.
##' 
##' @param n_sample_pairs Number of parameter pairs to be sampled. This will 
##'   determine how many trajectories are returned. Integer. Default = 10. This 
##'   will determine how many trajectories are returned. 
##'   
##' @param n_particles Number of particles. Positive Integer. Default = 100
##' 
##' @param forecast_days Number of days being forecast. Default = 0
##'   
##' @return \code{\link{list}}. First element (trajectories) is a 3 
##'   dimensional array of trajectories (time, state, tranjectories). Second 
##'   element (param_grid) is the parameters chosen when  sampling from the 
##'   \code{scan_results} grid and the third dimension (inputs) is a list of
##'   model inputs. 
##' 
##' @export
##' @import furrr
##' @importFrom utils tail
sample_grid_scan <- function(scan_results,
                             n_sample_pairs = 10, 
                             n_particles = 100, 
                             forecast_days = 0) {
  
  # checks on args
  assert_is(scan_results, "sircovid_scan")
  assert_pos_int(n_sample_pairs)
  assert_pos_int(n_particles)
  assert_pos_int(forecast_days)
  
  # grab the pobability matrix
  prob <- scan_results$renorm_mat_LL
  nr <- nrow(prob)
  nc <- ncol(prob)
  
  # construct what the grid of beta and start values that 
  # correspond to the z axis matrix
  x_grid <- matrix(scan_results$x, nrow = nr, ncol = nc)
  y_grid <- matrix(as.character(scan_results$y), nrow = nr, 
                   ncol = nc, byrow = TRUE)
  
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
  ## Sample one particle
  traces <- furrr::future_pmap(
    .l = param_grid, .f = beta_date_particle_filter,
    generate_beta_func = scan_results$generate_beta_func,
    model_params = model_params, data = data, 
    pars_obs = pars_obs, n_particles = n_particles,
    forecast_days = forecast_days, 
    save_particles = TRUE, return = "sample"
  )

  # collapse into an array of trajectories
  # the trajectories are different lengths in terms of dates
  # so we will fill the arrays with NAs where needed
  num_rows <- unlist(lapply(traces, nrow))
  max_rows <- max(num_rows)
  seq_max <- seq_len(max_rows)
  max_date_names <- rownames(traces[[which.max(unlist(lapply(traces, nrow)))]])
  
  trajectories <- array(NA, 
                        dim = c(max_rows, ncol(traces[[1]]), length(traces)),
                        dimnames = list(max_date_names, NULL, NULL))
  
  # fill the tail of the array slice
  # This is so that the end of the trajectories array is populated, 
  # and the start is padded with NA if it's shorter than the max. 
  for (i in seq_len(length(traces))){
    trajectories[tail(seq_max, nrow(traces[[i]])), , i] <- traces[[i]]
  }

  # combine and return
  res <- list("trajectories" = trajectories,
              "param_grid" = param_grid,
              inputs = list(
                model_params = model_params,
                pars_obs = pars_obs,
                data = data,
                model = sircovid(params = model_params)))
  
  class(res) <- "sample_grid_search"
  
  return(res)
  
}
