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
  sircovid_model <- scan_results$inputs$model
  model_params <- scan_results$inputs$model_params
  pars_obs <- scan_results$inputs$pars_obs
  data <- scan_results$inputs$data
  
  # Multi-core futures with furrr (parallel purrr)
  
  ## Particle filter outputs
  ## Sample one particle
  traces <- furrr::future_pmap(
    .l = param_grid, .f = run_grid_particle_filter,
    sircovid_model = sircovid_model, 
    model_params = model_params, data = data,
    pars_obs = pars_obs, n_particles = n_particles,
    forecast_days = forecast_days, 
    save_particles = TRUE, return = "sample")

  trajectories <- traces_to_trajectories(traces
  
  class(res) <- "sircovid_forecast"
  return(res)
  
}


##' Take a parameter search produced by \code{\link{pmcmc}} and 
##' sample \code{n_sample_pairs} from the parameter space
##' on their probability. For each parameter pair chosen, run particle
##' filter with \code{num_particles} and sample 1 trajectory
##' 
##' @title Sample pmcmc
##' 
##' @param mcmc_results Output of \code{\link{pmcmc}}.
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
sample_pmcmc <- function(mcmc_results,
                         burn_in = 101,
                         n_sample_pairs = 10, 
                         n_particles = 100, 
                         forecast_days = 0) {
  
  # checks on args
  assert_is(scan_results, "pmcmc_list")
  assert_pos_int(n_sample_pairs)
  assert_pos_int(n_particles)
  assert_pos_int(forecast_days)
  
  # discard 10% as burn-in
  pars_to_sample <- names(mcmc_results$inputs$pars$pars_init)
  chains <- create_master_chain(mcmc_results, burn_in)$results
  param_grid <- chains[sample.int(n = nrow(chains), 
                                  size  = n_sample_pairs, 
                                  replace = FALSE), pars_to_sample]
  
  forecasts <- function(sampled_pars) {
      pars <- as.list(sampled_pars)
      trace <- calc_loglikelihood(pars, 
              mcmc_results$inputs$data, 
              mcmc_results$inputs$sircovid_model, 
              mcmc_results$inputs$model_params,
              mcmc_results$inputs$steps_per_day, 
              mcmc_results$inputs$pars_obs, 
              n_particles,
              forecast_days,
              return = "full")
      trace
  }

  # Run forecasts in parallel
  traces <- furrr::future_pmap(
    .l = param_grid, 
    .f = forecasts)
  
  trajectories <- traces_to_trajectories(traces)
  
  class(sample_res) <- "sircovid_forecast"
  return(res)
  
}

##' Take a sampled grid search produced by \code{\link{sample_grid_scan}} and 
##' extract forecasts and quantiles in text format
##' 
##' @title Summarise grid forecast
##'
##' @param object Results from \code{\link{sample_grid_scan}}
##' 
##' @param ... other arguments to \code{summary()}
##'
##' @param what Output to summaries. "deaths", "icu" or "hosp"
##'
##' @param quantiles Quantiles to summarise forecast variance.
##' 
##' @export
summary.sircovid_forecast <- function(object, ...,
                                       what = "deaths",
                                       quantiles = seq(from=0.05, to=0.95, by=0.05)) {
  totals <- sum_over_compartments(object)[[what]]
  
  # Extract quantiles
  names(quantiles) <- sprintf("Quantile %s", quantiles)
  quantiles <- c(c(Value = 0.5), quantiles)
  
  i <- tail(seq_len(nrow(totals)), object$inputs$forecast_days)
  d <- round(t(apply(totals[i, ], 1, quantile, quantiles, names = FALSE)))
  colnames(d) <- names(quantiles)

  d
}

# Sum sampled model over compartments
sum_over_compartments <- function(sample_grid_res) {
  index <- odin_index(sample_grid_res$inputs$model$odin_model(
      user = sample_grid_res$inputs$model_params))

  ## filter out the ones that we actually have - this is unreasonably ugly:
  n <- vapply(index, min, numeric(1)) - 1
  keep <- names(which(n > 0 & n <= ncol(sample_grid_res$trajectories)))
  f <- function(k) {
    apply(sample_grid_res$trajectories[, k - 1, ], c(1, 3), sum)
  }
  totals <- lapply(index[keep], f)

  ## Compute deaths, icu and hosptialised:
  totals$deaths <- diff(totals$D)
  totals$icu <- totals$I_ICU_R + totals$I_ICU_D
  totals$hosp <- totals$I_triage + totals$I_hosp_R + totals$I_hosp_D +
    totals$I_ICU_R + totals$I_ICU_D + totals$R_stepdown

  totals
}

##' collapse into an array of trajectories
##' the trajectories are different lengths in terms of dates
##' so we will fill the arrays with NAs where needed
##' @importFrom utils tail
##' @noRd
traces_to_trajectories <- function(traces) {

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

  trajectories
}