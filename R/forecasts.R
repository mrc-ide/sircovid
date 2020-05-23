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
  y_grid <- matrix(scan_results$y, nrow = nr, 
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

  trajectories <- traces_to_trajectories(traces)
  
  # combine and return
  res <- list("trajectories" = trajectories,
              "param_grid" = param_grid,
              inputs = list(
                model_params = model_params,
                pars_obs = pars_obs,
                data = data,
                model = sircovid_model,
                forecast_days = forecast_days))
  
  class(res) <- "sircovid_forecast"
  return(res)
  
}


##' Take a parameter search produced by \code{\link{pmcmc}} and 
##' sample \code{n_sample} from the parameter space
##' on their probability. For each parameter pair chosen, run particle
##' filter with \code{num_particles} and sample 1 trajectory
##' 
##' @title Sample pmcmc
##' 
##' @param mcmc_results Output of \code{\link{pmcmc}}.
##' 
##' @param burn_in Number of burn-in samples to discard
##' 
##' @param n_sample Number of parameter pairs to be sampled. This will 
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
##' @importFrom furrr future_map
##' @importFrom purrr transpose
sample_pmcmc <- function(mcmc_results,
                         burn_in = 1,
                         n_sample = 10, 
                         n_particles = 100, 
                         forecast_days = 0) {
  
  # checks on args
  assert_is(mcmc_results, "pmcmc_list")
  assert_pos_int(n_sample)
  assert_pos_int(n_particles)
  assert_pos_int(forecast_days)
  
  # discard burn-in
  pars_to_sample <- names(mcmc_results$inputs$pars$pars_init)
  chains <- create_master_chain(mcmc_results, burn_in)
  param_grid <- chains[sample.int(n = nrow(chains), 
                                  size  = n_sample, 
                                  replace = FALSE), pars_to_sample]
  
  forecasts <- function(sampled_pars) {
    pars <- as.list(sampled_pars)
    pars$start_date <- sircovid_date(pars$start_date)
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
  traces <- furrr::future_map(.x = purrr::transpose(param_grid), .f = forecasts)
  trajectories <- traces_to_trajectories(traces)
  
  # combine and return
  res <- list("trajectories" = trajectories,
              "param_grid" = param_grid,
              inputs = list(
                model_params = mcmc_results$inputs$model_params,
                pars_obs = mcmc_results$inputs$pars_obs,
                data = mcmc_results$inputs$data,
                model = mcmc_results$inputs$sircovid_model,
                forecast_days = forecast_days))
  
  class(res) <- "sircovid_forecast"
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

##' Take a sampled grid search produced by \code{\link{sample_grid_scan}} summarise trajectories 
##' over compartmets
##' 
##' @title Summarise trajectories
##'
##' @param sample_grid_res Results from \code{\link{sample_grid_scan}}
##' 
##' @export
##' 
sum_over_compartments <- function(sample_grid_res) {
  index <- odin_index(sample_grid_res$inputs$model$odin_model(
      user = sample_grid_res$inputs$model_params))

  ## filter out the ones that we actually have - this is unreasonably ugly:
  n <- vapply(index, min, numeric(1)) - 1
  keep <- names(which(n > 0 & n <= ncol(sample_grid_res$trajectories)))
  f <- function(k) {
    apply(sample_grid_res$trajectories[, k - 1, , drop = FALSE], c(1, 3), sum)
  }
  totals <- lapply(index[keep], f)

  ## Compute deaths, icu and hosptialised:
  if ("sircovid_serology" %in% class(sample_grid_res$inputs$model)){
    totals$deaths <- diff(totals$D_hosp)
    totals$icu <- totals$I_ICU_R_conf + totals$I_ICU_D_conf
    totals$hosp <- totals$I_triage_R_conf + totals$I_triage_D_conf + totals$I_hosp_R_conf + totals$I_hosp_D_conf +
      totals$I_ICU_R_conf + totals$I_ICU_D_conf + totals$R_stepdown_conf
  } else {
    totals$deaths <- diff(totals$D)
    totals$icu <- totals$I_ICU_R + totals$I_ICU_D
    totals$hosp <- totals$I_triage + totals$I_hosp_R + totals$I_hosp_D +
      totals$I_ICU_R + totals$I_ICU_D + totals$R_stepdown
  } 

  totals
}

##' @export
plot.sircovid_forecast <- function(x, ..., what = "ICU", title = NULL, col = 'grey80') {
  idx <- odin_index(x$inputs$model$odin_model(
    user = x$inputs$model_params,
    unused_user_action = "message"
  ))
  
  # what are we plotting
  if (what == "ICU") {
    if ("sircovid_serology" %in% class(x$inputs$model)){
      index <- c(idx$I_ICU_D_conf, idx$I_ICU_R_conf) - 1L
      ylab <- "Confirmed covid-19 patients in ICU"
    } else {
      index <- c(idx$I_ICU_D, idx$I_ICU_R) - 1L
      ylab <- "ICU" 
    }
    particles <- vapply(seq_len(dim(x$trajectories)[3]), function(y) {
      rowSums(x$trajectories[, index, y], na.rm = TRUE)
    },
    FUN.VALUE = numeric(dim(x$trajectories)[1])
    )
    plot_particles(particles, ylab = ylab, title = title, col = col)
    points(as.Date(x$inputs$data$date), x$inputs$data$itu / x$inputs$pars_obs$phi_ICU, pch = 19)
  } else if (what == "general") {
    if ("sircovid_serology" %in% class(x$inputs$model)){
      index <- c(idx$I_triage_R_conf, idx$I_triage_D_conf, idx$I_hosp_R_conf, idx$I_hosp_D_conf, idx$R_stepdown_conf) - 1L
      ylab <- "Confirmed covid-19 patients in general beds"
    } else {
      index <- c(idx$I_triage, idx$I_hosp_R, idx$I_hosp_D, idx$R_stepdown) - 1L
      ylab <- "General beds" 
    }
    particles <- vapply(seq_len(dim(x$trajectories)[3]), function(y) {
      rowSums(x$trajectories[, index, y], na.rm = TRUE)
    },
    FUN.VALUE = numeric(dim(x$trajectories)[1])
    )
    plot_particles(particles, ylab = ylab, title = title, col = col)
    points(as.Date(x$inputs$data$date), x$inputs$data$general / x$inputs$pars_obs$phi_general, pch = 19)
  }
  
  else if (what == "deaths") {
    if ("sircovid_serology" %in% class(x$inputs$model)){
      index <- c(idx$D_hosp) - 1L
    } else {
      index <- c(idx$D) - 1L
    }
    
    ylab <- "Deaths"
    particles <- vapply(seq_len(dim(x$trajectories)[3]), function(y) {
      out <- c(0, diff(rowSums(x$trajectories[, index, y], na.rm = TRUE)))
      names(out)[1] <- rownames(x$trajectories)[1]
      out
    },
    FUN.VALUE = numeric(dim(x$trajectories)[1])
    )
    plot_particles(particles, ylab = ylab, title = title, col = col)
    points(as.Date(x$inputs$data$date),
           x$inputs$data$deaths / x$inputs$pars_obs$phi_death,
           pch = 19
    )
  } else if (what == "admitted") {
    index <- c(idx$cum_admit_conf) - 1L
    ylab <- "Patients admitted with covid-19"
    particles <- vapply(seq_len(dim(x$trajectories)[3]), function(y) {
      out <- c(0, diff(rowSums(x$trajectories[, index, y, drop = FALSE], na.rm = TRUE)))
      names(out)[1] <- rownames(x$trajectories)[1]
      out
    },
    FUN.VALUE = numeric(dim(x$trajectories)[1])
    )
    plot_particles(particles, ylab = ylab, title = title, col = col)
    points(as.Date(x$inputs$data$date),
           x$inputs$data$admitted / x$inputs$pars_obs$phi_admitted,
           pch = 19
    )
    
  } else if (what == "new") {
    index <- c(idx$cum_new_conf) - 1L
    ylab <- "Inpatients newly-diagnosed with covid-19"
    particles <- vapply(seq_len(dim(x$trajectories)[3]), function(y) {
      out <- c(0, diff(rowSums(x$trajectories[, index, y, drop = FALSE], na.rm = TRUE)))
      names(out)[1] <- rownames(x$trajectories)[1]
      out
    },
    FUN.VALUE = numeric(dim(x$trajectories)[1])
    )
    plot_particles(particles, ylab = ylab, title = title, col = col)
    points(as.Date(x$inputs$data$date),
           x$inputs$data$new / x$inputs$pars_obs$phi_new,
           pch = 19
    )
    
  } else {
    stop("Request what must be one of 'ICU', 'deaths', 'general', 'admitted' or 'new'")
  }
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