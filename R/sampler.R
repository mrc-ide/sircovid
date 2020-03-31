##' Run a particle filter
##'
##' @title Run a particle filter
##'
##' @param data Data to fit to.  We require a 'date' column here,
##'   everything else is up to you.
##'
##' @param model An odin model, used to generate stochastic samples
##'
##' @param compare A function to generate log-weights
##'
##' @param n_particles number of particles
##' @param start_date start date
##' @param time_steps_per_day Number of timesteps per day
##' @param forecast_days number of days to forecast forward from end states
##' @param output_states Logical, indicating if output of states is wanted
##' @param save_particles Logical, indicating if we save full particle histories
##' @export
particle_filter <- function(data, model, compare, n_particles,
                            start_date, time_steps_per_day, forecast_days = 0,
                            output_states = FALSE, save_particles = FALSE) {
  if (!("date" %in% names(data))) {
    stop("Expected a column 'date' within 'data'")
  }
  if (!inherits(model, "odin_model")) {
    stop("Expected 'model' to be an 'odin_model' object")
  }

  data2 <- sampler_prepare_data(data, start_date, time_steps_per_day)

  n_days <- nrow(data)
  log_ave_weight <- rep(0, n_days)

  state <- matrix(rep(mod$initial(), n_particles), ncol = n_particles)
  particles <- array(NA_real_, c(max(data2$day_end) + 1L, dim(state)))

  if (save_particles) {
    step <- seq(data2$step_start[[1]], data2$step_end[[1]], time_steps_per_day)
    tmp <- model$run(step, use_names = FALSE, replicate = n_particles)
    particles[seq_len(data2$day_end[[1]] + 1), , ] <-
      tmp[, seq_len(nrow(state)) + 1L, ]
  } else {
    step <- c(data2$step_start[[1]], data2$step_end[[1]])
    tmp <- model$run(step = step, use_names = FALSE, replicate = n_particles)
  }
  state <- tmp[nrow(tmp), seq_len(nrow(state)) + 1L, , drop = TRUE]

  ## Keep this in for ease of comparison later - unused here though.
  index <- model$transform_variables(seq_len(ncol(tmp)))

  for (t in seq_len(n_days)) {
    step <- c(data2$step_start[t + 1L], data2$step_end[t + 1L])
    prev_state <- state
    state <- particle_run_model(state, step, model)
    if (save_particles) {
      particles[data2$day_end[t + 1L] + 1L, , ] <- state
    }

    if (!is.na(data$itu[t]) && !is.na(data$deaths[t])) {
      log_weights <- compare(t, state, prev_state)

      tmp <- scale_log_weights(log_weights)
      log_ave_weight[t] <- tmp$average
      if (tmp$average == -Inf) {
        ## Everything is impossible, so stop here
        break
      }

      kappa <- resample(state, seq_len(n_particles), tmp$weights,
                        "systematic")
      state <- state[, kappa]
      if (save_particles) {
        particles <- particles[, , kappa]
      }
    }
  }

  if (forecast_days !=0){
    stop("This has not been checked")
    if (!save_particles){
      states <- array(states,dim=c(1,dim(states)))
    }
    for (t in seq_len(forecast_days)) {
      stop("compute step here")
      results <- particle_run_model(states[nrow(states),,], step, model)
      states <- abind::abind(states,results,along=1)
    }
  }

  ret <- list(log_likelihood = sum(log_ave_weight))
  if (save_particles) {
    ret$states <- particles
    ret$index <- index
  } else if (output_states) {
    ret$states <- state
    ret$index <- index
  }
  ret
}


particle_run_model <- function(y, step, model) {
  n_particles <- ncol(y)
  for (k in seq_len(n_particles)) {
    y[, k] <- model$run(step, y[, k],
                        use_names = FALSE, return_minimal = TRUE)
  }
  y
}


resample <- function(y, index, weights, method) {
  if (method == "multinomial") {
    sample(index, length(index), replace = TRUE, prob = weights)
  } else if (method == "systematic") {
    systematic_resample(index, length(index), weights)
  }
}


systematic_resample <- function(
  old_indexes,
  n_samples,
  weights){

  u = runif(1,0,1/n_samples)+seq(0,by=1/n_samples,length.out=n_samples)

  new_indexes <- integer(n_samples)
  weights <- weights/sum(weights)
  cum_weights <- cumsum(weights)
  k <- 1
  for (i in 1:n_samples) {
    found = FALSE
    while (!found){
      if (u[i] > cum_weights[k]) {
        k = k + 1
      } else {
        found = TRUE
      }
    }
    new_indexes[i] = old_indexes[k]
  }
  return(new_indexes)
}

ll_nbinom <- function(data, model,
                   phi, k, exp_noise) {
  out <- dnbinom(
    x = data,
    size = k,
    mu = phi*model+rexp(length(model),rate=exp_noise),
    log = TRUE
  )
  out
}

plot_particles <- function(particles, particle_dates, data, data_dates, ylab) {
  ## Need to set plot up first to get the dates to render on the axis
  ## (matplot does not cope with this)
  plot(particle_dates, particles[,1], type="n",
       ylab = ylab, ylim = range(particles))
  ## Individual traces
  matlines(particle_dates, particles, col="#ff444477", lty = 1)
  ## Observed data
  k <- !is.na(data)
  points(as.Date(data_dates[k]), data[k], col='black', pch=19)
  ## Quantiles
  quantiles <- t(apply(particles, 1, quantile, c(0.025, 0.5, 0.975)))
  matlines(particle_dates, quantiles, col = "black", lty = "dashed")
}


compare_icu <- function(index, pars_obs, data) {
  force(index)
  force(pars_obs)

  function(t, results, prev_states) {
    if (!is.na(data$itu[t]) && !is.na(data$deaths[t])) {
      ## calculate log weights from observation likelihood
      log_weights <- rep(0, ncol(results))

      ## contribution from itu/icu
      if (!is.na(data$itu[t])){
        ## sum model output across ages/infectivities
        model_icu <- colSums(results[c(index$I_ICU) - 1L, ])

        log_weights <- log_weights +
          ll_nbinom(data = data$itu[t], model = model_icu,
                    phi = pars_obs$phi_ICU, k = pars_obs$k_ICU,
                    exp_noise = pars_obs$exp_noise)
      }

      ## contribution from deaths
      if (!is.na(data$deaths[t])) {
        ## sum model output across ages/infectivities
        model_deaths <- colSums(results[c(index$D)-1L, ])-colSums(prev_states[c(index$D)-1L,])

        log_weights <- log_weights +
          ll_nbinom(data = data$deaths[t], model = model_deaths,
                    phi = pars_obs$phi_death, k = pars_obs$k_death,
                    exp_noise = pars_obs$exp_noise)
      }
    }
    log_weights
  }
}


scale_log_weights <- function(log_weights) {
  max_log_weights <- max(log_weights)
  if (max_log_weights == -Inf){
    ## if all log_weights at a time-step are -Inf, this should
    ## terminate the particle filter and output the marginal
    ## likelihood estimate as -Inf
    average <- -Inf
    weights <- rep(NaN, length(log_weights))
  } else {
    ## calculation of weights, there is some rescaling here to avoid
    ## issues where exp(log_weights) might give computationally zero
    ## values
    weights <- exp(log_weights - max_log_weights)
    average <- log(mean(weights)) + max_log_weights
  }
  list(weights = weights, average = average)
}

## This is very clumsy way of avoiding fencepost/off by one errors in
## translation.  This creates our set of periods to run from.
sampler_prepare_data <- function(data, start_date, steps_per_day) {
  if (!("date" %in% names(data))) {
    stop("Expected a column 'date' within 'data'")
  }
  data$date <- as.Date(data$date)
  if (any(diff(data$date) <= 0)) {
    stop("'date' must be strictly increasing")
  }
  start_date <- as.Date(start_date)
  if (start_date >= data$date[[1]]) {
    stop("'start_date' must be less than the first date in data")
  }

  ## Then for each timestep we work out the start and end date
  ret <- data
  ret$day_start <- as.integer(data$date - start_date - 1L)
  ret$day_end <- as.integer(data$date - start_date)

  d0 <- ret[1, ]
  d0[] <- NA
  d0$date <- start_date
  d0$day_start <- 0
  d0$day_end <- ret$day_start[[1]]
  ret <- rbind(d0, ret)
  rownames(ret) <- NULL

  ret$step_start <- ret$day_start * steps_per_day
  ret$step_end <- ret$day_end * steps_per_day

  ret
}
