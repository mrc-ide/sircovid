##' Run a particle filter
##'
##' @title Run a particle filter
##'
##' @param data Data to fit to.  This must be constructed with
##'   \code{particle_filter_data}
##'
##' @param model An odin model, used to generate stochastic samples
##'
##' @param compare A function to generate log-weights
##'
##' @param n_particles Number of particles
##'
##' @param forecast_days Number of days to forecast forward from end
##'   states.  Requires that \code{save_particles} is \code{TRUE}.
##'
##' @param save_particles Logical, indicating if we save full particle
##'   histories (this is slower).
##'
##' @export
particle_filter <- function(data, model, compare, n_particles,
                            forecast_days = 0, save_particles = FALSE) {
  if (!inherits(data, "particle_filter_data")) {
    stop("Expected a data set derived from particle_filter_data")
  }
  if (!inherits(model, "odin_model")) {
    stop("Expected 'model' to be an 'odin_model' object")
  }
  if (n_particles < 2) {
    stop("At least two particles required")
  }
  if (forecast_days > 0 && !save_particles) {
    stop("forecasting only possible if particles are saved")
  }
  if (forecast_days < 0) {
    stop("forecast_days must be positive")
  }

  i_state <- seq_along(model$initial()) + 1L

  ## Special treatment for the burn-in phase; later we might use this
  ## same approach for skipping steps though.
  if (save_particles) {
    step <- seq(data$step_start[[1L]], data$step_end[[1L]],
                attr(data, "steps_per_day"))
    tmp <- model$run(step, use_names = FALSE, replicate = n_particles)

    particles <- array(NA_real_,
                       c(max(data$day_end) + 1L + forecast_days,
                         length(i_state), n_particles))
    particles[seq_len(data$day_end[[1]] + 1), , ] <- tmp[, i_state, ]
  } else {
    particles <- NULL
    step <- c(data$step_start[[1L]], data$step_end[[1L]])
    tmp <- model$run(step = step, use_names = FALSE, replicate = n_particles)
  }

  ## Current state of the system after burn-in: we'll update this in-place
  state <- tmp[nrow(tmp), i_state, , drop = TRUE]

  log_likelihood <- 0
  for (t in seq_len(nrow(data))[-1L]) {
    step <- c(data$step_start[t], data$step_end[t])
    prev_state <- state
    state <- particle_run_model(state, step, model)
    if (save_particles) {
      particles[data$day_end[t] + 1L, , ] <- state
    }

    log_weights <- compare(t, state, prev_state)
    if (!is.null(log_weights)) {
      tmp <- scale_log_weights(log_weights)
      log_likelihood <- log_likelihood + tmp$average
      if (tmp$average == -Inf) {
        ## Everything is impossible, so stop here
        break
      }

      kappa <- resample(tmp$weights, "systematic")
      state <- state[, kappa]
      if (save_particles) {
        particles <- particles[, , kappa]
      }
    }
  }

  if (forecast_days > 0) {
    step <- seq(data$step_end[nrow(data)],
                length.out = forecast_days + 1L,
                by = attr(data, "steps_per_day"))
    tmp <- model$run(step, state, replicate = n_particles, use_names = FALSE)
    i <- seq(data$day_end[nrow(data)] + 1, length.out = forecast_days + 1L)
    particles[i, , ] <- tmp[, i_state, ]
  }

  ret <- list(log_likelihood = log_likelihood)
  if (save_particles) {
    date <- d$date[[1]] + seq_len(nrow(particles)) - 1L
    rownames(particles) <- as.character(date)
    attr(particles, "date") <- date
    ret$states <- particles
  }
  ret
}


##' Prepare data for use with the particle filter.  This function
##' exists to make explicit how time changes through the model
##' relative to the data and to odin's internal clock.
##' @title Prepare particle filter data
##'
##' @param data A data.frame of observed data.  There must be a column
##'   \code{date}, containing dates (or ISO-formatted strings for
##'   conversion with \code{\link{as.Date}}.
##'
##' @param start_date The date to start the simulation from.  Must be
##'   earlier than the first date in \code{data}.
##'
##' @param steps_per_day The number of steps per day
##'
##' @export
particle_filter_data <- function(data, start_date, steps_per_day) {
  if (!inherits(data, "data.frame")) {
    stop("Expected a data.frame for 'data'")
  }
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

  class(ret) <- c("particle_filter_data", "data.frame")
  attr(ret, "steps_per_day") <- steps_per_day

  ret
}


particle_run_model <- function(y, step, model) {
  model$run(step, y, use_names = FALSE, return_minimal = TRUE,
            replicate = ncol(y))[, 1, , drop = TRUE]
}


resample <- function(weights, method) {
  if (method == "multinomial") {
    sample.int(length(weights), replace = TRUE, prob = weights)
  } else if (method == "systematic") {
    systematic_resample(weights)
  }
}


##' @importFrom stats runif
systematic_resample <- function(weights) {
  n <- length(weights)
  u <- runif(1, 0, 1 / n) + seq(0, by = 1 / n, length.out = n)
  cum_weights <- cumsum(weights / sum(weights))
  findInterval(u, cum_weights) + 1L
}


##' @importFrom stats dnbinom rexp
ll_nbinom <- function(data, model, phi, k, exp_noise) {
  mu <- phi * model + rexp(length(model), rate = exp_noise)
  dnbinom(data, k, mu = mu, log = TRUE)
}

##' @importFrom graphics plot points matlines
##' @importFrom stats quantile
plot_particles <- function(particles, ylab) {
  ## Need to set plot up first to get the dates to render on the axis
  ## (matplot does not cope with this)
  dates <- as.Date(rownames(particles))
  plot(dates, particles[, 1], type = "n", ylab = ylab, ylim = range(particles))
  ## Individual traces
  matlines(dates, particles, col="#ff444477", lty = 1)
  ## Quantiles
  quantiles <- t(apply(particles, 1, quantile, c(0.025, 0.5, 0.975)))
  matlines(dates, quantiles, col = "black", lty = "dashed")
}


## Awful hack that will do for now:
odin_index <- function(model) {
  n_out <- environment(model$initialize)$private$n_out %||% 0
  n_state <- length(model$initial())
  model$transform_variables(seq_len(1L + n_state + n_out))
}


compare_icu <- function(index, pars_obs, data) {
  force(index)
  force(pars_obs)

  ## Unpack things that we will use repeatedly
  phi_ICU <- pars_obs$phi_ICU
  k_ICU <- pars_obs$k_ICU
  phi_death <- pars_obs$phi_death
  k_death <- pars_obs$k_death
  exp_noise <- pars_obs$exp_noise
  index_ICU <- c(index$I_ICU) - 1L
  index_D <- c(index$D) - 1L

  function(t, results, prev_states) {
    if (is.na(data$itu[t] || is.na(data$deaths[t]))) {
      return(NULL)
    }

    log_weights <- rep(0, ncol(results))

    if (!is.na(data$itu[t])) {
      ## sum model output across ages/infectivities
      model_icu <- colSums(results[index_ICU, ])
      log_weights <- log_weights +
        ll_nbinom(data$itu[t], model_icu, phi_ICU, k_ICU, exp_noise)
    }

    if (!is.na(data$deaths[t])) {
      ## new deaths summed across ages/infectivities
      model_deaths <- colSums(results[index_D, ]) -
        colSums(prev_states[index_D, ])
      log_weights <- log_weights +
        ll_nbinom(data$deaths[t], model_deaths, phi_death, k_death, exp_noise)
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
