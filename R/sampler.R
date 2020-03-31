##' Run a particle filter
##'
##' @title Run a particle filter
##'
##' @param data Data to fit to.  We require a 'date' column here,
##'   everything else is up to you.
##'
##' @param model An odin model, used to generate stochastic samples
##'
##' @param pars_obs parameters of the observations
##' @param n_particles number of particles
##' @param start_date start date
##' @param time_steps_per_day Number of timesteps per day
##' @param forecast_days number of days to forecast forward from end states
##' @param output_states Logical, indicating if output of states is wanted
##' @param save_particles Logical, indicating if we save full particle histories
##' @export
particle_filter <- function(data, model, pars_obs, n_particles,
                            start_date, time_steps_per_day, forecast_days = 0,
                            output_states = FALSE, save_particles = FALSE) {
  if (!("date" %in% names(data))) {
    stop("Expected a column 'date' within 'data'")
  }
  if (!inherits(model, "odin_model")) {
    stop("Expected 'model' to be an 'odin_model' object")
  }

  #setup
  n_days_before_data <- as.numeric(as.Date(data$date[1])-as.Date(start_date))
  n_days <- nrow(data)
  log_ave_weight <- rep(0,n_days)
  
  if (save_particles) {
    tmp <- model$run(step = seq(0, (n_days_before_data-1)*time_steps_per_day, time_steps_per_day),
                 replicate = n_particles)
    states <- unname(tmp[, seq_along(model$initial()) + 1L, ])
  } else {
    tmp <- model$run(step = c(0, (n_days_before_data-1)*time_steps_per_day),
                   replicate = n_particles)
    states <- unname(tmp[nrow(tmp), seq_along(model$initial()) + 1L, ])
  }
  index <- model$transform_variables(seq_len(ncol(tmp)))

  for (t in seq_len(n_days)) {
    if (save_particles) {
      prev_states <- states[nrow(states),,]
      results <- particle_run_model(states[nrow(states),,], time_steps_per_day, model)
      states <- abind::abind(states,results,along=1)
    } else {
      prev_states <- states
      results <- particle_run_model(states, time_steps_per_day, model)
      states <- results
    }
    
    if (!is.na(data$itu[t]) && !is.na(data$deaths[t])){
      
      #calculate log weights from observation likelihood
      log_weights <- rep(0, n_particles)
      
      #contribution from itu/icu
      if (!is.na(data$itu[t])){
        #sum model output across ages/infectivities
        model_icu <- colSums(results[c(index$I_ICU)-1L, ])
        
        log_weights <- log_weights+ll_nbinom(data = data$itu[t], model = model_icu,
                            phi = pars_obs$phi_ICU, k = pars_obs$k_ICU,
                            exp_noise = pars_obs$exp_noise)
      }
      
      #contribution from deaths
      if (!is.na(data$deaths[t])){
        #sum model output across ages/infectivities
        model_deaths <- colSums(results[c(index$D)-1L, ])-colSums(prev_states[c(index$D)-1L,])
        
        log_weights <- log_weights+ll_nbinom(data = data$deaths[t], model = model_deaths,
                              phi = pars_obs$phi_death, k = pars_obs$k_death,
                              exp_noise = pars_obs$exp_noise)
      }
      
      max_log_weights <- max(log_weights)
      if (max_log_weights == -Inf){
        #if all log_weights at a time-step are -Inf, this should terminate the particle filter
        #and output the marginal likelihood estimate as -Inf
        log_ave_weight[t] <- -Inf
        break
      } else {
        #calculation of weights, there is some rescaling here to avoid issues where exp(log_weights) might
        #give computationally zero values
        weights <- exp(log_weights - max_log_weights)
        log_ave_weight[t] <- log(mean(weights)) + max_log_weights
      }
      
      # resample
      states <- resample(states, seq_len(n_particles), weights, "systematic", save_particles)
      
    }
  }
  
  log_likelihood <- sum(log_ave_weight)
  
  if (forecast_days !=0){
    if (!save_particles){
      states <- array(states,dim=c(1,dim(states)))
    }
    for (t in seq_len(forecast_days)) {
      results <- particle_run_model(states[nrow(states),,], time_steps_per_day, model)
      states <- abind::abind(states,results,along=1)
    }
  }
  
  if (output_states){
    return(list(log_likelihood=log_likelihood,states=states,index=index))
  } else {
    return(list(log_likelihood=log_likelihood))
  }
}


particle_run_model <- function(y, step, model) {
  n_particles <- ncol(y)
  for (k in seq_len(n_particles)) {
    y[, k] <- model$run(c(0, step), y[, k],
                        use_names = FALSE, return_minimal = TRUE)
  }
  y
}


resample <- function(y, index, weights, method, save_particles) {
  if (method == "multinomial") {
    kappa <- sample(index, length(index), replace = TRUE, prob = weights)
  } else if (method == "systematic") {
    kappa <- systematic_resample(index, length(index), weights)
  }
  if (save_particles==TRUE){
    y[,, kappa]
  } else {
    y[, kappa]
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
