##' Run a pmcmc sampler
##'
##' @title Run a pmcmc sampler
##' 
##' @param data Data to fit to.  This must be constructed with
##'   \code{particle_filter_data}
##'   
##' @param model_params Model parameters, from a call to
##'   \code{generate_parameters()}. If NULL, uses defaults as
##'   in unit tests.
##'   
##' @param pars_obs list of parameters to use in comparison
##'   with \code{compare_icu}. If NULL, uses
##'   list(phi_ICU = 0.95, k_ICU = 2, phi_death = 926 / 1019, k_death = 2, exp_noise = 1e6)
##'   
##' @param n_mcmc number of mcmc mcmc iterations to perform
##' 
##' @param pars_init named list of initial inputs for parameters being sampled (currently beta and start_date)
##' 
##' @param pars_min named list of lower reflecting boundaries for parameter proposals
##' 
##' @param pars_max named list of upper reflecting boundaries for parameter proposals
##' 
##' @param pars_sd named list of proposal sds for parameters (currently normal dists, but can easily be made MVNorm)
##' 
##' @param pars_discrete named list of logicals, indicating if proposed jump should be discrete
##' 
##' @param log_likelihood function to calculate log likelihood, must take named parameter vector as input, 
##'                 allow passing of implicit arguments corresponding to the main function arguments. 
##'                 Returns a named list, with entries:
##'                   - $log_likelihood, a single numeric
##'                   - $sample_state, a numeric vector corresponding to the state of a single particle, chosen at random, 
##'                   at the final time point for which we have data.
##'                   If NULL, calculated using the function calc_loglikelihood.
##'                   
##' @param log_prior function to calculate log prior, must take named parameter vector as input, returns a single numeric.
##'                  If NULL, uses uninformative priors which do not affect the posterior
##'                   
##' @param n_particles Number of particles
##' 
##' @param steps_per_day Number of steps per day
##' 
##' @return list of length two containing
##' - List of inputs
##' - Matrix of accepted parameter samples, rows = iterations
##'   as well as log prior, (particle filter estimate of) log likelihood and log posterior
##'   
##' @description The user inputs initial parameter values for beta and sample_date
##' The log prior likelihood of these parameters is calculated based on the user-defined
##' prior distributions.
##' The log likelihood of the data given the initial parameters is estimated using a particle filter,
##' which has two functions:
##'      - Firstly, to generate a set of 'n_particles' samples of the model state space,
##'        at time points corresponding to the data, one of which is 
##'        selected randomly to serve as the proposed state sequence sample at the final
##'        data time point.
##'      - Secondly, to produce an unbiased estimate of the likelihood of the data given the proposed parameters. 
##' The log posterior of the initial parameters given the data is then estimated by adding the log prior and 
##' log likelihood estimate.
##' 
##' The pMCMC sampler then proceeds as follows, for n_mcmc iterations:
##' At each loop iteration the pMCMC sampler performs three steps:
##'   1. Propose new candidate samples for beta and start_date based on
##'     the current samples, using the proposal distribution 
##'     (currently i.i.d Gaussian with user-input sd pars_sd, and reflecting boundaries)
##'   2. Calculate the log prior of the proposed parameters, 
##'      Use the particle filter to estimate log likelihood of the data given the proposed parameters, as described above,
##'      as well as proposing a model state space.
##'      Add the log prior and log likelihood estimate to estimate the log posterior of the proposed parameters given the data.
##'   3. Metropolis-Hastings step: The joint canditate sample (consisting of the proposed parameters 
##'      and state space) is then accepted with probability min(1, a), where the acceptance ratio is
##'      simply the ratio of the posterior likelihood of the proposed parameters to the posterior likelihood
##'      of the current parameters. Note that by choosing symmetric proposal distributions by including
##'      reflecting boundaries, we avoid the the need to include the proposal likelihood in the MH ratio.
##'      
##'   If the proposed parameters and states are accepted then we update the current parameters and states
##'   to match the proposal, otherwise the previous parameters/states are retained for the next iteration.
##'   
##' @export
##' @import coda
pmcmc <- function(data,
                  n_mcmc, 
                  model_params = NULL,
                  pars_obs = list(phi_ICU = 0.95,
                                  k_ICU = 2,
                                  phi_death = 926 / 1019,
                                  k_death = 2,
                                  exp_noise = 1e6),
                  pars_init = list('beta' = 0.14, 
                                   'start_date' = as.Date("2020-02-07")),
                  pars_min = list('beta' = 0, 
                                  'start_date' = 0),
                  pars_max = list('beta' = 1, 
                                  'start_date' = 1e6),
                  pars_sd = list('beta' = 0.001, 
                                 'start_date' = 0.5),
                  pars_discrete = list('beta' = FALSE, 
                                       'start_date' = TRUE),
                  log_likelihood = NULL,
                  log_prior = NULL,
                  n_particles = 1e2,
                  steps_per_day = 4
                  ) {
  
  if (is.null(model_params)) {
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
      dt = 1/steps_per_day
    )
  } else {
    if (length(model_params$beta_y)) {
      stop("Time-varying beta not yet supported by pmcmc_sampler")
    }
  }
  
  
  inputs <- list(
    data = data,
    n_mcmc = n_mcmc,
    pars = list(pars_obs = pars_obs, 
                pars_init = pars_init, 
                pars_min = pars_min, 
                pars_max = pars_max, 
                pars_sd = pars_sd,
                pars_discrete = pars_discrete), 
    n_particles = n_particles, 
    steps_per_day = steps_per_day)
  
  
  par_names <- c('beta', 'start_date')
  
  
  
  correct_format <- function(pars) {
    is.list(pars) & setequal(names(pars), par_names)
  }

  if(!correct_format(pars_init)) {
    stop("pars_init must be a list of length two with names 'beta' and 'start_date'")
  }
  if(!correct_format(pars_min)) {
    stop("pars_min must be a list of length two with names 'beta' and 'start_date'")
  }
  if(!correct_format(pars_max)) {
    stop("pars_max must be a list of length two with names 'beta' and 'start_date'")
  }
  if(!correct_format(pars_sd)) {
    stop("pars_sd must be a list of length two with names 'beta' and 'start_date'")
  }
  if(!correct_format(pars_discrete)) {
    stop("pars_discrete must be a list of length two with names 'beta' and 'start_date'")
  }
  
  
  is.numeric.list <- function(x) all(vapply(X = x, is.numeric, logical(1)))
  is.logical.list <- function(x) all(vapply(X = x, is.logical, logical(1)))
  
  if(!is.numeric.list(pars_min)) {
    stop('pars_min entries must be numeric')
  }
  if(!is.numeric.list(pars_max)) {
    stop('pars_max entries must be numeric')
  }
  if(!is.numeric.list(pars_sd)) {
    stop('pars_sd entries must be numeric')
  }
  if(!is.logical.list(pars_discrete)) {
    stop('pars_discrete entries must be logical')
  }
  

  n_pars <- length(par_names)
  
  # needs to be a vector to pass to reflecting boundary function
  pars_min <- unlist(pars_min)  
  pars_max <- unlist(pars_max)
  pars_sd <- unlist(pars_sd)
  pars_discrete <- unlist(pars_discrete)
  
  
  
  # create prior and likelihood functions given the inputs
 
   if(is.null(log_prior)) {
     # set improper, uninformative prior
    log_prior <- function(pars) log(1e-10)
   }
  calc_lprior <- log_prior
  
  if(is.null(log_likelihood)) {
    log_likelihood <- calc_loglikelihood
  } else if (!('...' %in% names(formals(log_likelihood)))){
    stop('log_likelihood function must be able to take unnamed arguments')
  }
  
  
  
  # create shorthand function to calc ll given main inputs
  calc_ll <- function(pars) {  
    X <- log_likelihood(pars = pars, 
                            data = data, 
                            model_params = model_params,
                            steps_per_day = steps_per_day, 
                            pars_obs = pars_obs, 
                            n_particles = n_particles
    ) 
    X
  }
  
  
  
  # convert the current parameters into format easier for mcmc to deal with
  curr_pars <- unlist(pars_init)
  day_zero <- as.Date(data$date[1])
  curr_pars[['start_date']] <- as.numeric(day_zero - pars_init$start_date) # convert to numeric
  
  
  if(any(curr_pars < pars_min | curr_pars > pars_max)) {
    stop('initial parameters are outside of specified range')
  }
  

  ## calculate initial prior
  curr_lprior <- calc_lprior(pars = curr_pars)
  
  # run particle filter on initial parameters
  p_filter_est <- calc_ll(pars = curr_pars)
  
  
  # checks on log_prior and log_likelihood functions
  
  if(length(curr_lprior) > 1) {
    stop('log_prior must return a single numeric representing the log prior')
  }
  if(curr_lprior > 0 ) {
    stop('log_prior must be negative or zero')
  }
  if(is.infinite(curr_lprior)) {
    stop('initial parameters are not compatible with supplied prior')
  }
  
  if(length(p_filter_est) != 2) {
    stop('log_likelihood function must return a list containing elements log_likelihood and sample_state')
  }
  if(!setequal(names(p_filter_est), c('log_likelihood', 'sample_state'))) {
    stop('log_likelihood function must return a list containing elements log_likelihood and sample_state')
  }
  if(length(p_filter_est$log_likelihood) > 1) {
    stop('log_likelihood must be a single numeric representing the estimated log likelihood')
  }
  if(p_filter_est$log_likelihood > 0) {
    stop('log_likelihood must be negative or zero')
  }
  if(any(p_filter_est$sample_state < 0)) {
    stop('sample_state must be a vector of non-negative numbers')
  }
  
  
  # extract loglikelihhod estimate and sample state
  # calculate posterior
  curr_ll <- p_filter_est$log_likelihood
  curr_lpost <- curr_lprior + curr_ll
  
  curr_ss <- p_filter_est$sample_state
  

  res_init <- c(curr_pars, 
                'log_prior' = curr_lprior, 
                'log_likelihood' = curr_ll, 
                'log_posterior' = curr_lpost) 
  
  ## initialise output arrays
  
  res <- matrix(data = NA, 
                nrow = n_mcmc + 1L, 
                ncol = length(res_init), 
                dimnames = list(NULL, 
                                names(res_init)))
  
  states <- matrix(data = NA,
                   nrow = n_mcmc + 1L, 
                   ncol = length(curr_ss))
  
  
  
  ## record initial results
  res[1, ] <- res_init
  states[1, ] <- curr_ss
  
  
  
  ## create function to reflect proposal boundaries at pars_min and pars_max
  # this ensures the proposal is symetrical and we can simplify the MH step
  
  reflect_proposal <- function(x, floor, cap) {
    interval <- cap - floor
    abs((x + interval - floor) %% (2 * interval) - interval) + floor
  }
  
  
   # this can easily be adapted to being multivariate normal with covariance to improve mixing
  propose_jump <- function(pars) {
    
    ## proposed jumps are normal with mean pars and sd as input for parameter
    tmp <- rnorm(n = n_pars, mean = pars, sd = pars_sd)
    # discretise if necessary
    tmp[pars_discrete] <- round(tmp[pars_discrete])
    # reflect proposal if it exceeds upper or lower parameter boundary
    tmp <- reflect_proposal(x = tmp, 
                            floor = pars_min, 
                            cap = pars_max)
    tmp
  }


  # start progress bar  
  pb <- txtProgressBar(min = 0, max = n_mcmc, style = 3)
  
  # main pmcmc loop
  for(iter in seq_len(n_mcmc) + 1L) {
    
    setTxtProgressBar(pb, iter)
    
    # propose new parameters
    prop_pars <- propose_jump(curr_pars)
    
    ## calculate proposed prior / lhood / posterior 
    prop_lprior <- calc_lprior(pars = prop_pars)
    p_filter_est <- calc_ll(pars = curr_pars)
    prop_ll <- p_filter_est$log_likelihood
    prop_ss <- p_filter_est$sample_state
    prop_lpost <- prop_lprior + prop_ll
    
    
    # calculate probability of acceptance
    accept_prob <- exp(prop_lpost - curr_lpost)

    
    if(runif(1) < accept_prob) { # MH step
      # update parameters and calculated likelihoods
      curr_pars <- prop_pars
      curr_lprior <- prop_lprior
      curr_ll <- prop_ll
      curr_lpost <- prop_lpost
      curr_ss <- prop_ss
    }
    
    # record results
    res[iter, ] <- c(curr_pars, 
                     'log_prior' = curr_lprior, 
                     'log_likelihood' = curr_ll, 
                     'log_posterior' = curr_lpost) 
    states[iter, ] <- curr_ss

  }
  # end progress bar
  close(pb)
  
  
  res <- as.data.frame(res)
  
  tmp <- coda::as.mcmc(res)
  rejection_rate <- coda::rejectionRate(tmp)
  ess <- coda::effectiveSize(tmp)

  res$start_date <- as.Date(data$date[1]) - res$start_date
  
  
  
 list('inputs' = inputs, 
      'results' = as.data.frame(res),
      'states' = states, 
      'acceptance_rate' = 1-rejection_rate, 
      "ess" = ess
      )

}


calc_loglikelihood <- function(pars, data, model_params,  steps_per_day, pars_obs, n_particles) {
  
  start_date <- as.Date(data$date[1]) - pars[['start_date']]
  model_params$beta_y <-  pars[['beta']]
  
  # create sir_covid model
  mod <- sircovid(params = model_params)

  d <- particle_filter_data(data,
                            start_date = start_date, 
                            steps_per_day = steps_per_day)
  compare <- compare_icu(mod, pars_obs, d)
  X <- particle_filter(data = d, 
                       model = mod, 
                       compare = compare, 
                       n_particles = n_particles, 
                       save_sample_state = TRUE)

  X
}
