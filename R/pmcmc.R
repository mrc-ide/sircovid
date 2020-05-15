##' Run a pmcmc sampler
##'
##' @title Run a pmcmc sampler
##' 
##' @param data Data to fit to.  This must be constructed with
##'   \code{particle_filter_data}
##'   
##' @param model_params Model parameters, from a call to
##'   \code{generate_parameters()}. 
##'   
##' @param sircovid_model An odin model generator and comparison function.

##' @param pars_obs list of parameters to use in comparison
##'   with \code{compare_icu}. Must be a list containing, e.g.
##'   list(phi_general = 0.95, 
##'   k_general = 2,
##'   phi_ICU = 0.95, 
##'   k_ICU = 2, 
##'   phi_death = 926 / 1019, 
##'   k_death = 2, 
##'   exp_noise = 1e6)
##'   
##' @param n_mcmc number of mcmc mcmc iterations to perform
##' 
##' @param pars_to_sample \code{data.frame} detailing parameters to sample. Must contain columns
##'   'names' (parameter names), 'init' (initial values), 'min' (minimum values), 'max' (maximum values),
##'   'discrete (boolean indicating whether a discrete quantity)'
##'                   
##' @param pars_lprior functions to calculate log prior for each parameter. A named list for each parameter
##'   listed in \code{pars_to_sample}. Each value must be a function which takes named parameter vector as 
##'   input, returns a single numeric which is log of the prior probability.
##'       
##' @param proposal_kernel named matrix of proposal covariance for parameters 
##'
##' @param n_particles Number of particles
##' 
##' @param steps_per_day Number of steps per day
##' 
##' @param output_proposals Logical indicating whether proposed parameter jumps should be output along with results
##' 
##' @param n_chains Number of chains to run
##'
##' @return an mcmc object containing
##' - List of inputs
##' - Matrix of accepted parameter samples, rows = iterations
##'   as well as log prior, (particle filter estimate of) log likelihood and log posterior
##'   
##' @description The user inputs initial parameter values for beta_start and sample_date
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
##' At each loop iteration the pMCMC sampler perfsorms three steps:
##'   1. Propose new candidate samples for beta_start, beta_end and start_date based on
##'     the current samples, using the proposal distribution 
##'     (currently multivariate Gaussian with user-input covariance matrix (proposal_kernel), and reflecting boundaries defined by pars_min, pars_max)
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
##' @importFrom stats rnorm
##' @importFrom mvtnorm rmvnorm
pmcmc <- function(data,
                  n_mcmc, 
                  sircovid_model,
                  model_params,
                  pars_obs,
                  pars_to_sample = data.frame(
                                    names=c('beta_start',
                                     'beta_end', 
                                     'beta_pl',
                                     'start_date',  
                                     'gamma_triage', 
                                     'gamma_hosp_R', 
                                     'gamma_hosp_D', 
                                     'gamma_ICU_R', 
                                     'gamma_ICU_D', 
                                     'gamma_stepdown'),
                                    init=c(0.14, 
                                           0.14*0.238,
                                           0.14*0.238,
                                           as.Date("2020-02-07"),
                                           0.5099579,
                                           0.1092046,
                                           0.2911154,
                                           0.3541429,
                                           0.2913861,
                                           0.452381),
                                    min=c(0,
                                          0,
                                          0,
                                          0,
                                          0,
                                          0,
                                          0,
                                          0,
                                          0,
                                          0),
                                    max=c(1,
                                          1,
                                          1,
                                          1e6,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1),
                                    discrete=c(FALSE,
                                               FALSE,
                                               FALSE,
                                               TRUE,
                                               FALSE,
                                               FALSE,
                                               FALSE,
                                               FALSE,
                                               FALSE,
                                               FALSE),
                                    stringsAsFactors = FALSE),
                  pars_lprior = list('beta_start'     = function(pars) log(1e-10),
                                     'beta_end'       = function(pars) log(1e-10),
                                     'beta_pl'       = function(pars) log(1e-10),
                                     'start_date'     = function(pars) log(1e-10),
                                     'gamma_triage'   = function(pars) log(1e-10),
                                     'gamma_hosp_R'   = function(pars) log(1e-10),
                                     'gamma_hosp_D'   = function(pars) log(1e-10),
                                     'gamma_ICU_R'    = function(pars) log(1e-10),
                                     'gamma_ICU_D'    = function(pars) log(1e-10),
                                     'gamma_stepdown' = function(pars) log(1e-10)),
                  proposal_kernel,
                  n_particles = 1e2,
                  steps_per_day = 4, 
                  output_proposals = FALSE,
                  n_chains = 1) {

  if(length(output_proposals) != 1 || !is.logical(output_proposals)) {
    stop("output_proposals must be either TRUE or FALSE")
  }
  
  if (length(model_params$beta_y) > 1) {
    stop("Set beta variation through generate_beta_func in sircovid_model, not model_params")
  }
  
  #
  # Check pars_init input
  #
  if (!all(c("names", "init", "min", "max", "discrete") %in% colnames(pars_to_sample))) {
    stop("pars_to_samples must contain columns 'names', 'init', 'min', 'max' and 'discrete'")
  }
  par_names <- as.character(pars_to_sample$names)
  
  if (!all(c('beta_start', 'start_date') %in% par_names)) {
    stop("Turning off beta and start date sampling unsupported")
  }
  if(!all(names(pars_lprior) %in% par_names)) {
    stop("All sampled parameters must have a defined prior")
  }
  
  correct_format_mat <- function(pars, names) {
       setequal(rownames(pars),
                colnames(pars)) && 
         setequal(rownames(pars),
                   par_names) &&
      nrow(pars) == length(names) &&
      ncol(pars) == length(names)
  }

  if(!correct_format_mat(proposal_kernel, par_names)) {
    stop("proposal_kernel must be a matrix or vector with names corresponding to the parameters being sampled")
  }
  
  # Convert data frame into named lists
  df_col_to_list <- function(df, column) {
    list_out <- as.list(df[[column]])
    names(list_out) <- as.character(df$names)
    list_out
  }

  pars_init <- df_col_to_list(pars_to_sample, "init")
  pars_min <- df_col_to_list(pars_to_sample, "min")
  pars_max <- df_col_to_list(pars_to_sample, "max")
  pars_discrete <- df_col_to_list(pars_to_sample, "discrete")

  is.numeric.list <- function(x) all(vapply(X = x, is.numeric, logical(1)))
  is.logical.list <- function(x) all(vapply(X = x, is.logical, logical(1)))
  
  if(!is.numeric.list(pars_min)) {
    stop('pars_min entries must be numeric')
  }
  if(!is.numeric.list(pars_max)) {
    stop('pars_max entries must be numeric')
  }

  if(!is.logical.list(pars_discrete)) {
    stop('pars_discrete entries must be logical')
  }

  if(any(pars_to_sample$init < pars_to_sample$min | pars_to_sample$init > pars_to_sample$max)) {
    stop('initial parameters are outside of specified range')
  }
  if(pars_init['beta_start'] < 0) {
    stop('beta_start must not be negative')
  }
  if('beta_end' %in% par_names) {
    if(pars_init['beta_end'] < 0) {
      stop('beta_end must not be negative')
    }
  }
  if(start_date_to_offset(data$date[1], pars_init[['start_date']]) < 0) {
    stop('start date must not be before first date of supplied data')
  }

  #
  # Generate MCMC parameters
  #
  
  # This gets changed to numeric in data.frame
  pars_init[['start_date']] <- as.Date(pars_init[['start_date']], origin = "1970-01-01")
  
  inputs <- list(
    data = data,
    n_mcmc = n_mcmc,
    model_params = model_params,
    pars_obs = pars_obs,
    sircovid_model = sircovid_model,
    pars = list(pars_obs = pars_obs, 
                pars_init = pars_init, 
                pars_min = pars_min, 
                pars_max = pars_max, 
                proposal_kernel = proposal_kernel,
                pars_discrete = pars_discrete), 
    n_particles = n_particles, 
    steps_per_day = steps_per_day)
  
  # convert dates to be Date objects
  data$date <- as.Date(data$date) 
  
  # needs to be a vector to pass to reflecting boundary function
  pars_min <- unlist(pars_min)  
  pars_max <- unlist(pars_max)
  pars_discrete <- unlist(pars_discrete)
  curr_pars <- unlist(pars_init)

  # convert the current parameters into format easier for mcmc to deal with
  curr_pars['start_date'] <- start_date_to_offset(data$date[1], pars_init$start_date) # convert to numeric
  
  # create shorthand function to calc ll given main inputs
  # binds these input parameters
  calc_ll <- function(pars) {  
    X <- calc_loglikelihood(pars = pars, 
                            data = data,
                            sircovid_model = sircovid_model,
                            model_params = model_params,
                            steps_per_day = steps_per_day, 
                            pars_obs = pars_obs, 
                            n_particles = n_particles,
                            forecast_days = 0,
                            return = "ll"
    ) 
    X
  }

  # same thing for priors
  calc_lprior <- function(pars) {
    lprior <- 0
    for (par in names(pars)) {
      lprior <- lprior + pars_lprior[[par]](pars) 
    }
    lprior
  }

  # create shorthand function to propose new pars given main inputs
  propose_jump <- function(pars) {
    propose_parameters(pars = pars, 
                       proposal_kernel = proposal_kernel,
                       pars_discrete = pars_discrete,
                       pars_min = pars_min,
                       pars_max = pars_max)
  }

  chains <- furrr::future_pmap(
      .l =  list(n_mcmc = rep(n_mcmc, n_chains)), 
      .f = run_mcmc_chain,
      inputs = inputs,
      curr_pars = curr_pars,
      calc_lprior = calc_lprior,
      calc_ll = calc_ll,
      propose_jump = propose_jump,
      first_data_date = data$date[1],
      output_proposals = output_proposals,
      .progress = TRUE)
  
  if (n_chains > 1) {
    names(chains) <- paste0('chain', seq_len(n_chains))  
    
    # calculating rhat
    # convert parallel chains to a coda-friendly format
    chains_coda <- lapply(chains, function(x) {
        
        traces <- x$results
        if('start_date' %in% pars_to_sample$names) {
          traces$start_date <- start_date_to_offset(data$date[1], traces$start_date)
        }
        
      coda::as.mcmc(traces[, names(pars_init)])
    })

    rhat <- tryCatch(expr = {
      x <- coda::gelman.diag(chains_coda)
      x
    }, error = function(e) {
      print('unable to calculate rhat')
      })
    
    
    res <- list(inputs = chains[[1]]$inputs, 
                rhat = rhat,
                chains = lapply(chains, '[', -1))

    class(res) <- 'pmcmc_list'
  } else {
    res <- chains[[1]]
  }

  res

}

# Run a single pMCMC chain
run_mcmc_chain <- function(inputs,
                     curr_pars,
                     calc_lprior,
                     calc_ll,
                     n_mcmc,
                     propose_jump,
                     first_data_date,
                     output_proposals) {

  #
  # Set initial state
  #

  ## calculate initial prior
  curr_lprior <- calc_lprior(pars = curr_pars)
  
  # run particle filter on initial parameters
  p_filter_est <- calc_ll(pars = curr_pars)
  
  # checks on log_prior and log_likelihood functions
  
  if(length(curr_lprior) > 1) {
    stop('log_prior must return a single numeric representing the log prior')
  }

  if(is.infinite(curr_lprior)) {
    stop('initial parameters are not compatible with supplied prior')
  }
  
  # extract loglikelihood estimate and sample state
  # calculate posterior
  curr_ll <- p_filter_est$log_likelihood
  curr_lpost <- curr_lprior + curr_ll
  curr_ss <- p_filter_est$sample_state

  #
  # Create objects to store outputs
  #
  
  # initialise output arrays
  res_init <- c(curr_pars, 
                'log_prior' = curr_lprior, 
                'log_likelihood' = curr_ll, 
                'log_posterior' = curr_lpost) 
  res <- matrix(data = NA, 
                nrow = n_mcmc + 1L, 
                ncol = length(res_init), 
                dimnames = list(NULL, 
                                names(res_init)))
  
  states <- matrix(data = NA,
                   nrow = n_mcmc + 1L, 
                   ncol = length(curr_ss))
    
  
  if(output_proposals) {
    proposals <- matrix(data = NA, 
                        nrow = n_mcmc + 1L, 
                        ncol = length(res_init) + 1L, 
                        dimnames = list(NULL, 
                                        c(names(res_init), 
                                          'accept_prob')))
  }

  ## record initial results
  res[1, ] <- res_init
  states[1, ] <- curr_ss

  #
  # main pmcmc loop
  #

  for(iter in seq_len(n_mcmc) + 1L) {
    
    # propose new parameters
    prop_pars <- propose_jump(curr_pars)
    
    ## calculate proposed prior / lhood / posterior 
    prop_lprior <- calc_lprior(pars = prop_pars)
    p_filter_est <- calc_ll(pars = prop_pars)
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
                     curr_lprior, 
                     curr_ll, 
                     curr_lpost) 
    states[iter, ] <- curr_ss
    
    
    if(output_proposals) {
      proposals[iter, ] <- c(prop_pars, 
                             prop_lprior, 
                             prop_ll, 
                             prop_lpost, 
                             min(accept_prob, 1)) 
    }

  }
  
  res <- as.data.frame(res)
  
  coda_res <- coda::as.mcmc(res)
  rejection_rate <- coda::rejectionRate(coda_res)
  ess <- coda::effectiveSize(coda_res)

  res$start_date <- offset_to_start_date(first_data_date, res$start_date)
  
  out <- list('inputs' = inputs, 
              'results' = as.data.frame(res),
              'states' = states, 
              'acceptance_rate' = 1-rejection_rate, 
              "ess" = ess)
 
 if(output_proposals) {
   proposals <- as.data.frame(proposals)
   proposals$start_date <- offset_to_start_date(first_data_date, proposals$start_date)
   out$proposals <- proposals
 }
 
 class(out) <- 'pmcmc'
 out

}

# Run odin model to calculate log-likelihood
# 
# return: Set to 'll' to return the log-likelihood (for MCMC) or to
#
calc_loglikelihood <- function(pars, data, sircovid_model, model_params,
                               steps_per_day, pars_obs, n_particles,
                               forecast_days = 0, return = "ll") {
  if (return == "full") {
    save_particles <- TRUE
    pf_return <- "sample"
  } else if (return == "ll") {
    save_particles <- FALSE
    forecast_days <- 0
    pf_return <- "single"
  } else {
    stop("Unknown return type to calc_loglikelihood")
  }
  
  # defaults if not being sampled
  beta_start <- NULL
  beta_end <- NULL
  beta_pl <- NULL
  start_date <- data$date[1]
  
  # Update particle filter parameters from pars
  for (par in names(pars)) {
    if (par == "start_date") {
      start_date <- offset_to_start_date(data$date[1], pars[[par]])
    } else if (par == "beta_start") {
      beta_start <- pars[[par]]
    } else if (par == "beta_end") {
      beta_end <- pars[[par]]
    } else if (par == "beta_pl") {
      beta_pl <- pars[[par]]
    } else if (par %in% names(model_params)) {
      model_params[[par]] <- pars[[par]]
    } else if (par %in% names(pars_obs)) {
      pars_obs[[par]] <- pars[[par]]
    } else {
      stop(paste0("Don't know how to update parameter: ", par))
    }
  }

  # Beta needs a transform applied
  new_beta <- update_beta(sircovid_model, 
                          beta_start, 
                          beta_end, 
                          beta_pl,
                          start_date,
                          model_params$dt)
  model_params$beta_y <- new_beta$beta_y
  model_params$beta_t <- new_beta$beta_t

  pf_result <- run_particle_filter(data = data,
                                   sircovid_model = sircovid_model,
                                   model_params = model_params,
                                   model_start_date = start_date,
                                   obs_params = pars_obs,
                                   pars_seeding = NULL,
                                   n_particles = n_particles,
                                   forecast_days = forecast_days,
                                   save_particles = save_particles,
                                   return = pf_return)
  pf_result
}

propose_parameters <- function(pars, proposal_kernel, pars_discrete, pars_min, pars_max) {
  
  ## proposed jumps are normal with mean pars and sd as input for parameter
  jumps <- pars + drop(rmvnorm(n = 1,  sigma = proposal_kernel[names(pars), names(pars)]))


  # discretise if necessary
  jumps[pars_discrete] <- round(jumps[pars_discrete])
  # reflect proposal if it exceeds upper or lower parameter boundary
  jumps <- reflect_proposal(x = jumps, 
                          floor = pars_min, 
                          cap = pars_max)
  jumps
}


## create function to reflect proposal boundaries at pars_min and pars_max
# this ensures the proposal is symetrical and we can simplify the MH step

reflect_proposal <- function(x, floor, cap) {
  interval <- cap - floor
  abs((x + interval - floor) %% (2 * interval) - interval) + floor
}



##' @title create a master chain from a pmcmc_list object
##' @param x a pmcmc_list object
##' @param burn_in an integer denoting the number of samples to discard from each chain
##' @export
##' 
create_master_chain <- function(x, burn_in) {
  
  if(class(x) != 'pmcmc_list') {
    stop('x must be a pmcmc_list object')
  }
  if(!is.numeric(burn_in)) {
    stop('burn_in must be an integer')
  }
  if(burn_in < 0) {
    stop('burn_in must not be negative')
  }
  if(burn_in >= x$inputs$n_mcmc) {
    stop('burn_in is greater than chain length')
  }
  
  chains <- lapply(
    X = x$chains,
    FUN = function(z) z$results[-seq_len(burn_in), ]
  )
  
  do.call(what = rbind, args = chains)
}


##' @export
##' @importFrom stats cor sd 
summary.pmcmc <- function(object, ...) {
  
  par_names <- names(object$inputs$pars$pars_init)
  
  ## convert start_date to numeric to calculate stats
  data_start_date <- as.Date(object$inputs$data$date[1])
  traces <- object$results[,par_names] 
  traces$start_date <- start_date_to_offset(data_start_date, traces$start_date)
  
  # calculate correlation matrix
  corr_mat <- round(cor(traces),2)
  
  
  # add in reduction in beta parameter
  if('beta_end' %in% par_names) {
    traces <- cbind(traces, 
                    beta_red = traces$beta_end/traces$beta_start)
  }
  
  # compile summary
  summ <- rbind(mean = colMeans(traces),
                apply(traces, MARGIN = 2, quantile, c(0.025, 0.975)), 
                min = apply(traces, MARGIN = 2, min),
                max =  apply(traces, MARGIN = 2, max)
  )
  summ <- as.data.frame(summ)
  summ <- round(summ, 3)
  
  
  sds <- round(apply(traces, 2, sd), 3)
  # convert start_date back into dates
  summ$start_date <- as.Date(-summ$start_date, data_start_date)
  summ[c('2.5%', '97.5%', 'min', 'max'), 'start_date'] <- summ[c('97.5%', '2.5%', 'max', 'min'), 'start_date']

  
  out <- list('summary' = summ, 
              'corr_mat' = corr_mat, 
              'sd' = sds)
  out
  
}

##' @export
summary.pmcmc_list <- function(object, ..., burn_in = 101) {

  master_chain <- create_master_chain(x = object, 
                                      burn_in = burn_in)
  
  z <- list(inputs = object$inputs, 
            results = master_chain)
  summary.pmcmc(z)
  
}


##' @export
##' @importFrom viridis cividis
##' @importFrom graphics hist par plot.new text
plot.pmcmc <- function(x, ...) {
  
  summ <- summary(x)
  par_names <- names(x$inputs$pars$pars_init)
  
  traces <- x$results[, par_names]
  cols <- viridis::cividis(nrow(traces))
  cols <- cols[order(order(x$results$log_likelihood))]
  
  
  
  par_name <- 'beta_start'
  print_summ <- function(par_name) {
    x <- summ$summary
    paste0(x['mean', par_name], 
           '\n(', 
           x['2.5%', par_name], 
           ', ', 
           x['97.5%', par_name], ')')
  }
  
  
  
  n_pars <- length(par_names)
  
  
  par( bty = 'n',
       mfcol = c(n_pars, n_pars + 1L),
       mar = c(3,3,2,1),
       mgp = c(1.5, 0.5, 0), 
       oma = c(1,1,1,1))
  
  
  for (i in seq_len(n_pars)) {
    for(j in seq_len(n_pars)) {
      
      if (i == j) { # plot hists on diagonal
        par_name <- par_names[i]
        breaks = ifelse(par_name == 'start_date', 
                        yes = seq(as.Date('2019-12-01'), 
                                  as.Date(x$inputs$data$date[1]), 7),
                        no = 10)
        hist(traces[[i]], 
             main = print_summ(par_name), 
             xlab = par_name, 
             breaks = breaks,
             cex.main = 1,
             font.main = 1,
             freq = FALSE)
      } else if (i < j) {  # plot correlations on lower triangle
        plot(x = traces[[i]], 
             y = traces[[j]], 
             xlab = par_names[i], 
             ylab = par_names[j], 
             col = cols, 
             pch = 20)
      } else if (i > j) { # print rho on upper triangle
        plot.new()
        text(x = 0.5, 
             y=0.5, 
             labels = paste('r =', 
                            summ$corr_mat[i, j]))
      }
    }
  }
  
  # print traces in final column
  mapply(FUN = plot, traces, 
         type = 'l', 
         ylab = par_names, 
         xlab = "Iteration")
  
  
}

##' @export
##' @importFrom viridis cividis
##' @importFrom graphics hist par plot.new text lines legend
##' 
plot.pmcmc_list <- function(x, burn_in = 1, ...) {
  
  summ <- summary(x, burn_in = burn_in)
  par_names <- names(x$inputs$pars$pars_init)
  n_pars <- length(par_names)
  
  chains <- x$chains
  n_chains <- length(chains)
  cols_trace <- rev(viridis::viridis(n_chains))
  
  
  # compile master chain and order by log posterior for plotting
  master_chain <- create_master_chain(x, burn_in = burn_in)

  master_chain <- master_chain[order(master_chain$log_posterior), ]
  cols <- viridis::cividis(nrow(master_chain))
  cols <- cols[order(master_chain$log_posterior)]
  
  
  
  
  traces <- lapply(par_names, FUN = function(par_name) {
    lapply(X = chains, 
           FUN = function(z) z$results[-seq_len(burn_in), par_name])
  })
  names(traces) <- par_names
  
  plot_traces <- function(trace, col) {
    lines(x = seq_along(trace), 
          y = trace, 
          col = col)
  }
  
  
  
  breaks <- lapply(par_names, function(par_name){
    seq(from = min(master_chain[, par_name]), 
        to =  max(master_chain[, par_name]), 
        length.out = 20)
  })
  names(breaks) <- par_names
  
  hists <- lapply(par_names, FUN = function(par_name) {
    lapply(X = traces[[par_name]], 
           FUN = hist, 
           plot = FALSE, 
           breaks = breaks[[par_name]])
  })
  names(hists) <- par_names
  
  hist_ylim <- lapply(hists, function(h) {
    chain_max <- sapply(h, function(chain) max(chain$density) )
    c(0, max(chain_max))
  })
  
  plot_hists <- function(h, col, breaks) {
    with(h, lines(x =  breaks, 
                  y = c(density, 
                        density[length(density)]), 
                  type = 's',
                  col = col))
  }
  
  
  print_summ <- function(par_name) {
    x <- summ$summary
    paste0(x['mean', par_name], 
           '\n(', 
           x['2.5%', par_name], 
           ', ', 
           x['97.5%', par_name], ')')
  }
  
  
  par( bty = 'n',
       mfcol = c(n_pars, n_pars + 1L),
       mar = c(3,3,2,1),
       mgp = c(1.5, 0.5, 0), 
       oma = c(1,1,1,1))
  
  
  for (i in seq_len(n_pars)) {
    for(j in seq_len(n_pars)) {
      
      if (i == j) { # plot hists on diagonal
        par_name <- par_names[i]
        bs <- breaks[[par_name]]
        plot(x = bs[1] ,  # force date axis where needed
             y = 1, 
             type = 'n',
             xlim = c(bs[1], bs[length(bs)]),
             ylim = hist_ylim[[par_name]],
             xlab = par_name, 
             ylab = '',
             main = print_summ(par_name), 
             cex.main = 1,
             font.main = 1
        )
        
        mapply(FUN = plot_hists, 
               h = hists[[par_name]],
               breaks = bs,
               col = cols_trace)
        
        
      } else if (i < j) {  # plot correlations on lower triangle
        plot(x = master_chain[[i]], 
             y = master_chain[[j]], 
             xlab = par_names[i], 
             ylab = par_names[j], 
             col = cols, 
             pch = 20)
      } else if (i > j) { # print rho on upper triangle
        plot.new()
        text(x = 0.5, 
             y=0.5, 
             labels = paste('r =', 
                            summ$corr_mat[i, j]))
      }
    }
  }
  
  
  # print traces in final column
  n_iter <- nrow(master_chain) / n_chains
  
  mapply(FUN = function(par_name, leg) {
    plot(x = 1,
         y = breaks[[par_name]][1],
         type = 'n', 
         xlab = 'Iteration', 
         ylab = par_name, 
         xlim = c(0, n_iter), 
         ylim <- range(master_chain[, par_name]))
    
    mapply(FUN = plot_traces,
           trace = traces[[par_name]], 
           col = cols_trace)
    
    if(leg) {
      legend('top',
             ncol = n_chains, 
             legend = paste('Chain', seq_len(n_chains)),
             fill = cols_trace, 
             bty = 'n')
    }
  }, 
  par_name = par_names, 
  leg = c(TRUE, FALSE, FALSE))
  
}