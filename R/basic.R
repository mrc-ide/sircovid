##' @name basic
##' @title The basic sircovid model
##'
##' The "basic" sircovid model. This is a dust model.
##' @export basic
##' @examples
##' # Set up the basic model for England with default parameters and
##' # an initial seeding of early February:
##' p <- basic_parameters(sircovid_date("2020-02-07"), "england")
##' mod <- basic$new(p, 0, 10)
##'
##' # Set the initial state and index as we would us for a run
##' # (without setting an initial state there is no seeding)
##' initial <- basic_initial(mod$info(), 10, p)
##' mod$set_state(initial$state, initial$step)
##' mod$set_index(basic_index(mod$info())$run)
##'
##' # Run the model up to the end of march
##' step_end <- sircovid_date("2020-03-31") / p$dt
##'
##' # The filtered state is returned at the end of the run
##' mod$run(step_end)
##'
##' # More state can be retrieved using the "state" method
##' mod$state(1:6)
NULL

##' Parameters for the "[basic]" model.
##'
##' @title Parameters for the basic model
##'
##' @param start_date Thee start date, as a [sircovid_date()] (i.e.,
##'   the number of days into 2020)
##'
##' @param region The region to run the model for. This will bee used
##'   to get population data, which is currently fixed within the
##'   package and is limited to "uk", the four constituent nations
##'   ("england", "wales", "scotland", "northern_ireland") and the 7
##'   NHS regions (e.g., "midlands"). These names are case
##'   insensitive.
##'
##' @param beta_date A vector of dates for changes in beta (the
##'   contact rate parameter), or `NULL` if a single value is used for
##'   all times (see [sircovid_parameters_beta()], where this is
##'   passed as `date`).
##'
##' @param beta_value A vector of values for beta (the contact rate
##'   parameter). If not given, and if `beta_date` is `NULL` then a
##'   value of 0.08 will be used through the whole simulation,
##'   otherwise if `beta_date` is `NULL` this must be a scalar. If
##'   `beta_date` is given then `beta_date` and `beta_value` must have
##'   the same length (see [sircovid_parameters_beta()], where this is
##'   passed as `value`).
##'
##' @param severity Severity data, via Bob Verity's `markovid`
##'   package. This needs to be `NULL` (use the default bundled data
##'   version in the package), a [data.frame] object (for raw severity
##'   data) or a list (for data that has already been processed by
##'   `sircovid` for use).  New severity data comes from Bob Verity
##'   via the markovid package, and needs to be carefully calibrated
##'   with the progression parameters.
##'
##' @param exp_noise Rate of exponential noise used in the compare
##'   function - typically set to a large value so that noise is small
##'   but non-zero. If set to `Inf` then there is no noise in the
##'   observation process (not realistic but useful for testing).
##'
##' @return A list of inputs to the model, many of which are fixed and
##'   represent data. These correspond largely to `user()` calls
##'   within the odin code, though some are also used in processing
##'   just before the model is run.
##'
##' @export
##' @examples
##' basic_parameters(sircovid_date("2020-02-01"), "uk")
basic_parameters <- function(start_date, region,
                             beta_date = NULL, beta_value = NULL,
                             severity = NULL,
                             exp_noise = 1e6) {
  ret <- sircovid_parameters_shared(start_date, region,
                                    beta_date, beta_value)
  ret$m <- sircovid_transmission_matrix(region)
  ret$observation <- basic_parameters_observation(exp_noise)
  severity <- sircovid_parameters_severity(severity)

  ## Some additional processing of derived quantities used in the
  ## basic model; we could do these transformations in the odin code
  ## perhaps, which would reduce carrying redundant dependencies
  ## between parameters
  severity$p_recov_ICU <- 1 - severity[["p_ICU_D"]]
  severity$p_recov_hosp <-
    (1 - severity[["p_ICU"]]) *
    (1 - severity[["p_H_D"]])
  severity$p_death_hosp <-
    (1 - severity[["p_ICU"]]) *
    severity[["p_H_D"]]
  severity$p_recov_sympt <- 1 - severity[["p_H"]]

  c(ret,
    severity,
    basic_parameters_progression())
}


##' Index of "interesting" elements for the basic model. This function
##' conforms to the mcstate interface.
##' @title Index of basic model
##'
##' @param info The result of running the `$info()` method on an
##'   initialised [basic] model
##'
##' @return A list with element `run`, indicating the locations of the
##'   ICU and D compartments.
##'
##' @export
##' @examples
##' p <- basic_parameters(sircovid_date("2020-02-07"), "england")
##' mod <- basic$new(p, 0, 10)
##' basic_index(mod$info())
basic_index <- function(info) {
  index <- info$index
  list(run = c(icu = index[["I_ICU_tot"]],
               deaths = index[["D_tot"]]))
}


##' Compare observed and modelled data from the basic model. This
##' conforms to the mcstate interface.
##'
##' @title Compare observed and modelled data for the basic model
##'
##' @param state State vector for the end of the current day. This is
##'   assumed to be filtered following [basic_index()] so contains
##'   rows corresponding to ICU and deaths.
##'
##' @param prev_state State vector for the end of the previous day, as
##'   for `state`.
##'
##' @param observed Observed data. This will be a list with elements
##'   `icu` (number of ICU beds occupied) and `deaths` (number of
##'   deaths over this day).
##'
##' @param pars A list of parameters, as created by
##'   [basic_parameters()]
##'
##' @return A vector of log likelihoods, the same length as the number
##'   of particles (the number of columns in the modelled state)
##'
##' @export
##' @examples
##' state <- rbind(icu = 10:15, deaths = 1:6)
##' prev_state <- matrix(1, 2, 6, dimnames = dimnames(state))
##' observed <- list(icu = 13, deaths = 3)
##' pars <- basic_parameters(sircovid_date("2020-02-07"), "england")
##' basic_compare(state, prev_state, observed, pars)
##' basic_compare(state * 5, prev_state, observed, pars)
basic_compare <- function(state, prev_state, observed, pars) {
  if (is.na(observed$icu) && is.na(observed$deaths)) {
    return(NULL)
  }

  model_icu <- state["icu", ]
  model_deaths <- state["deaths", ] - prev_state["deaths", ]

  ## Noise parameter shared across both deaths and icu
  pars <- pars$observation
  exp_noise <- pars$exp_noise

  ll_icu <- ll_nbinom(observed$icu, pars$phi_ICU * model_icu,
                      pars$kappa_ICU, exp_noise)
  ll_deaths <- ll_nbinom(observed$deaths, pars$phi_death * model_deaths,
                         pars$kappa_death, exp_noise)

  ll_icu + ll_deaths
}


##' Create initial conditions for the basic model. This matches the
##' interface required for mcstate
##'
##' @title Initial conditions for the basic model
##'
##' @param info The result of running the `$info()` method on an
##'   initialised [basic] model
##'
##' @param n_particles The number of particles required. Currently
##'   only uniform initial seeding is implemented so this has no
##'   effect
##'
##' @param pars A parameter list created by [basic_parameters()]; from
##'   this list we will use the `population` and `initial_step`
##'   elements.
##'
##' @return A numeric vector of initial conditions
##' @export
##' @examples
##' p <- basic_parameters(sircovid_date("2020-02-07"), "england")
##' mod <- basic$new(p, 0, 10)
##' basic_initial(mod$info(), 10, p)
basic_initial <- function(info, n_particles, pars) {
  index <- info$index
  state <- numeric(info$len)

  ## Always start with 10, again for compatibility
  initial_I <- 10

  ## This corresponds to the 15-19y age bracket for compatibility with
  ## our first version, will be replaced by better seeding model, but
  ## probably has limited impact.
  seed_age_band <- 4L
  index_I <- index[["I_A"]][[1]] + seed_age_band - 1L

  ## ONS populations, subtracting the seed for pedantry.
  index_S <- index[["S"]]
  initial_S <- pars$population
  initial_S[seed_age_band] <- initial_S[seed_age_band] - initial_I

  index_N_tot <- index[["N_tot"]]

  state[index_S] <- initial_S
  state[index_I] <- initial_I
  state[index_N_tot] <- sum(pars$population)

  list(state = state,
       step = pars$initial_step)
}


basic_parameters_progression <- function() {
  ## These need to be aligned with Bob's severity outputs, and we will
  ## come up with a better way of correlating the two.

  ## The k_ parameters are the shape parameters for the Erlang
  ## distibution, while the gamma parameters are the gamma
  ## rate parameters of that distribution.
  list(k_E = 2,
       k_A = 1,
       k_C = 1,
       k_hosp = 2,
       k_ICU = 2,
       k_rec = 2,

       gamma_E = 1 / (4.59 / 2),
       gamma_A = 1 / 2.09,
       gamma_C = 1 / 4,
       gamma_hosp = 2,
       gamma_ICU = 2 / 5,
       gamma_rec = 2 / 5)
}


basic_parameters_observation <- function(exp_noise) {
  list(
    ## People currently in general beds
    phi_general = 0.95,
    kappa_general = 2,
    ## People currently in ICU
    phi_ICU = 0.95,
    kappa_ICU = 2,
    ## Daily deaths
    ##
    ## current proportion of England deaths over UK deaths (as of
    ## end of March 2020)
    phi_death = 926 / 1019,
    kappa_death = 2,
    ## rate for exponential noise, generally something big so noise is
    ## small (but non-zero))
    exp_noise = exp_noise)
}
