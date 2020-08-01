##' Parameters for the "basic" model.
##'
##' @title Parameters for the basic model
##'
##' @param start_date Thee start date, as a [sircovid_date()] (i.e.,
##'   the number of days into 2020)
##'
##' @param region The region to run the model for. This will bee used
##'   to get population data, which is currently fixed within the
##'   package and is limited to "uk", the four constituent nations
##'   ("england", "wales", "scotland", "northern_ireland") and thee 7
##'   NHS regions (e.g., "midlands"). These names are case
##'   insensitive.
##'
##' @param beta_date A vector of dates for changes in beta (the
##'   contact rate parameter), or `NULL` if a single value is used for
##'   all times (see [sircovid_parameters_beta()], where this is
##'   passed as `date`).
##'
##' @param beta_value A vector of values for beta (the contact rate
##'   parameter). If not given, and if `beta_date` is null then a
##'   value of 0.08 will be used through the whole simulation,
##'   otherwise if `beta_date` is `NULL` this must be a scalar. If
##'   `beta_date` is given then `beta_date` and `beta_value` must have
##'   the same length (see [sircovid_parameters_beta()], where this is
##'   passed as `value`).
##'
##' @param severity_data A `data.frame` of severity data, or `NULL` to
##'   use the default value within the package.  New severity data
##'   comes from Bob Verity via the markovid package, and needs to be
##'   carefully calibrated with the progression parameters.
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
                             severity_data = NULL) {
  ret <- sircovid_parameters_shared(start_date, region,
                                    beta_date, beta_value)
  ret$m <- sircovid_transmission_matrix()
  c(ret,
    sircovid_parameters_severity(severity_data),
    basic_parameters_progression())
}


basic_index <- function(info) {
  ## TODO: this will simplify once we get the index here, see
  ## odin.dust issue #24
  len <- vnapply(info, prod)
  start <- cumsum(len) - len + 1L
  list(run = c(start[["I_ICU_tot"]], start[["D_tot"]]))
}


basic_compare <- function(state, prev_state, observed, pars) {
  if (is.na(observed$itu) && is.na(observed$deaths)) {
    return(NULL)
  }

  ## TODO: tidy up in mcstate to pull index over - see mcstate issue #35
  model_icu <- state[1, ]
  model_deaths <- state[2, ] - prev_state[2, ]

  ## Noise parameter shared across both deaths and icu
  exp_noise <- pars$exp_noise

  ll_itu <- ll_nbinom(observed$itu, pars$phi_ICU * model_icu,
                      pars$k_ICU, exp_noise)
  ll_deaths <- ll_nbinom(observed$deaths, pars$phi_death * model_deaths,
                         pars$k_death, exp_noise)

  ll_itu + ll_deaths
}


basic_initial <- function(info, n_particles, pars) {
  ## TODO: this will simplify once we get the index here, see
  ## odin.dust issue #24
  len <- vnapply(info, prod)
  start <- cumsum(len) - len + 1L
  state <- numeric(sum(len))

  ## Always start with 10, again for compatibility
  initial_I <- 10

  ## This corresponds to the 15-19y age bracket for compatibility with
  ## our first version, will be replaced by better seeding model, but
  ## probably has limited impact.
  seed_age_band <- 4L
  index_I <- start[["I_asympt"]] + seed_age_band

  ## ONS populations, subtracting the seed for pedantry.
  index_S <- seq.int(start[["S"]], length.out = len[["S"]])
  initial_S <- pars$population
  initial_S[seed_age_band] <- initial_S[seed_age_band] - initial_I

  index_N_tot <- start[["N_tot"]]

  state[index_S] <- initial_S
  state[index_I] <- initial_I
  state[index_N_tot] <- sum(pars$population)

  list(state = state,
       step = pars$initial_step)
}


basic_parameters_progression <- function() {
  ## These need to be aligned with Bob's severity outputs, and we will
  ## come up with a better way of correlating the two.

  ## The s_ parameters are the scaling parameters for the Erlang
  ## distibution (a.k.a 'k'), while the gamma parameters are the gamma
  ## parameters of that distribution.
  list(s_E = 2,
       s_asympt = 1,
       s_mild = 1,
       s_ILI = 1,
       s_hosp = 2,
       s_ICU = 2,
       s_rec = 2,

       gamma_E = 1 / (4.59 / 2),
       gamma_asympt = 1 / 2.09,
       gamma_mild = 1 / 2.09,
       gamma_ILI = 1 / 4,
       gamma_hosp = 2,
       gamma_ICU = 2 / 5,
       gamma_rec = 2 / 5)
}


basic_parameters_observation <- function(exp_noise = 1e6) {
  list(
    ## People currently in general beds
    phi_general = 0.95,
    k_general = 2,
    ## People currently in ICU
    phi_ICU = 0.95,
    k_ICU = 2,
    ## Daily deaths
    ##
    ## current proportion of England deaths over UK deaths (as of
    ## end of March 2020)
    phi_death = 926 / 1019,
    k_death = 2,
    ## rate for exponential noise, generally something big so noise is
    ## small (but non-zero))
    exp_noise = exp_noise)
}
