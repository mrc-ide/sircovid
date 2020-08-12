##' @name carehomes
##' @title The carehomes sircovid model
##'
##' Our "current" sircovid model. This is a dust model.
##'
##' @export carehomes
NULL


##' Parameters for the "[carehomes]" model.
##'
##' @title Parameters for the carehomes model
##'
##' @inheritParams basic_parameters
##'
##' @return A list of inputs to the model, many of which are fixed and
##'   represent data. These correspond largely to `user()` calls
##'   within the odin code, though some are also used in processing
##'   just before the model is run.
##'
##' @export
##' @examples
##' carehomes_parameters(sircovid_date("2020-02-01"), "uk")
carehomes_parameters <- function(start_date, region,
                                 beta_date = NULL, beta_value = NULL,
                                 severity_data = NULL,
                                 exp_noise = 1e6) {
  ret <- sircovid_parameters_shared(start_date, region,
                                    beta_date, beta_value)

  ## TODO: These should be flexible and will be set in the pmcmc so
  ## will move up to the argument list of this function; these are
  ## used only here in the setup and are not used in the model itself,
  ## which means we need to think how this interacts with the pmcmc
  p_death_carehome <- 0.7
  eps <- 0.1
  C_1 <- 4e-5
  C_2 <- 5e-4

  ## These are only used here, and are fixed
  carehome_occupancy <- 0.742
  carehome_workers_per_resident <- 1

  ## These are used in constructing the initial population vectors (S0)
  carehome_beds <- sircovid_carehome_beds(region)
  carehome_residents <- round(carehome_beds * carehome_occupancy)
  carehome_workers <- round(carehome_residents * carehome_workers_per_resident)

  ## TODO: it's probably the case that having some tree structure here
  ## would make this nicer to work with, but we should do this
  ## consistently through the other parameters too. Keeping the
  ## progression and severity parameters together for example.
  ret$carehome_beds <- carehome_beds
  ret$carehome_residents <- carehome_residents
  ret$carehome_workers <- carehome_workers

  severity <- carehomes_parameters_severity(severity_data, ret$population,
                                            p_death_carehome)

  ret$m <- carehomes_transmission_matrix(eps, C_1, C_2, region, ret$population)

  ret$N_tot <- carehomes_population(ret$population, carehome_workers,
                                    carehome_residents)

  ## This is used to normalise the serology counts (converting them
  ## from number of positive/negative tests into a fraction). This is
  ## constant over the simulation, being the total population size of
  ## 15 to 64 year olds.
  N_tot_15_64 <- sum(ret$N_tot[4:13])

  ## All observation parameters:
  ret$observation <- carehomes_parameters_observation(exp_noise, N_tot_15_64)

  ## TODO: Adding this here, but better would be to pass N_age as-is,
  ## then update the leading dimension to something more accurate
  ## (e.g., N_groups, setting this as N_groups <- N_age + 2)
  ret$N_age <- ret$N_age + 2L

  c(ret,
    severity,
    carehomes_parameters_progression())
}


##' Index of "interesting" elements for the carehomes model. This function
##' conforms to the mcstate interface.
##'
##' @title Index of carehomes model
##'
##' @inheritParams basic_index
##'
##' @return A list with element `run`, indicating the locations of (in
##'   order) (1) ICU, (2) general, (3) deaths in community, (4) deaths
##'   in hospital, (5) total deaths, (6) cumulative confirmed
##'   admissions,(7) cumulative confirmed new admissions, (8)
##'   recovered pre seroconversion (15 - 64 year-olds only), (9)
##'   recovered seronegative and (10) recovered seropositive.
##'
##' @export
##' @examples
##' p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
##' mod <- carehomes$new(p, 0, 10)
##' carehomes_index(mod$info())
carehomes_index <- function(info) {
  index <- info$index
  list(run = c(icu = index[["I_ICU_tot"]],
               general = index[["general_tot"]],
               deaths_comm = index[["D_comm_tot"]],
               deaths_hosp = index[["D_hosp_tot"]],
               deaths_tot = index[["D_tot"]],
               admitted = index[["cum_admit_conf"]],
               new = index[["cum_new_conf"]],
               R_pre_15_64 = index[["R_pre_15_64"]],
               R_neg_15_64 = index[["R_neg_15_64"]],
               R_pos_15_64 = index[["R_pos_15_64"]]))
}


##' Compare observed and modelled data from the [carehomes] model. This
##' conforms to the mcstate interface.
##'
##' @title Compare observed and modelled data for the carehomes model
##'
##' @param state State vector for the end of the current day. This is
##'   assumed to be filtered following [carehomes_index()] so contains
##'   10 rows corresponding to itu, admissions, deaths and
##'   seroconversion compartments.
##'
##' @param prev_state State vector for the end of the previous day, as
##'   for `state`.
##'
##' @param observed Observed data. This will be a list with elements
##'   `itu` (number of itu/icu beds occupied), `general` (number of
##'   general beds occupied), `deaths_hosp` (cumulative deaths in
##'   hospital settings), `deaths_comm` (cumulative deaths in
##'   community settings), `deaths` (combined deaths - used by some
##'   nations - if given then `deaths_hosp` and `deaths_comm` must be
##'   `NA`), `admitted` (hospital admissions), `new` (new cases),
##'   `npos_15_64` (number of people seropositive in ages 15-64) and
##'   `ntot_15_64` (number of people tested in ages 15-64).
##'
##'
##' @param pars A list of parameters, as created by
##'   [carehomes_parameters()]
##'
##' @return A vector of log likelihoods, the same length as the number
##'   of particles (the number of columns in the modelled state)
##'
##' @export
##' @importFrom stats dbinom
carehomes_compare <- function(state, prev_state, observed, pars) {
  ## TODO: we might refactor this to produce a subset of comparisons
  ## (and indices above) to suit either the SPI-M or paper fits as
  ## we're using different streams; that will make the comparisons a
  ## touch faster and data copying smaller. This requires a bit of a
  ## rethink though as it affects how we do index functions too.

  model_icu <- state["icu", ]
  model_general <- state["general", ]
  model_deaths_comm <- state["deaths_comm", ] - prev_state["deaths_comm", ]
  model_deaths_hosp <- state["deaths_hosp", ] - prev_state["deaths_hosp", ]
  model_deaths_tot <- state["deaths_tot", ] - prev_state["deaths_tot", ]
  model_admitted <- state["admitted", ] - prev_state["admitted", ]
  model_new <- state["new", ] - prev_state["new", ]
  model_R_pre_15_64 <- state["R_pre_15_64", ]
  model_R_neg_15_64 <- state["R_neg_15_64", ]
  model_R_pos_15_64 <- state["R_pos_15_64", ]

  pars <- pars$observation
  exp_noise <- pars$exp_noise

  ll_itu <- ll_nbinom(observed$itu, pars$phi_ICU * model_icu,
                      pars$k_ICU, exp_noise)
  ll_general <- ll_nbinom(observed$general, pars$phi_general * model_general,
                          pars$k_general, exp_noise)
  ll_deaths_hosp <- ll_nbinom(observed$deaths_hosp,
                              pars$phi_death_hosp * model_deaths_hosp,
                              pars$k_death_hosp, exp_noise)
  ll_deaths_comm <- ll_nbinom(observed$deaths_comm,
                              pars$phi_death_comm * model_deaths_comm,
                              pars$k_death_comm, exp_noise)
  ll_deaths <- ll_nbinom(observed$deaths,
                         pars$phi_death_hosp * model_deaths_hosp +
                         pars$phi_death_comm * model_deaths_comm,
                         pars$k_death, exp_noise)
  ll_admitted <- ll_nbinom(observed$admitted,
                           pars$phi_admitted * model_admitted,
                           pars$k_admitted, exp_noise)
  ll_new <- ll_nbinom(observed$new, pars$phi_new * model_new,
                      pars$k_new, exp_noise)

  ## TODO: it would be easy to return the true_pos and positive tests
  ## as two numbers rather than these three from the odin code
  true_pos <- model_R_pos_15_64 + model_R_neg_15_64 + model_R_pre_15_64
  prob_true_pos <- model_R_pos_15_64 / pars$N_tot_15_64
  prob_false_pos <- (1 - pars$p_specificity) * (1 - true_pos / pars$N_tot_15_64)

  ll_serology <- ll_binom(observed$npos_15_64, 
                          observed$ntot_15_64,
                          prob_true_pos + prob_false_pos)
  
  ll_itu + ll_general + ll_deaths_hosp + ll_deaths_comm + ll_deaths +
    ll_admitted + ll_new + ll_serology
}


## We store within the severity parameters information on severity for
## carehome workers and residents. The vector ends up structured as
##
##   [1..N_age, workers, residents]
##
## so we have length of N_age + 2
##' @importFrom stats weighted.mean
carehomes_severity <- function(p, population) {
  index_workers <- carehomes_index_workers()
  p_workers <- weighted.mean(p[index_workers], population[index_workers])
  p_residents <- p[length(p)]
  c(p, p_workers, p_residents)
}


carehomes_parameters_severity <- function(severity_data, population,
                                          p_death_carehome) {
  severity <- sircovid_parameters_severity(severity_data)
  severity <- lapply(severity, carehomes_severity, population)
  severity$p_death_comm[length(severity$p_death_comm)] <- p_death_carehome
  severity
}


carehomes_index_workers <- function() {
  age_bins <- sircovid_age_bins()
  which(age_bins$start >= 25 & age_bins$start < 65)
}


carehomes_transmission_matrix <- function(eps, C_1, C_2, region, population) {
  index_workers <- carehomes_index_workers()
  m <- sircovid_transmission_matrix(region)
  N_age <- nrow(m)

  m_chw <- apply(m[seq_len(N_age), index_workers], 1, weighted.mean,
                 population[index_workers])
  m_chr <- eps * m[N_age, seq_len(N_age)]

  ## Construct a block matrix:
  ##
  ##   M     m_chw m_chr
  ##   m_chw C_1   C_1
  ##   m_chr C_1   C_2

  i <- seq_len(N_age)
  i_chw <- N_age + 1L
  i_chr <- N_age + 2L

  ret <- matrix(0.0, N_age + 2, N_age + 2)
  ret[i, i] <- m
  ret[i, i_chw] <- ret[i_chw, i] <- m_chw
  ret[i, i_chr] <- ret[i_chr, i] <- m_chr
  ret[i_chw:i_chr, i_chw:i_chr] <- c(C_1, C_1, C_1, C_2)

  nms <- c(rownames(m), "CHW", "CHR")
  dimnames(ret) <- list(nms, nms)

  ret
}


##' Create initial conditions for the carehomes model. This matches the
##' interface required for mcstate
##'
##' @title Initial conditions for the carehomes model
##'
##' @inheritParams basic_initial
##'
##' @return A numeric vector of initial conditions
##' @export
##' @examples
##' p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
##' mod <- carehomes$new(p, 0, 10)
##' carehomes_initial(mod$info(), 10, p)
carehomes_initial <- function(info, n_particles, pars) {
  index <- info$index
  state <- numeric(info$len)

  ## Always start with 10, again for compatibility
  initial_I <- 10

  ## This corresponds to the 15-19y age bracket for compatibility with
  ## our first version, will be replaced by better seeding model, but
  ## probably has limited impact.
  seed_age_band <- 4L
  index_I <- index[["I_asympt"]][[1L]] + seed_age_band - 1L
  index_R_pre <- index[["R_pre"]][[1L]] + seed_age_band - 1L
  index_PCR_pos <- index[["PCR_pos"]][[1L]] + seed_age_band - 1L
  index_N_tot2 <- index[["N_tot2"]][[1L]]

  index_S <- index[["S"]]
  index_N_tot <- index[["N_tot"]]

  ## S0 is the population totals, minus the seeded infected
  ## individuals
  initial_S <- pars$N_tot
  initial_S[seed_age_band] <- initial_S[seed_age_band] - initial_I

  state[index_S] <- initial_S
  state[index_I] <- initial_I
  state[index_R_pre] <- initial_I
  state[index_PCR_pos] <- initial_I
  state[index_N_tot] <- pars$N_tot
  state[index_N_tot2] <- sum(pars$N_tot)

  list(state = state,
       step = pars$initial_step)
}

carehomes_parameters_progression <- function() {
  ## These need to be aligned with Bob's severity outputs, and we will
  ## come up with a better way of correlating the two.

  ## The s_ parameters are the scaling parameters for the Erlang
  ## distibution (a.k.a 'k'), while the gamma parameters are the gamma
  ## parameters of that distribution.
  list(s_E = 2,
       s_asympt = 1,
       s_mild = 1,
       s_ILI = 1,
       s_comm_D = 2,
       s_hosp_D = 2,
       s_hosp_R = 2,
       s_ICU_D = 2,
       s_ICU_R = 2,
       s_triage = 2,
       s_stepdown = 2,
       s_PCR_pos = 2,

       gamma_E = 1 / (4.59 / 2),
       gamma_asympt = 1 / 2.09,
       gamma_mild = 1 / 2.09,
       gamma_ILI = 1 / 4,
       gamma_comm_D = 2 / 5,
       gamma_hosp_D = 2 / 5,
       gamma_hosp_R = 2 / 10,
       gamma_ICU_D = 2 / 5,
       gamma_ICU_R = 2 / 10,
       gamma_triage = 2,
       gamma_stepdown = 2 / 5,
       gamma_R_pre_1 = 1 / 5,
       gamma_R_pre_2 = 1 / 10,
       gamma_test = 3 / 10,
       gamma_PCR_pos = 1 / 5)
}


sircovid_carehome_beds <- function(region) {
  if (is.null(region)) {
    stop("'region' must not be NULL")
  }

  if (is.null(cache$carehomes)) {
    cache$carehomes <- read_csv(sircovid_file("extdata/carehomes.csv"))
  }

  i <- match(tolower(region), cache$carehomes$region)
  if (is.na(i)) {
    valid <- paste(squote(cache$carehomes$region), collapse = ", ")
    stop(sprintf("Carehome beds not found for '%s': must be one of %s",
                 region, valid))
  }

  cache$carehomes$carehome_beds[[i]]
}


carehomes_parameters_observation <- function(exp_noise, N_tot_15_64) {
  list(
    ## People currently in ICU
    phi_ICU = 0.95,
    k_ICU = 2,
    ## People currently in general beds
    phi_general = 0.95,
    k_general = 2,
    ## Daily hospital deaths
    phi_death_hosp = 1.15,
    k_death_hosp = 2,
    ## Daily community deaths
    phi_death_comm = 1,
    k_death_comm = 2,
    ## Daily total deaths (if not split)
    k_death = 2,
    ## Daily new confirmed admissions
    phi_admitted = 0.95,
    k_admitted = 2,
    ## Daily new inpatient diagnoses
    phi_new = 0.95,
    k_new = 2,
    ## Specificity for serology tests
    ##
    ## TODO: p_specificity here needs to be tuneable, as that will be
    ## fit within the mcmc
    p_specificity = 0.9,
    ## Population size of eligible 15-64 year olds for serology testing
    N_tot_15_64 = N_tot_15_64,
    ## rate for exponential noise, generally something big so noise is
    ## small (but non-zero))
    exp_noise = exp_noise)
}


carehomes_population <- function(population, carehome_workers,
                                 carehome_residents) {
  ## Our core S0 calculation is more complicated than the basic model
  ## because we have to add the carehome workers and residents, *and*
  ## remove them from the core population. This extracts carehome
  ## residents from the older groups of the population, weighted
  ## towards the oldest, and extracts carehome workers from most
  ## working ages, evenly across the population.
  N_tot <- c(population, carehome_workers, carehome_residents)

  index_workers <- carehomes_index_workers()
  weights_workers <- N_tot[index_workers] / sum(N_tot[index_workers])
  index_residents <- which(sircovid_age_bins()$start >= 65)
  weights_residents <- c(0.05, 0.05, 0.15, 0.75)

  N_tot[index_residents] <-
    round(N_tot[index_residents] - carehome_residents * weights_residents)
  N_tot[index_workers] <-
    round(N_tot[index_workers] - carehome_workers * weights_workers)

  if (any(N_tot[index_residents] < 0)) {
    stop("Not enough population to meet care home occupancy")
  }
  if (any(N_tot[index_workers] < 0)) {
    stop("Not enough population to be care workers")
  }
  N_tot
}
