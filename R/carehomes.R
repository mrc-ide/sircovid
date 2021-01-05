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
##' @param severity Severity data, via Bob Verity's `markovid`
##'   package. This needs to be `NULL` (use the default bundled data
##'   version in the package), a [data.frame] object (for raw severity
##'   data) or a list (for data that has already been processed by
##'   `sircovid` for use).  New severity data comes from Bob Verity
##'   via the markovid package, and needs to be carefully calibrated
##'   with the progression parameters.
##'
##' @param p_death_carehome Probability of death within carehomes
##'   (conditional on having an "inflenza-like-illness")
##'
##' @param sero_specificity Specificity of the serology test
##'
##' @param sero_sensitivity Sensitivity of the serology test
##'
##' @param progression Progression data
##'
##' @param eps Change in contact rate for carehome residents
##'
##' @param C_1 Contact rate between carehome workers and either
##'   residents or workers
##'
##' @param C_2 Contact rate between carehome residents
##'
##' @param pillar2_specificity Specificity of the Pillar 2 test
##'
##' @param pillar2_sensitivity Sensitivity of the Pillar 2 test
##'
##' @param react_specificity Specificity of the REACT test
##'
##' @param react_sensitivity Sensitivity of the REACT test
##'
##' @param prop_noncovid_sympt Proportion of population who do not have
##'   covid but have covid-like symptoms
##'
##' @param strain_transmission Vector of relative transmissibility of each
##'   strain modelled. First element should be 1. Length will define the
##'   number of strains used in the model
##'
##' @param strain_seed_date Either `NULL` (no seeding) or a vector of
##'   exactly two [sircovid::sircovid_date] values corresponding to
##'   the start and stop dates of seeding (inclusive). For example, to
##'   seed for a single day, these two values should have the same
##'   date
##'
##' @param strain_seed_value Either `NULL` (no seeding) or a single
##'   integer value represending the *daily* rate of seeding. Because
##'   each day is dividied into `1 / dt` steps, this value will be
##'   spread out fairly evenly across the steps that occur within each
##'   date. Seeding is not stochastic; this many individuals *will*
##'   become infected with the new strain, unless the pool of
##'   susceptibles has been exhausted.
##'
##' @param rel_susceptibility A vector or matrix of values representing the
##'   relative susceptibility of individuals in different vaccination groups.
##'   If a vector, the first value should be 1 (for the non-vaccinated group)
##'   and subsequent values be between 0 and 1. In that case relative
##'   susceptibility will be the same across all age groups within one
##'   vaccination category. Specifying a matrix instead of a vector allows
##'   different relative susceptibilities by age (rows of the matrix) and
##'   vaccination group (columns of the matrix); in that case, in each row of
##'   the matrix, the first value should be 1 (for the non-vaccinated group)
##'   and subsequent values be between 0 and 1
##'
##' @param rel_p_sympt A vector or matrix of values of same dimension as
##'   rel_susceptibility representing the
##'   relative probability of symptomatic infection in different
##'   vaccination groups. If a vector, the first value should be 1 (for the
##'   non-vaccinated group) and subsequent values be between 0 and 1.
##'   In that case the relative reduction in probability of symptomatic
##'   infection will be the same across all age groups within one vaccination
##'   category.
##'   Specifying a matrix instead of a vector allows different relative
##'   reductions in probability of symptomatic infection by age (rows of the
##'   matrix) and vaccination group (columns of the matrix); in that case,
##'   in each row of the matrix, the first value should be 1 (for the
##'   non-vaccinated group) and subsequent values be between 0 and 1
##'
##' @param rel_p_hosp_if_sympt A vector or matrix of values of same dimension as
##'   rel_susceptibility representing the
##'   relative probability of hospitalisation for symptomatic cases in different
##'   vaccination groups. If a vector, the first value should be 1 (for the
##'   non-vaccinated group) and subsequent values be between 0 and 1.
##'   In that case the relative reduction in probability of hospitalisation for
##'   symptomatic cases will be the same across all age groups within one
##'   vaccination category.
##'   Specifying a matrix instead of a vector allows different relative
##'   reductions in probability of hospitalisation for symptomatic cases by age
##'   (rows of the matrix) and vaccination group (columns of the matrix);
##'   in that case, in each row of the matrix, the first value should be 1
##'   (for the non-vaccinated group) and subsequent values be between 0 and 1
##'
##' @param vaccine_progression_rate A vector or matrix of values of same
##'   dimension as rel_susceptibility representing
##'   the rate of movement between different vaccination classes. If a vector,
##'   it should have as many values as vaccination classes, and the same rates
##'   of progression will be used for all age
##'   groups (the first rate is the vaccination rate and the last of the rates
##'   is the rate of returning to the initial
##'   vacination class); if a matrix, the element on row i and column j is the
##'   rate of progression from the jth vaccination class to the (j+1)th for age
##'   group i.
##'
##' @param vaccine_uptake A vector of length 19 with the proportion of
##'   the population who are able to be vaccinated.
##'
##' @param vaccine_daily_doses A single value indicating the number of
##'   (first) vaccine doses per day to distribute. The actual number
##'   distributed will be stochastically distributed around this,
##'   provided there are sufficient eligible candidates.
##'
##' @param waning_rate A single value or a vector of values representing the
##'   rates of waning of immunity after infection; if a single value the same
##'   rate is used for all age groups; if a vector of values if used it should
##'   have one value per age group.
##'
##' @param model_pcr_and_serology_user A value of 1 or 0 so switch on or off the
##'   flows out of PCR_neg and R_neg and the corresponding cap on the number of
##'   individuals leaving the R compartments
##'
##' @return A list of inputs to the model, many of which are fixed and
##'   represent data. These correspond largely to `user()` calls
##'   within the odin code, though some are also used in processing
##'   just before the model is run.
##'
##' @export
##' @examples
##' carehomes_parameters(sircovid_date("2020-02-01"), "uk")
##'
##' # example set up of vaccination parameters independent of age
##' # 3 groups: 1) unvaccinated, 2) vaccinated with partial immunity
##' # 3) fully vaccinated (but with an imperfect vaccine). People return
##' # to group 1 after a period of time in group 3 to mirror waning immunity
##'
##' # Assumption: immediately after vaccination susceptibility is reduced by
##' # 20%, and then by 50% when you reach full effect of the vaccination,
##' # then susceptibility returns to 100% upon waning of vaccine=induced
##' # immunity
##' # effect of vaccination similar across all age groups
##' rel_susceptibility <- c(1, 0.8, 0.5)
##'
##' # The vaccine also reduces the risk of symptoms
##' rel_p_sympt <- c(1, 0.6, 0.3)
##' # and the risk of hospitalisation for those with symptoms
##' rel_p_hosp_if_sympt <- c(1, 0.95, 0.95)
##'
##' # Vaccination occurs at a constant rate of 0.03 per day,
##' # (i.e. average time to vaccination is 33 days)
##' # similar across all age groups.
##' # the period of build-up of immunity following vaccination is exponentially
##' # distributed and lasts on average two weeks;
##' # vaccine-induced immunity wanes after a period which is exponentially
##' # distributed and lasts on average 26 weeks (half a year);
##' # they are similar across all age groups
##' vaccine_progression_rate <- c(0, 1/(2*7), 1/(26*7))
##'
##' # generate model parameters
##' p <- carehomes_parameters(
##'        sircovid_date("2020-02-01"), "uk",
##'        rel_susceptibility = rel_susceptibility,
##'        rel_p_sympt = rel_p_sympt,
##'        rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
##'        vaccine_progression_rate = vaccine_progression_rate,
##'        vaccine_daily_doses = 10000)
##'
##' # vaccination parameters are automatically copied across all age groups
##' p$rel_susceptibility
##' p$rel_p_sympt
##' p$rel_p_hosp_if_sympt
##' # Note that this is only the "base" rate as we fill in the first
##' # column dynamically based on vaccine_daily_doses
##' p$vaccine_progression_rate
##'
##' ### same example as above BUT assume a different effect of vaccine in the
##' ### first age group
##' n_groups <- 19
##'
##' # Assumption: vaccine is twice more effective at reducing susceptibility
##' # in the first age group
##' rel_susceptibility_agegp1 <- c(1, 0.4, 0.25)
##' rel_susceptibility_other_agegp <- c(1, 0.8, 0.5)
##' rel_susceptibility <- matrix(NA, nrow = n_groups, ncol = 3)
##' rel_susceptibility[1, ] <- rel_susceptibility_agegp1
##' for (i in seq(2, n_groups)) {
##'   rel_susceptibility[i, ] <- rel_susceptibility_other_agegp
##' }
##' rel_susceptibility
##'
##' # But vaccine has the same impact on probability of symptoms and
##' # hospitalisation for the symptomatic across all age groups
##' rel_p_sympt <- matrix(rep(rel_p_sympt, n_groups), nrow = n_groups,
##'   byrow = TRUE)
##' rel_p_hosp_if_sympt <-
##'   matrix(rep(rel_p_hosp_if_sympt, n_groups), nrow = n_groups, byrow = TRUE)
##'
##' # the period of build-up of immunity is the same for all age groups,
##' # lasting on average 2 weeks,
##' # but the first age group loses immunity more quickly
##' # (on average after 3 months) than the other age groups
##' # (on average after 6 months)
##' vaccine_progression_rate <- cbind(0,
##'                                   rep(1 / (2 * 7), n_groups),
##'                                   c(1 / (13 * 7),
##'                                   rep( 1 / (26 * 7), n_groups - 1)))
##'
##' # generate model parameters
##' p <- carehomes_parameters(
##'        sircovid_date("2020-02-01"), "uk",
##'        rel_susceptibility = rel_susceptibility,
##'        rel_p_sympt = rel_p_sympt,
##'        rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
##'        vaccine_progression_rate = vaccine_progression_rate,
##'        vaccine_daily_doses = 10000)
##'
carehomes_parameters <- function(start_date, region,
                                 beta_date = NULL, beta_value = NULL,
                                 severity = NULL,
                                 p_death_carehome = 0.7,
                                 sero_specificity = 0.9,
                                 sero_sensitivity = 0.99,
                                 progression = NULL,
                                 eps = 0.1,
                                 C_1 = 4e-6,
                                 C_2 = 5e-5,
                                 pillar2_specificity = 0.99,
                                 pillar2_sensitivity = 0.99,
                                 react_specificity = 0.99,
                                 react_sensitivity = 0.99,
                                 prop_noncovid_sympt = 0.01,
                                 strain_transmission = 1,
                                 strain_seed_date = NULL,
                                 strain_seed_value = NULL,
                                 rel_susceptibility = 1,
                                 rel_p_sympt = 1,
                                 rel_p_hosp_if_sympt = 1,
                                 vaccine_progression_rate = NULL,
                                 vaccine_uptake = NULL,
                                 vaccine_daily_doses = 0,
                                 waning_rate = 0,
                                 model_pcr_and_serology_user = 1,
                                 exp_noise = 1e6) {
  ret <- sircovid_parameters_shared(start_date, region,
                                    beta_date, beta_value)

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

  severity <- carehomes_parameters_severity(severity, p_death_carehome)

  ## TODO Rich, these parameters are now time-varying. We may want to rethink
  ## implementation of severity parameters
  ## probability of symptomatic individual requiring hospital treatment
  severity$psi_hosp_sympt <- severity$p_hosp_sympt / max(severity$p_hosp_sympt)
  severity$p_hosp_sympt_step <- max(severity$p_hosp_sympt)
  ## probability of hospitalised patient going to ICU
  severity$psi_ICU_hosp <- severity$p_ICU_hosp / max(severity$p_ICU_hosp)
  severity$p_ICU_hosp_step <- max(severity$p_ICU_hosp)
  ## probability of ICU patient dying
  severity$psi_death_ICU <- severity$p_death_ICU / max(severity$p_death_ICU)
  severity$p_death_ICU_step <- max(severity$p_death_ICU)
  ## probability of non-ICU hospital patient dying
  severity$psi_death_hosp_D <- severity$p_death_hosp_D /
    max(severity$p_death_hosp_D)
  severity$p_death_hosp_D_step <- max(severity$p_death_hosp_D)
  ## probability of stepdown hospital patient dying
  severity$psi_death_stepdown <- severity$p_death_stepdown /
    max(severity$p_death_stepdown)
  severity$p_death_stepdown_step <- max(severity$p_death_stepdown)
  ## probability of patient requiring hospital treatment dying in community
  severity$psi_death_comm <- severity$p_death_comm / max(severity$p_death_comm)
  severity$p_death_comm_step <- max(severity$p_death_comm)
  ## probability of an admission already being confirmed covid
  severity$psi_admit_conf <- severity$p_admit_conf / max(severity$p_admit_conf)
  severity$p_admit_conf_step <- max(severity$p_admit_conf)

  progression <- progression %||% carehomes_parameters_progression()

  waning <- carehomes_parameters_waning(waning_rate)

  ret$m <- carehomes_transmission_matrix(eps, C_1, C_2, region)

  ret$N_tot <- carehomes_population(ret$population, carehome_workers,
                                    carehome_residents)

  ## This is used to normalise the serology counts (converting them
  ## from number of positive/negative tests into a fraction). This is
  ## constant over the simulation, being the total population size of
  ## 15 to 64 year olds.
  ret$N_tot_15_64 <- sum(ret$N_tot[4:13])

  ## Specificity for serology tests
  ret$sero_specificity <- sero_specificity
  ret$sero_sensitivity <- sero_sensitivity

  ## Specificity and sensitivity for Pillar 2 testing
  ret$pillar2_specificity <- pillar2_specificity
  ret$pillar2_sensitivity <- pillar2_sensitivity

  ## Specificity and sensitivity for REACT testing
  ret$react_specificity <- react_specificity
  ret$react_sensitivity <- react_sensitivity

  ## Proportion of population with covid-like symptoms without covid
  ret$prop_noncovid_sympt <- prop_noncovid_sympt

  ## All observation parameters:
  ret$observation <- carehomes_parameters_observation(exp_noise)

  ret$n_groups <- ret$n_age_groups + 2L

  ## number of strains and relative transmissibility
  strain <- carehomes_parameters_strain(
    strain_transmission, strain_seed_date, strain_seed_value, ret$dt)

  ## vaccination
  vaccination <- carehomes_parameters_vaccination(ret$N_tot,
                                                  rel_susceptibility,
                                                  rel_p_sympt,
                                                  rel_p_hosp_if_sympt,
                                                  vaccine_progression_rate,
                                                  vaccine_uptake,
                                                  vaccine_daily_doses)

  model_pcr_and_serology_user <-
    list(model_pcr_and_serology_user = model_pcr_and_serology_user)

  c(ret, severity, progression, strain, vaccination, waning,
    model_pcr_and_serology_user)

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
##'   probability of a positive sero test, (9) probability of a positive
##'   pillar 2 test, and with element `state` containing the same values
##'   followed by 17 S compartments (one per age group, then one for
##'   carehome workers and carehome residents respectively) and 17
##'   "cumulative admission" compartments.
##'
##' @export
##' @examples
##' p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
##' mod <- carehomes$new(p, 0, 10)
##' carehomes_index(mod$info())
carehomes_index <- function(info) {
  index <- info$index

  ## Variables required for the particle filter to run:
  index_run <- c(icu = index[["I_ICU_tot"]],
                 general = index[["general_tot"]],
                 deaths_comm = index[["D_comm_tot"]],
                 deaths_hosp = index[["D_hosp_tot"]],
                 admitted = index[["cum_admit_conf"]],
                 new = index[["cum_new_conf"]],
                 sero_pos = index[["sero_pos"]],
                 sympt_cases = index[["cum_sympt_cases"]],
                 sympt_cases_over25 = index[["cum_sympt_cases_over25"]],
                 react_pos = index[["react_pos"]])

  ## Variables that we want to save for post-processing
  index_save <- c(hosp = index[["hosp_tot"]],
                  deaths = index[["D_tot"]],
                  infections = index[["cum_infections"]])
  suffix <- paste0("_", c(sircovid_age_bins()$start, "CHW", "CHR"))
  ## NOTE: We do use the S category for the Rt calculation in some
  ## downstream work, so this is going to require some work to get
  ## right.

  n_vacc_classes <- info$dim$S[[2]]

  ## To name our S categories following age and vaccine classes, we
  ## use two suffixes. The first vaccination class is special and so
  ## has an empty suffix so that we retain our model without
  ## vaccination if needed.
  ## S_0, S_5, ..., S_CHW, S_CHR, S_0_1, S_5_1, ..., S_CHW_1, S_CHR_1, ...
  s_type <- rep(c("", sprintf("_%s", seq_len(n_vacc_classes - 1L))),
                each = length(suffix))

  index_S <- set_names(index[["S"]],
                       paste0("S", suffix, s_type))
  index_cum_admit <- set_names(index[["cum_admit_by_age"]],
                               paste0("cum_admit", suffix))

  list(run = index_run,
       state = c(index_run, index_save, index_S, index_cum_admit))
}


##' Compare observed and modelled data from the [carehomes] model. This
##' conforms to the mcstate interface.
##'
##' @title Compare observed and modelled data for the carehomes model
##'
##' @param state State vector for the end of the current day. This is
##'   assumed to be filtered following [carehomes_index()] so contains
##'   10 rows corresponding to ICU, general beds, admissions, deaths and
##'   seroconversion compartments.
##'
##' @param prev_state State vector for the end of the previous day, as
##'   for `state`.
##'
##' @param observed Observed data. This will be a list with elements
##'   `icu` (number of ICU beds occupied), `general` (number of
##'   general beds occupied), `deaths_hosp` (cumulative deaths in
##'   hospital settings), `deaths_comm` (cumulative deaths in
##'   community settings), `deaths` (combined deaths - used by some
##'   nations - if given then `deaths_hosp` and `deaths_comm` must be
##'   `NA`), `admitted` (hospital admissions), `new` (new cases),
##'   `npos_15_64` (number of people seropositive in ages 15-64) and
##'   `ntot_15_64` (number of people tested in ages 15-64), `pillar2_pos`
##'   (number of pillar 2 positives), `pillar2_tot` (number of pillar 2
##'   tests)
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
  model_hosp <- model_icu + model_general
  model_deaths_comm <- state["deaths_comm", ] - prev_state["deaths_comm", ]
  model_deaths_hosp <- state["deaths_hosp", ] - prev_state["deaths_hosp", ]
  model_admitted <- state["admitted", ] - prev_state["admitted", ]
  model_new <- state["new", ] - prev_state["new", ]
  model_new_admitted <- model_admitted + model_new
  model_sero_pos <- state["sero_pos", ]
  model_sympt_cases <- state["sympt_cases", ] - prev_state["sympt_cases", ]
  model_sympt_cases_over25 <- state["sympt_cases_over25", ] -
    prev_state["sympt_cases_over25", ]
  model_react_pos <- state["react_pos", ]

  ## calculate test positive probabilities for the various test data streams
  ## Pillar 2
  pillar2_negs <- pars$prop_noncovid_sympt * (sum(pars$N_tot)
                                              - model_sympt_cases)
  model_pillar2_prob_pos <- test_prob_pos(model_sympt_cases,
                                          pillar2_negs,
                                          pars$pillar2_sensitivity,
                                          pars$pillar2_specificity,
                                          pars$observation$exp_noise)

  ## Pillar 2 over 25s
  pillar2_over25_negs <- pars$prop_noncovid_sympt * (sum(pars$N_tot[6:19])
                                              - model_sympt_cases_over25)
  model_pillar2_over25_prob_pos <- test_prob_pos(model_sympt_cases_over25,
                                                 pillar2_over25_negs,
                                                 pars$pillar2_sensitivity,
                                                 pars$pillar2_specificity,
                                                 pars$observation$exp_noise)

  ## REACT (Note that for REACT we exclude group 1 (0-4) and 19 (CHR))
  model_react_prob_pos <- test_prob_pos(model_react_pos,
                                        sum(pars$N_tot[2:18]) - model_react_pos,
                                        pars$react_sensitivity,
                                        pars$react_specificity,
                                        pars$observation$exp_noise)

  ## serology
  model_sero_prob_pos <- test_prob_pos(model_sero_pos,
                                       pars$N_tot_15_64 - model_sero_pos,
                                       pars$sero_sensitivity,
                                       pars$sero_specificity,
                                       pars$observation$exp_noise)

  pars <- pars$observation
  exp_noise <- pars$exp_noise

  ## Note that in ll_nbinom, the purpose of exp_noise is to allow a
  ## non-zero probability when the model value is 0 and the observed
  ## value is non-zero (i.e. there is overreporting)
  ll_icu <- ll_nbinom(observed$icu, pars$phi_ICU * model_icu,
                      pars$k_ICU, exp_noise)
  ll_general <- ll_nbinom(observed$general, pars$phi_general * model_general,
                          pars$k_general, exp_noise)
  ll_hosp <- ll_nbinom(observed$hosp, pars$phi_hosp * model_hosp,
                       pars$k_hosp, exp_noise)
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
  ll_new_admitted <- ll_nbinom(observed$new_admitted,
                               pars$phi_new_admitted * model_new_admitted,
                               pars$k_new_admitted, exp_noise)

  ll_serology <- ll_binom(observed$npos_15_64,
                          observed$ntot_15_64,
                          model_sero_prob_pos)

  ll_pillar2_tests <- ll_betabinom(observed$pillar2_pos,
                         observed$pillar2_tot,
                         model_pillar2_prob_pos,
                         pars$rho_pillar2_tests)

  ll_pillar2_cases <- ll_nbinom(observed$pillar2_cases,
                                pars$phi_pillar2_cases * model_sympt_cases,
                                pars$k_pillar2_cases, exp_noise)

  ll_pillar2_over25_tests <- ll_betabinom(observed$pillar2_over25_pos,
                                          observed$pillar2_over25_tot,
                                          model_pillar2_over25_prob_pos,
                                          pars$rho_pillar2_tests)

  ll_pillar2_over25_cases <- ll_nbinom(observed$pillar2_over25_cases,
                                       pars$phi_pillar2_cases *
                                         model_sympt_cases_over25,
                                       pars$k_pillar2_cases, exp_noise)

  ll_react <- ll_binom(observed$react_pos,
                       observed$react_tot,
                       model_react_prob_pos)

  ll_icu + ll_general + ll_hosp + ll_deaths_hosp + ll_deaths_comm + ll_deaths +
    ll_admitted + ll_new + ll_new_admitted + ll_serology + ll_pillar2_tests +
    ll_pillar2_cases + ll_pillar2_over25_tests + ll_pillar2_over25_cases +
    ll_react
}


## We store within the severity parameters information on severity for
## carehome workers and residents. The vector ends up structured as
##
##   [1..n_age_groups, workers, residents]
##
## so we have length of n_groups = n_age_groups + 2
##
carehomes_severity <- function(p) {
  index_workers <- carehomes_index_workers()
  p_workers <- mean(p[index_workers])
  p_residents <- p[length(p)]
  c(p, p_workers, p_residents)
}


carehomes_parameters_severity <- function(severity, p_death_carehome) {
  severity <- sircovid_parameters_severity(severity)
  severity <- lapply(severity, carehomes_severity)
  severity$p_death_comm[length(severity$p_death_comm)] <- p_death_carehome
  severity
}


carehomes_index_workers <- function() {
  age_bins <- sircovid_age_bins()
  which(age_bins$start >= 25 & age_bins$start < 65)
}


carehomes_transmission_matrix <- function(eps, C_1, C_2, region) {
  index_workers <- carehomes_index_workers()
  m <- sircovid_transmission_matrix(region)
  n_age_groups <- nrow(m)

  m_chw <- apply(m[seq_len(n_age_groups), index_workers], 1, mean)
  m_chr <- eps * m[n_age_groups, seq_len(n_age_groups)]

  ## Construct a block matrix:
  ##
  ##   M     m_chw m_chr
  ##   m_chw C_1   C_1
  ##   m_chr C_1   C_2

  i <- seq_len(n_age_groups)
  i_chw <- n_age_groups + 1L
  i_chr <- n_age_groups + 2L

  ret <- matrix(0.0, n_age_groups + 2, n_age_groups + 2)
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
  index_react_pos <- index[["react_pos"]][[1L]]
  index_N_tot2 <- index[["N_tot2"]][[1L]]
  index_N_tot3 <- index[["N_tot3"]][[1L]]

  index_S <- index[["S"]]
  index_S_no_vacc <- index_S[seq_len(length(pars$N_tot))]
  index_N_tot <- index[["N_tot"]]

  ## S0 is the population totals, minus the seeded infected
  ## individuals
  initial_S <- pars$N_tot
  initial_S[seed_age_band] <- initial_S[seed_age_band] - initial_I

  state[index_S_no_vacc] <- initial_S
  state[index_I] <- initial_I
  state[index_R_pre] <- initial_I
  state[index_PCR_pos] <- initial_I
  state[index_react_pos] <- initial_I
  state[index_N_tot] <- pars$N_tot
  state[index_N_tot2] <- sum(pars$N_tot)
  state[index_N_tot3] <- sum(pars$N_tot)

  list(state = state,
       step = pars$initial_step)
}


carehomes_parameters_vaccination <- function(N_tot,
                                             rel_susceptibility = 1,
                                             rel_p_sympt = 1,
                                             rel_p_hosp_if_sympt = 1,
                                             vaccine_progression_rate = NULL,
                                             vaccine_uptake = NULL,
                                             vaccine_daily_doses = 0) {
  stopifnot(length(N_tot) == carehomes_n_groups())
  calc_n_vacc_classes <- function(x) {
    if (is.matrix(x)) ncol(x) else length(x)
  }
  rel_params <- list(rel_susceptibility = rel_susceptibility,
                     rel_p_sympt = rel_p_sympt,
                     rel_p_hosp_if_sympt = rel_p_hosp_if_sympt)

  n <- vnapply(rel_params, calc_n_vacc_classes)

  if (any(n > 1) && length(unique(n[n > 1])) != 1) {
    msg1 <- paste(names(rel_params), collapse = ", ")
    msg2 <- "should have the same dimension"
    stop(paste(msg1, msg2))
  }

  ret <- Map(function(value, name) build_rel_param(value, max(n), name),
             rel_params, names(rel_params))

  ret$vaccine_progression_rate_base <- build_vaccine_progression_rate(
    vaccine_progression_rate, max(n))

  if (is.null(vaccine_uptake)) {
    vaccine_uptake <- rep(1, carehomes_n_groups())
  } else if (length(vaccine_uptake) == 1L) {
    vaccine_uptake <- rep(vaccine_uptake, carehomes_n_groups())
  } else if (length(vaccine_uptake) != carehomes_n_groups()) {
    stop(sprintf("Invalid length %d for 'vaccine_uptake', must be 1 or %d",
                 length(vaccine_uptake), carehomes_n_groups()))
  }

  ret$vaccine_population_reluctant <- (1 - vaccine_uptake) * N_tot
  ret$vaccine_daily_doses <- vaccine_daily_doses

  ret
}

carehomes_parameters_strain <- function(strain_transmission, strain_seed_date,
                                        strain_seed_value, dt) {
  if (length(strain_transmission) == 0) {
    stop("At least one value required for 'strain_transmission'")
  }
  if (length(strain_transmission) > 2) {
    stop(paste(
      "Only 1 or 2 strains valid ('strain_transmission' too long)'.",
      "See 'n_S_progress' in the odin code to fix this"))
  }
  if (any(strain_transmission < 0)) {
    stop("'strain_transmission' must have only non-negative values")
  }
  if (strain_transmission[[1]] != 1) {
    stop("'strain_transmission[1]' must be 1")
  }

  if (is.null(strain_seed_date)) {
    if (!is.null(strain_seed_value)) {
      stop(paste("As 'strain_seed_date' is NULL, expected 'strain_seed_value'",
                 "to be NULL"))
    }
    strain_seed_step <- 0
  } else {
    if (length(strain_transmission) == 1L) {
      stop("Can't use 'strain_seed_date' if only using one strain")
    }
    if (length(strain_seed_date) != 2L) {
      stop("'strain_seed_date', if given, must be exactly two elements")
    }
    assert_sircovid_date(strain_seed_date)
    assert_increasing(strain_seed_date, strict = FALSE)

    if (length(strain_seed_value) != 1L) {
      stop("'strain_seed_value' must be a scalar if 'strain_seed_date' is used")
    }
    assert_integer(strain_seed_value)

    ## The + 1 here prevents the start of the next day having the
    ## seeding value
    strain_seed_step <- numeric((strain_seed_date[[2]] + 1) / dt)
    i <- (strain_seed_date[[1]] / dt):((strain_seed_date[[2]] + 1) / dt - 1)
    strain_seed_step[i] <- spread_integer(strain_seed_value, 1 / dt)
  }

  list(n_strains = length(strain_transmission),
       strain_transmission = strain_transmission,
       strain_seed_step = strain_seed_step)
}


carehomes_parameters_waning <- function(waning_rate) {
  waning_rate <- build_waning_rate(waning_rate)
  list(
    waning_rate = waning_rate
  )
}


##' Carehomes progression parameters.  The `s_` parameters are the
##' scaling parameters for the Erlang distibution (a.k.a 'k'), while
##' the `gamma_` parameters are the gamma parameters of that
##' distribution.  These need to be aligned with Bob's severity
##' outputs, and we will come up with a better way of coordinating the
##' two.
##'
##' @title Carehomes progression parameters
##'
##' @return A list of parameter values
##'
##' @export
carehomes_parameters_progression <- function() {

  ## The s_ parameters are the scaling parameters for the Erlang
  ## distibution (a.k.a 'k'), while the gamma parameters are the gamma
  ## parameters of that distribution.
  list(s_E = 2,
       s_asympt = 1,
       s_sympt = 1,
       s_comm_D = 2,
       s_hosp_D = 2,
       s_hosp_R = 2,
       s_ICU_D = 2,
       s_ICU_S_R = 2,
       s_ICU_S_D = 2,
       s_triage = 2,
       s_stepdown_R = 2,
       s_stepdown_D = 2,
       s_R_pos = 2,
       s_PCR_pre = 2,
       s_PCR_pos = 2,

       gamma_E = 1 / (4.59 / 2),
       gamma_asympt = 1 / 2.09,
       gamma_sympt = 1 / 4,
       gamma_comm_D = 2 / 5,
       gamma_hosp_D = 2 / 5,
       gamma_hosp_R = 2 / 10,
       gamma_ICU_D = 2 / 5,
       gamma_ICU_S_R = 2 / 10,
       gamma_ICU_S_D = 2 / 10,
       gamma_triage = 2,
       gamma_stepdown_R = 2 / 5,
       gamma_stepdown_D = 2 / 5,
       gamma_R_pre_1 = 1 / 5,
       gamma_R_pre_2 = 1 / 10,
       gamma_R_pos = 1 / 25,
       gamma_test = 3 / 10,
       gamma_PCR_pre = 2 / 3,
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


carehomes_parameters_observation <- function(exp_noise) {
  list(
    ## People currently in ICU
    phi_ICU = 1,
    k_ICU = 2,
    ## People currently in general beds
    phi_general = 1,
    k_general = 2,
    ## People currently in all hospital beds
    phi_hosp = 1,
    k_hosp = 2,
    ## Daily hospital deaths
    phi_death_hosp = 1,
    k_death_hosp = 2,
    ## Daily community deaths
    phi_death_comm = 1,
    k_death_comm = 2,
    ## Daily total deaths (if not split)
    k_death = 2,
    ## Daily new confirmed admissions
    phi_admitted = 1,
    k_admitted = 2,
    ## Daily new inpatient diagnoses
    phi_new = 1,
    k_new = 2,
    ## Daily combined new confirmed admissions and new inpatient diagnoses
    phi_new_admitted = 1,
    k_new_admitted = 2,
    ## Pillar 2 testing
    phi_pillar2_cases = 1,
    k_pillar2_cases = 2,
    ##
    rho_pillar2_tests = 0.1,
    ##
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
  ## working ages, extracting equal amounts from those age groups.
  N_tot <- c(population, carehome_workers, carehome_residents)

  index_workers <- carehomes_index_workers()
  index_residents <- which(sircovid_age_bins()$start >= 65)
  weights_residents <- c(0.05, 0.05, 0.15, 0.75)

  N_tot[index_residents] <-
    round(N_tot[index_residents] - carehome_residents * weights_residents)
  N_tot[index_workers] <-
    round(N_tot[index_workers] - carehome_workers / length(index_workers))

  if (any(N_tot[index_residents] < 0)) {
    stop("Not enough population to meet care home occupancy")
  }
  if (any(N_tot[index_workers] < 0)) {
    stop("Not enough population to be care workers")
  }
  N_tot
}


##' Construct a particle filter using the [`carehomes`] model. This is
##' a convenience function, ensuring that compatible functions are
##' used together.
##'
##' @title Carehomes particle filter
##'
##' @param data Data suitable for use with with the
##'   [`mcstate::particle_filter`], created by created by
##'   [mcstate::particle_filter_data()]. We require columns "icu",
##'   "general", "deaths_hosp", "deaths_comm", "deaths", "admitted",
##'   "new", "npos_15_64", "ntot_15_64", though thse may be entirely
##'   `NA` if no data are present.
##'
##' @param n_particles Number of particles to use
##'
##' @param n_threads Number of threads to use
##'
##' @param seed Random seed to use
##'
##' @return A [`mcstate::particle_filter`] object
##' @export
carehomes_particle_filter <- function(data, n_particles,
                                      n_threads = 1L, seed = NULL) {
  mcstate::particle_filter$new(
    carehomes_particle_filter_data(data),
    carehomes,
    n_particles,
    carehomes_compare,
    carehomes_index,
    carehomes_initial,
    n_threads,
    seed)
}


carehomes_particle_filter_data <- function(data) {
  required <- c("icu", "general", "hosp", "deaths_hosp", "deaths_comm",
                "deaths", "admitted", "new", "new_admitted", "npos_15_64",
                "ntot_15_64", "pillar2_pos", "pillar2_tot", "pillar2_cases",
                "pillar2_over25_pos", "pillar2_over25_tot",
                "pillar2_over25_cases", "react_pos", "react_tot")

  verify_names(data, required, allow_extra = TRUE)

  if (any(!is.na(data$deaths) &
           (!is.na(data$deaths_comm) | !is.na(data$deaths_hosp)))) {
    stop("Deaths are not consistently split into total vs community/hospital")
  }

  pillar2_streams <- sum(c(any(!is.na(data$pillar2_pos)) |
                             any(!is.na(data$pillar2_tot)),
                           any(!is.na(data$pillar2_cases)),
                           any(!is.na(data$pillar2_over25_pos)) |
                             any(!is.na(data$pillar2_over25_tot)),
                           any(!is.na(data$pillar2_over25_cases))))
  if (pillar2_streams > 1) {
    stop("Cannot fit to more than one pillar 2 data stream")
  }

  data
}

carehomes_n_groups <- function() {
  length(sircovid_age_bins()$start) + 2L
}


##' Forecast from the carehomes model; this provides a wrapper around
##' [mcstate::pmcmc_sample] and [mcstate::pmcmc_predict] that samples
##' the trajectories then creates samples, setting the sircovid dates
##' and adding trajectories of incidence.
##'
##' @title Forecast the carehomes model
##'
##' @inheritParams mcstate::pmcmc_sample
##' @inheritParams mcstate::pmcmc_predict
##'
##' @param samples Results of running [mcstate::pmcmc()]
##'
##' @param n_sample Number of samples to take. If you provide a value
##'   of 0 then sampling is skipped (i.e., the entire set of
##'   parameters in `samples` will be forecasted). Otherwise this is
##'   passed to [mcstate::pmcmc_sample]
##'
##' @param forecast_days The number of days to create a forecast for
##'
##' @param incidence_states A character vector of states for which
##'   incidnce should be computed (from cumulative compartments, such
##'   as deaths). These will end up in the final trajectories object
##'   with the sufix `_inc` (e.g., `deaths` becomes `deaths_inc`).
##'
##' @export
carehomes_forecast <- function(samples, n_sample, burnin, forecast_days,
                               incidence_states,
                               prepend_trajectories = TRUE) {
  if (n_sample == 0) {
    ret <- samples
  } else {
    ret <- mcstate::pmcmc_sample(samples, n_sample, burnin)
  }
  steps_predict <- seq(ret$predict$step,
                       length.out = forecast_days + 1L,
                       by = ret$predict$rate)
  ret$trajectories <- mcstate::pmcmc_predict(
    ret, steps_predict,
    prepend_trajectories = prepend_trajectories)

  ret$trajectories$date <- ret$trajectories$step / ret$trajectories$rate
  ret$trajectories <- add_trajectory_incidence(
    ret$trajectories, incidence_states)
  ret
}
