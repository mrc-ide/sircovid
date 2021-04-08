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
##' @param initial_I Initial number of infected indidviduals; these
##'   will enter the model as asymptomatic 15-19 year olds at
##'   `start_date`. The default is 10 individuals.
##'
##' @param eps Change in contact rate for carehome residents
##'
##' @param m_CHW Contact rate between carehome workers and either
##'   residents or workers
##'
##' @param m_CHR Contact rate between carehome residents
##'
##' @param pillar2_specificity Specificity of the Pillar 2 test
##'
##' @param pillar2_sensitivity Sensitivity of the Pillar 2 test
##'
##' @param react_specificity Specificity of the REACT test
##'
##' @param react_sensitivity Sensitivity of the REACT test
##'
##' @param p_NC Proportion of population who do not have
##'   covid but have covid-like symptoms
##'
##' @param strain_transmission Vector of relative transmissibility of each
##'   strain modelled. First element should be 1. Length will define the
##'   number of strains used in the model
##'
##' @param strain_seed_date Either `NULL` (no seeding) or a vector of
##'   [sircovid::sircovid_date] values corresponding to
##'   the dates that `strain_seed_rate` should change. For example, to
##'   seed a constant rate from a given date, provide one value for each of
##'   `strain_seed_date` and `strain_seed_rate`; or to seed with a constant rate
##'   for a set period, provide two dates and two
##'   rates with the second rate equal to 0.
##'
##' @param strain_seed_rate Either `NULL` (no seeding) or a vector of
##'   values representing the *daily* rate of seeding, starting on the dates set
##'   in `strain_seed_date`.
##'   Seeding is drawn from Poisson(strain_seed_rate * dt)
##'   at each **day** and so the rate is spread evenly across the steps that
##'   occur within each date. For example, to
##'   seed a constant rate from a given date, provide one value for each of
##'   `strain_seed_date` and `strain_seed_rate`; or to seed with a constant rate
##'   for a set period, provide two dates and two
##'   rates with the second rate equal to 0.
##'
##' @param strain_rel_gamma_A Vector of relative rates of progression out of
##' I_A (gamma_A) for each
##'   strain modelled. If `1` all strains have same rates. Otherwise vector of
##'   same length as `strain_transmission`, with entries that determines the
##'   relative scaling of the defaults for each strain.
##'
##' @param strain_rel_gamma_P Vector of relative rates of progression out of
##' I_P (gamma_P) for each
##'   strain modelled. If `1` all strains have same rates. Otherwise vector of
##'   same length as `strain_transmission`, with entries that determines the
##'   relative scaling of the defaults for each strain.
##'
##' @param strain_rel_gamma_C_1 Vector of relative rates of progression out of
##' I_C_1 (gamma_C_1) for each
##'   strain modelled. If `1` all strains have same rates. Otherwise vector of
##'   same length as `strain_transmission`, with entries that determines the
##'   relative scaling of the defaults for each strain.
##'
##' @param strain_rel_gamma_C_2 Vector of relative rates of progression out of
##' I_C_2 (gamma_C_2) for each
##'   strain modelled. If `1` all strains have same rates. Otherwise vector of
##'   same length as `strain_transmission`, with entries that determines the
##'   relative scaling of the defaults for each strain.
##'
##' @param strain_rel_severity Vector of relative probabilities of death for
##'   each strain modelled. If `1` all strains have same
##'   probabilities of death. Otherwise vector of same length as
##'   `strain_transmission`, where the first value should be 1 (for the first
##'   strain) and subsequent values between 0 and 1. To ensure valid
##'   probabilities, severity is upper-truncated at 1 after scaling.
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
##' @param rel_infectivity A vector or matrix of values representing the
##'   relative infectivity of individuals in different vaccination groups,
##'   if they are infected.
##'   If a vector, the first value should be 1 (for the non-vaccinated group)
##'   and subsequent values be between 0 and 1. In that case relative
##'   infectivity will be the same across all age groups within one
##'   vaccination category. Specifying a matrix instead of a vector allows
##'   different relative infectivities by age (rows of the matrix) and
##'   vaccination group (columns of the matrix); in that case, in each row of
##'   the matrix, the first value should be 1 (for the non-vaccinated group)
##'   and subsequent values be between 0 and 1
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
##' @param vaccine_schedule A [vaccine_schedule] object indicating the
##'   people to be vaccinated by group over time
##'
##' @param vaccine_index_dose2 The index to use for the second dose
##'
##' @param vaccine_catchup_fraction A value between 0 and 1 indicating the
##'   proportion of doses not distributed according to schedule (e.g. because
##'   too many people were in the I or H compartments and could not be
##'   vaccinated at the scheduled time) that we postpone to a later date.
##'   A value of 0 means we do not catch up at all on any missed doses; a
##'   value of 1 means we try to catch up for all missed doses. This is set
##'   to 1 by default
##'
##' @param waning_rate A single value or a vector of values representing the
##'   rates of waning of immunity after infection; if a single value the same
##'   rate is used for all age groups; if a vector of values if used it should
##'   have one value per age group.
##'
##' @return A list of inputs to the model, many of which are fixed and
##'   represent data. These correspond largely to `user()` calls
##'   within the odin code, though some are also used in processing
##'   just before the model is run.
##'
##' @export
##' @examples
##'
##' region <- "london"
##' carehomes_parameters(sircovid_date("2020-02-01"), region)
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
##' # The vaccine also reduces infectivity of infected individuals by half
##' rel_infectivity <- c(1, 0.5, 0.5)
##'
##' # On average 10000 first doses of vaccine are given every day
##' # Second doses are given on average 12 weeks after the first dose
##' # vaccine-induced immunity wanes after a period which is exponentially
##' # distributed and lasts on average 26 weeks (half a year);
##' # they are similar across all age groups
##' vaccine_progression_rate <- c(0, 0, 1/(26*7))
##'
##' daily_doses <- rep(10000, 365)
##' mean_days_between_doses <- 12 * 7
##' n <- vaccine_priority_population(region, uptake = 1)
##' schedule <- vaccine_schedule_future(
##'   0, daily_doses, mean_days_between_doses, n)
##'
##' # generate model parameters
##' p <- carehomes_parameters(
##'        sircovid_date("2020-02-01"), region,
##'        rel_susceptibility = rel_susceptibility,
##'        rel_p_sympt = rel_p_sympt,
##'        rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
##'        rel_infectivity = rel_infectivity,
##'        vaccine_progression_rate = vaccine_progression_rate,
##'        vaccine_schedule = schedule,
##'        vaccine_index_dose2 = 2)
##'
##' # vaccination parameters are automatically copied across all age groups
##' p$rel_susceptibility
##' p$rel_p_sympt
##' p$rel_p_hosp_if_sympt
##' p$rel_infectivity
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
##' # And vaccine has the same impact on onwards infectivity across age groups
##' rel_infectivity <- matrix(rep(rel_infectivity, n_groups), nrow = n_groups,
##'   byrow = TRUE)
##'
##' # the period of build-up of immunity is the same for all age groups,
##' # lasting on average 2 weeks,
##' # but the first age group loses immunity more quickly
##' # (on average after 3 months) than the other age groups
##' # (on average after 6 months)
##' vaccine_progression_rate <- cbind(0, 0,
##'                                   c(1 / (13 * 7),
##'                                   rep( 1 / (26 * 7), n_groups - 1)))
##'
##' # generate model parameters
##' p <- carehomes_parameters(
##'        sircovid_date("2020-02-01"), region,
##'        rel_susceptibility = rel_susceptibility,
##'        rel_p_sympt = rel_p_sympt,
##'        rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
##'        rel_infectivity = rel_infectivity,
##'        vaccine_progression_rate = vaccine_progression_rate,
##'        vaccine_schedule = schedule,
##'        vaccine_index_dose2 = 2)
##'
##' # TODO: add an example of manually set up vaccine schedule
##'
carehomes_parameters <- function(start_date, region,
                                 beta_date = NULL, beta_value = NULL,
                                 severity = NULL,
                                 p_death_carehome = 0.7,
                                 sero_specificity = 0.9,
                                 sero_sensitivity = 0.99,
                                 progression = NULL,
                                 initial_I = 10,
                                 eps = 0.1,
                                 m_CHW = 4e-6,
                                 m_CHR = 5e-5,
                                 pillar2_specificity = 0.99,
                                 pillar2_sensitivity = 0.99,
                                 react_specificity = 0.99,
                                 react_sensitivity = 0.99,
                                 p_NC = 0.01,
                                 strain_transmission = 1,
                                 strain_seed_date = NULL,
                                 strain_seed_rate = NULL,
                                 strain_rel_gamma_A = 1,
                                 strain_rel_gamma_P = 1,
                                 strain_rel_gamma_C_1 = 1,
                                 strain_rel_gamma_C_2 = 1,
                                 strain_rel_severity = 1,
                                 rel_susceptibility = 1,
                                 rel_p_sympt = 1,
                                 rel_p_hosp_if_sympt = 1,
                                 rel_infectivity = 1,
                                 vaccine_progression_rate = NULL,
                                 vaccine_schedule = NULL,
                                 vaccine_index_dose2 = NULL,
                                 vaccine_catchup_fraction = 1,
                                 waning_rate = 0,
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
  strain_rel_severity <- mcstate:::recycle(
                                 assert_relatives(strain_rel_severity),
                                 length(strain_transmission))
  severity <- scale_severity(severity, strain_rel_severity)

  ## TODO Rich, these parameters are now time-varying. We may want to rethink
  ## implementation of severity parameters
  ## probability of symptomatic individual requiring hospital treatment

  get_psi <- function(p) {
    if (all(p == 0)) {
      res <- p
    } else {
      res <- p / max(p)
    }
  }

  severity$psi_H <- get_psi(severity$p_H)
  severity$p_H_step <- max(severity$p_H)
  ## probability of hospitalised patient going to ICU
  severity$psi_ICU <- get_psi(severity$p_ICU)
  severity$p_ICU_step <- max(severity$p_ICU)
  ## probability of ICU patient dying
  severity$psi_ICU_D <- apply(severity$p_ICU_D, 2, get_psi)
  severity$p_ICU_D_step <- matrix(apply(severity$p_ICU_D, 2, max), nrow = 1)
  ## probability of non-ICU hospital patient dying
  severity$psi_H_D <- apply(severity$p_H_D, 2, get_psi)
  severity$p_H_D_step <- matrix(apply(severity$p_H_D, 2, max), nrow = 1)
  ## probability of stepdown hospital patient dying
  severity$psi_W_D <- apply(severity$p_W_D, 2, get_psi)
  severity$p_W_D_step <- matrix(apply(severity$p_W_D, 2, max), nrow = 1)
  ## probability of patient requiring hospital treatment dying in community
  severity$psi_G_D <- apply(severity$p_G_D, 2, get_psi)
  severity$p_G_D_step <- matrix(apply(severity$p_G_D, 2, max), nrow = 1)
  ## probability of an admission already being confirmed covid
  severity$psi_star <- get_psi(severity$p_star)
  severity$p_star_step <- max(severity$p_star)

  strain_rel_gamma_A <- mcstate:::recycle(assert_relatives(strain_rel_gamma_A),
                                          length(strain_transmission))
  strain_rel_gamma_P <- mcstate:::recycle(assert_relatives(strain_rel_gamma_P),
                                          length(strain_transmission))
  strain_rel_gamma_C_1 <-
    mcstate:::recycle(assert_relatives(strain_rel_gamma_C_1),
                                       length(strain_transmission))
  strain_rel_gamma_C_2 <-
    mcstate:::recycle(assert_relatives(strain_rel_gamma_C_2),
                                       length(strain_transmission))

  progression <- progression %||%
                  carehomes_parameters_progression(strain_rel_gamma_A,
                                                   strain_rel_gamma_P,
                                                   strain_rel_gamma_C_1,
                                                   strain_rel_gamma_C_2)

  ## implementation of time-varying progression gammas
  progression$gamma_H_R_step <- progression$gamma_H_R
  progression$gamma_W_R_step <- progression$gamma_W_R
  progression$gamma_ICU_W_R_step <- progression$gamma_ICU_W_R
  progression$gamma_H_D_step <- progression$gamma_H_D
  progression$gamma_W_D_step <- progression$gamma_W_D
  progression$gamma_ICU_W_D_step <- progression$gamma_ICU_W_D
  progression$gamma_ICU_D_step <- progression$gamma_ICU_D
  progression$gamma_ICU_pre_step <- progression$gamma_ICU_pre

  waning <- carehomes_parameters_waning(waning_rate)

  ret$m <- carehomes_transmission_matrix(eps, m_CHW, m_CHR, region)

  ret$N_tot <- carehomes_population(ret$population, carehome_workers,
                                    carehome_residents)

  ret$initial_I <- initial_I

  ## This is used to normalise the serology counts (converting them
  ## from number of positive/negative tests into a fraction). This is
  ## constant over the simulation, being the total population size of
  ## 15 to 64 year olds. Similarly, we need the total population over
  ## all ranges, the population corresponding to REACT (all but CHR
  ## and <5) and the over 25
  ret$N_tot_15_64 <- sum(ret$N_tot[4:13])
  ret$N_tot_all <- sum(ret$N_tot)
  ret$N_tot_over25 <- sum(ret$N_tot[6:19])
  ret$N_tot_react <- sum(ret$N_tot[2:18])

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
  ret$p_NC <- p_NC

  ## relative transmissibility of various I compartments
  ret$I_A_transmission <- 0.363
  ret$I_P_transmission <- 1
  ret$I_C_1_transmission <- 1
  ret$I_C_2_transmission <- 0

  ## All observation parameters:
  observation <- carehomes_parameters_observation(exp_noise)

  ret$n_groups <- ret$n_age_groups + 2L

  ## number of strains and relative transmissibility
  strain <- carehomes_parameters_strain(
    strain_transmission, strain_seed_date, strain_seed_rate, ret$dt)

  ## vaccination
  vaccination <- carehomes_parameters_vaccination(ret$N_tot,
                                                  ret$dt,
                                                  rel_susceptibility,
                                                  rel_p_sympt,
                                                  rel_p_hosp_if_sympt,
                                                  rel_infectivity,
                                                  vaccine_progression_rate,
                                                  vaccine_schedule,
                                                  vaccine_index_dose2,
                                                  vaccine_catchup_fraction)

  c(ret, severity, progression, strain, vaccination, waning, observation)
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
  index_core <- c(icu = index[["ICU_tot"]],
                 general = index[["general_tot"]],
                 deaths_comm = index[["D_comm_tot"]],
                 deaths_carehomes = index[["D_carehomes_tot"]],
                 deaths_hosp = index[["D_hosp_tot"]],
                 admitted = index[["cum_admit_conf"]],
                 diagnoses = index[["cum_new_conf"]],
                 deaths_carehomes_inc = index[["D_carehomes_inc"]],
                 deaths_comm_inc = index[["D_comm_inc"]],
                 deaths_hosp_inc = index[["D_hosp_inc"]],
                 admitted_inc = index[["admit_conf_inc"]],
                 diagnoses_inc = index[["new_conf_inc"]],
                 sero_pos = index[["sero_pos"]],
                 sympt_cases = index[["cum_sympt_cases"]],
                 sympt_cases_over25 = index[["cum_sympt_cases_over25"]],
                 sympt_cases_non_variant_over25 =
                 index[["cum_sympt_cases_non_variant_over25"]],
                 sympt_cases_inc = index[["sympt_cases_inc"]],
                 sympt_cases_over25_inc = index[["sympt_cases_over25_inc"]],
                 sympt_cases_non_variant_over25_inc =
                   index[["sympt_cases_non_variant_over25_inc"]],
                 react_pos = index[["react_pos"]])

  ## Only incidence versions for the likelihood now:
  index_run <- index_core[c("icu", "general", "deaths_carehomes_inc",
                            "deaths_comm_inc", "deaths_hosp_inc",
                            "admitted_inc", "diagnoses_inc",
                            "sero_pos", "sympt_cases_inc",
                            "sympt_cases_over25_inc",
                            "sympt_cases_non_variant_over25_inc",
                            "react_pos")]
  ## But cumulative versions for everything else:
  index_state_core <- index_core[sub("_inc$", "", names(index_run))]

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
  index_I_weighted <- set_names(index[["I_weighted"]],
                                paste0("I_weighted", suffix, s_type))
  index_cum_admit <- set_names(index[["cum_admit_by_age"]],
                               paste0("cum_admit", suffix))
  index_D_hosp <- set_names(index[["D_hosp"]],
                            paste0("D_hosp", suffix))
  index_cum_n_vaccinated <- set_names(index[["cum_n_vaccinated"]],
                                    paste0("cum_n_vaccinated", suffix, s_type))

  ## prob_strain is named similarly to S, with the second suffix representing
  ## strain instead of vacc_class

  n_strains <- info$dim$prob_strain[[2]]
  strain_type <- rep(c("", sprintf("_%s", seq_len(n_strains - 1L))),
                     each = length(suffix))
  index_prob_strain <- set_names(index[["prob_strain"]],
                                 paste0("prob_strain", suffix, strain_type))

  list(run = index_run,
       state = c(index_state_core, index_save, index_S, index_cum_admit,
                 index_D_hosp, index_I_weighted, index_prob_strain,
                 index_cum_n_vaccinated))
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
##' @param observed Observed data. At the moment please see the tests
##'   for a full list as this changes frequently (and this function
##'   may be removed in future).
##'
##' @param pars A list of parameters, as created by
##'   [carehomes_parameters()]
##'
##' @return A vector of log likelihoods, the same length as the number
##'   of particles (the number of columns in the modelled state)
##'
##' @export
##' @importFrom stats dbinom
carehomes_compare <- function(state, observed, pars) {
  ## TODO: we might refactor this to produce a subset of comparisons
  ## (and indices above) to suit either the SPI-M or paper fits as
  ## we're using different streams; that will make the comparisons a
  ## touch faster and data copying smaller. This requires a bit of a
  ## rethink though as it affects how we do index functions too.

  model_icu <- state["icu", ]
  model_general <- state["general", ]
  model_hosp <- model_icu + model_general
  model_deaths_carehomes <- state["deaths_carehomes_inc", ]
  model_deaths_comm <- state["deaths_comm_inc", ]
  model_deaths_hosp <- state["deaths_hosp_inc", ]
  model_admitted <- state["admitted_inc", ]
  model_diagnoses <- state["diagnoses_inc", ]
  model_all_admission <- model_admitted + model_diagnoses
  model_sero_pos <- state["sero_pos", ]
  model_sympt_cases <- state["sympt_cases_inc", ]
  model_sympt_cases_over25 <- state["sympt_cases_over25_inc", ]
  model_sympt_cases_non_variant_over25 <-
    state["sympt_cases_non_variant_over25_inc", ]
  model_react_pos <- state["react_pos", ]

  ## calculate test positive probabilities for the various test data streams
  ## Pillar 2
  pillar2_negs <- pars$p_NC * (pars$N_tot_all - model_sympt_cases)
  model_pillar2_prob_pos <- test_prob_pos(model_sympt_cases,
                                          pillar2_negs,
                                          pars$pillar2_sensitivity,
                                          pars$pillar2_specificity,
                                          pars$exp_noise)

  ## Pillar 2 over 25s
  pillar2_over25_negs <- pars$p_NC * (pars$N_tot_over25 -
                                      model_sympt_cases_over25)
  model_pillar2_over25_prob_pos <- test_prob_pos(model_sympt_cases_over25,
                                                 pillar2_over25_negs,
                                                 pars$pillar2_sensitivity,
                                                 pars$pillar2_specificity,
                                                 pars$exp_noise)

  ## REACT (Note that for REACT we exclude group 1 (0-4) and 19 (CHR))
  model_react_prob_pos <- test_prob_pos(model_react_pos,
                                        pars$N_tot_react - model_react_pos,
                                        pars$react_sensitivity,
                                        pars$react_specificity,
                                        pars$exp_noise)

  ## serology
  model_sero_prob_pos <- test_prob_pos(model_sero_pos,
                                       pars$N_tot_15_64 - model_sero_pos,
                                       pars$sero_sensitivity,
                                       pars$sero_specificity,
                                       pars$exp_noise)

  ## Strain
  model_strain_over25_prob_pos <- test_prob_pos(
    model_sympt_cases_non_variant_over25,
    model_sympt_cases_over25 - model_sympt_cases_non_variant_over25,
    1, 1, pars$exp_noise)

  exp_noise <- pars$exp_noise

  ## Note that in ll_nbinom, the purpose of exp_noise is to allow a
  ## non-zero probability when the model value is 0 and the observed
  ## value is non-zero (i.e. there is overreporting)
  ll_icu <- ll_nbinom(observed$icu, pars$phi_ICU * model_icu,
                      pars$kappa_ICU, exp_noise)
  ll_general <- ll_nbinom(observed$general, pars$phi_general * model_general,
                          pars$kappa_general, exp_noise)
  ll_hosp <- ll_nbinom(observed$hosp, pars$phi_hosp * model_hosp,
                       pars$kappa_hosp, exp_noise)
  ll_deaths_hosp <- ll_nbinom(observed$deaths_hosp,
                              pars$phi_death_hosp * model_deaths_hosp,
                              pars$kappa_death_hosp, exp_noise)
  ll_deaths_carehomes <- ll_nbinom(observed$deaths_carehomes,
                                   pars$phi_death_carehomes *
                                     model_deaths_carehomes,
                                   pars$kappa_death_carehomes, exp_noise)
  ll_deaths_comm <- ll_nbinom(observed$deaths_comm,
                              pars$phi_death_comm * model_deaths_comm,
                              pars$kappa_death_comm, exp_noise)
  ll_deaths_non_hosp <- ll_nbinom(observed$deaths_non_hosp,
                                  pars$phi_death_comm * model_deaths_comm +
                                    pars$phi_death_carehomes *
                                    model_deaths_carehomes,
                                  pars$kappa_death_non_hosp, exp_noise)
  ll_deaths <- ll_nbinom(observed$deaths,
                         pars$phi_death_hosp * model_deaths_hosp +
                           pars$phi_death_carehomes * model_deaths_carehomes +
                           pars$phi_death_comm * model_deaths_comm,
                         pars$kappa_death, exp_noise)
  ll_admitted <- ll_nbinom(observed$admitted,
                           pars$phi_admitted * model_admitted,
                           pars$kappa_admitted, exp_noise)
  ll_diagnoses <- ll_nbinom(observed$diagnoses,
                            pars$phi_diagnoses * model_diagnoses,
                            pars$kappa_diagnoses, exp_noise)
  ll_all_admission <- ll_nbinom(observed$all_admission,
                               pars$phi_all_admission * model_all_admission,
                               pars$kappa_all_admission, exp_noise)

  ll_serology <- ll_binom(observed$npos_15_64,
                          observed$ntot_15_64,
                          model_sero_prob_pos)

  ll_pillar2_tests <- ll_betabinom(observed$pillar2_pos,
                         observed$pillar2_tot,
                         model_pillar2_prob_pos,
                         pars$rho_pillar2_tests)

  ll_pillar2_cases <- ll_nbinom(observed$pillar2_cases,
                                pars$phi_pillar2_cases * model_sympt_cases,
                                pars$kappa_pillar2_cases, exp_noise)

  ll_pillar2_over25_tests <- ll_betabinom(observed$pillar2_over25_pos,
                                          observed$pillar2_over25_tot,
                                          model_pillar2_over25_prob_pos,
                                          pars$rho_pillar2_tests)

  ll_pillar2_over25_cases <- ll_nbinom(observed$pillar2_over25_cases,
                                       pars$phi_pillar2_cases *
                                         model_sympt_cases_over25,
                                       pars$kappa_pillar2_cases, exp_noise)

  ll_react <- ll_binom(observed$react_pos,
                       observed$react_tot,
                       model_react_prob_pos)

  ll_strain_over25 <- ll_binom(observed$strain_non_variant,
                               observed$strain_tot,
                               model_strain_over25_prob_pos)

  ll_icu + ll_general + ll_hosp + ll_deaths_hosp + ll_deaths_carehomes +
    ll_deaths_comm + ll_deaths_non_hosp + ll_deaths + ll_admitted +
    ll_diagnoses + ll_all_admission + ll_serology + ll_pillar2_tests +
    ll_pillar2_cases + ll_pillar2_over25_tests + ll_pillar2_over25_cases +
    ll_react + ll_strain_over25
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
  severity$p_G_D[length(severity$p_G_D)] <- p_death_carehome
  severity
}


carehomes_index_workers <- function() {
  age_bins <- sircovid_age_bins()
  which(age_bins$start >= 25 & age_bins$start < 65)
}


carehomes_transmission_matrix <- function(eps, m_CHW, m_CHR, region) {
  index_workers <- carehomes_index_workers()
  m <- sircovid_transmission_matrix(region)
  n_age_groups <- nrow(m)

  m_gen_chw <- apply(m[seq_len(n_age_groups), index_workers], 1, mean)
  m_gen_chr <- eps * m[n_age_groups, seq_len(n_age_groups)]

  ## Construct a block matrix:
  ##
  ##   M          m_gen_chw  m_gen_chr
  ##   m_gen_chw  m_CHW      m_CHW
  ##   m_gen_chr  m_CHW      m_CHR

  i <- seq_len(n_age_groups)
  i_chw <- n_age_groups + 1L
  i_chr <- n_age_groups + 2L

  ret <- matrix(0.0, n_age_groups + 2, n_age_groups + 2)
  ret[i, i] <- m
  ret[i, i_chw] <- ret[i_chw, i] <- m_gen_chw
  ret[i, i_chr] <- ret[i_chr, i] <- m_gen_chr
  ret[i_chw:i_chr, i_chw:i_chr] <- c(m_CHW, m_CHW, m_CHW, m_CHR)

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

  ## Default will be to start with 10 individuals, but this is tuneable
  initial_I <- pars$initial_I

  ## This corresponds to the 15-19y age bracket for compatibility with
  ## our first version, will be replaced by better seeding model, but
  ## probably has limited impact.
  seed_age_band <- 4L
  index_I <- index[["I_A"]][[1L]] + seed_age_band - 1L
  index_I_weighted <- index[["I_weighted"]][[1L]] + seed_age_band - 1L
  index_T_sero_pre <- index[["T_sero_pre"]][[1L]] + seed_age_band - 1L
  index_T_PCR_pos <- index[["T_PCR_pos"]][[1L]] + seed_age_band - 1L
  index_react_pos <- index[["react_pos"]][[1L]]
  index_N_tot2 <- index[["N_tot2"]][[1L]]
  index_N_tot3 <- index[["N_tot3"]][[1L]]

  index_S <- index[["S"]]
  index_S_no_vacc <- index_S[seq_len(length(pars$N_tot))]
  index_N_tot <- index[["N_tot"]]

  index_prob_strain <- index[["prob_strain"]]
  index_prob_strain_original <- index_prob_strain[seq_len(pars$n_groups)]

  ## S0 is the population totals, minus the seeded infected
  ## individuals
  initial_S <- pars$N_tot
  initial_S[seed_age_band] <- initial_S[seed_age_band] - initial_I

  state[index_S_no_vacc] <- initial_S
  state[index_I] <- initial_I
  state[index_I_weighted] <- pars$I_A_transmission * initial_I
  state[index_T_sero_pre] <- initial_I
  state[index_T_PCR_pos] <- initial_I
  state[index_react_pos] <- initial_I
  state[index_N_tot] <- pars$N_tot
  state[index_N_tot2] <- sum(pars$N_tot)
  state[index_N_tot3] <- sum(pars$N_tot)
  state[index_prob_strain_original] <- 1

  list(state = state,
       step = pars$initial_step)
}


carehomes_parameters_vaccination <- function(N_tot,
                                             dt,
                                             rel_susceptibility = 1,
                                             rel_p_sympt = 1,
                                             rel_p_hosp_if_sympt = 1,
                                             rel_infectivity = 1,
                                             vaccine_progression_rate = NULL,
                                             vaccine_schedule = NULL,
                                             vaccine_index_dose2 = NULL,
                                             vaccine_catchup_fraction = 1) {
  n_groups <- carehomes_n_groups()
  stopifnot(length(N_tot) == n_groups)
  calc_n_vacc_classes <- function(x) {
    if (is.matrix(x)) ncol(x) else length(x)
  }
  rel_params <- list(rel_susceptibility = rel_susceptibility,
                     rel_p_sympt = rel_p_sympt,
                     rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
                     rel_infectivity = rel_infectivity)

  n <- vnapply(rel_params, calc_n_vacc_classes)

  if (any(n > 1) && length(unique(n[n > 1])) != 1) {
    msg1 <- paste(names(rel_params), collapse = ", ")
    msg2 <- "should have the same dimension"
    stop(paste(msg1, msg2))
  }
  n_vacc_classes <- max(n)

  ret <- Map(function(value, name) build_rel_param(value, n_vacc_classes, name),
             rel_params, names(rel_params))

  n_doses <- 2L # fixed; see the odin code

  if (is.null(vaccine_schedule)) {
    if (!is.null(vaccine_index_dose2) && vaccine_index_dose2 != 1L) {
      stop("'vaccine_index_dose2' set without schedule")
    }
    ret$vaccine_dose_step <- array(0, c(n_groups, n_doses, 1))
    ret$index_dose <- c(1L, 1L)
  } else {
    assert_is(vaccine_schedule, "vaccine_schedule")
    vaccine_index_dose2 <- vaccine_index_dose2 %||% 1L
    if (vaccine_index_dose2 > n_vacc_classes) {
      stop(sprintf(
        "Invalid value for 'vaccine_index_dose2', must be in [1, %d]",
        n_vacc_classes))
    }

    n_days <- dim(vaccine_schedule$doses)[[3]]
    i <- rep(seq_len(n_days), each = 1 / dt)
    len <- vaccine_schedule$date / dt
    ret$index_dose <- c(1L, vaccine_index_dose2)

    ret$vaccine_dose_step <- mcstate::array_bind(
      array(0, c(n_groups, n_doses, len)),
      (vaccine_schedule$doses * dt)[, , i])
  }

  ret$n_vacc_classes <- n_vacc_classes
  ret$vaccine_progression_rate_base <- build_vaccine_progression_rate(
    vaccine_progression_rate, n_vacc_classes, ret$index_dose)


  assert_scalar(vaccine_catchup_fraction)
  assert_proportion(vaccine_catchup_fraction)
  ret$vaccine_catchup_fraction <- vaccine_catchup_fraction

  ret
}

carehomes_parameters_strain <- function(strain_transmission, strain_seed_date,
                                        strain_seed_rate, dt) {
  if (length(strain_transmission) == 0) {
    stop("At least one value required for 'strain_transmission'")
  }
  if (length(strain_transmission) > 2) {
    stop(paste(
      "Only 1 or 2 strains valid ('strain_transmission' too long)'.",
      "See 'n_S_progress' in the odin code to fix this"))
  }

  assert_relatives(strain_transmission)

  if (is.null(strain_seed_date)) {
    if (!is.null(strain_seed_rate)) {
      stop(paste("As 'strain_seed_date' is NULL, expected 'strain_seed_rate'",
                 "to be NULL"))
    }
    strain_seed_step <- 0
  } else {
    if (length(strain_transmission) == 1L) {
      stop("Can't use 'strain_seed_date' if only using one strain")
    }
    if (length(strain_seed_date) != length(strain_seed_rate)) {
      stop("'strain_seed_date' and 'strain_seed_rate' must be the same length")
    }
    assert_sircovid_date(strain_seed_date)
    assert_increasing(strain_seed_date, strict = FALSE)
    assert_non_negative(strain_seed_rate)

    ## The + 1 here prevents the start of the next day having the
    ## same seeding value
    strain_seed_step <-
      numeric(strain_seed_date[[length(strain_seed_date)]] / dt)
    for (j in seq_along(strain_seed_date)) {
      if (j == length(strain_seed_date)) {
        i <- length(strain_seed_step)
      } else {
        i <- seq.int(
          strain_seed_date[[j]] / dt,
          strain_seed_date[[j + 1]] / dt - 1
        )
      }
      strain_seed_step[i] <- strain_seed_rate[[j]] * dt
    }

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
##' @param rel_gamma_A Vector of relative rates of gamma_A for each
##'   strain modelled. If `1` all strains have same rates. Otherwise vector of
##'   same length as `strain_transmission`, with entries that determines the
##'   relative scaling of the defaults for each strain.
##'
##' @param rel_gamma_P Vector of relative rates of gamma_P for each
##'   strain modelled. If `1` all strains have same rates. Otherwise vector of
##'   same length as `strain_transmission`, with entries that determines the
##'   relative scaling of the defaults for each strain.
##'
##' @param rel_gamma_C_1 Vector of relative rates of gamma_C_1 for each
##'   strain modelled. If `1` all strains have same rates. Otherwise vector of
##'   same length as `strain_transmission`, with entries that determines the
##'   relative scaling of the defaults for each strain.
##'
##' @param rel_gamma_C_2 Vector of relative rates of gamma_C_2 for each
##'   strain modelled. If `1` all strains have same rates. Otherwise vector of
##'   same length as `strain_transmission`, with entries that determines the
##'   relative scaling of the defaults for each strain.
##'
##' @export
carehomes_parameters_progression <- function(rel_gamma_A = 1,
                                             rel_gamma_P = 1,
                                             rel_gamma_C_1 = 1,
                                             rel_gamma_C_2 = 1) {

  ## The k_ parameters are the shape parameters for the Erlang
  ## distribution, while the gamma parameters are the rate
  ## parameters of that distribution.
  list(k_E = 2,
       k_A = 1,
       k_P = 1,
       k_C_1 = 1,
       k_C_2 = 1,
       k_G_D = 2,
       k_H_D = 2,
       k_H_R = 2,
       k_ICU_D = 2,
       k_ICU_W_R = 2,
       k_ICU_W_D = 2,
       k_ICU_pre = 2,
       k_W_R = 2,
       k_W_D = 2,
       k_sero_pos = 2,
       k_PCR_pre = 2,
       k_PCR_pos = 2,

       gamma_E = 1 / (3.42 / 2),
       gamma_A = 1 / 2.88 * rel_gamma_A,
       gamma_P = 1 / 1.68 * rel_gamma_P,
       gamma_C_1 = 1 / 2.14 * rel_gamma_C_1,
       gamma_C_2 = 1 / 1.86 * rel_gamma_C_2,
       gamma_G_D = 1 / (3 / 2),
       gamma_H_D = 2 / 5,
       gamma_H_R = 2 / 10,
       gamma_ICU_D = 2 / 5,
       gamma_ICU_W_R = 2 / 10,
       gamma_ICU_W_D = 2 / 10,
       gamma_ICU_pre = 2,
       gamma_W_R = 2 / 5,
       gamma_W_D = 2 / 5,
       gamma_sero_pre_1 = 1 / 5,
       gamma_sero_pre_2 = 1 / 10,
       gamma_sero_pos = 1 / 25,
       gamma_U = 3 / 10,
       gamma_PCR_pre = 2 / 3,
       gamma_PCR_pos = 1 / 5
       )
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
    kappa_ICU = 2,
    ## People currently in general beds
    phi_general = 1,
    kappa_general = 2,
    ## People currently in all hospital beds
    phi_hosp = 1,
    kappa_hosp = 2,
    ## Daily hospital deaths
    phi_death_hosp = 1,
    kappa_death_hosp = 2,
    ## Daily care home deaths
    phi_death_carehomes = 1,
    kappa_death_carehomes = 2,
    ## Daily community deaths
    phi_death_comm = 1,
    kappa_death_comm = 2,
    ## Daily total non-hospital deaths (if not split)
    kappa_death_non_hosp = 2,
    ## Daily total deaths (if not split)
    kappa_death = 2,
    ## Daily new confirmed admissions
    phi_admitted = 1,
    kappa_admitted = 2,
    ## Daily new inpatient diagnoses
    phi_diagnoses = 1,
    kappa_diagnoses = 2,
    ## Daily combined new confirmed admissions and new inpatient diagnoses
    phi_all_admission = 1,
    kappa_all_admission = 2,
    ## Pillar 2 testing
    phi_pillar2_cases = 1,
    kappa_pillar2_cases = 2,
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
##'   "general", "hosp", "deaths_hosp", "deaths_carehomes",
##'   "deaths_non_hosp", "deaths_comm", "deaths",
##'   "admitted", "diagnoses", "all_admission", "npos_15_64",
##'   "ntot_15_64", "pillar2_pos", "pillar2_tot", "pillar2_cases",
##'   "pillar2_over25_pos", "pillar2_over25_tot", "pillar2_over25_cases",
##'   "react_pos", "react_tot", though thse may be entirely `NA`
##'   if no data are present.
##'
##' @param n_particles Number of particles to use
##'
##' @param n_threads Number of threads to use
##'
##' @param seed Random seed to use
##'
##' @param compiled_compare Logical, indicating if we should use the
##'   new compiled compare function (this will shortly become the
##'   default).
##'
##' @return A [`mcstate::particle_filter`] object
##' @export
carehomes_particle_filter <- function(data, n_particles,
                                      n_threads = 1L, seed = NULL,
                                      compiled_compare = FALSE) {
  mcstate::particle_filter$new(
    carehomes_particle_filter_data(data),
    carehomes,
    n_particles,
    if (compiled_compare) NULL else carehomes_compare,
    carehomes_index,
    carehomes_initial,
    n_threads,
    seed)
}


carehomes_particle_filter_data <- function(data) {
  required <- c("icu", "general", "hosp", "deaths_hosp", "deaths_carehomes",
                "deaths_comm", "deaths_non_hosp", "deaths", "admitted",
                "diagnoses", "all_admission", "npos_15_64", "ntot_15_64",
                "pillar2_pos", "pillar2_tot", "pillar2_cases",
                "pillar2_over25_pos", "pillar2_over25_tot",
                "pillar2_over25_cases", "react_pos", "react_tot")

  verify_names(data, required, allow_extra = TRUE)

  if (any(!is.na(data$deaths) &
           (!is.na(data$deaths_comm) | !is.na(data$deaths_hosp) |
            !is.na(data$deaths_carehomes) | !is.na(data$deaths_non_hosp)))) {
    stop("Deaths are not consistently split into total vs hospital/non-hospital
          or hospital/care homes/community")
  }
  if (any(!is.na(data$deaths_non_hosp) &
          (!is.na(data$deaths_comm) | !is.na(data$deaths_carehomes)))) {
    stop("Non-hospital deaths are not consistently split into total vs care
         homes/community")
  }

  check_pillar2_streams <- function(df) sum(c(any(!is.na(df$pillar2_pos)) |
                             any(!is.na(df$pillar2_tot)),
                           any(!is.na(df$pillar2_cases)),
                           any(!is.na(df$pillar2_over25_pos)) |
                             any(!is.na(df$pillar2_over25_tot)),
                           any(!is.na(df$pillar2_over25_cases))))

  if (is.null(data$population)) {
    if (check_pillar2_streams(data) > 1) {
      stop("Cannot fit to more than one pillar 2 data stream")
    }
  } else {
    lapply(split(data, data$population), function(x) {
      if (check_pillar2_streams(x) > 1) {
        stop(sprintf("Cannot fit to more than one pillar 2 data stream for
                      region %s", x$population[[1]]))
      }
    })
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


carehomes_data <- function(data, start_date, dt) {
  expected <- c(icu = NA_real_, general = NA_real_, hosp = NA_real_,
                deaths_hosp = NA_real_, deaths_non_hosp = NA_real_,
                deaths_comm = NA_real_, deaths_carehomes = NA_real_,
                deaths = NA_real_, admitted = NA_real_, diagnoses = NA_real_,
                all_admission = NA_real_, npos_15_64 = NA_real_,
                ntot_15_64 = NA_real_, pillar2_pos = NA_real_,
                pillar2_tot = NA_real_, pillar2_cases = NA_real_,
                pillar2_over25_pos = NA_real_, pillar2_over25_tot = NA_real_,
                pillar2_over25_cases = NA_real_,  react_pos = NA_real_,
                react_tot = NA_real_, strain_non_variant = NA_real_,
                strain_tot = NA_real_)
  data <- sircovid_data(data, start_date, dt, expected)
  carehomes_particle_filter_data(data)
}


scale_severity <- function(severity, strain_rel_severity,
                           which = c("p_G_D", "p_H_D", "p_W_D", "p_ICU_D")) {
  severity[which] <- lapply(severity[which], function(x) {
    x <- matrix(x, nrow = length(x), ncol = length(strain_rel_severity))
    prob <- matrix(strain_rel_severity, nrow = nrow(x),
                  ncol = length(strain_rel_severity), byrow = TRUE)
    x <- x * prob
    x[x > 1] <- 1
    x
  })
  severity
}
