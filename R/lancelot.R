##' @name lancelot
##' @title The lancelot sircovid model
##'
##' Our "current" sircovid model. This is a dust model.
##'
##' @export lancelot
NULL


##' Parameters for the "[lancelot]" model.
##'
##' @title Parameters for the lancelot model
##'
##' @inheritParams basic_parameters
##'
##' @param population Population data. A vector of length 17 for the population
##'   size for age groups 0-4, 5-9, ..., 75-79, 80+. If `NULL`, the population
##'   data will be sourced within the package for the specified `region` (only
##'   if available).
##'
##' @param carehome_beds The number of care home beds in the region. If `NULL`,
##'   this will be sourced within the package for the specified `region` (only
##'   if available).
##'
##' @param severity Severity data, via Bob Verity's `markovid`
##'   package. This needs to be `NULL` (use the default bundled data
##'   version in the package), a [data.frame] object (for raw severity
##'   data) or a list (for data that has already been processed by
##'   `sircovid` for use).  New severity data comes from Bob Verity
##'   via the markovid package, and needs to be carefully calibrated
##'   with the progression parameters.
##'
##' @param progression Progression data
##'
##' @param observation Either `NULL` or a list of observation parameters. If
##'   `NULL`, then a list of observation parameters will be generated using
##'   `lancelot_parameters_observation(exp_noise)`
##'
##' @param sens_and_spec Either `NULL` or a list of diagnostic test sensitivity
##'   and specificity parameters. If `NULL`, then a list of sensitivity and
##'   specificity parameters will be generated using
##'   `lancelot_parameters_sens_and_spec()`
##'
##' @param eps Change in contact rate for carehome residents
##'
##' @param m_CHW Contact rate between carehome workers and either
##'   residents or workers
##'
##' @param m_CHR Contact rate between carehome residents
##'
##' @param strain_transmission Vector of length two for relative
##'   transmissibility of each strain modelled. Length will define the number of
##'   strains used in the model, either 1 or 2.
##'
##' @param strain_seed_date Either `NULL` (no seeding) or a
##'   [sircovid::sircovid_date] corresponding to the date the seeding of strain
##'   2 begins.
##'
##' @param strain_seed_size Either `NULL` (no seeding) or the size of strain 2
##'   seeding from the S to E compartment; all seeding is in the 15-19 year old
##'   group from `strain_seed_date` according to `strain_seed_pattern`.
##'
##' @param strain_seed_pattern Either `NULL` (no seeding) or a vector of seeding
##'   weights for the initial seeding. The length represents the number of steps
##'   to seed over from the `start_date`, and the `strain_seed_size` is split
##'   over these steps according to those weights. If `start_seed_date` is not a
##'   multiple of the step size (and thus falls between two steps) then we
##'   weight over an additional step and adjust the weights according to how far
##'   the `strain_seed_date` is from the previous full step.
##'
##' @param strain_rel_gamma_E Vector of relative rates of progression out of
##' E (gamma_E) for each
##'   strain modelled. If `1` all strains have same rates. Otherwise vector of
##'   same length as `strain_transmission`, with entries that determines the
##'   relative scaling of the defaults for each strain.
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
##' @param strain_rel_p_sympt Vector of relative probabilities of
##'   symptoms for
##'   each strain modelled. If `1` all strains have same
##'   probabilities of symptoms. Otherwise vector of same length as
##'   `strain_transmission`, where the first value should be 1 (for the first
##'   strain) and subsequent values between 0 and 1. In this case parameters
##'   will be "mirrored" for pseudostrains i.e. the relative probability of
##'   symptoms will be assume the same irrespective of previous infection
##'   with another strain. Alternatively, a vector of twice the length of
##'   `strain_transmission` can be provided to allow specifying directly
##'   relative probability of symptoms for each pseudostrain
##'   (with strain 3 - 1.2 and strain 4 = 2.1). To ensure valid
##'   probabilities, p_sympt is upper-truncated at 1 after scaling.
##'
##' @param strain_rel_p_hosp_if_sympt Vector of relative probabilities of
##'   hospitalisation given symptoms for
##'   each strain modelled. If `1` all strains have same
##'   probabilities of hospitalisation. Otherwise vector of same length as
##'   `strain_transmission`, where the first value should be 1 (for the first
##'   strain) and subsequent values between 0 and 1. In this case parameters
##'   will be "mirrored" for pseudostrains i.e. the relative probability of
##'   hospitalisation will be assume the same irrespective of previous infection
##'   with another strain. Alternatively, a vector of twice the length of
##'   `strain_transmission` can be provided to allow specifying directly
##'   relative probability of hospitalisation for each pseudostrain
##'   (with strain 3 - 1.2 and strain 4 = 2.1). To ensure valid
##'   probabilities, p_hosp_if_sympt is upper-truncated at 1 after scaling.
##'
##' @param strain_rel_p_icu Vector of relative probabilities of
##'   icu given hospitalised for
##'   each strain modelled. If `1` all strains have same
##'   probabilities of icu admission. Otherwise vector of same length as
##'   `strain_transmission`, where the first value should be 1 (for the first
##'   strain) and subsequent values between 0 and 1. In this case parameters
##'   will be "mirrored" for pseudostrains i.e. the relative probability of
##'   icu admission will be assume the same irrespective of previous infection
##'   with another strain. Alternatively, a vector of twice the length of
##'   `strain_transmission` can be provided to allow specifying directly
##'   relative probability of icu admission for each pseudostrain
##'   (with strain 3 - 1.2 and strain 4 = 2.1). To ensure valid
##'   probabilities, p_icu is upper-truncated at 1 after scaling.
##'
##' @param strain_rel_p_death Vector of relative probabilities of death for
##'   each strain modelled. If `1` all strains have same
##'   probabilities of death. Otherwise vector of same length as
##'   `strain_transmission`, where the first value should be 1 (for the first
##'   strain) and subsequent values between 0 and 1. In this case parameters
##'   will be "mirrored" for pseudostrains i.e. the relative probability of
##'   death will be assume the same irrespective of previous infection with
##'   another strain. Alternatively, a vector of twice the length of
##'   `strain_transmission` can be provided to allow specifying directly
##'   relative probability of death for each pseudostrain (with strain 3 - 1.2
##'   and strain 4 = 2.1). To ensure valid
##'   probabilities, p_death is upper-truncated at 1 after scaling.
##'
##' @param rel_susceptibility A vector or array of values representing the
##'   relative susceptibility of individuals in different vaccination groups.
##'   If a vector, the first value should be 1 (for the non-vaccinated group)
##'   and subsequent values be between 0 and 1. In that case relative
##'   susceptibility will be the same across all age groups within one
##'   vaccination category, and will be the same for all pathogen strains.
##'   Specifying an array instead of a vector allows
##'   different relative susceptibilities by age (first dimension of the array),
##'   pathogen strain (second dimension) and vaccination group
##'   (third dimension); in that case, the first layer (3rd dimension) of
##'   rel_susceptibility should be 1 (for the non-vaccinated group)
##'   for the first column (first infection with first strain) and other
##'   values between 0 and 1
##'
##' @param rel_p_sympt A vector or matrix of values of same dimension as
##'   rel_susceptibility representing the
##'   relative probability of symptomatic infection in different
##'   vaccination groups. If a vector, the first value should be 1 (for the
##'   non-vaccinated group) and subsequent values be between 0 and 1.
##'   In that case the relative reduction in probability of symptomatic
##'   infection will be the same across all age groups and all strains
##'   within one vaccination category.
##'   Specifying an array instead of a vector allows different relative
##'   reductions in probability of symptomatic infection by age (first
##'   dimension of the array), pathogen strain (second dimension) and
##'   vaccination group (third dimension); in that case,
##'   the first layer of rel_p_sympt should be 1 (for the non-vaccinated group)
##'   for the first column (first infection with first strain) and other
##'   values between 0 and 1
##'
##' @param rel_p_hosp_if_sympt A vector or array of values of same dimension as
##'   rel_susceptibility representing the
##'   relative probability of hospitalisation for symptomatic cases in different
##'   vaccination groups. If a vector, the first value should be 1 (for the
##'   non-vaccinated group) and subsequent values be between 0 and 1.
##'   In that case the relative reduction in probability of hospitalisation for
##'   symptomatic cases will be the same across all age groups and all strains
##'   within one vaccination category.
##'   Specifying an array instead of a vector allows different relative
##'   reductions in probability of hospitalisation for symptomatic cases by age
##'   (first dimension of the array), pathogen strain (second dimension) and
##'   vaccination group (third dimension); in that case,
##'   the first layer of rel_p_hosp_if_sympt should be 1 (for the
##'   non-vaccinated group) for the first column
##'   (first infection with first strain) and other values between 0 and 1
##'
##' @param rel_p_death A vector or array of values of same dimension as
##'   rel_susceptibility representing the relative probability of death for
##'    severe cases (either in hospital or the community) in different
##'   vaccination groups. If a vector, the first value should be 1 (for the
##'   non-vaccinated group) and subsequent values be between 0 and 1.
##'   In that case the relative reduction in probability of death for
##'   severe cases will be the same across all age groups and all strains
##'   within one vaccination category.
##'   Specifying an array instead of a vector allows different relative
##'   reductions in probability of death for severe cases by age
##'   (first dimension of the array), pathogen strain (second dimension) and
##'   vaccination group (third dimension); in that case,
##'   the first layer of rel_p_death should be 1 (for the non-vaccinated group)
##'   for the first column (first infection with first strain) and other
##'   values between 0 and 1
##'
##' @param rel_infectivity A vector or array of values representing the
##'   relative infectivity of individuals in different vaccination groups,
##'   if they are infected.
##'   If a vector, the first value should be 1 (for the non-vaccinated group)
##'   and subsequent values be between 0 and 1. In that case relative
##'   infectivity will be the same across all age groups and all strains
##'   within one vaccination category.
##'   Specifying an array instead of a vector allows
##'   different relative infectivities by age
##'   (first dimension of the array), pathogen strain (second dimension) and
##'   vaccination group (third dimension); in that case,
##'   the first layer of rel_infectivity should be 1 (for the
##'   non-vaccinated group)
##'   for the first column (first infection with first strain) and other
##'   values between 0 and 1
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
##' @param vaccine_index_booster The index to use for the booster dose
##'
##' @param vaccine_catchup_fraction A value between 0 and 1 indicating the
##'   proportion of doses not distributed according to schedule (e.g. because
##'   too many people were in the I or H compartments and could not be
##'   vaccinated at the scheduled time) that we postpone to a later date.
##'   A value of 0 means we do not catch up at all on any missed doses; a
##'   value of 1 means we try to catch up for all missed doses. This is set
##'   to 1 by default
##'
##' @param vacc_skip_from An integer, the vaccine stratum index from which we
##'   allow a "vaccine skip move", which enables individuals to progress through
##'   more than one vaccine class in a single step. Default is 1
##'
##' @param vacc_skip_to An integer, the vaccine stratum index representing the
##'   vaccine stratum an individual in stratum `vacc_skip_from` can move to in a
##'   "vaccine skip move". Must be either equal to `vacc_skip_from` (to allow
##'   the vaccine skip move to be switched "off"), or greater than or equal to
##'   `vacc_skip_from + 2` (we already account for moves between successive
##'   vaccine strata). Default is 1
##'
##' @param vacc_skip_progression_rate Either a vector of length 19 representing
##'   the base progression rate for "vaccine skip moves" for each group, or a
##'   scalar representing the base progression rate for all groups. Must be 0
##'   (which is the default) if `vacc_skip_to == vacc_skip_from`
##'
##' @param vacc_skip_weight A scalar weight. If movement into `vacc_skip_to`
##'   from `vacc_skip_to - 1` is controlled by doses, then the "vaccine skip
##'   move" is too, and the weight represents how much those in `vacc_skip_from`
##'   are weighted in dose distribution relative to those in `vacc_skip_to - 1`.
##'   Must be between 0 and 1. Must be 0 (which is the default) if
##'   `vacc_skip_to == vacc_skip_from`
##'
##' @param n_doses Number of doses given out, including boosters. Default is
##'   2.
##'
##' @param waning_rate A single value or a vector of values representing the
##'   rates of waning of immunity after infection; if a single value the same
##'   rate is used for all age groups; if a vector of values if used it should
##'   have one value per age group.
##'
##' @param cross_immunity A value or vector of same length as
##'   `strain_transmission` that controls the amount of immunity conferred by
##'   previous infection with one strain. If a scalar is given
##'   then same level of cross immunity is assumed between both strains.
##'   Otherwise a vector of length two should be provided where the first value
##'   is the relative protection against infection with strain 2
##'   following infection with strain 1 (i.e. while in the R1 compartment),
##'   and vice versa for the second value. Values between 0 and 1 are allowed
##'   with values of 1 (default) indicating complete cross-immunity,
##'   and values of 0 mean no cross-immunity. Modelling
##'  'superinfections' (being exposed to one strain after
##'   recovering from another) can be turned off by setting
##'   `cross_immunity = 1`.
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
##' lancelot_parameters(sircovid_date("2020-02-01"), region)
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
##' # and the risk of death for those with severe disease
##' rel_p_death <- c(1, 0.9, 0.9)
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
##' p <- lancelot_parameters(
##'        sircovid_date("2020-02-01"), region,
##'        rel_susceptibility = rel_susceptibility,
##'        rel_p_sympt = rel_p_sympt,
##'        rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
##'        rel_p_death = rel_p_death,
##'        rel_infectivity = rel_infectivity,
##'        vaccine_progression_rate = vaccine_progression_rate,
##'        vaccine_schedule = schedule,
##'        vaccine_index_dose2 = 2)
##'
##' # vaccination parameters are automatically copied across all age groups#
##' # (and across strains but here we only have 1 strain which is the 2nd
##' # dimension here)
##' p$rel_susceptibility
##' p$rel_p_sympt
##' p$rel_p_hosp_if_sympt
##' p$rel_p_death
##' p$rel_infectivity
##' # Note that this is only the "base" rate as we fill in the first
##' # column dynamically based on vaccine_daily_doses
##' p$vaccine_progression_rate
##'
##' ### same example as above BUT assume a different effect of vaccine in the
##' ### first age group
##' n_groups <- 19
##' n_strains <- 1
##'
##' # Assumption: vaccine is twice more effective at reducing susceptibility
##' # in the first age group
##' rel_susceptibility_agegp1 <- c(1, 0.4, 0.25)
##' rel_susceptibility_other_agegp <- c(1, 0.8, 0.5)
##' rel_susceptibility <- array(NA, dim = c(n_groups, n_strains, 3))
##' rel_susceptibility[1, , ] <- rel_susceptibility_agegp1
##' for (i in seq(2, n_groups)) {
##'   rel_susceptibility[i, , ] <- rel_susceptibility_other_agegp
##' }
##' rel_susceptibility
##'
##' # But vaccine has the same impact on probability of symptoms and
##' # hospitalisation for the symptomatic across all age groups
##' rel_p_sympt <- array(rep(rel_p_sympt, each = n_groups),
##'   dim = c(n_groups, n_strains, 3))
##' rel_p_hosp_if_sympt <-
##'   array(rep(rel_p_hosp_if_sympt, each = n_groups),
##'   dim = c(n_groups, n_strains, 3))
##' rel_p_death <-
##'   array(rep(rel_p_hosp_if_sympt, each = n_groups),
##'   dim = c(n_groups, n_strains, 3))
##'
##' # And vaccine has the same impact on onwards infectivity across age groups
##' rel_infectivity <- array(rep(rel_infectivity, each = n_groups),
##'   dim = c(n_groups, n_strains, 3))
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
##' p <- lancelot_parameters(
##'        sircovid_date("2020-02-01"), region,
##'        rel_susceptibility = rel_susceptibility,
##'        rel_p_sympt = rel_p_sympt,
##'        rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
##'        rel_p_death = rel_p_death,
##'        rel_infectivity = rel_infectivity,
##'        vaccine_progression_rate = vaccine_progression_rate,
##'        vaccine_schedule = schedule,
##'        vaccine_index_dose2 = 2)
##'
##' # TODO: add an example of manually set up vaccine schedule
##'
lancelot_parameters <- function(start_date, region,
                                beta_date = NULL, beta_value = NULL,
                                beta_type = "piecewise-linear",
                                population = NULL,
                                carehome_beds = NULL,
                                severity = NULL,
                                progression = NULL,
                                observation = NULL,
                                sens_and_spec = NULL,
                                initial_seed_size = 30,
                                initial_seed_pattern = 1,
                                eps = 0.1,
                                m_CHW = 4e-6,
                                m_CHR = 5e-5,
                                strain_transmission = 1,
                                strain_seed_date = NULL,
                                strain_seed_size = NULL,
                                strain_seed_pattern = NULL,
                                strain_rel_gamma_E = 1,
                                strain_rel_gamma_A = 1,
                                strain_rel_gamma_P = 1,
                                strain_rel_gamma_C_1 = 1,
                                strain_rel_gamma_C_2 = 1,
                                strain_rel_p_sympt = 1,
                                strain_rel_p_hosp_if_sympt = 1,
                                strain_rel_p_icu = 1,
                                strain_rel_p_death = 1,
                                rel_susceptibility = 1,
                                rel_p_sympt = 1,
                                rel_p_hosp_if_sympt = 1,
                                rel_p_death = 1,
                                rel_infectivity = 1,
                                vaccine_progression_rate = NULL,
                                vaccine_schedule = NULL,
                                vaccine_index_dose2 = NULL,
                                vaccine_index_booster = NULL,
                                vaccine_catchup_fraction = 1,
                                n_doses = 2L,
                                vacc_skip_progression_rate = 0,
                                vacc_skip_from = 1L,
                                vacc_skip_to = 1L,
                                vacc_skip_weight = 0,
                                waning_rate = 0,
                                exp_noise = 1e6,
                                cross_immunity = 1) {

  n_real_strains <- length(strain_transmission)

  if (!is.null(population)) {
    if (!is.null(dim(population)) || length(population) != 17L) {
      stop("If population is specified it must be a vector of length 17")
    }
    population <- assert_integer(population)
    population <- assert_non_negative(population)
  }

  ret <- sircovid_parameters_shared(start_date, region,
                                    beta_date, beta_value, beta_type,
                                    population,
                                    initial_seed_pattern, initial_seed_size)

  ## These are only used here, and are fixed
  carehome_occupancy <- 0.742
  carehome_workers_per_resident <- 1

  ## These are used in constructing the initial population vectors (S0)
  carehome_beds <- carehome_beds %||% sircovid_carehome_beds(region)
  carehome_residents <- round(carehome_beds * carehome_occupancy)
  carehome_workers <- round(carehome_residents * carehome_workers_per_resident)

  ## TODO: it's probably the case that having some tree structure here
  ## would make this nicer to work with, but we should do this
  ## consistently through the other parameters too. Keeping the
  ## progression and severity parameters together for example.
  ret$carehome_beds <- carehome_beds
  ret$carehome_residents <- carehome_residents
  ret$carehome_workers <- carehome_workers

  if (n_real_strains > 2) {
    stop("Only 1 or 2 strains valid ('strain_transmission' too long)'.")
  }

  severity <- severity %||% lancelot_parameters_severity(ret$dt, severity)

  progression <- progression %||% lancelot_parameters_progression(ret$dt)

  waning <- lancelot_parameters_waning(waning_rate)

  ret$m <- lancelot_transmission_matrix(eps, m_CHW, m_CHR, region, population)

  ret$N_tot <- lancelot_population(ret$population, carehome_workers,
                                   carehome_residents)

  ret$initial_seed_size <- initial_seed_size

  ## control cross-immunity
  ret$cross_immunity <- assert_proportion(
    recycle(cross_immunity, n_real_strains)
  )

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
  ret$N_tot_under15 <- sum(ret$N_tot[1:3])
  ret$N_tot_15_24 <- sum(ret$N_tot[4:5])
  ## assume CHW [18] are equally distributed amongst 25-64 age bands
  ret$N_tot_25_49 <- sum(ret$N_tot[6:10]) + (sum(ret$N_tot[18]) / 8) * 5
  ret$N_tot_50_64 <- sum(ret$N_tot[11:13]) + (sum(ret$N_tot[18]) / 8) * 3
  ## assume CHR [19] are 1/4 aged 65-79 and 3/4 80 plus
  ret$N_tot_65_79 <- sum(ret$N_tot[14:16]) + sum(ret$N_tot[19]) * 0.25
  ret$N_tot_80_plus <- sum(ret$N_tot[17]) + sum(ret$N_tot[19]) * 0.75

  ## relative transmissibility of various I compartments
  ret$I_A_transmission <- 0.223
  ret$I_P_transmission <- 1
  ret$I_C_1_transmission <- 1
  ret$I_C_2_transmission <- 0

  ## All observation parameters:
  observation <- observation %||% lancelot_parameters_observation(exp_noise)

  sens_and_spec <- sens_and_spec %||% lancelot_parameters_sens_and_spec()

  ret$n_groups <- ret$n_age_groups + 2L

  ## number of strains and relative transmissibility
  strain <- lancelot_parameters_strain(
    strain_transmission, strain_seed_date, strain_seed_size,
    strain_seed_pattern, ret$dt)

  ## vaccination
  vaccination <- lancelot_parameters_vaccination(ret$N_tot,
                                                 ret$dt,
                                                 rel_susceptibility,
                                                 rel_p_sympt,
                                                 rel_p_hosp_if_sympt,
                                                 rel_p_death,
                                                 rel_infectivity,
                                                 vaccine_progression_rate,
                                                 vaccine_schedule,
                                                 vaccine_index_dose2,
                                                 vaccine_index_booster,
                                                 strain$n_strains,
                                                 vaccine_catchup_fraction,
                                                 n_doses)

  ## vacc_skip parameters
  vacc_classes <- seq_len(vaccination$n_vacc_classes)
  if (!vacc_skip_from %in% vacc_classes || !vacc_skip_to %in% vacc_classes) {
    stop(sprintf("There are %s vaccine classes so 'vacc_skip_from' and
                 'vacc_skip_to' must each be one of: %s",
                 vaccination$n_vacc_classes,
                 paste(vacc_classes, collapse = ", ")))
  }
  if (vacc_skip_to != vacc_skip_from && vacc_skip_to < vacc_skip_from + 2) {
    stop("Require vacc_skip_to = vacc_skip_from or vacc_skip_to >=
         vacc_skip_from + 2")
  }
  if (vacc_skip_to == vacc_skip_from && vacc_skip_weight != 0) {
    stop("Require vacc_skip_weight = 0 as vacc_skip_to = vacc_skip_from")
  }
  if (vacc_skip_to == vacc_skip_from && any(vacc_skip_progression_rate != 0)) {
    stop("Require vacc_skip_progression_rate = 0 as vacc_skip_to =
         vacc_skip_from")
  }
  if (length(vacc_skip_progression_rate) == 1) {
    ret$vacc_skip_progression_rate_base <-
      rep(vacc_skip_progression_rate, ret$n_groups)
  } else {
    if (length(vacc_skip_progression_rate) != ret$n_groups) {
      stop(sprintf("'vacc_skip_progression_rate' must be a scalar or a vector of
                   length %s", ret$n_groups))
    } else {
      ret$vacc_skip_progression_rate_base <- vacc_skip_progression_rate
    }
  }
  ret$vacc_skip_to <- vacc_skip_to
  ret$vacc_skip_from <- vacc_skip_from
  ret$vacc_skip_weight <- assert_proportion(vacc_skip_weight)
  if (vacc_skip_to == 1) {
    ret$vacc_skip_dose <- 0
  } else {
    ret$vacc_skip_dose <- vaccination$index_dose_inverse[vacc_skip_to - 1]
  }

  strain_rel_p_death <- process_strain_rel_p(strain_rel_p_death,
                                             strain$n_strains,
                                             n_real_strains)
  ret$strain_rel_p_ICU_D <- strain_rel_p_death
  ret$strain_rel_p_H_D <- strain_rel_p_death
  ret$strain_rel_p_W_D <- strain_rel_p_death
  ret$strain_rel_p_G_D <- strain_rel_p_death

  strain_rel_p_icu <- process_strain_rel_p(strain_rel_p_icu,
                                           strain$n_strains,
                                           n_real_strains)
  ret$strain_rel_p_icu <- strain_rel_p_icu

  strain_rel_p_hosp_if_sympt <- process_strain_rel_p(strain_rel_p_hosp_if_sympt,
                                                     strain$n_strains,
                                                     n_real_strains)
  ret$strain_rel_p_hosp_if_sympt <- strain_rel_p_hosp_if_sympt

  strain_rel_p_sympt <- process_strain_rel_p(strain_rel_p_sympt,
                                                     strain$n_strains,
                                                     n_real_strains)
  ret$strain_rel_p_sympt <- strain_rel_p_sympt

  # combine strain-specific relative severity with vaccination reduction in
  # probability of death

  rel_p_death <- build_rel_param(rel_p_death, strain$n_strains,
                                 vaccination$n_vacc_classes, "rel_p_death")

  ret$rel_p_ICU <- array(1, c(ret$n_groups, strain$n_strains,
                              vaccination$n_vacc_classes))

  ret$rel_p_R <- array(1, c(ret$n_groups, strain$n_strains,
                            vaccination$n_vacc_classes))

  ret$rel_p_ICU_D <- rel_p_death
  ret$rel_p_H_D <- rel_p_death
  ret$rel_p_W_D <- rel_p_death
  ret$rel_p_G_D <- rel_p_death

  strain_rel_gammas <- list(E = strain_rel_gamma_E,
                            A = strain_rel_gamma_A,
                            P = strain_rel_gamma_P,
                            C_1 = strain_rel_gamma_C_1,
                            C_2 = strain_rel_gamma_C_2,
                            ICU_pre = NULL,
                            H_D = NULL,
                            H_R = NULL,
                            ICU_D = NULL,
                            ICU_W_D = NULL,
                            ICU_W_R = NULL,
                            W_D = NULL,
                            W_R = NULL,
                            G_D = NULL)
  for (name in names(strain_rel_gammas)) {
    rel_gamma <- strain_rel_gammas[[name]]
    rel_gamma_name <- paste0("rel_gamma_", name)
    if (is.null(rel_gamma)) {
      ret[[rel_gamma_name]] <- rep(1, strain$n_strains)
    } else {
      rel_gamma <- recycle(assert_relatives(rel_gamma),
                           n_real_strains)
      if (length(rel_gamma) == 2) {
        ret[[rel_gamma_name]] <- mirror_strain(rel_gamma)
      } else {
        ret[[rel_gamma_name]] <- rep(1, strain$n_strains)
      }
    }
  }

  out <- c(ret, severity, progression, strain, vaccination, waning,
           observation, sens_and_spec)

  lancelot_check_severity(out)
}

process_strain_rel_p <- function(p, n_strains, n_real_strains) {
  if (length(p) < n_strains) {
    p <- recycle(p, n_real_strains)
    if (n_real_strains > 1) {
      p <- mirror_strain(p)
    }
  }
  p
}

##' Index of "interesting" elements for the lancelot model. This function
##' conforms to the mcstate interface.
##'
##' @title Index of lancelot model
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
##' p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")
##' mod <- lancelot$new(p, 0, 10)
##' lancelot_index(mod$info())
lancelot_index <- function(info) {
  index <- info$index
  dim <- info$dim

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
                  deaths_hosp_0_49_inc = index[["D_hosp_0_49_inc"]],
                  deaths_hosp_50_54_inc = index[["D_hosp_50_54_inc"]],
                  deaths_hosp_55_59_inc = index[["D_hosp_55_59_inc"]],
                  deaths_hosp_60_64_inc = index[["D_hosp_60_64_inc"]],
                  deaths_hosp_65_69_inc = index[["D_hosp_65_69_inc"]],
                  deaths_hosp_70_74_inc = index[["D_hosp_70_74_inc"]],
                  deaths_hosp_75_79_inc = index[["D_hosp_75_79_inc"]],
                  deaths_hosp_80_plus_inc = index[["D_hosp_80_plus_inc"]],
                  admitted_inc = index[["admit_conf_inc"]],
                  diagnoses_inc = index[["new_conf_inc"]],
                  sero_pos_1 = index[["sero_pos_1"]],
                  sero_pos_2 = index[["sero_pos_2"]],
                  sympt_cases = index[["cum_sympt_cases"]],
                  sympt_cases_non_variant =
                    index[["cum_sympt_cases_non_variant"]],
                  sympt_cases_over25 = index[["cum_sympt_cases_over25"]],
                  sympt_cases_non_variant_over25 =
                    index[["cum_sympt_cases_non_variant_over25"]],
                  sympt_cases_inc = index[["sympt_cases_inc"]],
                  sympt_cases_non_variant_inc =
                    index[["sympt_cases_non_variant_inc"]],
                  sympt_cases_over25_inc = index[["sympt_cases_over25_inc"]],
                  sympt_cases_under15_inc = index[["sympt_cases_under15_inc"]],
                  sympt_cases_15_24_inc = index[["sympt_cases_15_24_inc"]],
                  sympt_cases_25_49_inc = index[["sympt_cases_25_49_inc"]],
                  sympt_cases_50_64_inc = index[["sympt_cases_50_64_inc"]],
                  sympt_cases_65_79_inc = index[["sympt_cases_65_79_inc"]],
                  sympt_cases_80_plus_inc = index[["sympt_cases_80_plus_inc"]],
                  sympt_cases_non_variant_over25_inc =
                    index[["sympt_cases_non_variant_over25_inc"]],
                  react_pos = index[["react_pos"]])

  ## Only incidence versions for the likelihood now. We add time here so it
  ## can be used in the compare, without having to save it
  index_run <- c(time = index[["time"]],
                 index_core[c("icu", "general", "deaths_carehomes_inc",
                              "deaths_comm_inc", "deaths_hosp_inc",
                              "deaths_hosp_0_49_inc", "deaths_hosp_50_54_inc",
                              "deaths_hosp_55_59_inc", "deaths_hosp_60_64_inc",
                              "deaths_hosp_65_69_inc", "deaths_hosp_70_74_inc",
                              "deaths_hosp_75_79_inc",
                              "deaths_hosp_80_plus_inc", "admitted_inc",
                              "diagnoses_inc", "sero_pos_1", "sero_pos_2",
                              "sympt_cases_inc", "sympt_cases_non_variant_inc",
                              "sympt_cases_over25_inc",
                              "sympt_cases_under15_inc",
                              "sympt_cases_15_24_inc", "sympt_cases_25_49_inc",
                              "sympt_cases_50_64_inc", "sympt_cases_65_79_inc",
                              "sympt_cases_80_plus_inc",
                              "sympt_cases_non_variant_over25_inc",
                              "react_pos")])

  ## Variables that we want to save for post-processing
  index_save <- c(hosp = index[["hosp_tot"]],
                  deaths = index[["D_tot"]],
                  deaths_inc = index[["D_inc"]],
                  infections = index[["cum_infections"]],
                  infections_inc = index[["infections_inc"]])
  suffix <- paste0("_", c(sircovid_age_bins()$start, "CHW", "CHR"))
  ## NOTE: We do use the S category for the Rt calculation in some
  ## downstream work, so this is going to require some work to get
  ## right.

  n_vacc_classes <- info$dim$S[[2]]
  n_strains <- info$dim$prob_strain
  if (n_strains == 2) {
    n_tot_strains <- 4
  } else {
    n_tot_strains <- 1
  }

  ## age varying only
  index_cum_admit <- calculate_index(index$cum_admit_by_age,
                                     list(), suffix, "cum_admit")

  ## age x vacc class
  index_S <- calculate_index(index$S, list(n_vacc_classes), suffix, "S")
  index_diagnoses_admitted <- calculate_index(index$diagnoses_admitted,
                                              list(n_vacc_classes), suffix,
                                              "diagnoses_admitted")
  index_cum_infections_disag <- calculate_index(index$cum_infections_disag,
                                                list(n_vacc_classes), suffix,
                                                "cum_infections_disag")
  index_cum_n_vaccinated <- calculate_index(index$cum_n_vaccinated,
                                            list(n_vacc_classes), suffix,
                                            "cum_n_vaccinated")
  index_D <- calculate_index(index$D, list(n_vacc_classes), suffix, "D_all")

  ## (real) strain only
  index_prob_strain <- calculate_index(index$prob_strain, list(n_strains),
                                       NULL, "prob_strain")

  ## age x (total) strain x vacc class
  if (n_strains == 1) {
    index_R <- index$R
  } else {
    index_R <- array(index$R, dim$R)[, 1:2, , drop = FALSE]
  }
  index_R <- calculate_index(index_R,
                             list(S = n_strains, V = n_vacc_classes),
                             suffix, "R")

  list(run = index_run,
       state = c(index_core, index_save, index_S, index_R,
                 index_cum_admit, index_D,
                 index_diagnoses_admitted, index_cum_infections_disag,
                 index_prob_strain, index_cum_n_vaccinated
       ))
}


##' Compare observed and modelled data from the [lancelot] model. This
##' conforms to the mcstate interface.
##'
##' @title Compare observed and modelled data for the lancelot model
##'
##' @param state State vector for the end of the current day. This is
##'   assumed to be filtered following [lancelot_index()] so contains
##'   10 rows corresponding to ICU, general beds, admissions, deaths and
##'   seroconversion compartments.
##'
##' @param observed Observed data. At the moment please see the tests
##'   for a full list as this changes frequently (and this function
##'   may be removed in future).
##'
##' @param pars A list of parameters, as created by
##'   [lancelot_parameters()]
##'
##' @return A vector of log likelihoods, the same length as the number
##'   of particles (the number of columns in the modelled state)
##'
##' @export
##' @importFrom stats dbinom
lancelot_compare <- function(state, observed, pars) {
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
  model_deaths_hosp_0_49 <- state["deaths_hosp_0_49_inc", ]
  model_deaths_hosp_50_54 <- state["deaths_hosp_50_54_inc", ]
  model_deaths_hosp_55_59 <- state["deaths_hosp_55_59_inc", ]
  model_deaths_hosp_60_64 <- state["deaths_hosp_60_64_inc", ]
  model_deaths_hosp_65_69 <- state["deaths_hosp_65_69_inc", ]
  model_deaths_hosp_70_74 <- state["deaths_hosp_70_74_inc", ]
  model_deaths_hosp_75_79 <- state["deaths_hosp_75_79_inc", ]
  model_deaths_hosp_80_plus <- state["deaths_hosp_80_plus_inc", ]
  model_admitted <- state["admitted_inc", ]
  model_diagnoses <- state["diagnoses_inc", ]
  model_all_admission <- model_admitted + model_diagnoses
  model_sero_pos_1 <- state["sero_pos_1", ]
  model_sero_pos_2 <- state["sero_pos_2", ]
  model_sympt_cases <- state["sympt_cases_inc", ]
  model_sympt_cases_non_variant <- state["sympt_cases_non_variant_inc", ]
  model_sympt_cases_over25 <- state["sympt_cases_over25_inc", ]
  model_sympt_cases_under15 <- state["sympt_cases_under15_inc", ]
  model_sympt_cases_15_24 <- state["sympt_cases_15_24_inc", ]
  model_sympt_cases_25_49 <- state["sympt_cases_25_49_inc", ]
  model_sympt_cases_50_64 <- state["sympt_cases_50_64_inc", ]
  model_sympt_cases_65_79 <- state["sympt_cases_65_79_inc", ]
  model_sympt_cases_80_plus <- state["sympt_cases_80_plus_inc", ]
  model_sympt_cases_non_variant_over25 <-
    state["sympt_cases_non_variant_over25_inc", ]
  model_react_pos <- state["react_pos", ]

  ## calculate test positive probabilities for the various test data streams

  ## Pillar 2
  ## First determine which value is used for p_NC or phi_pillar2, based on
  ## whether it is a weekday or a weekend and age band. Note all values of the
  ## time state will be the same so we can just use the first value
  time <- state["time", 1L]

  p_NC_under15 <- if ((time + 3) %% 7 < 2) pars$p_NC_weekend_under15 else
    pars$p_NC_under15
  p_NC_15_24 <- if ((time + 3) %% 7 < 2) pars$p_NC_weekend_15_24 else
    pars$p_NC_15_24
  p_NC_25_49 <- if ((time + 3) %% 7 < 2) pars$p_NC_weekend_25_49 else
    pars$p_NC_25_49
  p_NC_50_64 <- if ((time + 3) %% 7 < 2) pars$p_NC_weekend_50_64 else
    pars$p_NC_50_64
  p_NC_65_79 <- if ((time + 3) %% 7 < 2) pars$p_NC_weekend_65_79 else
    pars$p_NC_65_79
  p_NC_80_plus <- if ((time + 3) %% 7 < 2) pars$p_NC_weekend_80_plus else
    pars$p_NC_80_plus

  phi_pillar2_cases_under15 <- if ((time + 3) %% 7 < 2)
    pars$phi_pillar2_cases_weekend_under15 else pars$phi_pillar2_cases_under15
  phi_pillar2_cases_15_24 <- if ((time + 3) %% 7 < 2)
    pars$phi_pillar2_cases_weekend_15_24 else pars$phi_pillar2_cases_15_24
  phi_pillar2_cases_25_49 <- if ((time + 3) %% 7 < 2)
    pars$phi_pillar2_cases_weekend_25_49 else pars$phi_pillar2_cases_25_49
  phi_pillar2_cases_50_64 <- if ((time + 3) %% 7 < 2)
    pars$phi_pillar2_cases_weekend_50_64 else pars$phi_pillar2_cases_50_64
  phi_pillar2_cases_65_79 <- if ((time + 3) %% 7 < 2)
    pars$phi_pillar2_cases_weekend_65_79 else pars$phi_pillar2_cases_65_79
  phi_pillar2_cases_80_plus <- if ((time + 3) %% 7 < 2)
    pars$phi_pillar2_cases_weekend_80_plus else pars$phi_pillar2_cases_80_plus

  ## Calculate Pillar 2 probability of positive test
  pillar2_under15_negs <- p_NC_under15 * (pars$N_tot_under15 -
                                            model_sympt_cases_under15)
  model_pillar2_under15_prob_pos <- test_prob_pos(model_sympt_cases_under15,
                                                  pillar2_under15_negs,
                                                  pars$pillar2_sensitivity,
                                                  pars$pillar2_specificity,
                                                  pars$exp_noise)

  pillar2_15_24_negs <- p_NC_15_24 * (pars$N_tot_15_24 -
                                        model_sympt_cases_15_24)
  model_pillar2_15_24_prob_pos <- test_prob_pos(model_sympt_cases_15_24,
                                                pillar2_15_24_negs,
                                                pars$pillar2_sensitivity,
                                                pars$pillar2_specificity,
                                                pars$exp_noise)

  pillar2_25_49_negs <- p_NC_25_49 * (pars$N_tot_25_49 -
                                        model_sympt_cases_25_49)
  model_pillar2_25_49_prob_pos <- test_prob_pos(model_sympt_cases_25_49,
                                                pillar2_25_49_negs,
                                                pars$pillar2_sensitivity,
                                                pars$pillar2_specificity,
                                                pars$exp_noise)

  pillar2_50_64_negs <- p_NC_50_64 * (pars$N_tot_50_64 -
                                        model_sympt_cases_50_64)
  model_pillar2_50_64_prob_pos <- test_prob_pos(model_sympt_cases_50_64,
                                                pillar2_50_64_negs,
                                                pars$pillar2_sensitivity,
                                                pars$pillar2_specificity,
                                                pars$exp_noise)

  pillar2_65_79_negs <- p_NC_65_79 * (pars$N_tot_65_79 -
                                        model_sympt_cases_65_79)
  model_pillar2_65_79_prob_pos <- test_prob_pos(model_sympt_cases_65_79,
                                                pillar2_65_79_negs,
                                                pars$pillar2_sensitivity,
                                                pars$pillar2_specificity,
                                                pars$exp_noise)

  pillar2_80_plus_negs <- p_NC_80_plus * (pars$N_tot_80_plus -
                                            model_sympt_cases_80_plus)
  model_pillar2_80_plus_prob_pos <- test_prob_pos(model_sympt_cases_80_plus,
                                                  pillar2_80_plus_negs,
                                                  pars$pillar2_sensitivity,
                                                  pars$pillar2_specificity,
                                                  pars$exp_noise)

  ## Pillar 2 over 25s
  pillar2_over25_negs <- pillar2_25_49_negs + pillar2_50_64_negs +
    pillar2_65_79_negs + pillar2_80_plus_negs
  model_pillar2_over25_prob_pos <- test_prob_pos(model_sympt_cases_over25,
                                                 pillar2_over25_negs,
                                                 pars$pillar2_sensitivity,
                                                 pars$pillar2_specificity,
                                                 pars$exp_noise)

  pillar2_negs <- pillar2_under15_negs + pillar2_15_24_negs +
    pillar2_over25_negs
  model_pillar2_prob_pos <- test_prob_pos(model_sympt_cases,
                                          pillar2_negs,
                                          pars$pillar2_sensitivity,
                                          pars$pillar2_specificity,
                                          pars$exp_noise)

  model_pillar2_under15_cases <-
    phi_pillar2_cases_under15 * model_sympt_cases_under15
  model_pillar2_15_24_cases <- phi_pillar2_cases_15_24 * model_sympt_cases_15_24
  model_pillar2_25_49_cases <- phi_pillar2_cases_25_49 * model_sympt_cases_25_49
  model_pillar2_50_64_cases <- phi_pillar2_cases_50_64 * model_sympt_cases_50_64
  model_pillar2_65_79_cases <- phi_pillar2_cases_65_79 * model_sympt_cases_65_79
  model_pillar2_80_plus_cases <-
    phi_pillar2_cases_80_plus * model_sympt_cases_80_plus
  model_pillar2_over25_cases <-
    model_pillar2_25_49_cases + model_pillar2_50_64_cases +
    model_pillar2_65_79_cases + model_pillar2_80_plus_cases
  model_pillar2_cases <- model_pillar2_under15_cases +
    model_pillar2_15_24_cases + model_pillar2_over25_cases


  ## REACT (Note that for REACT we exclude group 1 (0-4) and 19 (CHR))
  ## It is possible that model_react_pos > pars$N_tot_react, so we cap it to
  ## avoid probabilities > 1 here
  model_react_pos_capped <- pmin(model_react_pos, pars$N_tot_react)
  model_react_prob_pos <- test_prob_pos(model_react_pos_capped,
                                        pars$N_tot_react -
                                          model_react_pos_capped,
                                        pars$react_sensitivity,
                                        pars$react_specificity,
                                        pars$exp_noise)

  ## serology assay 1
  ## It is possible that model_sero_pos_1 > pars$N_tot_15_64, so we cap it to
  ## avoid probabilities > 1 here
  model_sero_pos_1_capped <- pmin(model_sero_pos_1, pars$N_tot_15_64)
  model_sero_prob_pos_1 <- test_prob_pos(model_sero_pos_1_capped,
                                         pars$N_tot_15_64 -
                                           model_sero_pos_1_capped,
                                         pars$sero_sensitivity_1,
                                         pars$sero_specificity_1,
                                         pars$exp_noise)

  ## serology assay 2
  ## It is possible that model_sero_pos_2 > pars$N_tot_15_64, so we cap it to
  ## avoid probabilities > 1 here
  model_sero_pos_2_capped <- pmin(model_sero_pos_2, pars$N_tot_15_64)
  model_sero_prob_pos_2 <- test_prob_pos(model_sero_pos_2_capped,
                                         pars$N_tot_15_64 -
                                           model_sero_pos_2_capped,
                                         pars$sero_sensitivity_2,
                                         pars$sero_specificity_2,
                                         pars$exp_noise)

  ## Strain
  model_strain_prob_pos <- test_prob_pos(
    model_sympt_cases_non_variant,
    model_sympt_cases - model_sympt_cases_non_variant,
    1, 1, pars$exp_noise)

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
  ll_deaths_hosp_0_49 <- ll_nbinom(observed$deaths_hosp_0_49,
                                   pars$phi_death_hosp * model_deaths_hosp_0_49,
                                   pars$kappa_death_hosp, exp_noise)
  ll_deaths_hosp_50_54 <-
    ll_nbinom(observed$deaths_hosp_50_54,
              pars$phi_death_hosp * model_deaths_hosp_50_54,
              pars$kappa_death_hosp, exp_noise)
  ll_deaths_hosp_55_59 <-
    ll_nbinom(observed$deaths_hosp_55_59,
              pars$phi_death_hosp * model_deaths_hosp_55_59,
              pars$kappa_death_hosp, exp_noise)
  ll_deaths_hosp_60_64 <-
    ll_nbinom(observed$deaths_hosp_60_64,
              pars$phi_death_hosp * model_deaths_hosp_60_64,
              pars$kappa_death_hosp, exp_noise)
  ll_deaths_hosp_65_69 <-
    ll_nbinom(observed$deaths_hosp_65_69,
              pars$phi_death_hosp * model_deaths_hosp_65_69,
              pars$kappa_death_hosp, exp_noise)
  ll_deaths_hosp_70_74 <-
    ll_nbinom(observed$deaths_hosp_70_74,
              pars$phi_death_hosp * model_deaths_hosp_70_74,
              pars$kappa_death_hosp, exp_noise)
  ll_deaths_hosp_75_79 <-
    ll_nbinom(observed$deaths_hosp_75_79,
              pars$phi_death_hosp * model_deaths_hosp_75_79,
              pars$kappa_death_hosp, exp_noise)
  ll_deaths_hosp_80_plus <-
    ll_nbinom(observed$deaths_hosp_80_plus,
              pars$phi_death_hosp * model_deaths_hosp_80_plus,
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
                           pars$phi_death_hosp * model_deaths_hosp_0_49 +
                           pars$phi_death_hosp * model_deaths_hosp_50_54 +
                           pars$phi_death_hosp * model_deaths_hosp_55_59 +
                           pars$phi_death_hosp * model_deaths_hosp_60_64 +
                           pars$phi_death_hosp * model_deaths_hosp_65_69 +
                           pars$phi_death_hosp * model_deaths_hosp_70_74 +
                           pars$phi_death_hosp * model_deaths_hosp_75_79 +
                           pars$phi_death_hosp * model_deaths_hosp_80_plus +
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

  ll_serology_1 <- ll_binom(observed$sero_pos_15_64_1,
                            observed$sero_tot_15_64_1,
                            model_sero_prob_pos_1)

  ll_serology_2 <- ll_binom(observed$sero_pos_15_64_2,
                            observed$sero_tot_15_64_2,
                            model_sero_prob_pos_2)

  ll_pillar2_tests <- ll_betabinom(observed$pillar2_pos,
                                   observed$pillar2_tot,
                                   model_pillar2_prob_pos,
                                   pars$rho_pillar2_tests)

  ll_pillar2_over25_tests <- ll_betabinom(observed$pillar2_over25_pos,
                                          observed$pillar2_over25_tot,
                                          model_pillar2_over25_prob_pos,
                                          pars$rho_pillar2_tests)

  ll_pillar2_under15_tests <- ll_betabinom(observed$pillar2_under15_pos,
                                           observed$pillar2_under15_tot,
                                           model_pillar2_under15_prob_pos,
                                           pars$rho_pillar2_tests)

  ll_pillar2_15_24_tests <- ll_betabinom(observed$pillar2_15_24_pos,
                                         observed$pillar2_15_24_tot,
                                         model_pillar2_15_24_prob_pos,
                                         pars$rho_pillar2_tests)

  ll_pillar2_25_49_tests <- ll_betabinom(observed$pillar2_25_49_pos,
                                         observed$pillar2_25_49_tot,
                                         model_pillar2_25_49_prob_pos,
                                         pars$rho_pillar2_tests)

  ll_pillar2_50_64_tests <- ll_betabinom(observed$pillar2_50_64_pos,
                                         observed$pillar2_50_64_tot,
                                         model_pillar2_50_64_prob_pos,
                                         pars$rho_pillar2_tests)

  ll_pillar2_65_79_tests <- ll_betabinom(observed$pillar2_65_79_pos,
                                         observed$pillar2_65_79_tot,
                                         model_pillar2_65_79_prob_pos,
                                         pars$rho_pillar2_tests)

  ll_pillar2_80_plus_tests <- ll_betabinom(observed$pillar2_80_plus_pos,
                                           observed$pillar2_80_plus_tot,
                                           model_pillar2_80_plus_prob_pos,
                                           pars$rho_pillar2_tests)

  ll_pillar2_under15_cases <- ll_nbinom(observed$pillar2_under15_cases,
                                        model_pillar2_under15_cases,
                                        pars$kappa_pillar2_cases, exp_noise)

  ll_pillar2_15_24_cases <- ll_nbinom(observed$pillar2_15_24_cases,
                                      model_pillar2_15_24_cases,
                                      pars$kappa_pillar2_cases, exp_noise)

  ll_pillar2_25_49_cases <- ll_nbinom(observed$pillar2_25_49_cases,
                                      model_pillar2_25_49_cases,
                                      pars$kappa_pillar2_cases, exp_noise)

  ll_pillar2_50_64_cases <- ll_nbinom(observed$pillar2_50_64_cases,
                                      model_pillar2_50_64_cases,
                                      pars$kappa_pillar2_cases, exp_noise)

  ll_pillar2_65_79_cases <- ll_nbinom(observed$pillar2_65_79_cases,
                                      model_pillar2_65_79_cases,
                                      pars$kappa_pillar2_cases, exp_noise)


  ll_pillar2_80_plus_cases <- ll_nbinom(observed$pillar2_80_plus_cases,
                                        model_pillar2_80_plus_cases,
                                        pars$kappa_pillar2_cases, exp_noise)

  ll_pillar2_over25_cases <- ll_nbinom(observed$pillar2_over25_cases,
                                       model_pillar2_over25_cases,
                                       pars$kappa_pillar2_cases, exp_noise)

  ll_pillar2_cases <- ll_nbinom(observed$pillar2_cases,
                                model_pillar2_cases,
                                pars$kappa_pillar2_cases, exp_noise)

  ll_react <- ll_binom(observed$react_pos,
                       observed$react_tot,
                       model_react_prob_pos)

  ll_strain <- ll_binom(observed$strain_non_variant,
                        observed$strain_tot,
                        model_strain_prob_pos)

  ll_strain_over25 <- ll_binom(observed$strain_over25_non_variant,
                               observed$strain_over25_tot,
                               model_strain_over25_prob_pos)

  ll_icu + ll_general + ll_hosp + ll_deaths_hosp + ll_deaths_carehomes +
    ll_deaths_comm + ll_deaths_non_hosp + ll_deaths + ll_deaths_hosp_0_49 +
    ll_deaths_hosp_50_54 + ll_deaths_hosp_55_59 + ll_deaths_hosp_60_64 +
    ll_deaths_hosp_65_69 + ll_deaths_hosp_70_74 + ll_deaths_hosp_75_79 +
    ll_deaths_hosp_80_plus + ll_admitted + ll_diagnoses + ll_all_admission +
    ll_serology_1 + ll_serology_2 + ll_pillar2_tests + ll_pillar2_cases +
    ll_pillar2_over25_tests + ll_pillar2_under15_tests +
    ll_pillar2_15_24_tests + ll_pillar2_25_49_tests + ll_pillar2_50_64_tests +
    ll_pillar2_65_79_tests + ll_pillar2_80_plus_tests +
    ll_pillar2_over25_cases + ll_pillar2_under15_cases +
    ll_pillar2_15_24_cases + ll_pillar2_25_49_cases + ll_pillar2_50_64_cases +
    ll_pillar2_65_79_cases + ll_pillar2_80_plus_cases + ll_react + ll_strain +
    ll_strain_over25
}


## We store within the severity parameters information on severity for
## carehome workers and residents. The vector ends up structured as
##
##   [1..n_age_groups, workers, residents]
##
## so we have length of n_groups = n_age_groups + 2
##
lancelot_severity <- function(p) {
  index_workers <- lancelot_index_workers()
  p_workers <- mean(p[index_workers])
  p_residents <- p[length(p)]
  c(p, p_workers, p_residents)
}


##' Lancelot severity parameters
##'
##' @title Lancelot severity parameters
##'
##' @param dt The step size
##'
##'
##' @param severity Severity data, used to determine default severity parameter
##'   age-scalings and to provide default severity parameter values. Can be
##'   `NULL` (use the default bundled data version in the package), or a
##'   [data.frame] object (for raw severity data).
##'
##' @section Time-varying parameters:
##' Every time varying parameter has the same format, which can be `NULL` (in
##'   which case the value from `severity` is used) or a list with `date` and
##'   `value` for changes in the parameter. If `value` is scalar then
##'   `date` can be `NULL` or missing. If `value` is a vector then `date`
##'   must be a vector of sircovid dates of the same length as `value`.
##'
##' @param p_C Time-varying parameters for p_C (the probability of an infected
##'   individual becoming symptomatic). See Details.
##'
##' @param p_H Time-varying parameters for p_H (the probability of a symptomatic
##'   individual requiring hospitalisation). See Details.
##'
##' @param p_H_CHR Time-varying parameters for p_H (the probability of a
##'   symptomatic individual requiring hospitalisation) for care home residents.
##'   If `NULL` then the value for the oldest age group is used. See Details.
##'
##' @param p_ICU Time-varying parameters for p_ICU (the probability of a
##'   hospitalised individual going to ICU). See Details.
##'
##' @param p_H_D Time-varying parameters for p_H_D (the probability of death in
##'   general beds).See Details.
##'
##' @param p_ICU_D Time-varying parameters for p_ICU_D (the probability of death
##'   in ICU). See Details.
##'
##' @param p_W_D Time-varying parameters for p_W_D (the probability of death in
##'   stepdown). See Details.
##'
##' @param p_G_D Time-varying parameters for p_G_D (the probability of
##'   individuals requiring hospitalisation dying in the community or a care
##'   home). See Details.
##'
##' @param p_G_D_CHR Time-varying parameters for p_G_D (the probability of
##'   individuals requiring hospitalisation dying in a care home) for care home
##'   residents. If `NULL` then the value for the oldest age group is used. See
##'   Details.
##'
##' @param p_R Time-varying parameters for p_R (the probability of an non-
##'   fatally infected individual having immunity post-infection). See Details.
##'
##'
##' @param p_star Time-varying parameters for p_star (the probability of
##'   patients being confirmed as covid on admission to hospital). See Details.
##'
##' @return A list of severity parameters
##'
##' @export
lancelot_parameters_severity <- function(dt,
                                         severity = NULL,
                                         p_C = NULL,
                                         p_H = NULL,
                                         p_H_CHR = NULL,
                                         p_ICU = NULL,
                                         p_H_D = NULL,
                                         p_ICU_D = NULL,
                                         p_W_D = NULL,
                                         p_G_D = NULL,
                                         p_G_D_CHR = NULL,
                                         p_R = NULL,
                                         p_star = NULL) {

  severity <- sircovid_parameters_severity(severity)
  severity <- lapply(severity, lancelot_severity)

  time_varying_severity <- list(C = p_C,
                                H = p_H,
                                ICU = p_ICU,
                                H_D = p_H_D,
                                ICU_D = p_ICU_D,
                                W_D = p_W_D,
                                G_D = p_G_D,
                                R = p_R,
                                star = p_star)

  time_varying_severity_CHR <- list(C = NULL,
                                    H = p_H_CHR,
                                    ICU = NULL,
                                    H_D = NULL,
                                    ICU_D = NULL,
                                    W_D = NULL,
                                    G_D = p_G_D_CHR,
                                    R = NULL,
                                    star = NULL)

  get_p_step <- function(x, name) {

    p_name <- paste0("p_", name)
    p <- x[[p_name]]
    time_vary <- time_varying_severity[[name]]
    time_vary_CHR <- time_varying_severity_CHR[[name]]

    if (is.null(time_vary)) {
      p_value <- NULL
      p_date <- NULL
    } else {
      p_value <- time_vary$value
      if ("date" %in% names(time_vary)) {
        p_date <- time_vary$date
      } else {
        p_date <- NULL
      }
    }

    if (is.null(time_vary_CHR)) {
      p_CHR_value <- NULL
      p_CHR_date <- NULL
    } else {
      p_CHR_value <- time_vary_CHR$value
      if ("date" %in% names(time_vary_CHR)) {
        p_CHR_date <- time_vary_CHR$date
      } else {
        p_CHR_date <- NULL
      }
    }

    if (!is.null(p_value)) {
      assert_proportion(p_value, p_name)
      if (length(p_value) == 1L) {
        if (length(p_date) != 0) {
          stop(sprintf(
            "As '%s' has a single 'value', expected NULL or missing 'date'",
            p_name))
        }
      } else if (length(p_date) != length(p_value)) {
        stop(sprintf("'date' and 'value' for '%s' must have the same length",
                     p_name))
      }
    }


    CHR <- FALSE
    if (!is.null(p_CHR_value)) {
      assert_proportion(p_CHR_value, paste0(p_name, "_CHR"))
      CHR <- TRUE
      p <- p[1:18]
      if (length(p_CHR_value) == 1L) {
        if (length(p_CHR_date) != 0) {
          stop(sprintf(
            "As '%s' has a single 'value', expected NULL or missing 'date'",
            paste0(p_name, "_CHR")))
        }
      } else if (length(p_CHR_date) != length(p_CHR_value)) {
        stop(sprintf("'date' and 'value' for '%s' must have the same length",
                     paste0(p_name, "_CHR")))
      }
    }

    if (all(p == 0)) {
      psi <- p
    } else {
      psi <- p / max(p)
    }

    if (is.null(p_value)) {
      p_step <- max(p)
    } else {
      p_step <- sircovid_parameters_piecewise_linear(p_date, p_value, dt)
    }

    p_step <- outer(p_step, psi)

    if (CHR) {
      p_CHR_step <- sircovid_parameters_piecewise_linear(p_CHR_date,
                                                         p_CHR_value, dt)
      n_p_steps <- dim(p_step)[1]
      n_p_CHR_steps <- length(p_CHR_step)[1]
      if (n_p_steps < n_p_CHR_steps) {
        p_step <- vapply(seq_len(18),
                         function(i) sircovid_parameters_expand_step(
                           seq_len(n_p_CHR_steps), p_step[, i]),
                         rep(0, n_p_CHR_steps))
      } else if (n_p_steps > n_p_CHR_steps) {
        p_CHR_step <- sircovid_parameters_expand_step(
          seq_len(n_p_steps),
          p_CHR_step)
      }
      p_step <- cbind(p_step, p_CHR_step, deparse.level = 0)
    }

    x[[paste0(p_name, "_step")]] <- p_step
    x[[paste0(p_name)]] <- NULL
    x[[paste0("n_", p_name, "_steps")]] <- dim(p_step)[1]

    x
  }

  for (name in names(time_varying_severity)) {
    severity <- get_p_step(severity, name)
  }

  severity
}


lancelot_index_workers <- function() {
  age_bins <- sircovid_age_bins()
  which(age_bins$start >= 25 & age_bins$start < 65)
}


lancelot_transmission_matrix <- function(eps, m_CHW, m_CHR, region,
                                         population) {
  index_workers <- lancelot_index_workers()
  m <- sircovid_transmission_matrix(region, population)
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


##' Create initial conditions for the lancelot model. This matches the
##' interface required for mcstate
##'
##' @title Initial conditions for the lancelot model
##'
##' @inheritParams basic_initial
##'
##' @return A numeric vector of initial conditions
##' @export
##' @examples
##' p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")
##' mod <- lancelot$new(p, 0, 10)
##' lancelot_initial(mod$info(), 10, p)
lancelot_initial <- function(info, n_particles, pars) {
  index <- info$index
  state <- numeric(info$len)

  index_S <- index[["S"]]
  index_S_no_vacc <- index_S[seq_len(length(pars$N_tot))]
  index_N_tot <- index[["N_tot"]]
  index_N_tot_sero_1 <- index[["N_tot_sero_1"]][[1L]]
  index_N_tot_sero_2 <- index[["N_tot_sero_2"]][[1L]]
  index_N_tot_PCR <- index[["N_tot_PCR"]][[1L]]

  index_prob_strain <- index[["prob_strain"]]

  seed_age_band <- 4L
  index_I_weighted <- index[["I_weighted"]][[1L]] + seed_age_band - 1L

  ## S0 is the population totals
  initial_S <- pars$N_tot

  state[index_S_no_vacc] <- initial_S
  state[index_N_tot] <- pars$N_tot
  state[index_N_tot_sero_1] <- sum(pars$N_tot)
  state[index_N_tot_sero_2] <- sum(pars$N_tot)
  state[index_N_tot_PCR] <- sum(pars$N_tot)
  state[index_prob_strain] <- c(1L, numeric(length(index_prob_strain) - 1L))
  state[index_I_weighted] <- 1

  state
}


lancelot_parameters_vaccination <- function(N_tot,
                                            dt,
                                            rel_susceptibility = 1,
                                            rel_p_sympt = 1,
                                            rel_p_hosp_if_sympt = 1,
                                            rel_p_death = 1,
                                            rel_infectivity = 1,
                                            vaccine_progression_rate = NULL,
                                            vaccine_schedule = NULL,
                                            vaccine_index_dose2 = NULL,
                                            vaccine_index_booster = NULL,
                                            n_strains = 1,
                                            vaccine_catchup_fraction = 1,
                                            n_doses = 2L) {
  n_groups <- lancelot_n_groups()
  stopifnot(length(N_tot) == n_groups)
  calc_n_vacc_classes <- function(x) {
    if (is.array(x)) nlayer(x) else length(x)
  }

  assert_proportion(rel_susceptibility)

  rel_params <- list(rel_susceptibility = rel_susceptibility,
                     rel_p_sympt = rel_p_sympt,
                     rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
                     rel_p_death = rel_p_death,
                     rel_infectivity = rel_infectivity)

  n <- vnapply(rel_params, calc_n_vacc_classes)

  if (any(n > 1) && length(unique(n[n > 1])) != 1) {
    msg1 <- paste(names(rel_params), collapse = ", ")
    msg2 <- "should have the same dimension"
    stop(paste(msg1, msg2))
  }
  n_vacc_classes <- max(n)

  ret <- Map(function(value, name)
    build_rel_param(value, n_strains, n_vacc_classes, name),
    rel_params, names(rel_params))

  if (is.null(vaccine_schedule)) {
    if (!is.null(vaccine_index_dose2) && vaccine_index_dose2 != 1L) {
      stop("'vaccine_index_dose2' set without schedule")
    }
    ret$vaccine_dose_step <- array(0, c(n_groups, n_doses, 1))
    ret$index_dose <- rep(1L, n_doses)
  } else {
    assert_is(vaccine_schedule, "vaccine_schedule")
    vaccine_index_dose2 <- vaccine_index_dose2 %||% 1L
    if (vaccine_index_dose2 > n_vacc_classes) {
      stop(sprintf(
        "Invalid value for 'vaccine_index_dose2', must be in [1, %d]",
        n_vacc_classes))
    }
    stopifnot(vaccine_schedule$n_doses == n_doses)

    if (is.null(vaccine_index_booster)) {
      if (n_doses != 2L) {
        stop("'n_doses' must be 2 as boosters not used")
      }
    } else {
      if (n_doses != 3L) {
        stop("'n_doses' must be 3 as boosters are used")
      }
      if (vaccine_index_booster > n_vacc_classes) {
        stop(sprintf(
          "Invalid value for 'vaccine_index_booster', must be in [1, %d]",
          n_vacc_classes))
      }
    }

    n_days <- dim(vaccine_schedule$doses)[[3]]
    i <- rep(seq_len(n_days), each = 1 / dt)
    len <- vaccine_schedule$date / dt
    ret$index_dose <- c(1L, vaccine_index_dose2, vaccine_index_booster)

    ret$vaccine_dose_step <- mcstate::array_bind(
      array(0, c(n_groups, n_doses, len)),
      (vaccine_schedule$doses * dt)[, , i])
  }

  ret$index_dose_inverse <- create_index_dose_inverse(n_vacc_classes,
                                                      ret$index_dose)

  ret$n_vacc_classes <- n_vacc_classes
  ret$vaccine_progression_rate_base <- build_vaccine_progression_rate(
    vaccine_progression_rate, n_vacc_classes, ret$index_dose)


  assert_scalar(vaccine_catchup_fraction)
  assert_proportion(vaccine_catchup_fraction)
  ret$vaccine_catchup_fraction <- vaccine_catchup_fraction

  ret$n_doses <- n_doses

  ret
}

lancelot_parameters_strain <- function(strain_transmission, strain_seed_date,
                                       strain_seed_size, strain_seed_pattern,
                                       dt) {
  if (length(strain_transmission) == 0) {
    stop("At least one value required for 'strain_transmission'")
  }
  if (length(strain_transmission) > 2) {
    stop(paste(
      "Only 1 or 2 strains valid ('strain_transmission' too long)'.",
      "See 'n_S_progress' in the odin code to fix this"))
  }

  assert_non_negative(strain_transmission)

  if (is.null(strain_seed_date)) {
    if (!is.null(strain_seed_size)) {
      stop(paste("As 'strain_seed_date' is NULL, expected 'strain_seed_size'",
                 "to be NULL"))
    }
    if (!is.null(strain_seed_pattern)) {
      stop(paste("As 'strain_seed_date' is NULL, expected",
                  "'strain_seed_pattern' to be NULL"))
    }
    strain_seed_step_start <- 0
    strain_seed_value <- 0
  } else {
    if (length(strain_transmission) == 1L) {
      stop("Can't use 'strain_seed_date' if only using one strain")
    }
    if (length(strain_seed_date) != 1L) {
      stop("'strain_seed_date' must be a single date")
    }
    if (length(strain_seed_size) != 1L) {
      stop("'strain_seed_size' must be a single value")
    }
    assert_sircovid_date(strain_seed_date)
    assert_non_negative(strain_seed_size)
    assert_positive(strain_seed_pattern)

    strain_seed_step <- strain_seed_date / dt

    strain_seed_value <- strain_seed_size *
      seed_over_steps(strain_seed_step, strain_seed_pattern)

    strain_seed_step_start <- floor(strain_seed_step)
  }

  if (length(strain_transmission) == 2) {
    strain_transmission <- mirror_strain(strain_transmission)
  }

  list(n_strains = length(strain_transmission),
       strain_transmission = strain_transmission,
       strain_seed_step_start = strain_seed_step_start,
       strain_seed_value = strain_seed_value)
}


lancelot_parameters_waning <- function(waning_rate) {
  waning_rate <- build_waning_rate(waning_rate)
  list(
    waning_rate = waning_rate
  )
}


##' Lancelot progression parameters.  The `s_` parameters are the
##' scaling parameters for the Erlang distibution (a.k.a 'k'), while
##' the `gamma_` parameters are the gamma parameters of that
##' distribution.  These need to be aligned with Bob's severity
##' outputs, and we will come up with a better way of coordinating the
##' two.
##'
##' @title Lancelot progression parameters
##'
##' @param dt The step size
##'
##' @section Time-varying parameters:
##' Every time varying parameter has the same format, which can be `NULL` (in
##'   which case a default single value is used) or a list with `date` and
##'   `value` for changes in the parameter. If `value` is scalar then
##'   `date` can be `NULL` or missing. If `value` is a vector then `date`
##'   must be a vector of sircovid dates of the same length as `value`.
##'
##' @param gamma_E Time-varying parameters for the Erlang rate parameter of the
##'   duration in the E (exposed) compartment. See Details.
##'
##' @param gamma_A Time-varying parameters for the Erlang rate parameter of the
##'   duration in the I_A (asymptomatic) compartment. See Details.
##'
##' @param gamma_P Time-varying parameters for the Erlang rate parameter of the
##'   duration in the I_P (presymptomatic) compartment. See Details.
##'
##' @param gamma_C_1 Time-varying parameters for the Erlang rate parameter of
##'   the duration in the I_C_1 (first stage symptomatic) compartment. See
##'   Details.
##'
##' @param gamma_C_2 Time-varying parameters for the Erlang rate parameter of
##'   the duration in the I_C_2 (second stage symptomatic) compartment. See
##'   Details.
##'
##' @param gamma_ICU_pre Time-varying parameters for the Erlang rate parameter
##'   of the duration in the ICU_pre (general beds stay before ICU) compartment.
##'   See Details.
##'
##' @param gamma_ICU_D Time-varying parameters for the Erlang rate parameter of
##'   the duration in the ICU_D (ICU stay for individuals who die in ICU)
##'   compartment. See Details.
##'
##' @param gamma_ICU_W_D Time-varying parameters for the Erlang rate parameter
##'   of the duration in the ICU_W_D (ICU stay for individuals who go on to die
##'   in stepdown) compartment. See Details.
##'
##' @param gamma_ICU_W_R Time-varying parameters for the Erlang rate parameter
##'   of the duration in the ICU_W_R (ICU stay for individuals who go on to
##'   recover in stepdown) compartment. See Details.
##'
##' @param gamma_H_D Time-varying parameters for the Erlang rate parameter of
##'   the duration in the H_D (general beds stay for individuals who die in
##'   general beds) compartment. See Details.
##'
##' @param gamma_H_R Time-varying parameters for the Erlang rate parameter of
##'   the duration in the H_R (general beds stay for individuals who recover in
##'   general beds) compartment. See Details.
##'
##' @param gamma_W_D Time-varying parameters for the Erlang rate parameter of
##'   the duration in the W_D (stepdown for individuals who die in stepdown)
##'   compartment. See Details.
##'
##' @param gamma_W_R Time-varying parameters for the Erlang rate parameter of
##'   the duration in the W_R (stepdown for individuals who recover in stepdown)
##'   compartment. See Details.
##'
##' @param gamma_G_D Time-varying parameters for the Erlang rate parameter of
##'   the duration in the G_D (delay to death in the community/care homes)
##'   compartment. See Details.
##'
##' @param gamma_U Time-varying parameters for the Erlang rate parameter of
##'   the duration of the delay from hospital admission to diagnosis for those
##'   not confirmed on admission. Note this duration is a single-stage Erlang.
##'   See Details.
##'
##' @return A list of parameter values
##'
##' @export
lancelot_parameters_progression <- function(dt,
                                            gamma_E = NULL,
                                            gamma_A = NULL,
                                            gamma_P = NULL,
                                            gamma_C_1 = NULL,
                                            gamma_C_2 = NULL,
                                            gamma_ICU_pre = NULL,
                                            gamma_ICU_D = NULL,
                                            gamma_ICU_W_D = NULL,
                                            gamma_ICU_W_R = NULL,
                                            gamma_H_D = NULL,
                                            gamma_H_R = NULL,
                                            gamma_W_D = NULL,
                                            gamma_W_R = NULL,
                                            gamma_G_D = NULL,
                                            gamma_U = NULL) {

  ## The k_ parameters are the shape parameters for the Erlang
  ## distribution, while the gamma parameters are the rate
  ## parameters of that distribution.
  ret <- list(k_E = 2,
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
              k_sero_pre_1 = 2,
              k_sero_pos_1 = 2,
              k_sero_pre_2 = 2,
              k_sero_pos_2 = 2,
              k_PCR_pre = 2,
              k_PCR_pos = 2,

              gamma_E = 1 / (3.42 / 2),
              gamma_A = 1 / 2.88,
              gamma_P = 1 / 1.68,
              gamma_C_1 = 1 / 2.14,
              gamma_C_2 = 1 / 1.86,
              gamma_G_D = 1 / (3 / 2),
              gamma_H_D = 2 / 5,
              gamma_H_R = 2 / 10,
              gamma_ICU_D = 2 / 5,
              gamma_ICU_W_R = 2 / 10,
              gamma_ICU_W_D = 2 / 10,
              gamma_ICU_pre = 2,
              gamma_W_R = 2 / 5,
              gamma_W_D = 2 / 5,
              gamma_sero_pre_1 = 1 / 10,
              gamma_sero_pos_1 = 1 / 25,
              gamma_sero_pre_2 = 1 / 10,
              gamma_sero_pos_2 = 1 / 25,
              gamma_U = 3 / 10,
              gamma_PCR_pre = 2 / 3,
              gamma_PCR_pos = 1 / 5
  )

  time_varying_gammas <- list(E = gamma_E,
                              A = gamma_A,
                              P = gamma_P,
                              C_1 = gamma_C_1,
                              C_2 = gamma_C_2,
                              ICU_pre = gamma_ICU_pre,
                              ICU_D = gamma_ICU_D,
                              ICU_W_D = gamma_ICU_W_D,
                              ICU_W_R = gamma_ICU_W_R,
                              H_D = gamma_H_D,
                              H_R = gamma_H_R,
                              W_D = gamma_W_D,
                              W_R = gamma_W_R,
                              G_D = gamma_G_D,
                              U = gamma_U)

  get_gamma_step <- function(x, name) {

    gamma_name <- paste0("gamma_", name)
    gamma <- x[[gamma_name]]
    time_vary <- time_varying_gammas[[name]]

    if (is.null(time_vary)) {
      gamma_value <- NULL
      gamma_date <- NULL
    } else {
      gamma_value <- time_vary$value
      if ("date" %in% names(time_vary)) {
        gamma_date <- time_vary$date
      } else {
        gamma_date <- NULL
      }
    }

    if (!is.null(gamma_value)) {
      assert_non_negative(gamma_value, gamma_name)
      if (length(gamma_value) == 1L) {
        if (length(gamma_date) != 0) {
          stop(sprintf(
            "As '%s' has a single 'value', expected NULL or missing 'date'",
            gamma_name))
        }
      } else if (length(gamma_date) != length(gamma_value)) {
        stop(sprintf("'date' and 'value' for '%s' must have the same length",
                     gamma_name))
      }
    }

    if (is.null(gamma_value)) {
      gamma_step <- gamma
    } else {
      gamma_step <- sircovid_parameters_piecewise_linear(gamma_date,
                                                         gamma_value, dt)
    }

    x[[paste0("gamma_", name, "_step")]] <- gamma_step
    x[[paste0("gamma_", name)]] <- NULL
    x[[paste0("n_gamma_", name, "_steps")]] <- length(gamma_step)

    x
  }

  for (name in names(time_varying_gammas)) {
    ret <- get_gamma_step(ret, name)
  }

  ret

}


##' Lancelot observation parameters
##'
##' @title Lancelot observation parameters
##'
##' @return A list of parameter values
##'
##' @param exp_noise Rate parameter for the exponentially-distributed noise used
##'   in the likelihood calculation. Typically a large value so noise is small.
##'
##' @export
lancelot_parameters_observation <- function(exp_noise = 1e6) {
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
    kappa_pillar2_cases = 2,
    phi_pillar2_cases_under15 = 1,
    phi_pillar2_cases_weekend_under15 = 0.8,
    phi_pillar2_cases_15_24 = 1,
    phi_pillar2_cases_weekend_15_24 = 0.8,
    phi_pillar2_cases_25_49 = 1,
    phi_pillar2_cases_weekend_25_49 = 0.8,
    phi_pillar2_cases_50_64 = 1,
    phi_pillar2_cases_weekend_50_64 = 0.8,
    phi_pillar2_cases_65_79 = 1,
    phi_pillar2_cases_weekend_65_79 = 0.8,
    phi_pillar2_cases_80_plus = 1,
    phi_pillar2_cases_weekend_80_plus = 0.8,
    ##
    rho_pillar2_tests = 0.1,
    p_NC_under15 = 0.01,
    p_NC_weekend_under15 = 0.005,
    p_NC_15_24 = 0.01,
    p_NC_weekend_15_24 = 0.005,
    p_NC_25_49 = 0.01,
    p_NC_weekend_25_49 = 0.005,
    p_NC_50_64 = 0.01,
    p_NC_weekend_50_64 = 0.005,
    p_NC_65_79 = 0.01,
    p_NC_weekend_65_79 = 0.005,
    p_NC_80_plus = 0.01,
    p_NC_weekend_80_plus = 0.005,
    ##
    ## rate for exponential noise, generally something big so noise is
    ## small (but non-zero))
    exp_noise = exp_noise)
}


##' Lancelot observation parameters
##'
##' @title Lancelot sensitivity and specificity parameters
##'
##' @return A list of parameter values
##'
##' @param sero_specificity_1 Specificity of the first serology test assay
##'
##' @param sero_sensitivity_1 Sensitivity of the first serology test assay
##'
##' @param sero_specificity_2 Specificity of the second serology test assay
##'
##' @param sero_sensitivity_2 Sensitivity of the second serology test assay
##'
##' @param pillar2_specificity Specificity of the Pillar 2 test
##'
##' @param pillar2_sensitivity Sensitivity of the Pillar 2 test
##'
##' @param react_specificity Specificity of the REACT test
##'
##' @param react_sensitivity Sensitivity of the REACT test
##'
##' @export
lancelot_parameters_sens_and_spec <- function(sero_specificity_1 = 0.9,
                                              sero_sensitivity_1 = 0.99,
                                              sero_specificity_2 = 0.9,
                                              sero_sensitivity_2 = 0.99,
                                              pillar2_specificity = 0.99,
                                              pillar2_sensitivity = 0.99,
                                              react_specificity = 0.99,
                                              react_sensitivity = 0.99) {
  ret <- list(
    ## Specificity and sensitivity for serology tests
    sero_specificity_1 = sero_specificity_1,
    sero_sensitivity_1 = sero_sensitivity_1,
    sero_specificity_2 = sero_specificity_2,
    sero_sensitivity_2 = sero_sensitivity_2,
    ## Specificity and sensitivity for Pillar 2 testing
    pillar2_specificity = pillar2_specificity,
    pillar2_sensitivity = pillar2_sensitivity,
    ## Specificity and sensitivity for REACT testing
    react_specificity = react_specificity,
    react_sensitivity = react_sensitivity)

  lapply(seq_len(length(ret)),
         function(i) assert_proportion(ret[i], names(ret)[i]))

  ret
}


lancelot_population <- function(population, carehome_workers,
                                carehome_residents) {
  ## Our core S0 calculation is more complicated than the basic model
  ## because we have to add the carehome workers and residents, *and*
  ## remove them from the core population. This extracts carehome
  ## residents from the older groups of the population, weighted
  ## towards the oldest, and extracts carehome workers from most
  ## working ages, extracting equal amounts from those age groups.
  N_tot <- c(population, carehome_workers, carehome_residents)

  index_workers <- lancelot_index_workers()
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


##' Check that data for the particle filter has required columns and
##' constraints among columns.  This function does not alter the data
##' at present (though it does return it invisibly).
##'
##' @title Check data for particle filter
##' @param data A data.frame of data
##' @return Invisibly, a data frame, identical to `data`
##' @export
lancelot_check_data <- function(data) {
  ## Column names apply at the level of the whole data set, and would
  ## be confusing if we reported it for each region together
  required <- c("icu", "general", "hosp", "deaths_hosp", "deaths_carehomes",
                "deaths_comm", "deaths_non_hosp", "deaths", "admitted",
                "diagnoses", "all_admission", "sero_pos_15_64_1",
                "sero_tot_15_64_1",  "sero_pos_15_64_2", "sero_tot_15_64_2",
                "pillar2_pos", "pillar2_tot", "pillar2_cases",
                "pillar2_over25_pos", "pillar2_over25_tot",
                "pillar2_over25_cases", "react_pos", "react_tot",
                "strain_non_variant", "strain_tot", "strain_over25_non_variant",
                "strain_over25_tot",
                "pillar2_under15_cases", "pillar2_15_24_cases",
                "pillar2_25_49_cases", "pillar2_50_64_cases",
                "pillar2_65_79_cases", "pillar2_80_plus_cases",
                "pillar2_under15_tot", "pillar2_15_24_tot", "pillar2_25_49_tot",
                "pillar2_50_64_tot", "pillar2_65_79_tot", "pillar2_80_plus_tot",
                "pillar2_under15_pos", "pillar2_15_24_pos", "pillar2_25_49_pos",
                "pillar2_50_64_pos", "pillar2_65_79_pos", "pillar2_80_plus_pos",
                "deaths_hosp_0_49", "deaths_hosp_50_54", "deaths_hosp_55_59",
                "deaths_hosp_60_64", "deaths_hosp_65_69", "deaths_hosp_70_74",
                "deaths_hosp_75_79", "deaths_hosp_80_plus")
  assert_is(data, "data.frame")
  verify_names(data, required, allow_extra = TRUE)

  if (length(unique(data$region)) > 1) {
    res <- lapply(split(as.data.frame(data), data$region), function(d)
      tryCatch(lancelot_check_data(d), error = identity))
    err <- vlapply(res, inherits, "error")
    if (any(err)) {
      msg <- sprintf("  %s: %s",
                     names(which(err)), vcapply(res[err], "[[", "message"))
      stop("Validation failures:\n", paste(msg, collapse = "\n"))
    }
    return(invisible(data))
  }

  has <- as.data.frame(!is.na(data[required]))
  has_any <- function(nms) {
    apply(has[nms], 1, any)
  }

  nms_deaths_split <- c("deaths_comm", "deaths_hosp", "deaths_carehomes",
                        "deaths_non_hosp")
  err_deaths <- has$deaths & has_any(nms_deaths_split)
  if (any(has$deaths & has_any(nms_deaths_split))) {
    stop(paste("Deaths are not consistently split into total vs",
               "hospital/non-hospital or hospital/care homes/community"))
  }

  nms_deaths_non_hosp <- c("deaths_comm", "deaths_carehomes")
  if (any(has$deaths_non_hosp & has_any(nms_deaths_non_hosp))) {
    stop(paste("Non-hospital deaths are not consistently split into total vs",
               "care homes/community"))
  }

  deaths_ages <- c("0_49", "50_54", "55_59", "60_64", "65_69", "70_74",
                   "75_79", "80_plus")
  P2_over25_ages <- c("_25_49", "_50_64", "_65_79", "_80_plus")
  P2_all_ages <- c("_under15", "_15_24", "_over25", P2_over25_ages)
  P2_all <- c("", P2_all_ages)

  ## NOTE: I (RGF) think that nms_deaths_aggr might really be deaths
  ## and nms_deaths_split more completely, but practically I expect
  ## this is equivalent.
  ##
  ## NOTE: This is asserted at the level of the _entire_ time series,
  ## not per day; that might be overly strict?
  nms_deaths_ages <- paste0("deaths_hosp_", deaths_ages)
  nms_deaths_aggr <- c("deaths", "deaths_hosp")
  err_deaths <- any(has_any(nms_deaths_ages)) && any(has_any(nms_deaths_aggr))
  if (err_deaths) {
    ## Perhaps it's just me, but I am not sure that would find this
    ## actionable; we probably need to indicate the appropriate
    ## columns.
    stop(paste("Cannot fit to all ages aggregated for deaths if fitting",
               "to any sub-groups"))
  }

  nms_pillar2_pos <- sprintf("pillar2%s_pos", P2_all)
  nms_pillar2_tot <- sprintf("pillar2%s_tot", P2_all)
  nms_pillar2_cases <- sprintf("pillar2%s_cases", P2_all)
  err_pillar2 <- any(has_any(nms_pillar2_cases)) &&
    any(has_any(c(nms_pillar2_pos, nms_pillar2_tot)))
  if (err_pillar2) {
    stop("Cannot fit to pillar 2 cases and positivity together")
  }

  err_strain_stream <-
    any(has_any(c("strain_non_variant", "strain_tot"))) &&
    any(has_any(c("strain_over25_non_variant", "strain_over25_tot")))
  if (err_strain_stream) {
    stop("Cannot fit to more than one strain data stream")
  }

  nms_pillar2_over25_pos <- sprintf("pillar2%s_pos", P2_over25_ages)
  nms_pillar2_over25_tot <- sprintf("pillar2%s_tot", P2_over25_ages)
  nms_pillar2_over25_cases <- sprintf("pillar2%s_cases", P2_over25_ages)
  nms_pillar2_over25_aggregated <- c("pillar2_over25_pos",
                                     "pillar2_over25_tot",
                                     "pillar2_over25_cases")
  err_pillar2_over25 <-
    any(has_any(nms_pillar2_over25_aggregated)) &&
    any(has_any(c(nms_pillar2_over25_pos,
                  nms_pillar2_over25_tot,
                  nms_pillar2_over25_cases)))
  if (err_pillar2_over25) {
    stop(paste("Cannot fit to over 25s for pillar 2 if fitting to any over 25",
               "sub-groups"))
  }

  nms_pillar2_all_ages_pos <- sprintf("pillar2%s_pos", P2_all_ages)
  nms_pillar2_all_ages_tot <- sprintf("pillar2%s_tot", P2_all_ages)
  nms_pillar2_all_ages_cases <- sprintf("pillar2%s_cases", P2_all_ages)
  nms_pillar2_all_ages_aggregated <- c("pillar2_pos",
                                     "pillar2_tot",
                                     "pillar2_cases")
  err_pillar2_all_ages <-
    any(has_any(nms_pillar2_all_ages_aggregated)) &&
    any(has_any(c(nms_pillar2_all_ages_pos,
                  nms_pillar2_all_ages_tot,
                  nms_pillar2_all_ages_cases)))
  if (err_pillar2_all_ages) {
    stop(paste("Cannot fit to all ages aggregated for pillar 2",
               "if fitting to any sub-groups"))
  }

  invisible(data)
}


lancelot_n_groups <- function() {
  length(sircovid_age_bins()$start) + 2L
}


##' Checks severity probabilities are valid
##'
##' @title Check severity parameters
##'
##' @param pars A list of probabilities from [lancelot_parameters].
##'
##' @return Returns `pars` invisibly if all severity probabilities lie in
##   (0, 1), otherwise errors.
##'
##' @export
lancelot_check_severity <- function(pars) {
  check_parameters <- function(p_step, rel_p) {
    vapply(
      seq_len(pars$n_groups),
      function(i) {
        if (is.null(pars[[rel_p]])) {
          stop(sprintf("Parameter '%s' is missing", rel_p))
        } else if (is.null(pars[[p_step]])) {
          stop(sprintf("Parameter '%s' is missing", p_step))
        }

        ## simple but crude checks
        min_p <- min(pars[[rel_p]][i, , ]) * min(pars[[p_step]][, i])

        if (any(min_p < 0)) {
          stop(sprintf("%s * %s is not in [0, 1] for group %d",
                       rel_p, p_step, i))
        } else {
          if (!(ncol(pars[[rel_p]]) %in% c(1, 4))) {
            stop(sprintf("%s should have 1 or 4 columns", rel_p))
          }
        }
        TRUE
      }, logical(1)
    )
  }

  step_pars <- c("p_C_step", "p_H_step", "p_ICU_step", "p_ICU_D_step",
                 "p_H_D_step", "p_W_D_step", "p_G_D_step", "p_R_step")
  rel_pars <- c("rel_p_sympt", "rel_p_hosp_if_sympt", "rel_p_ICU",
                "rel_p_ICU_D", "rel_p_H_D", "rel_p_W_D", "rel_p_G_D", "rel_p_R")

  Map(check_parameters,
      p_step = step_pars,
      rel_p = rel_pars
  )



  invisible(pars)
}


create_index_dose_inverse <- function(n_vacc_classes, index_dose) {
  index_dose_inverse <- integer(n_vacc_classes)
  index_dose_inverse[index_dose] <- seq_along(index_dose)
  index_dose_inverse
}


calculate_index <- function(index, suffix_list, suffix0 = NULL,
                            state_name = state) {
  if (is.null(suffix0)) {
    suffixes <- list()
  } else {
    suffixes <- list(suffix0)
  }
  for (i in seq_along(suffix_list)) {
    nm <- names(suffix_list)[[i]]
    if (length(nm) == 0) {
      nm <- ""
    }
    suffixes <- c(suffixes,
                  list(c("", sprintf("_%s%s", nm,
                                     seq_len(suffix_list[[i]] - 1L)))))
  }
  suffixes <- expand.grid(suffixes)
  nms <- apply(suffixes, 1,
               function(x) sprintf("%s%s",
                                   state_name, paste0(x, collapse = "")))
  set_names(index, nms)
}
