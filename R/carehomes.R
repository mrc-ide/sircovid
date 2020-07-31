## To add: C_1, C_2 (defaults 1e-7)
carehomes_parameters <- function(start_date, region,
                                 beta_date = NULL, beta_value = NULL,
                                 severity_path = NULL) {
  ret <- sircovid_parameters_shared(start_date, region,
                                    beta_date, beta_value)

  ## These should be flexible and will be set in the pmcmc so will
  ## move up to the argument list of this function; these are used
  ## only here in the setup and are not used in the model itself.
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

  severity <- carehomes_parameters_severity(severity_path, ret$population,
                                            p_death_carehome)

  ret$m <- carehomes_transmission_matrix(eps, C_1, C_2, ret$population)

  ## Our core S0 calculation is more complicated than the basic model
  ## because we have to add the carehome workers and residents, *and*
  ## remove them from the core population. This extracts carehome
  ## residents from the older groups of the population, weighted
  ## towards the oldest, and extracts carehome workers from most
  ## working ages, evenly across the population.
  N_tot <- c(ret$population, carehome_workers, carehome_residents)

  index_workers <- carehomes_index_workers()
  weights_workers <- N_tot[index_workers] / sum(N_tot[index_workers])
  index_residents <- which(sircovid_age_bins()$start >= 65)
  weights_residents <- c(0.05, 0.05, 0.15, 0.75)

  N_tot[index_residents] <-
    round(N_tot[index_residents] - carehome_residents * weights_residents)
  N_tot[index_workers] <-
    round(N_tot[index_workers] - carehome_workers * weights_workers)

  ## This is used to normalise the serology counts (converting them
  ## from number of positive/negative tests into a fraction). This is
  ## constant over the simulation, being the total population size of
  ## 15 to 64 year olds.
  N_tot_15_64 <- sum(N_tot[4:13])

  if (any(N_tot[index_residents] < 0)) {
    stop("Not enough population to meet care home occupancy")
  }
  if (any(N_tot[index_workers] < 0)) {
    stop("Not enough population to be care workers")
  }

  ret$N_tot <- N_tot

  ## Adding this here, but better would be to pass N_age as-is, then
  ## update the leading dimension to something more accurate (e.g.,
  ## N_groups, setting this as N_groups <- N_age + 2)
  ret$N_age <- ret$N_age + 2L

  c(ret,
    severity,
    carehomes_parameters_progression())
}


carehomes_index <- function(info) {
  len <- vnapply(info, prod)
  start <- cumsum(len) - len + 1L
  list(run = c(start[["I_ICU_tot"]],
               start[["D_comm_tot"]],
               start[["D_hosp_tot"]],
               start[["D_tot"]],
               start[["cum_admit_conf"]],
               start[["cum_new_conf"]],
               start[["R_pre_15_64"]],
               start[["R_neg_15_64"]],
               start[["R_pos_15_64"]]))
}


## We store within the severity parameters information on severity for
## carehome workers and residents. The vector ends up structured as
##
##   [1..N_age, workers, residents]
##
## so we have length of N_age + 2
carehomes_severity <- function(p, population) {
  index_workers <- carehomes_index_workers()
  p_workers <- weighted.mean(p[index_workers], population[index_workers])
  p_residents <- p[length(p)]
  c(p, p_workers, p_residents)
}


carehomes_parameters_severity <- function(severity_path, population,
                                          p_death_carehome) {
  severity <- sircovid_parameters_severity(severity_path)
  severity <- lapply(severity, carehomes_severity, population)
  severity$p_death_comm[length(severity$p_death_comm)] <- p_death_carehome
  severity
}


carehomes_index_workers <- function() {
  age_bins <- sircovid_age_bins()
  which(age_bins$start >= 25 & age_bins$start < 65)
}


carehomes_transmission_matrix <- function(eps, C_1, C_2, population) {
  index_workers <- carehomes_index_workers()
  m <- sircovid_transmission_matrix()
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


carehomes_initial <- function(info, n_particles, pars) {
  ## TODO: this will simplify once we get the index here, see
  ## odin.dust issue #24
  len <- vnapply(info, prod)
  start <- cumsum(len) - len + 1L
  state <- numeric(sum(len))

  ## This corresponds to the 15-19y age bracket for compatibility with
  ## our first version, will be replaced by better seeding model, but
  ## probably has limited impact.
  seed_age_band <- 3L
  index_I <- start[["I_asympt"]] + seed_age_band
  index_R_pre <- start[["R_pre"]] + seed_age_band
  index_PCR_pos <- start[["PCR_pos"]] + seed_age_band
  index_S <- seq.int(start[["S"]], length.out = len[["S"]])
  index_N_tot <- seq.int(start[["N_tot"]], length.out = len[["N_tot"]])
  index_N_tot2 <- start[["N_tot2"]]

  ## Always start with 10, again for compatibility
  initial_I <- 10

  ## S0 is the population totals, minus the seeded infected
  ## individuals
  initial_S <- pars$N_tot
  initial_S[seed_age_band] <- initial_S[seed_age_band] - initial_I

  state[index_S] <- initial_S
  state[index_I] <- initial_I
  state[index_R_pre] <- initial_I
  state[index_PCR_pos] <- initial_I
  state[index_N_tot] <- pars$N_tot
  state[index_N_tot2] <- sum(pars$population)

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
       s_hosp_D = 2 ,
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

  ## TODO: cache this file read as it's constant within a session
  data <- read_csv(sircovid_file("extdata/carehomes.csv"))
  i <- match(tolower(region), data$region)
  if (is.na(i)) {
    valid <- paste(squote(data$region), collapse = ", ")
    stop(sprintf("Carehome beds not found for '%s': must be one of %s",
                 region, valid))
  }

  data$carehome_beds[[i]]
}
