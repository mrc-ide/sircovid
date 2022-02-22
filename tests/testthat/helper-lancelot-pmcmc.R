reference_data_lancelot_mcmc <- function() {
  load_reference("data/lancelot_pmcmc.rds", {
    start_date <- sircovid_date("2020-02-02")
    pars <- lancelot_parameters(
      start_date, "england",
      beta_date = sircovid_date(c("2020-03-10", "2020-03-20")),
      beta_value = c(0.08, 0.04))
    data <- helper_sircovid_data(read_csv(sircovid_file("extdata/example.csv")),
                                 0, pars$dt)
    v <- c("deaths_comm", "deaths", "deaths_carehomes", "deaths_non_hosp",
           "deaths_hosp_0_49", "deaths_hosp_50_54", "deaths_hosp_55_59",
           "deaths_hosp_60_64", "deaths_hosp_65_69", "deaths_hosp_70_74",
           "deaths_hosp_75_79", "deaths_hosp_80_plus",
           "general", "hosp", "admitted", "diagnoses", "all_admission",
           "all_admission_0_9", "all_admission_10_19", "all_admission_20_29",
           "all_admission_30_39", "all_admission_40_49", "all_admission_50_59",
           "all_admission_60_69", "all_admission_70_79",
           "all_admission_80_plus", "sero_pos_15_64_1", "sero_tot_15_64_1",
           "sero_pos_15_64_2", "sero_tot_15_64_2", "pillar2_pos", "pillar2_tot",
           "pillar2_cases", "pillar2_over25_pos", "pillar2_over25_tot",
           "pillar2_over25_cases", "pillar2_under15_pos", "pillar2_under15_tot",
           "pillar2_under15_cases", "pillar2_15_24_pos", "pillar2_15_24_tot",
           "pillar2_15_24_cases", "pillar2_25_49_pos", "pillar2_25_49_tot",
           "pillar2_25_49_cases", "pillar2_50_64_pos", "pillar2_50_64_tot",
           "pillar2_50_64_cases", "pillar2_65_79_pos", "pillar2_65_79_tot",
           "pillar2_65_79_cases", "pillar2_80_plus_pos", "pillar2_80_plus_tot",
           "pillar2_80_plus_cases", "react_pos", "react_tot",
           "strain_non_variant", "strain_tot",
           "strain_over25_non_variant", "strain_over25_tot")
    for (i in v) {
      data[[i]] <- NA
    }
    data$deaths_hosp <- data$deaths

    filter <- helper_lancelot_particle_filter(data, 10)

    ## Completely dummy mcmc parameters object
    pars_mcmc <- mcstate::pmcmc_parameters$new(
      list(
        mcstate::pmcmc_parameter("a", 1),
        mcstate::pmcmc_parameter("b", 1)),
      diag(2),
      function(p) pars)

    control <- mcstate::pmcmc_control(10, n_chains = 1,
                                      save_trajectories = TRUE,
                                      progress = FALSE)
    mcstate::pmcmc(pars_mcmc, filter, control = control)
  })
}


## This is the remains of lancelot_forecast, since removed from the
## package; it is used in some integration tests elsewhere in the
## package, but the objects these would be used on in practice are
## build by spimalot.
reference_data_lancelot_trajectories <- function() {
  load_reference("data/lancelot_trajectories.rds", {
    dat <- reference_data_lancelot_mcmc()
    incidence <- "deaths"
    ret <- mcstate::pmcmc_thin(dat)
    steps_predict <- seq(ret$predict$step, length.out = 11,
                         by = ret$predict$rate)
    ret$trajectories <- mcstate::pmcmc_predict(
      ret, steps_predict, prepend_trajectories = TRUE)
    ret$trajectories$date <- ret$trajectories$step / ret$trajectories$rate
    ret$trajectories <- add_trajectory_incidence(ret$trajectories, "deaths")
    ret
  })
}


load_reference <- function(path, code) {
  version <- packageVersion("sircovid")
  if (file.exists(path)) {
    prev <- readRDS(path)
    if (identical(prev$version, version)) {
      return(prev$value)
    }
    message(sprintf("Reference data '%s' is out of date - regenerating",
                    path))
  } else {
    message(sprintf("Reference data '%s' does not exist - generating",
                    path))
  }

  value <- force(code)
  dir.create(dirname(path), FALSE, TRUE)
  ## Ignore a warning about package:sircovid not being available
  suppressWarnings(
    saveRDS(list(version = version, value = value), path))
  value
}
