reference_data_mcmc <- function() {
  path <- "data/pmcmc.rds"
  if (!file.exists(path)) {
    message(sprintf("Reference data '%s' does not exist - generating",
                    path))

    start_date <- sircovid_date("2020-02-02")
    pars <- carehomes_parameters(start_date, "england")
    data <- sircovid_data(read_csv(sircovid_file("extdata/example.csv")),
                          start_date, pars$dt)
    v <- c("deaths_comm", "deaths", "general", "hosp", "admitted",
           "new", "new_admitted", "npos_15_64", "ntot_15_64",
           "pillar2_pos", "pillar2_tot", "pillar2_cases",
           "pillar2_over25_pos", "pillar2_over25_tot",
           "pillar2_over25_cases", "react_pos", "react_tot")
    for (i in v) {
      data[[i]] <- NA
    }
    data$deaths_hosp <- data$deaths

    filter <- carehomes_particle_filter(data, 10)

    ## Completely dummy mcmc parameters object
    pars_mcmc <- mcstate::pmcmc_parameters$new(
      list(
        mcstate::pmcmc_parameter("a", 1),
        mcstate::pmcmc_parameter("b", 1)),
      diag(2),
      function(p) pars)

    samples <- mcstate::pmcmc(pars_mcmc, filter, n_steps = 10,
                              n_chains = 1,
                              save_trajectories = TRUE, progress = FALSE)

    dir.create(dirname(path), FALSE, TRUE)
    ## Ignore a warning about package:sircovid not being available
    suppressWarnings(saveRDS(samples, path))
  }
  readRDS(path)
}


reference_data_trajectories <- function() {
  path <- "data/trajectories.rds"
  if (!file.exists(path)) {
    message(sprintf("Reference data '%s' does not exist - generating",
                    path))

    dat <- reference_data_mcmc()
    incidence <- c("deaths", "deaths_hosp", "infections")
    res <- carehomes_forecast(dat, 3, 5, 10, incidence, TRUE)

    dir.create(dirname(path), FALSE, TRUE)
    ## Ignore a warning about package:sircovid not being available
    suppressWarnings(saveRDS(res, path))
  }
  readRDS(path)
}
