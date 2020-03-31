context("sampler")

## This is not a real test but simply tries to run the model
test_that("sampler runs without error", {
  set.seed(1)
  time_steps_per_day <- 4

  pars_model <- generate_parameters(
    transmission_model = "POLYMOD",
    beta = rep(0.125, 3),
    age_limits = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80),
    infection_seeding = c(0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    #severity_data_file = "extdata/severity.csv", ## the format of this file is rather different
    progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 1, hosp = 2, ICU = 2, rec = 2),
    gammas = list(E = 1/2.5, asympt = 1/2.09, mild = 1/2.09, ILI = 1/4, hosp = 2/1, ICU = 2/5, rec = 2/5),
    hosp_transmission = 0,
    ICU_transmission = 0,
    dt = 1/time_steps_per_day
  )
  
  # Manually set this parameter, to match previous results
  ## Carried over from the initial NHS meeting
  prop_symp_seek_HC <- c(0.3570377550, 0.3570377550, 0.3712946230,0.3712946230,	0.420792849,0.420792849,
                         0.459552523,0.459552523,	0.488704572,0.488704572,	0.578769171,0.578769171,	0.65754772,0.65754772,	0.73278164,0.73278164,0.76501082)
  #Proportion seeking healthcare
  pars_model$p_sympt_ILI <- rep(0.66,length(prop_symp_seek_HC))*prop_symp_seek_HC

  pars_obs <- list(
    ## what should this be?
    phi_ICU = 0.95,
    ## what should this be?
    k_ICU = 2,
    ## current proportion of England deaths over UK deaths
    phi_death = 926/1019,
    ## what should this be?
    k_death = 2,
    #rate for exponential noise, something big so noise is small (but non-zero))
    exp_noise=1e6)

  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)
  data$date <- as.Date(data$date)

  start_date <- as.Date("2020-02-02")
  set.seed(1)
  X <- particle_filter(data = data, pars_model, pars_obs, n_particles = 100,
                       start_date = start_date,
                       time_steps_per_day = time_steps_per_day)

  expect_is(X, "list")
  expect_equal(names(X), "log_likelihood")

  set.seed(1)
  Y <- particle_filter(data = data, pars_model, pars_obs, n_particles = 100,
                       start_date = start_date,
                       time_steps_per_day = time_steps_per_day,
                       output_states = TRUE,
                       save_particles = TRUE)
  expect_equal(X$log_likelihood, Y$log_likelihood)
  expect_setequal(names(Y), c("log_likelihood", "states", "index"))
  ##                            t   state  particles
  expect_equal(dim(Y$states), c(58, 238,   100))

  ## saveRDS(Y, "reference.rds")
  expect_equal(Y, readRDS("reference.rds"))

  ## Testing plotting is always a nightmare
  if (FALSE) {
    date <- as.Date(start_date) + seq_len(nrow(Y$states)) - 1L
    particles <- apply(Y$states[, c(Y$index$I_ICU) - 1L, ], c(1, 3), sum)
    plot_particles(particles = particles,
                   particle_dates = date,
                   data = data$itu / pars_obs$phi_ICU,
                   data_dates = data$date,
                   ylab = "ICU")
  }
})
