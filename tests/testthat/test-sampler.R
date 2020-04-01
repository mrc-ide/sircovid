context("sampler")

## This is not a real test but simply tries to run the model
test_that("sampler runs without error", {
  set.seed(1)

  time_steps_per_day <- 4
  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)
  d <- particle_filter_data(data, "2020-02-02", time_steps_per_day)

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
    phi_death = 926 / 1019,
    ## what should this be?
    k_death = 2,
    #rate for exponential noise, something big so noise is small (but non-zero))
    exp_noise = 1e6)

  mod <- sircovid(params = pars_model)
  compare <- compare_icu(mod, pars_obs, d)

  cmp <- readRDS("reference.rds")
  dimnames(cmp$states) <- NULL
  cmp$index <- NULL

  set.seed(1)
  X <- particle_filter(d, mod, compare, n_particles = 100)
  expect_is(X, "list")
  expect_equal(names(X), "log_likelihood")
  expect_equal(X$log_likelihood, cmp$log_likelihood)

  set.seed(1)
  Y <- particle_filter(d, mod, compare, n_particles = 100,
                       save_particles = TRUE)
  expect_equal(X$log_likelihood, Y$log_likelihood)
  expect_setequal(names(Y), c("log_likelihood", "states"))
  ##                            t   state  particles
  expect_equal(dim(Y$states), c(58, 238,   100))

  date <- as.Date("2020-02-02") + seq_len(nrow(Y$states)) - 1L
  expect_equal(attr(Y$states, "date"), date)
  expect_equal(rownames(Y$states), as.character(date))
  ## Reset for comparison
  attr(Y$states, "date") <- NULL
  dimnames(Y$states) <- NULL
  expect_equal(Y, cmp)

  ## Testing plotting is always a nightmare
  if (FALSE) {
    set.seed(1)
    Y <- particle_filter(d, mod, compare, n_particles = 100,
                         save_particles = TRUE)
    index <- c(odin_index(mod)$I_ICU) - 1L
    particles <- apply(Y$states[, index, ], c(1, 3), sum)
    plot_particles(particles, ylab = "ICU")
    points(as.Date(data$date), data$itu / pars_obs$phi_ICU, pch = 19)
  }
})


test_that("systematic resample", {
  ## Reference version, from 9f59938d9e5f7a241f5e27c9a71da2f8f643079c
  ## with minimal changes for style.
  ref <- function(weights) {
    n <- length(weights)
    old_indexes <- seq_len(n)
    u <- runif(1, 0, 1 / n) + seq(0, by = 1 / n,length.out = n)
    new_indexes <- integer(n)
    weights <- weights / sum(weights)
    cum_weights <- cumsum(weights)
    k <- 1
    for (i in seq_len(n)) {
      found <- FALSE
      while (!found){
        if (u[i] > cum_weights[k]) {
          k <- k + 1
        } else {
          found <- TRUE
        }
      }
      new_indexes[i] <- old_indexes[k]
    }
    new_indexes
  }

  w <- lapply(rpois(400, 30), runif)

  set.seed(1)
  cmp <- lapply(w, ref)

  set.seed(1)
  ans <- lapply(w, systematic_resample)
  expect_identical(cmp, ans)
})


test_that("particle_filter error cases", {
  set.seed(1)

  time_steps_per_day <- 4
  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)
  d <- particle_filter_data(data, "2020-02-02", time_steps_per_day)

  pars_model <- generate_parameters(
    dt = 1 / time_steps_per_day,
    beta = rep(0.125, 3))
  pars_obs <- list(phi_ICU = 0.95, k_ICU = 2, phi_death = 926 / 1019,
                   k_death = 2, exp_noise = 1e6)

  mod <- sircovid(params = pars_model)
  compare <- compare_icu(mod, pars_obs, d)

  expect_error(
    particle_filter(NULL, mod, compare, 100),
    "Expected a data set derived from particle_filter_data")
  expect_error(
    particle_filter(data, mod, compare, 100),
    "Expected a data set derived from particle_filter_data")
  expect_error(
    particle_filter(d, NULL, compare, 100),
    "Expected 'model' to be an 'odin_model' object")
  expect_error(
    particle_filter(d, mod, compare, 1),
    "At least two particles required")
  expect_error(
    particle_filter(d, mod, compare, 100, forecast_days = 1),
    "forecasting only possible if particles are saved")
  expect_error(
    particle_filter(d, mod, compare, 100, forecast_days = -1),
    "forecast_days must be positive")
})


test_that("particle filter data; error validation", {
  start <- Sys.Date()
  date <- start + 1:4
  steps_per_day <- 4
  expect_error(
    particle_filter_data(1, start, steps_per_day),
    "Expected a data.frame")
  expect_error(
    particle_filter_data(data.frame(a = 1:10), start, steps_per_day),
    "Expected a column 'date' within 'data'")
  expect_error(
    particle_filter_data(data.frame(date = start - 0:4), start,
                         steps_per_day),
    "'date' must be strictly increasing")
  expect_error(
    particle_filter_data(data.frame(date = start + c(0, 1, 1, 2)),
                         start, steps_per_day),
    "'date' must be strictly increasing")
  expect_error(
    particle_filter_data(data.frame(date = date - 1), start, steps_per_day),
    "'start_date' must be less than the first date in data")
})


test_that("particle filter data", {
  start <- as.Date("2020-02-02")
  data <- data.frame(date = as.Date("2020-03-01") + 0:30,
                     a = 0:30,
                     b = 0:30 * 2)
  steps_per_day <- 4

  d <- particle_filter_data(data, start, steps_per_day)
  expect_equal(attr(d, "steps_per_day"), steps_per_day)
  attr(d, "steps_per_day") <-  NULL

  expect_equal(as.list(d[1, ]),
               list(date = start, a = NA_integer_, b = NA_real_,
                    day_start = 0, day_end = 27,
                    step_start = 0, step_end = 108))
  expect_equal(as.list(d[2, ]),
               list(date = data$date[[1]], a = 0L, b = 0,
                    day_start = 27, day_end = 28,
                    step_start = 108, step_end = 112))
  expect_equal(nrow(d), 32)
  expect_equal(as.list(d[32, ]),
               list(date = data$date[[31]], a = 30L, b = 60,
                    day_start = 57, day_end = 58,
                    step_start = 228, step_end = 232))
})
