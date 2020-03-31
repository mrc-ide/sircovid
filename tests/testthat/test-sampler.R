context("sampler")

## This is not a real test but simply tries to run the model
test_that("sampler runs without error", {
  set.seed(1)

  time_steps_per_day <- 4
  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)
  d <- particle_filter_data(data, "2020-02-02", time_steps_per_day)

  pars_model <- generate_parameter_rtm(
    seed_SSP=10,
    dt = 1 / time_steps_per_day,
    beta = 0.125)

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

  mod <- sircovid(params = pars_model)
  compare <- compare_icu(odin_index(mod), pars_obs, d)

  set.seed(1)
  X <- particle_filter(d, mod, compare, n_particles = 100)
  expect_is(X, "list")
  expect_equal(names(X), "log_likelihood")

  set.seed(1)
  Y <- particle_filter(d, mod, compare, n_particles = 100,
                       save_particles = TRUE)
  expect_equal(X$log_likelihood, Y$log_likelihood)
  expect_setequal(names(Y), c("log_likelihood", "states"))
  ##                            t   state  particles
  expect_equal(dim(Y$states), c(58, 238,   100))

  expect_equal(X$log_likelihood, Y$log_likelihood)

  cmp <- readRDS("reference.rds")
  dimnames(cmp$states) <- NULL
  cmp$index <- NULL
  expect_equal(Y, cmp)

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
