context("vaccine schedule")

test_that("validate a vaccine schedule", {
  doses <- array(1, c(19, 2, 10))
  expect_error(
    vaccine_schedule(1:10, doses),
    "'date' must be a scalar")
  expect_error(
    vaccine_schedule("2020-01-01", doses),
    "'date' must be numeric - did you forget sircovid_date")
  expect_error(
    vaccine_schedule(10, doses[-1, , ]),
    "'doses' must have 19 rows")
  expect_error(
    vaccine_schedule(10, doses[, rep(1, 3), ]),
    "'doses' must have 2 columns")
  expect_error(
    vaccine_schedule(10, doses[, , integer(0)]),
    "'doses' must have at least one element in the 3rd dimension")
  expect_error(
    vaccine_schedule(10, doses[, 1, , drop = TRUE]),
    "Expected a 3d array for 'doses'")
  expect_error(
    vaccine_schedule(10, -doses),
    "'doses' must all be non-negative")
  doses[30] <- NA_real_
  expect_error(
    vaccine_schedule(10, doses),
    "'doses' must all be non-NA")
})


test_that("vaccine_schedule_future produces consistent schedule", {
  region <- "london"
  uptake_by_age <- test_example_uptake()
  daily_doses <- rep(20000, 365) # a vector of number of doses to give each day
  mean_days_between_doses <- 12 * 7

  n <- vaccine_priority_population(region, uptake_by_age)
  dose_schedule <- vaccine_schedule_future(
    0, daily_doses, mean_days_between_doses, n)

  doses <- dose_schedule$doses
  n_to_vaccinate1 <- dose_schedule$doses[, 1, ]
  n_to_vaccinate2 <- dose_schedule$doses[, 2, ]

  ## initially only first doses are delivered
  ## and they add up to the target daily_doses
  ## not exactly equal because of some rounding
  phase1 <- seq_len(mean_days_between_doses)
  expect_true(all(abs(
    colSums(n_to_vaccinate1[, phase1]) - daily_doses[phase1]) < 5))
  expect_true(all(colSums(n_to_vaccinate2[, phase1]) == 0))

  ## during phase 2 only second doses are delivered
  ## as doses are constant over time
  ## and they add up to the target daily_doses
  ## not exactly equal because of some rounding
  phase2 <- mean_days_between_doses + seq_len(mean_days_between_doses)
  expect_true(all(abs(
    colSums(n_to_vaccinate2[, phase2]) - daily_doses[phase2]) < 5))
  expect_true(all(colSums(n_to_vaccinate1[, phase2]) == 0))

  ## check that n_to_vaccinate1 + n_to_vaccinate2 = daily_doses
  ## not exactly equal because of some rounding
  expect_true(all(abs(
      colSums(n_to_vaccinate1 + n_to_vaccinate2) - daily_doses) < 5))
})


test_that("vaccine_priority_proportion adds up to uptake", {
  uptake_by_age <- test_example_uptake()
  p <- vaccine_priority_proportion(uptake_by_age)
  expect_true(all(signif(rowSums(p), 5) == signif(uptake_by_age, 5)))
})


test_that("vaccine_priority_proportion with null uptake is 100%", {
  p <- vaccine_priority_proportion(NULL)
  expect_equal(p, vaccine_priority_proportion(rep(1, 19)))
  expect_equal(rowSums(p), rep(1, 19))
})


test_that("vaccine_priority_population adds to correct population", {
  ## Expecting rows to equal population size * uptake
  region <- "london"
  uptake_by_age <- test_example_uptake()
  n <- vaccine_priority_population(region, uptake_by_age)

  ## check that n to vaccinate adds up to uptake * population
  pop_by_age <- carehomes_parameters(1, region)$N_tot

  expect_true(all(rowSums(n) - uptake_by_age * pop_by_age <= 2))
})


test_that("Can append a schedule", {
  region <- "london"
  uptake_by_age <- test_example_uptake()
  daily_doses <- rep(20000, 365) # a vector of number of doses to give each day
  mean_days_between_doses <- 12 * 7

  n <- vaccine_priority_population(region, uptake_by_age)
  expected <- vaccine_schedule_future(
    0, daily_doses, mean_days_between_doses, n)

  ## We'll break this schedule in a few places and put it back
  ## together, aiming to do this in times where there are 2nd doses
  ## and not.
  i <- 50
  j <- seq_len(i)
  pre <- expected
  pre$doses <- pre$doses[, , j, drop = FALSE]

  check <- vaccine_schedule_future(
    pre, daily_doses[-j], mean_days_between_doses, n)

  ## This is subject to lots of rounding error. So we need to do some
  ## tricks to make sure it's broadly correct.
  expect_equal(check, expected, tolerance = 1 / 1000)

  ## The majority of cases will be identical, as rounding error only
  ## is a real issue at transitions.
  i <- expected$doses != 0
  expect_equal(median(check$doses[i] - expected$doses[i]), 0)
})


test_that("Cope with decreasing vaccination schedule", {
  region <- "london"
  uptake_by_age <- test_example_uptake()
  daily_doses <- seq(20000, length.out = 365, by = -50)
  mean_days_between_doses <- 12 * 7

  n <- vaccine_priority_population(region, uptake_by_age)
  res <- vaccine_schedule_future(
    0, daily_doses, mean_days_between_doses, n)

  d <- t(apply(res$doses, 2:3, sum))
  ## Daily dose schedule is followed within rounding error:
  expect_true(all(abs(rowSums(d) - daily_doses) < 10))
  tmp <- rle(apply(d, 1, which.max))
  expect_equal(tmp$values, c(1:2, 1:2))
  ## Our second window is quite a bit longer than the first
  r1 <- tmp$lengths[[2]] / tmp$lengths[[1]]
  expect_true(r1 > 1.3 && r1 < 1.5)
  ## But the second round of first doses is the same as the first
  ## because the decrease in dose rate is the same
  r2 <- tmp$lengths[[3]] / tmp$lengths[[1]]
  expect_true(r2 > 0.95 && r2 < 1.05)
})


test_that("Can add carehome residents to vaccine data", {
  data <- test_vaccine_data()

  uptake_by_age <- test_example_uptake()
  region <- "london"
  n_carehomes <-
    vaccine_priority_population(region, uptake_by_age)[18:19, 1]

  sched1 <- vaccine_schedule_from_data(data, c(0, 0))
  sched2 <- vaccine_schedule_from_data(data, c(1, 1))
  sched3 <- vaccine_schedule_from_data(data, n_carehomes)

  expect_equal(dim(sched1$doses), c(19, 2, 25))
  expect_equal(apply(sched1$doses, 2, sum),
               unname(colSums(data[c("dose1", "dose2")])))

  expect_equal(apply(sched1$doses, 3, sum),
               apply(sched2$doses, 3, sum))
  expect_equal(apply(sched1$doses, 3, sum),
               apply(sched3$doses, 3, sum))
  expect_equal(sched1$date, sched2$date)
  expect_equal(sched1$date, sched3$date)
})


test_that("can create a schedule that covers past and future", {
  set.seed(1)
  data <- test_vaccine_data()
  uptake <- test_example_uptake()
  n_carehomes <-
    vaccine_priority_population("london", uptake)[18:19, 1]
  set.seed(1)
  cmp <- vaccine_schedule_from_data(data, n_carehomes)
  set.seed(1) # we use rmultinom so be exactly the same
  end_date <- cmp$date + 25 + 14
  schedule <- vaccine_schedule_data_future(data, "london", uptake, end_date, 90)

  i <- seq_len(dim(cmp$doses)[[3]])
  expect_equal(schedule$doses[, , i], cmp$doses)
  doses_future <- schedule$doses[, , -i]
  expect_equal(dim(doses_future), c(19, 2, 14))

  mean_doses <- round(sum(cmp$doses[, , 19:25]) / 7)
  expect_equal(apply(doses_future, 3, sum), rep(mean_doses, 14))
})


test_that("Validate inputs in vaccine_schedule_from_data", {
  data <- test_vaccine_data()
  expect_error(
    vaccine_schedule_from_data(data[-1], c(0, 0)),
    "Required columns missing from 'data': 'date'")
  expect_error(
    vaccine_schedule_from_data(data[-seq_len(2)], c(0, 0)),
    "Required columns missing from 'data': 'age_band_min', 'date'")
  expect_error(
    vaccine_schedule_from_data(data, NULL),
    "Expected a vector of length 2 for n_carehomes")
  expect_error(
    vaccine_schedule_from_data(data, 1),
    "Expected a vector of length 2 for n_carehomes")

  data$age_band_min <- data$age_band_min + 1
  expect_error(
    vaccine_schedule_from_data(data, c(0, 0)),
    "Invalid values for data$age_band_min: 16, 21, 26, 31, 36, 41,",
    fixed = TRUE)
})
