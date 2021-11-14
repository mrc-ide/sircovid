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
  expect_vector_equal(
    colSums(n_to_vaccinate1[, phase1]), daily_doses[phase1], tol = 5)
  expect_vector_equal(colSums(n_to_vaccinate2[, phase1]), 0)

  ## during phase 2 only second doses are delivered
  ## as doses are constant over time
  ## and they add up to the target daily_doses
  ## not exactly equal because of some rounding
  phase2 <- mean_days_between_doses + seq_len(mean_days_between_doses)
  expect_vector_equal(
    colSums(n_to_vaccinate2[, phase2]), daily_doses[phase2], tol = 5)
  expect_vector_equal(colSums(n_to_vaccinate1[, phase2]), 0)

  ## check that n_to_vaccinate1 + n_to_vaccinate2 = daily_doses
  ## not exactly equal because of some rounding
  expect_vector_equal(
    colSums(n_to_vaccinate1 + n_to_vaccinate2), daily_doses, tol = 5)
})


test_that("vaccine_priority_proportion adds up to uptake", {
  uptake_by_age <- test_example_uptake()
  p <- vaccine_priority_proportion(uptake_by_age)
  expect_vector_equal(rowSums(p), uptake_by_age, digits = 5)
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
  pop_by_age <- lancelot_parameters(1, region)$N_tot

  expect_vector_equal(rowSums(n), uptake_by_age * pop_by_age, tol = 2)
})


test_that("vaccine_priority_population adds to correct pop - full uptake", {
  region <- "london"
  n <- vaccine_priority_population(region, NULL)
  pop_by_age <- lancelot_parameters(1, region)$N_tot
  expect_vector_equal(rowSums(n), pop_by_age, tol = 2)
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
  expect_vector_equal(rowSums(d), daily_doses, tol = 10)
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
  uptake_by_age <- array(uptake_by_age, c(length(uptake_by_age), 2))

  uptake_by_age1 <- uptake_by_age
  uptake_by_age1[18:19, ] <- 0
  uptake_by_age2 <- uptake_by_age
  uptake_by_age2[18:19, ] <- 1

  region <- "london"

  sched1 <- vaccine_schedule_from_data(data, region, uptake_by_age1)
  sched2 <- vaccine_schedule_from_data(data, region, uptake_by_age2)
  sched3 <- vaccine_schedule_from_data(data, region, uptake_by_age)

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


test_that("can input data with 3 doses", {
  data <- test_vaccine_data()

  data$dose3 <- rpois(nrow(data), 50)

  region <- "london"
  uptake_by_age <- test_example_uptake()
  uptake_by_age <- array(uptake_by_age, c(length(uptake_by_age), 3))

  data_schedule <- vaccine_schedule_from_data(data, region, uptake_by_age)

  expect_equal(data_schedule$n_doses, 3)
  expect_equal(dim(data_schedule$doses), c(19, 3, 25))
})


test_that("can create data with age-aggregated doses", {
  data <- test_vaccine_data()

  region <- "london"
  uptake_by_age <- test_example_uptake()
  uptake_by_age <- array(uptake_by_age, c(length(uptake_by_age), 2))

  data_schedule <- vaccine_schedule_from_data(data, region, uptake_by_age)

  ## Only age-aggregated data
  agg_data <- data
  agg_data$age_band_min <- 0
  agg_data <- stats::aggregate(agg_data[c("dose1", "dose2")],
                               agg_data[c("age_band_min", "date")],
                               sum)
  agg_data$age_band_min <- NA

  data_schedule2 <- vaccine_schedule_from_data(agg_data, region, uptake_by_age)


  ## Mixed age-specific and age-aggregated
  ## Aggregate half of all doses
  agg_data2 <- data.frame(
    date = data$date,
    age_band_min = 0,
    dose1 = round(0.5 * data$dose1),
    dose2 = round(0.5 * data$dose2)
  )
  mixed_data <- data
  mixed_data[, c("dose1", "dose2")] <-
    mixed_data[, c("dose1", "dose2")] - agg_data2[, c("dose1", "dose2")]
  agg_data2 <- stats::aggregate(agg_data2[c("dose1", "dose2")],
                                agg_data2[c("age_band_min", "date")],
                                sum)
  agg_data2$age_band_min <- NA
  mixed_data <- rbind(mixed_data, agg_data2)

  data_schedule3 <-
    vaccine_schedule_from_data(mixed_data, region, uptake_by_age)

  ## Total doses will not be exact but should be close
  expect_equal(dim(data_schedule$doses), dim(data_schedule2$doses))
  expect_equal(dim(data_schedule$doses), dim(data_schedule3$doses))
  expect_approx_equal(apply(data_schedule$doses, 2, sum),
                      apply(data_schedule2$doses, 2, sum), 10)
  expect_approx_equal(apply(data_schedule$doses, 2, sum),
                      apply(data_schedule3$doses, 2, sum), 10)

})


test_that("Validate inputs in vaccine_schedule_from_data", {
  data <- test_vaccine_data()
  uptake_by_age <- test_example_uptake()
  uptake_by_age <- array(uptake_by_age, c(length(uptake_by_age), 2))
  region <- "london"

  expect_error(
    vaccine_schedule_from_data(data[-1], region, uptake_by_age),
    "Required columns missing from 'data': 'date'")
  expect_error(
    vaccine_schedule_from_data(data[-seq_len(2)], region, uptake_by_age),
    "Required columns missing from 'data': 'age_band_min', 'date'")

  data$age_band_min <- data$age_band_min + 1
  expect_error(
    vaccine_schedule_from_data(data, region, uptake_by_age),
    "Invalid values for data$age_band_min: 16, 21, 26, 31, 36, 41,",
    fixed = TRUE)

  data$age_band_min <- data$age_band_min - 1
  names(data)[names(data) == "dose2"] <- "dose3"
  expect_error(
    vaccine_schedule_from_data(data, region, uptake_by_age),
    "There are 2 dose columns so expected dose column names: 'dose1', 'dose2'")

  data$dose2 <- data$dose3
  expect_error(
    vaccine_schedule_from_data(data, region, uptake_by_age),
    "Data has 3 dose columns so expected uptake to have 3")

  expect_error(
    vaccine_schedule_from_data(data, region, uptake_by_age[-19, ]),
    "Expected uptake to have 19 rows as there are 19 groups")
})


test_that("create schedule scenario", {
  data <- test_vaccine_data()

  region <- "london"
  uptake_by_age <- test_example_uptake()
  n <- vaccine_priority_population(region, uptake_by_age)
  uptake_by_age <- array(uptake_by_age, c(length(uptake_by_age), 2))
  past <- vaccine_schedule_from_data(data, region, uptake_by_age)

  mean_days_between_doses <- 30
  doses_future <- c(
    "2021-04-10" = 6000,
    "2021-04-20" = 7000,
    "2021-04-30" = 9000)
  end_date <- "2021-08-01"

  res <- vaccine_schedule_scenario(past, doses_future, end_date,
                                   mean_days_between_doses, n)

  i <- seq_len(dim(past$doses)[[3]])
  expect_equal(res$doses[, , i], past$doses)
  doses_future <- res$doses[, , -i]
  expect_equal(dim(doses_future), c(19, 2, 125))

  n <- apply(doses_future, 3, sum)
  expect_equal(
    signif(n, 1),
    rep(c(5000, 6000, 7000, 9000), c(12, 10, 10, 93)))
  expect_lt(max(abs(n - signif(n, 1))[-seq_len(12)]), 10)
})


test_that("create schedule scenario where future doses drift into past", {
  data <- test_vaccine_data()

  region <- "london"
  uptake_by_age <- test_example_uptake()
  n <- vaccine_priority_population(region, uptake_by_age)
  uptake_by_age <- array(uptake_by_age, c(length(uptake_by_age), 2))
  past <- vaccine_schedule_from_data(data, region, uptake_by_age)

  ## Our schedule runs from 2021-03-05..2021-03-29
  mean_days_between_doses <- 30
  doses_future <- c(
    "2021-03-20" = 6000,
    "2021-04-20" = 7000,
    "2021-04-30" = 9000)
  end_date <- "2021-08-01"

  res <- vaccine_schedule_scenario(past, doses_future, end_date,
                                   mean_days_between_doses, n)

  i <- seq_len(dim(past$doses)[[3]])
  expect_equal(res$doses[, , i], past$doses)
  doses_future <- res$doses[, , -i]
  expect_equal(dim(doses_future), c(19, 2, 125))

  n <- apply(doses_future, 3, sum)
  expect_equal(
    signif(n, 1),
    rep(c(6000, 7000, 9000), c(22, 10, 93)))
  expect_lt(max(abs(n - signif(n, 1))), 10)
})


test_that("create schedule scenario with no extra dose info", {
  set.seed(10)
  data <- test_vaccine_data()

  region <- "london"
  uptake_by_age <- test_example_uptake()
  n <- vaccine_priority_population(region, uptake_by_age)
  set.seed(1)
  past <- vaccine_schedule_from_data(data, region,
                                     array(uptake_by_age,
                                           c(length(uptake_by_age), 2)))

  mean_days_between_doses <- 30
  end_date <- "2021-06-01"

  res <- vaccine_schedule_scenario(past, NULL, end_date,
                                   mean_days_between_doses, n)
  set.seed(1)
  ## Commenting this out for the moment
  ## cmp <- vaccine_schedule_data_future(data, region, uptake_by_age,
  ##                                     end_date, mean_days_between_doses)

  ## expect_equal(dim(res$doses), dim(cmp$doses))
  # ## rounding error nightmare, even with making the above deterministic
  ## expect_lt(max(abs(res$doses - cmp$doses)), 5)

  d <- seq(res$date, length.out = dim(res$doses)[[3]])
  expect_equal(last(d), sircovid_date(end_date))
})


test_that("prevent impossible scenarios", {
  data <- test_vaccine_data()

  region <- "london"
  uptake_by_age <- test_example_uptake()
  n <- vaccine_priority_population(region, uptake_by_age)
  uptake_by_age <- array(uptake_by_age, c(length(uptake_by_age), 2))
  past <- vaccine_schedule_from_data(data, region, uptake_by_age)

  mean_days_between_doses <- 30
  doses_future <- c(
    "2021-04-10" = 6000,
    "2021-04-20" = 7000,
    "2021-04-30" = 9000)
  end_date <- "2021-08-01"

  expect_error(
    vaccine_schedule_scenario(past, unname(doses_future), end_date,
                              mean_days_between_doses, n),
    "'doses_future' must be named")
  expect_error(
    vaccine_schedule_scenario(past, set_names(doses_future, c("a", "b", "c")),
                              end_date, mean_days_between_doses, n),
    "Expected ISO dates or R dates for 'names(doses_future)'",
    fixed = TRUE)
  expect_error(
    vaccine_schedule_scenario(past, rev(doses_future), end_date,
                              mean_days_between_doses, n),
    "'names(doses_future)' must be strictly increasing",
    fixed = TRUE)
  expect_error(
    vaccine_schedule_scenario(past, doses_future[c(1, 2, 2, 3)], end_date,
                              mean_days_between_doses, n),
    "'names(doses_future)' must be strictly increasing",
    fixed = TRUE)
  expect_error(
    vaccine_schedule_scenario(past, doses_future, "2021-04-29",
                              mean_days_between_doses, n),
    paste("'end_date' must be at least 2021-04-30 (last doses_future date)",
          "but was 2021-04-29"),
    fixed = TRUE)
  expect_error(
    vaccine_schedule_scenario(past, NULL, "2021-03-01",
                              mean_days_between_doses, n),
    paste("'end_date' must be at least 2021-03-29 (previous end date)",
          "but was 2021-03-01"),
    fixed = TRUE)
})


test_that("vaccine_schedule_future functions with boosters", {
  region <- "london"
  uptake_by_age <- test_example_uptake()
  daily_doses <- rep(20000, 100) # a vector of number of doses to give each day
  mean_days_between_doses <- 12 * 7

  booster_doses <- c(rep(0, 150), rep(10000, 50))

  n <- vaccine_priority_population(region, uptake_by_age)
  dose_schedule <- vaccine_schedule_future(
    0, daily_doses, mean_days_between_doses, n,
    booster_daily_doses_value = booster_doses)

  ## errors when trying to add a lag
  expect_error(vaccine_schedule_future(
    0, daily_doses, mean_days_between_doses, n,
    booster_daily_doses_value = booster_doses, lag_days = 1))

  doses <- dose_schedule$doses
  n_to_vaccinate1 <- dose_schedule$doses[, 1, ]
  n_to_vaccinate2 <- dose_schedule$doses[, 2, ]
  n_to_vaccinate3 <- dose_schedule$doses[, 3, ]

  ## initially only first doses are delivered
  ## and they add up to the target daily_doses
  ## not exactly equal because of some rounding
  phase1 <- seq_len(mean_days_between_doses)
  expect_vector_equal(
    colSums(n_to_vaccinate1[, phase1]), daily_doses[phase1], tol = 5)
  expect_vector_equal(colSums(n_to_vaccinate2[, phase1]), 0)
  expect_vector_equal(colSums(n_to_vaccinate3[, phase1]), 0)

  ## during phase 2 only second doses are delivered
  ## as doses are constant over time
  ## and they add up to the target daily_doses
  ## not exactly equal because of some rounding
  phase2 <- 85:100
  expect_vector_equal(
    colSums(n_to_vaccinate2[, phase2]), daily_doses[phase2], tol = 5)
  expect_vector_equal(colSums(n_to_vaccinate1[, phase2]), 0)
  expect_vector_equal(colSums(n_to_vaccinate3[, phase2]), 0)

  ## during phase 3, no doses or boosters
  phase3 <- 101:150
  expect_vector_equal(colSums(n_to_vaccinate1[, phase3]), 0)
  expect_vector_equal(colSums(n_to_vaccinate2[, phase3]), 0)
  expect_vector_equal(colSums(n_to_vaccinate3[, phase3]), 0)

  ## during phase 4 only boosters are delivered
  phase4 <- 151:200
  expect_vector_equal(
    colSums(n_to_vaccinate3[, phase4]), booster_doses[phase4], tol = 5)
  expect_vector_equal(colSums(n_to_vaccinate1[, phase4]), 0)
  expect_vector_equal(colSums(n_to_vaccinate2[, phase4]), 0)

  ## check that n_to_vaccinate1 + n_to_vaccinate2 = daily_doses
  ## not exactly equal because of some rounding
  expect_vector_equal(
    colSums(n_to_vaccinate1 + n_to_vaccinate2)[1:100],
    daily_doses[1:100], tol = 5)

  ## check that n_to_vaccinate3 = booster_doses
  ## not exactly equal because of some rounding
  expect_vector_equal(
    colSums(n_to_vaccinate3), booster_doses, tol = 5)
})


test_that("Can add boosters to schedule", {
  region <- "london"
  uptake_by_age <- test_example_uptake()
  daily_doses <- rep(20000, 100) # a vector of number of doses to give each day
  mean_days_between_doses <- 12 * 7

  booster_doses <- c(rep(0, 150), rep(10000, 50))

  n <- vaccine_priority_population(region, uptake_by_age)

  ## schedule with doses and boosters
  expected_schedule <- vaccine_schedule_future(
    0, daily_doses, mean_days_between_doses, n,
    booster_daily_doses_value = booster_doses)

  ## schedule with doses only
  dose_schedule <- vaccine_schedule_future(
    0, daily_doses, mean_days_between_doses, n)

  ## schedule with doses only
  daily_doses <- numeric(100)
  booster_doses <- c(rep(0, 50), rep(10000, 50))

  complete_schedule <- vaccine_schedule_future(
    dose_schedule, daily_doses, mean_days_between_doses, n,
    booster_daily_doses_value = booster_doses)

  expect_equal(expected_schedule, complete_schedule)
})



test_that("check schedule scenario prepends zeros when needed", {
  data <- test_vaccine_data()

  region <- "london"
  uptake_by_age <- test_example_uptake()
  n <- vaccine_priority_population(region, uptake_by_age)
  uptake_by_age <- array(uptake_by_age, c(length(uptake_by_age), 2))
  past <- vaccine_schedule_from_data(data, region, uptake_by_age)

  sircovid_date_as_date(past$date + dim(past$doses)[[3]] - 1L)

  mean_days_between_doses <- 30
  doses_future <- c(
    "2021-04-10" = 6000,
    "2021-04-20" = 7000,
    "2021-04-30" = 9000,
    "2021-05-10" = 0)

  end_date <- "2021-08-01"

  expect_zeros <- rep(c(5000, 6000, 7000, 9000, 0, 600, 700, 900),
                      c(12, 10, 10, 10, 31, 10, 10, 32))
  i <- seq_len(dim(past$doses)[[3]])

  ## 1. manually and automatically add zeros
  boosters_future <- c(
    "2021-03-29" = 0,
    "2021-06-10" = 600,
    "2021-06-20" = 700,
    "2021-06-30" = 900)


  res <- vaccine_schedule_scenario(past, doses_future, end_date,
                                   mean_days_between_doses, n,
                                   boosters_future = boosters_future,
                                   boosters_prepend_zero = TRUE)

  expect_equal(signif(apply(res$doses[, , -i], 3, sum), 1), expect_zeros)

  ## 2. automatically add zeros
  boosters_future <- c(
    "2021-06-10" = 600,
    "2021-06-20" = 700,
    "2021-06-30" = 900)

  res <- vaccine_schedule_scenario(past, doses_future, end_date,
                                   mean_days_between_doses, n,
                                   boosters_future = boosters_future,
                                   boosters_prepend_zero = TRUE)

  expect_equal(signif(apply(res$doses[, , -i], 3, sum), 1), expect_zeros)

  ## 3. manually add zeros
  boosters_future <- c(
    "2021-03-29" = 0,
    "2021-06-10" = 600,
    "2021-06-20" = 700,
    "2021-06-30" = 900)

  res <- vaccine_schedule_scenario(past, doses_future, end_date,
                                   mean_days_between_doses, n,
                                   boosters_future = boosters_future,
                                   boosters_prepend_zero = FALSE)

  expect_equal(signif(apply(res$doses[, , -i], 3, sum), 1), expect_zeros)

  ## 4. no zeros
  boosters_future <- c(
    "2021-06-10" = 600,
    "2021-06-20" = 700,
    "2021-06-30" = 900)

  res <- vaccine_schedule_scenario(past, doses_future, end_date,
                                   mean_days_between_doses, n,
                                   boosters_future = boosters_future,
                                   boosters_prepend_zero = FALSE)

  ## 5000 booster doses are being added to everything before 2021-06-10
  expect_no_zeros <- rep(c(10000, 5000, 600, 700, 900),
                      c(42, 31, 10, 10, 32))
  expect_equal(signif(apply(res$doses[, , -i], 3, sum), 1), expect_no_zeros)
})


test_that("create schedule scenario with doses and boosters", {
  data <- test_vaccine_data()

  region <- "london"
  uptake_by_age <- test_example_uptake()
  n <- vaccine_priority_population(region, uptake_by_age)
  uptake_by_age <- array(uptake_by_age, c(length(uptake_by_age), 2))
  past <- vaccine_schedule_from_data(data, region, uptake_by_age)

  sircovid_date_as_date(past$date + dim(past$doses)[[3]] - 1L)

  mean_days_between_doses <- 30
  doses_future <- c(
    "2021-04-10" = 6000,
    "2021-04-20" = 7000,
    "2021-04-30" = 9000,
    "2021-05-10" = 0)

  boosters_future <- c(
    "2021-06-10" = 600,
    "2021-06-20" = 700,
    "2021-06-30" = 900)
  end_date <- "2021-08-01"

  res <- vaccine_schedule_scenario(past, doses_future, end_date,
                                   mean_days_between_doses, n,
                                   boosters_future = boosters_future)

  i <- seq_len(dim(past$doses)[[3]])
  expect_equal(res$doses[, 1:2, i], past$doses)
  doses_future <- res$doses[, , -i]
  expect_equal(dim(doses_future), c(19, 3, 125))

  n <- apply(doses_future, 3, sum)
  expect_equal(
    signif(n, 1),
    rep(c(5000, 6000, 7000, 9000, 0, 600, 700, 900),
        c(12, 10, 10, 10, 31, 10, 10, 32)))
  expect_lt(max(abs(n - signif(n, 1))[-seq_len(12)]), 10)

  expect_vector_equal(doses_future[, 3, 1:73], 0)
  expect_false(all(doses_future[, 3, 74] == 0))
})


test_that("create schedule scenario with boosters only", {
  data <- test_vaccine_data()

  region <- "london"
  uptake_by_age <- test_example_uptake()
  n <- vaccine_priority_population(region, uptake_by_age)
  uptake_by_age <- array(uptake_by_age, c(length(uptake_by_age), 2))
  past <- vaccine_schedule_from_data(data, region, uptake_by_age)

  sircovid_date_as_date(past$date + dim(past$doses)[[3]] - 1L)

  mean_days_between_doses <- 30
  boosters_future <- c(
    "2021-06-10" = 600,
    "2021-06-20" = 700,
    "2021-06-30" = 900
  )
  end_date <- "2021-08-01"

  res <- vaccine_schedule_scenario(past, doses_future = NULL, end_date,
                                   mean_days_between_doses, n,
                                   boosters_future = boosters_future
  )

  i <- seq_len(dim(past$doses)[[3]])
  expect_equal(res$doses[, 1:2, i], past$doses)
  doses_future <- res$doses[, , -i]
  expect_equal(dim(doses_future), c(19, 3, 125))

  n <- apply(doses_future, 3, sum)
  expect_equal(
    signif(n, 1),
    rep(c(5000, 6000), c(73, 52))
  )

  expect_vector_equal(doses_future[, 3, 1:73], 0)
  expect_false(all(doses_future[, 3, 74] == 0))
})


test_that("can exclude groups from vaccine_schedule_future", {
  region <- "london"
  uptake_by_age <- test_example_uptake()
  daily_doses <- rep(20000, 100) # a vector of number of doses to give each day
  mean_days_between_doses <- 12 * 7

  booster_doses <- c(rep(0, 150), rep(10000, 50))

  n <- vaccine_priority_population(region, uptake_by_age)
  dose_schedule_all <- vaccine_schedule_future(
    0, daily_doses, mean_days_between_doses, n,
    booster_daily_doses_value = booster_doses)

  dose_schedule_ch <- vaccine_schedule_future(
    0, daily_doses, mean_days_between_doses, n,
    booster_daily_doses_value = booster_doses,
    booster_proportion = c(numeric(17), 1, 1))

  expect_equal(dose_schedule_all$date, dose_schedule_ch$date)
  expect_equal(dose_schedule_all$n_doses, dose_schedule_ch$n_doses)

  doses_all <- dose_schedule_all$doses
  doses_ch <- dose_schedule_ch$doses

  ## check 0 when expected
  expect_vector_gte(doses_all, doses_ch)
  expect_equal(doses_all[18:19, , ], doses_ch[18:19, , ])
  expect_vector_equal(doses_ch[1:17, 3, ], 0)
  expect_false(all(doses_ch[18:19, 3, ] == 0))

  ## check totals still as expected
  n_to_vaccinate1 <- doses_ch[, 1, ]
  n_to_vaccinate2 <- doses_ch[, 2, ]
  n_to_vaccinate3 <- doses_ch[, 3, ]

  ## initially only first doses are delivered
  ## and they add up to the target daily_doses
  ## not exactly equal because of some rounding
  phase1 <- seq_len(mean_days_between_doses)
  expect_vector_equal(
    colSums(n_to_vaccinate1[, phase1]), daily_doses[phase1], tol = 5)
  expect_vector_equal(colSums(n_to_vaccinate2[, phase1]), 0)
  expect_vector_equal(colSums(n_to_vaccinate3[, phase1]), 0)

  ## during phase 2 only second doses are delivered
  ## as doses are constant over time
  ## and they add up to the target daily_doses
  ## not exactly equal because of some rounding
  phase2 <- 85:100
  expect_vector_equal(
    colSums(n_to_vaccinate2[, phase2]), daily_doses[phase2], tol = 5)
  expect_vector_equal(colSums(n_to_vaccinate1[, phase2]), 0)
  expect_vector_equal(colSums(n_to_vaccinate3[, phase2]), 0)

  ## during phase 3, no doses or boosters
  phase3 <- 101:150
  expect_vector_equal(colSums(n_to_vaccinate1[, phase3]), 0)
  expect_vector_equal(colSums(n_to_vaccinate2[, phase3]), 0)
  expect_vector_equal(colSums(n_to_vaccinate3[, phase3]), 0)

  ## during phase 4 only boosters are delivered
  phase4 <- 151:200
  ## expect not to be equal because not all doses given once 18/19 filled
  expect_vector_nequal(
    colSums(n_to_vaccinate3[, phase4]), booster_doses[phase4], tol = 5)
  expect_equal(sum(n_to_vaccinate3[, phase4]), sum(n[18:19, ]))
  expect_vector_equal(colSums(n_to_vaccinate1[, phase4]), 0)
  expect_vector_equal(colSums(n_to_vaccinate2[, phase4]), 0)

  ## check that n_to_vaccinate1 + n_to_vaccinate2 = daily_doses
  ## not exactly equal because of some rounding
  expect_vector_equal(
    colSums(n_to_vaccinate1 + n_to_vaccinate2)[1:100],
    daily_doses[1:100], tol = 5)

  ## check that n_to_vaccinate3 = booster_doses
  ## not exactly equal because of some rounding
  expect_vector_nequal(
    colSums(n_to_vaccinate3), booster_doses, tol = 5)
  expect_equal(
    sum(n_to_vaccinate3), sum(n[18:19, ]), tolerance = 5)
})


test_that("can exclude booster groups from schedule scenario", {
  data <- test_vaccine_data()

  region <- "london"
  uptake_by_age <- test_example_uptake()
  n <- vaccine_priority_population(region, uptake_by_age)
  uptake_by_age <- array(uptake_by_age, c(length(uptake_by_age), 2))
  past <- vaccine_schedule_from_data(data, region, uptake_by_age)

  sircovid_date_as_date(past$date + dim(past$doses)[[3]] - 1L)

  mean_days_between_doses <- 30
  boosters_future <- c(
    "2021-06-10" = 6000,
    "2021-06-20" = 7000,
    "2021-06-30" = 9000
  )
  end_date <- "2021-08-01"

  doses_all <- vaccine_schedule_scenario(past, doses_future = NULL, end_date,
                                   mean_days_between_doses, n,
                                   boosters_future = boosters_future)$doses

  doses_ch <- vaccine_schedule_scenario(
    past, doses_future = NULL, end_date,
    mean_days_between_doses, n,
    boosters_future = boosters_future,
    booster_proportion = c(numeric(17), 1, 1))$doses

  expect_vector_nequal(doses_all, doses_ch)

  i <- seq_len(dim(past$doses)[[3]])
  expect_equal(doses_ch[, 1:2, i], past$doses)
  doses_future <- doses_ch[, , -i]
  expect_equal(dim(doses_future), c(19, 3, 125))

  expect_vector_equal(doses_future[1:17, 3, ], 0)
  expect_vector_nequal(doses_future[18:19, 3, ], 0)
})



test_that("can partially exclude booster groups from schedule scenario", {
  data <- test_vaccine_data()

  region <- "london"
  uptake_by_age <- test_example_uptake()
  n <- vaccine_priority_population(region, uptake_by_age)
  uptake_by_age <- array(uptake_by_age, c(length(uptake_by_age), 2))
  past <- vaccine_schedule_from_data(data, region, uptake_by_age)

  sircovid_date_as_date(past$date + dim(past$doses)[[3]] - 1L)

  mean_days_between_doses <- 30
  boosters_future <- c(
    "2021-06-10" = 6000,
    "2021-06-20" = 7000,
    "2021-06-30" = 9000
  )
  end_date <- "2021-08-01"

  doses_half <- vaccine_schedule_scenario(
    past, doses_future = NULL, end_date, mean_days_between_doses, n,
    boosters_future = boosters_future,
    booster_proportion = c(numeric(17), 0.5, 0.5))$doses

  doses_full <- vaccine_schedule_scenario(
    past, doses_future = NULL, end_date, mean_days_between_doses, n,
    boosters_future = boosters_future,
    booster_proportion = c(numeric(17), 1, 1))$doses

  expect_vector_equal(doses_half[1:17, 3, ], 0)
  expect_vector_equal(doses_full[1:17, 3, ], 0)

  ## Check that double the number of 18:19 are vaccinated in doses_full than
  ##  doses_half
  expect_equal(sum(colSums(doses_full[18:19, 3, ]) > 0),
               sum(colSums(doses_half[18:19, 3, ]) > 0) * 2)
})
