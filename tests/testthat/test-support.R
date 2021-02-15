context("date")

test_that("sircovid_date can convert to dates into 2020", {
  expect_equal(sircovid_date("2020-01-01"), 1)
  expect_equal(sircovid_date("2020-12-31"), 366) # leap year!
  r <- seq(as_date("2020-01-01"), as_date("2020-12-31"), by = 1)
  expect_equal(sircovid_date(r), 1:366)

  expect_equal(sircovid_date_as_date(sircovid_date(r)), r)

  expect_error(sircovid_date("2019-01-01"),
               "Negative dates, sircovid_date likely applied twice")
  expect_error(sircovid_date(c("2020-01-01", "2019-01-01")),
               "Negative dates, sircovid_date likely applied twice")
})


test_that("assert sircovid date throws on non sircovid dates", {
  expect_silent(assert_sircovid_date(1))
  expect_error(assert_sircovid_date(as_date("2020-01-01")),
               "'date' must be numeric - did you forget sircovid_date()?",
               fixed = TRUE)
})


test_that("helper function can convert somewhat helpfully", {
  expect_equal(as_sircovid_date("2020-02-01"), 32)
  expect_equal(as_sircovid_date(as_date("2020-02-01")), 32)
  expect_equal(as_sircovid_date(32), 32)
})


test_that("negative numbers are not allowed", {
  expect_error(assert_sircovid_date(-1),
               "Negative dates, sircovid_date likely applied twice")
  expect_error(assert_sircovid_date(c(10, -1, 29)),
               "Negative dates, sircovid_date likely applied twice")
})


test_that("read population", {
  clear_cache()
  n <- sircovid_population("south_west")
  expect_s3_class(cache$population, "data.frame")
  expected <- c(296855, 325583, 307586, 304574, 339548, 337183, 327525,
                327910, 311774, 374176, 403024, 380109, 338931, 332910,
                333351, 225646, 339312)
  expect_equal(n, expected)
})


test_that("reject invalid regions", {
  expect_error(sircovid_population("oxfordshire"),
               "Population not found for 'oxfordshire': must be one of 'uk',")
  expect_error(sircovid_population(NULL),
               "'region' must not be NULL")
})


test_that("Downcase region name", {
  expect_identical(sircovid_population("UK"), sircovid_population("uk"))
})


test_that("Use constant age bins", {
  expect_equal(
    sircovid_age_bins(),
    list(start = seq(0, 80, by = 5),
         end = c(seq(4, 79, by = 5), 100)))
})


test_that("validate data against age bins", {
  good <- sprintf("%d to %d", seq(0, 80, by = 5), c(seq(4, 79, by = 5), 100))
  expect_equal(check_age_bins(good), sircovid_age_bins())

  expect_error(check_age_bins(good[-1]),
               "Incorrect age bands:")
  expect_error(check_age_bins(good[-1]),
               "expected: '0 to 4', '5 to 9'")
  expect_error(check_age_bins(good[-1]),
               "expected: '0 to 4', '5 to 9'")
  expect_error(check_age_bins(good[-1]),
               "given: '5 to 9', '10 to 14'")
})


test_that("ll_nbinom", {
  f <- function(model) {
    set.seed(1)
    dnbinom(10, 0.1, mu = model + rexp(length(model), rate = 1e6), log = TRUE)
  }

  set.seed(1)
  expect_equal(
    ll_nbinom(10, 10, 0.1, 1e6),
    f(10))

  x <- 1:10
  set.seed(1)
  expect_equal(
    ll_nbinom(10, x, 0.1, 1e6),
    f(x))

  x <- rep(10, 10)
  expect_lt(
    diff(range(ll_nbinom(10, x, 0.1, 1e6))),
    diff(range(ll_nbinom(10, x, 0.1, 1e2))))
})


test_that("ll_nbinom returns a vector of zeros if data missing", {
  expect_equal(
    ll_nbinom(NA, 10, 0.1, 1e6),
    0)
  expect_equal(
    ll_nbinom(NA, rep(10, 5), 0.1, 1e6),
    rep(0, 5))
})

test_that("ll_binom", {
  f <- function(model_prob) {
    dbinom(5, 10, prob = model_prob, log = TRUE)
  }

  set.seed(1)
  expect_equal(
    ll_binom(5, 10, 0.2),
    f(0.2))

  x <- seq(0, 1, 0.1)
  expect_equal(
    ll_binom(5, 10, x),
    f(x))
})


test_that("ll_binom returns a vector of zeros if data missing", {
  expect_equal(
    ll_binom(NA, NA, 0.5),
    0)
  expect_equal(
    ll_binom(5, NA, 0.5),
    0)
  expect_equal(
    ll_binom(NA, 10, 0.5),
    0)
  expect_equal(
    ll_binom(NA, NA, rep(0.5, 5)),
    rep(0, 5))
  expect_equal(
    ll_binom(5, NA, rep(0.5, 5)),
    rep(0, 5))
  expect_equal(
    ll_binom(NA, 10, rep(0.5, 5)),
    rep(0, 5))
})

test_that("ll_betabinom", {
  f <- function(model_prob, log = TRUE) {
    dbetabinom(5, 10, prob = model_prob, 0.01, log)
  }

  expect_equal(
    ll_betabinom(5, 10, 0.2, 0.01),
    f(0.2))

  x <- seq(0, 1, 0.1)
  expect_equal(
    ll_betabinom(5, 10, x, 0.01),
    f(x))
})


test_that("ll_betabinom returns a vector of zeros if data missing", {
  expect_equal(
    ll_betabinom(NA, NA, 0.5, 0.01),
    0)
  expect_equal(
    ll_betabinom(5, NA, 0.5, 0.01),
    0)
  expect_equal(
    ll_betabinom(NA, 10, 0.5, 0.01),
    0)
  expect_equal(
    ll_betabinom(NA, NA, rep(0.5, 5), 0.01),
    rep(0, 5))
  expect_equal(
    ll_betabinom(5, NA, rep(0.5, 5), 0.01),
    rep(0, 5))
  expect_equal(
    ll_betabinom(NA, 10, rep(0.5, 5), 0.01),
    rep(0, 5))
})


test_that("dbetabinom", {

  ## when prob = 1/2, rho = 1/3 (equivalently a = b = 1),
  ## equivalent to discrete uniform from 0 to size
  expect_equal(dbetabinom(4, 35, 1 / 2, 1 / 3), 1 / 36)
  expect_equal(dbetabinom(4, 35, 1 / 2, 1 / 3, log = TRUE), log(1 / 36))


  ## compare with extraDistr::dbbinom (uses a and b as parameters - to match
  ## use prob = a / (a + b) and rho = 1 / (a + b + 1))
  f <- function(x, size, a, b, log = FALSE) {
   prob  <- a / (a + b)
   rho <- 1 / (a + b + 1)

   dbetabinom(x, size, prob, rho, log)
  }

  ## compare to extraDistr::dbbinom(15, 63, 2, 5) ~ 0.03613356
  expect_equal(f(15, 63, 2, 5), 0.03613356, tolerance = 1e-8)
  ## compare to extraDistr::dbbinom(15, 63, 2, 5, log = TRUE) ~ -3.320533
  expect_equal(f(15, 63, 2, 5, log = TRUE), -3.320533, tolerance = 1e-6)

  ## compare to extraDistr::dbbinom(672, 50454, 3, 2) ~ 4.18089e-08
  expect_equal(f(672, 50454, 3, 2), 4.18089e-08, tolerance = 1e-13)
  ## compare to extraDistr::dbbinom(15, 63, 2, 5, log = TRUE) ~ -16.99016
  expect_equal(f(672, 50454, 3, 2, log = TRUE), -16.99016, tolerance = 1e-5)
})

test_that("test_prob_pos", {
  ## test_prob_pos outputs between 0 and 1 even when specificity and sensitivity
  ## are 100% and there are no positive individuals
  p <- test_prob_pos(0, 5e7, 1, 1, 1e6)
  expect_true(
    p > 0 & p < 1
  )

  ## test_prob_pos outputs between 0 and 1 even when specificity and sensitivity
  ## are 100% and there are no negative individuals
  p <- test_prob_pos(5e7, 0, 1, 1, 1e6)
  expect_true(
    p > 0 & p < 1
  )
})


test_that("transmission matrix", {
  clear_cache()
  m <- sircovid_transmission_matrix("uk")
  expect_identical(m, cache$transmission_matrix$uk)
  expect_equal(dim(m), c(17, 17))

  nms <- c("[0,4)", "[5,9)", "[10,14)", "[15,19)", "[20,24)", "[25,29)",
           "[30,34)", "[35,39)", "[40,44)", "[45,49)", "[50,54)", "[55,59)",
           "[60,64)", "[65,69)", "[70,74)", "[75,79)", "[80,100)")
  expect_equal(dimnames(m), list(nms, nms))

  ## Be notified when this changes
  skip_on_cran()
  expect_equal(sum(m), 4.4980269728090259e-05)
})


test_that("transmission matrices are region specific", {
  m1 <- sircovid_transmission_matrix("uk")
  m2 <- sircovid_transmission_matrix("south_west")
  expect_false(all(abs(m1 - m2) < 1e-5))
})


test_that("read default severity", {
  clear_cache()
  d <- severity_default()
  expect_identical(d, cache$severity_default)
})


test_that("Can add new betas", {
  dat <- reference_data_mcmc()

  ## In order to be able to update future beta values we need to
  ## compute Rt for this simulation. This should probably be
  ## automated, as it's still too fiddly.
  rt <- local({
    p <- lapply(seq_len(nrow(dat$pars)), function(i)
      dat$predict$transform(dat$pars[i, ]))
    i <- grep("S_", rownames(dat$trajectories$state))
    S <- dat$trajectories$state[i, , ]
    carehomes_Rt_trajectories(dat$trajectories$step, S, p)
  })
  
  n_par <- ncol(rt$eff_Rt_all)
  set.seed(1)

  future <- list(
    "2020-04-01" = future_Rt(1.5),
    "2020-05-01" = future_Rt(0.5, "2020-03-27"),
    "2020-06-01" = future_Rt(runif(n_par, 1.4, 1.6)))
  res <- add_future_betas(dat, rt, future)
  expect_is(res, "mcstate_pmcmc")

  p_base <- dat$predict$transform(dat$pars[1, ])
  p_new <- res$predict$transform(res$pars[1, ])

  ## Nothing except beta_step wil have changed:
  expect_mapequal(p_base[names(p_base) != "beta_step"],
                  p_new[names(p_new) != "beta_step"])

  beta_base <- p_base$beta_step
  beta_new <- p_new$beta_step

  ## Unchanged up to the point where we start changing
  expect_identical(beta_new[seq_along(beta_base)], beta_base)

  ## We have a few change points to consider here:
  i <- c(sircovid_date("2020-04-01") - 1, sircovid_date("2020-04-01"),
         sircovid_date("2020-05-01") - 1, sircovid_date("2020-05-01"),
         sircovid_date("2020-06-01") - 1, sircovid_date("2020-06-01"))
  j <- i / p_new$dt

  ## Unchanged until the first change point
  expect_true(all(beta_new[seq(length(beta_base), j[[1]])] == last(beta_base)))

  ## Our parameters are expanded as expected:
  v <- beta_new[j[c(1, 2, 2, 4, 4, 6)] + 1L]
  cmp <- sircovid_parameters_beta(i, v, p_new$dt)
  expect_equal(beta_new[seq(j[[1]], length(beta_new))],
               cmp[seq(j[[1]], length(cmp))])

  tmp <- future_relative_beta(future, rt$date[, 1], rt$Rt_general)
  expect_equal(tmp$value[1, ] * last(beta_base), v)
})

test_that("Compute relative betas", {
  dat <- reference_data_mcmc()
  rt <- local({
    p <- lapply(seq_len(nrow(dat$pars)), function(i)
      dat$predict$transform(dat$pars[i, ]))
    i <- grep("S_", rownames(dat$trajectories$state))
    S <- dat$trajectories$state[i, , ]
    carehomes_Rt_trajectories(dat$trajectories$step, S, p)
  })

  ## NOTE: inject a bit of wobble to the numbers here so that it's
  ## easier to check for exact correctness, otherwise because the toy
  ## simulated data don't vary across parameter sets we don't trigger
  ## all corner cases.
  rt$Rt_general <- rt$Rt_general * rlnorm(length(rt$Rt_general), 0, 0.05)
  future <- list(
    "2020-05-01" = future_Rt(1.5, "2020-03-20"),
    "2020-06-01" = future_Rt(1.2, "2020-03-19"),
    "2020-07-01" = future_Rt(1.1, "2020-03-18"))
  res <- future_relative_beta(future, rt$date[, 1], rt$Rt_general)

  expect_equal(res$date, c(121, 122, 152, 153, 182, 183))

  np <- ncol(rt$Rt_general)

  expect_equal(res$value[, 1], rep(1, np))
  expect_equal(dim(res$value), c(np, 6))
  expect_equal(res$value[, c(2, 4)], res$value[, c(3, 5)])

  ## which leaves us with our actual comparisons to make; do we manage
  ## to set things to the correct values? We do this by an equivalent,
  ## though different, calculation to that used in the package.
  i <- match(sircovid_date(c("2020-03-20", "2020-03-19", "2020-03-18")),
             rt$date[, 1])
  r <- t(rt$Rt_general[i, ] * c(1.5, 1.2, 1.1)) /
    rt$Rt_general[nrow(rt$Rt_general), ]
  expect_equal(res$value[, c(2, 4, 6)], r)
})


test_that("strip projections", {
  dat1 <- reference_data_mcmc()
  rt1 <- calculate_rt_simple(dat1)

  dat2 <- carehomes_forecast(dat1, 0, 0, 10, NULL)
  rt2 <- calculate_rt_simple(dat2)

  future <- list(
    "2020-04-01" = future_Rt(1.5),
    "2020-05-01" = future_Rt(0.5, "2020-03-27"))

  cmp <- add_future_betas(dat1, rt1, future)
  expect_equal(add_future_betas(dat2, rt1, future)$pars, cmp$pars)
  expect_equal(add_future_betas(dat2, rt1, future)$pars, cmp$pars)
  expect_equal(add_future_betas(dat2, rt2, future)$pars, cmp$pars)
})


test_that("validate future beta values", {
  dat <- reference_data_mcmc()
  rt <- calculate_rt_simple(dat)

  f1 <- future_Rt(1.5)
  f2 <- future_Rt(1.5)

  expect_error(
    add_future_betas(dat, rt, NULL),
    "Expected at least one element in 'future'")
  expect_error(
    add_future_betas(dat, rt, list()),
    "Expected at least one element in 'future'")
  expect_error(
    add_future_betas(dat, rt, list(f1, f2)),
    "Expected 'future' to be named")
  expect_error(
    add_future_betas(dat, rt, list("2020-04-01" = f1, "2020-05-01" = 2)),
    "Expected all elements of 'future' to be 'future_Rt' objects")
  expect_error(
    add_future_betas(dat, rt, list("2020-05-01" = f1, f2)),
    "Expected ISO dates or R dates - please convert")
  expect_error(
    add_future_betas(dat, rt, list("2020-05-01" = f1, "2020-04-01" = f2)),
    "Future change dates must be increasing")
  expect_error(
    add_future_betas(dat, rt, list("2020-04-01" = f1, "2020-04-02" = f2)),
    "At least one full date required between all future change dates")
  expect_error(
    add_future_betas(dat, rt, list("2020-03-01" = f1, "2020-04-02" = f2)),
    "The first future date must be at least 2020-04-01 (but was 2020-03-01)",
    fixed = TRUE)
  expect_error(
    add_future_betas(dat, rt, list("2020-05-01" = future_Rt(1, "2020-01-01"))),
    "Relative date not found in rt set: 2020-01-01")
  expect_error(
    add_future_betas(dat, rt, list("2020-05-01" = future_Rt(1:2))),
    "Future value must be of length 1 or match rt_value")
})


test_that("Can add new betas with incomplete Rt calculation", {
  ## Due to staggered start dates, we may end up with Rt calculations
  ## that miss off the first few dates. This test checks that we can
  ## handle these appropriately.
  dat <- reference_data_mcmc()
  rt1 <- calculate_rt_simple(dat)
  rt2 <- rt1

  ## Then we remove some rt calculations:
  j <- 1:5
  for (i in names(rt2)) {
    rt2[[i]] <- rt1[[i]][-j, ]
  }

  ## Calculations work the same as the non-trimmed version
  future <- list("2020-05-01" = future_Rt(1, "2020-03-24"))
  res1 <- add_future_betas(dat, rt1, future)
  res2 <- add_future_betas(dat, rt2, future)

  expect_equal(
    res1$predict$transform(res1$pars[1, ])$beta_step,
    res2$predict$transform(res1$pars[1, ])$beta_step)

  ## Need to zero these out or the comparison fails on R-devel
  res1$predict$transform <- NULL
  res2$predict$transform <- NULL
  expect_equal(res1, res2)

  ## Error when setting to a date that has been trimmed
  date_error <- sircovid_date_as_date(rt1$date[[2]])
  expect_error(
    add_future_betas(dat, rt2, list("2020-05-01" = future_Rt(1, date_error))),
    "Relative date not found in rt set: 2020-02-29")
  expect_silent(
    add_future_betas(dat, rt1, list("2020-05-01" = future_Rt(1, date_error))))
})
