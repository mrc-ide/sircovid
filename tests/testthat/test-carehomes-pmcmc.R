context("pmcmc")

test_that("adding incidence adds appropriate states", {
  dat <- reference_data_mcmc()
  res <- add_trajectory_incidence(dat$trajectories, c("deaths", "deaths_hosp"))
  expect_true(all(c("deaths_inc", "deaths_hosp_inc") %in% rownames(res$state)))

  tmp <- res$state["deaths_inc", , ]
  expect_true(all(is.na(tmp[, 1:2])))
  deaths <- t(apply(tmp[, -(1:2)], 1, cumsum))
  expect_equal(
    deaths,
    res$state["deaths", , -(1:2)] - res$state["deaths", , 2])
})


test_that("can compute incidence for a single variable", {
  dat <- reference_data_mcmc()
  cmp <- add_trajectory_incidence(dat$trajectories, c("deaths", "deaths_hosp"))
  res <- add_trajectory_incidence(dat$trajectories, "deaths")
  expect_identical(res$state["deaths_inc", , ],
                   cmp$state["deaths_inc", , ])
})


test_that("Can compute forecasts from mcmc output", {
  dat <- reference_data_mcmc()
  res <- carehomes_forecast(dat, 3, 5, 10, c("deaths", "deaths_hosp"))

  expect_equal(dim(res$pars), c(3, 2))
  expect_equal(dim(res$probabilities), c(3, 3))
  expect_equal(dim(res$state), c(nrow(dat$state), 3))
  expect_equal(dim(res$trajectories$state),
               dim(dat$trajectories$state) + c(2, -8, 10))

  expect_true(all(c("deaths_inc", "deaths_hosp_inc") %in%
                  rownames(res$trajectories$state)))
})


test_that("Can compute forecasts from mcmc output without prepending", {
  dat <- reference_data_mcmc()
  res <- carehomes_forecast(dat, 3, 5, 10, c("deaths", "deaths_hosp"),
                            FALSE)

  expect_equal(dim(res$pars), c(3, 2))
  expect_equal(dim(res$probabilities), c(3, 3))
  expect_equal(dim(res$state), c(nrow(dat$state), 3))
  expect_equal(dim(res$trajectories$state),
               c(nrow(dat$trajectories$state) + 2, 3, 11))
  expect_true(all(c("deaths_inc", "deaths_hosp_inc") %in%
                  rownames(res$trajectories$state)))
})
