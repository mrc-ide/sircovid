context("pmcmc")

test_that("adding incidence adds appropriate states", {
  dat <- reference_data_mcmc()
  res <- calculate_incidence(dat$trajectories, c("deaths", "deaths_hosp"))
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
  cmp <- calculate_incidence(dat$trajectories, c("deaths", "deaths_hosp"))
  res <- calculate_incidence(dat$trajectories, "deaths")
  expect_identical(res$state["deaths_inc", , ],
                   cmp$state["deaths_inc", , ])
})
