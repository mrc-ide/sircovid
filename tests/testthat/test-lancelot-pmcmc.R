context("pmcmc")

test_that("adding no incidence leaves object the same", {
  dat <- reference_data_lancelot_mcmc()
  res <- add_trajectory_incidence(dat$trajectories, NULL)
  expect_identical(res, dat$trajectories)
})


test_that("adding incidence adds appropriate states", {
  dat <- reference_data_lancelot_mcmc()
  res <- add_trajectory_incidence(dat$trajectories, "icu")
  expect_true("icu_inc" %in% rownames(res$state))

  tmp <- res$state["icu_inc", , ]
  expect_true(all(is.na(tmp[, 1:2])))
  icu <- t(apply(tmp[, -c(1, 2)], 1, cumsum))
  expect_equal(
    icu,
    res$state["icu", , -c(1, 2)] - res$state["icu", , 2])

  expect_equal(drop_trajectory_incidence(res),
               drop_trajectory_incidence(dat$trajectories))
})


test_that("can add and remove trajectories from mcstate_pmcmc objects", {
  dat <- reference_data_lancelot_mcmc()
  v <- c("deaths", "icu")
  res <- add_trajectory_incidence(dat, v)
  expect_identical(res$trajectories,
                   add_trajectory_incidence(dat$trajectories, v))
  expect_identical(drop_trajectory_incidence(res),
                   drop_trajectory_incidence(dat))
})


test_that("can compute incidence for a single variable", {
  dat <- reference_data_lancelot_mcmc()
  cmp <- add_trajectory_incidence(dat$trajectories, c("deaths", "icu"))
  res <- add_trajectory_incidence(dat$trajectories, "deaths")
  expect_identical(res$state["deaths_inc", , ],
                   cmp$state["deaths_inc", , ])
})

test_that("Can combine trajectories of equal size - rank = FALSE", {
  dat <- reference_data_lancelot_trajectories()
  res <- combine_trajectories(list(dat, dat), rank = FALSE)
  expect_equal(res$step, dat$trajectories$step)
  expect_equal(res$rate, dat$trajectories$rate)
  expect_equal(res$predicted, dat$trajectories$predicted)
  expect_equal(res$date, dat$trajectories$date)
  expect_equal(res$state, dat$trajectories$state * 2)
})


test_that("Can combine trajectories of equal size", {
  dat <- reference_data_lancelot_trajectories()
  res <- combine_trajectories(list(dat, dat))
  expect_equal(res$step, dat$trajectories$step)
  expect_equal(res$rate, dat$trajectories$rate)
  expect_equal(res$predicted, dat$trajectories$predicted)
  expect_equal(res$date, dat$trajectories$date)

  expect_equal(sum(res$state), sum(dat$trajectories$state * 2))
  expect_equal(apply(res$state, c(1, 3), sum),
               apply(dat$trajectories$state, c(1, 3), sum) * 2)
  ## TODO: Lilith to check that the trajectories increase over the
  ## particle index.
})


test_that("Can combine trajectories with missing times", {
  dat <- reference_data_lancelot_trajectories()

  ## Create a set where the first "real" data point is missing,
  ## shifting predictions back one day but keeping the same total
  ## number of forecast days.
  err <- dat
  n <- length(err$trajectories$date)
  err$trajectories$predicted <- err$trajectories$predicted[-1]
  err$trajectories$date <- err$trajectories$date[-n]
  err$trajectories$state <- err$trajectories$state[, , -n] * 2

  res <- combine_trajectories(list(dat, err))

  ## Predictions have shifted, so this is the same total sum
  expect_equal(sum(res$predicted), sum(dat$trajectories$predicted))
  ## But one less real data point
  expect_equal(sum(!res$predicted), sum(!dat$trajectories$predicted) - 1L)
  expect_equal(res$predicted, err$trajectories$predicted)
  expect_equal(res$date, err$trajectories$date)
  expect_equal(dim(res$state), dim(err$trajectories$state))

  f <- function(x) apply(x, c(1, 3), sum)
  tmp <- dat$trajectories$state[, , -n] + err$trajectories$state
  expect_equal(f(res$state), f(tmp))
})


test_that("can combine rt calculations over trajectories", {
  dat <- reference_data_lancelot_trajectories()

  index_S <- grep("^S_", names(dat$predict$index))
  S <- dat$trajectories$state[index_S, , , drop = FALSE]
  pars <- lapply(seq_len(nrow(dat$pars)), function(i)
    dat$predict$transform(dat$pars[i, ]))
  rt <- lancelot_Rt_trajectories(
    dat$trajectories$step, S, pars,
    initial_step_from_parameters = TRUE,
    shared_parameters = FALSE)

  res <- combine_rt(list(rt, rt), list(dat, dat))
  cmp <- rt
  for (i in setdiff(names(cmp), c("step", "date"))) {
    cmp[[i]][1, ] <- NA
  }
  ## This should pass for everything, but the effective Rt
  ## calculations are different after aggregation for reasons as-yet
  ## unexplained. Possibly this is reasonable based on how Rt is
  ## calculated.
  expect_equal(res$Rt_general, cmp$Rt_general)
  expect_equal(res$Rt_all, cmp$Rt_all)
})


test_that("when combining rt calculations the output has the expected order", {
  dat <- reference_data_lancelot_trajectories()

  index_S <- grep("^S_", names(dat$predict$index))
  S <- dat$trajectories$state[index_S, , , drop = FALSE]
  pars <- lapply(seq_len(nrow(dat$pars)), function(i)
    dat$predict$transform(dat$pars[i, ]))
  rt <- lancelot_Rt_trajectories(
    dat$trajectories$step, S, pars,
    initial_step_from_parameters = TRUE,
    shared_parameters = FALSE)

  res <- combine_rt(list(rt, rt), list(dat, dat))
  cmp <- rt
  for (i in setdiff(names(cmp), c("step", "date"))) {
    cmp[[i]][1, ] <- NA
  }

  for (what in c("eff_Rt_all", "eff_Rt_general", "Rt_all", "Rt_general")) {
    cum_rt <- colSums(rt[[what]])
    expected_order <- order(cum_rt, decreasing = FALSE)
    expect_equal(res[[what]], cmp[[what]][, expected_order])
  }

})


test_that("can combine rt calculations over trajectories without reordering", {
  dat <- reference_data_lancelot_trajectories()

  index_S <- grep("^S_", names(dat$predict$index))
  S <- dat$trajectories$state[index_S, , , drop = FALSE]
  pars <- lapply(seq_len(nrow(dat$pars)), function(i)
    dat$predict$transform(dat$pars[i, ]))
  rt <- lancelot_Rt_trajectories(
    dat$trajectories$step, S, pars,
    initial_step_from_parameters = TRUE,
    shared_parameters = FALSE)

  res <- combine_rt(list(rt, rt), list(dat, dat), rank = FALSE)
  cmp <- rt
  for (i in setdiff(names(cmp), c("step", "date"))) {
    cmp[[i]][1, ] <- NA
  }

  expect_equal(res$Rt_general, cmp$Rt_general)
  expect_equal(res$Rt_all, cmp$Rt_all)
  expect_equal(res$eff_Rt_general, cmp$eff_Rt_general)
  expect_equal(res$eff_Rt_all, cmp$eff_Rt_all)
})


test_that("adding incidence adds appropriate states - nested", {
  dat <- reference_data_lancelot_mcmc()
  dat$trajectories$state <- array(
    dat$trajectories$state, c(158, 11, 2, 32),
    dimnames = c(list(dimnames(dat$trajectories$state)[[1]], NULL,
                      letters[1:2], NULL)))
  res <- add_trajectory_incidence(dat$trajectories, "icu")
  expect_true("icu_inc" %in% rownames(res$state))

  tmp <- res$state["icu_inc", , , ]
  expect_true(all(is.na(tmp[, 2, 1:2])))
  icu <- t(apply(tmp[, 2, -c(1, 2)], 1, cumsum))
  expect_equal(
    icu,
    res$state["icu", , 2,  -c(1, 2)] -
      res$state["icu", , 2,  2])

  expect_equal(drop_trajectory_incidence(res),
               drop_trajectory_incidence(dat$trajectories))
})



test_that("add and remove trajectories from nested mcstate_pmcmc objects", {
  dat <- reference_data_lancelot_mcmc()
  dat$trajectories$state <- array(
    dat$trajectories$state, c(158, 11, 2, 32),
    dimnames = c(list(dimnames(dat$trajectories$state)[[1]], NULL,
                      letters[1:2], NULL)))
  v <- "icu"
  res <- add_trajectory_incidence(dat, v)
  expect_identical(res$trajectories,
                   add_trajectory_incidence(dat$trajectories, v))
  expect_identical(drop_trajectory_incidence(res),
                   drop_trajectory_incidence(dat))
})


test_that("can compute incidence for a single variable - nested", {
  dat <- reference_data_lancelot_mcmc()
  dat$trajectories$state <- array(
    dat$trajectories$state, c(158, 11, 2, 32),
    dimnames = c(list(dimnames(dat$trajectories$state)[[1]], NULL,
                      letters[1:2], NULL)))
  cmp <- add_trajectory_incidence(dat$trajectories, "icu")
  res <- add_trajectory_incidence(dat$trajectories, "icu")
  expect_identical(res$state["icu_inc", , , ],
                   cmp$state["icu_inc", , , ])
})


test_that("can combine EpiEstim rt calculations over trajectories", {
  dat <- reference_data_lancelot_trajectories()

  index_inc <- grep("infections_inc", names(dat$predict$index))
  inc <- dat$trajectories$state[index_inc, , ]

  p <- dat$predict$transform(dat$pars[1, ])

  ## Fix p_C across age groups for the rest of the test
  p$p_C <- rep(0.6, 19)

  rt <- lancelot_rt_trajectories_epiestim(
    dat$trajectories$step, inc, p)

  res <- combine_rt_epiestim(list(rt, rt), list(dat, dat))
  cmp <- rt
  cmp$Rt[, 1] <- NA
  cmp$Rt_summary[, 1] <- NA

  ## Rts are ordered differently
  expect_equal(sort(as.vector(res$Rt)), sort(as.vector(cmp$Rt)))
  expect_equal(res$Rt_summary, cmp$Rt_summary)
})


test_that("can combine EpiEstim rt over trajectories without reordering", {
  dat <- reference_data_lancelot_trajectories()

  index_inc <- grep("infections_inc", names(dat$predict$index))
  inc <- dat$trajectories$state[index_inc, , ]

  p <- dat$predict$transform(dat$pars[1, ])

  ## Fix p_C across age groups for the rest of the test
  p$p_C <- rep(0.6, 19)

  rt <- lancelot_rt_trajectories_epiestim(
    dat$trajectories$step, inc, p)

  res <- combine_rt_epiestim(list(rt, rt), list(dat, dat), rank = FALSE)
  cmp <- rt
  cmp$Rt[, 1] <- NA
  cmp$Rt_summary[, 1] <- NA

  ## Rts are ordered differently
  expect_equal(res$Rt, cmp$Rt)
  expect_equal(res$Rt_summary, cmp$Rt_summary)
})


test_that("Combining EpiEstim rt reject invalid inputs", {
  dat <- reference_data_lancelot_trajectories()

  index_inc <- grep("infections_inc", names(dat$predict$index))
  inc <- dat$trajectories$state[index_inc, , ]

  p <- dat$predict$transform(dat$pars[1, ])

  ## Fix p_C across age groups for the rest of the test
  p$p_C <- rep(0.6, 19)

  rt <- lancelot_rt_trajectories_epiestim(
    dat$trajectories$step, inc, p, save_all_Rt_sample = FALSE)

  error_msg <-
    paste("rt$Rt missing. Did you forget 'save_all_Rt_sample = TRUE'",
          "in 'lancelot_EpiEstim_Rt_trajectories'?")
  expect_error(combine_rt_epiestim(list(rt, rt), list(dat, dat)),
               error_msg, fixed = TRUE)

})


test_that("get_sample_rank rejects invalid inputs", {
  dat <- reference_data_lancelot_trajectories()

  expect_error(get_sample_rank(dat$state, by = "not_the_right_thing"),
               "'sample' should be an 'mcstate_pmcmc' object",
               fixed = TRUE)

  expect_error(get_sample_rank(dat, by = "not_the_right_thing"),
               "Unkwnown 'by' argument. Should be one of: ",
               fixed = TRUE)
})


test_that("get_sample_rank returns expected output", {
  dat <- reference_data_lancelot_trajectories()

  dim3 <- dim(dat$trajectories$state)[3]
  expected <- order(dat$trajectories$state["infections", , dim3])
  res1 <- get_sample_rank(dat, by = "infections")

  expect_equal(expected, res1)
})


test_that("reorder_sample rejects invalid inputs", {
  dat <- reference_data_lancelot_trajectories()

  expect_error(reorder_sample(dat$state, 1:3),
               "'sample' should be an 'mcstate_pmcmc' object",
               fixed = TRUE)

  expect_error(reorder_sample(dat, 1:50),
               "Unexpected length for 'rank': 50 ; should have length 10",
               fixed = TRUE)
})


test_that("reorder_sample returns expected output", {
  dat <- reference_data_lancelot_trajectories()

  ## maintaining the initial order returns the same as input
  expect_equal(reorder_sample(dat, 1:10), dat)

  ## ordering and then ordering back returns the same as input
  # order by increasing cumulative incidence
  rnk1 <- get_sample_rank(dat, by = "infections")
  dat2 <- reorder_sample(dat, rnk1)
  rnk2 <- get_sample_rank(dat2, by = "infections")
  # check this is now ordered by increasing incidence
  expect_equal(rnk2, 1:10)
  # apply revert ordering
  dat3 <- reorder_sample(dat2, match(1:10, rnk1))
  # check we are back to initial object
  expect_equal(dat3, dat)

})


test_that("reorder_sample tolerates multiple chains", {
  dat <- reference_data_lancelot_trajectories()
  dat$chain <- rep(1:2, each = 5)
  dat$iteration <- rep(1:5, 2)

  ## maintaining the initial order returns the same as input
  expect_equal(reorder_sample(dat, 1:10), dat)

  ## ordering and then ordering back returns the same as input
  # order by increasing cumulative incidence
  rnk1 <- get_sample_rank(dat, by = "infections")
  dat2 <- reorder_sample(dat, rnk1)
  expect_equal(get_sample_rank(dat2, by = "infections"), 1:10)
  expect_equal(dat2$chain, dat$chain[rnk1])
})


test_that("reorder_rt_ifr rejects invalid inputs", {
  dat <- reference_data_lancelot_trajectories()

  index_S <- grep("^S_", names(dat$predict$index))
  S <- dat$trajectories$state[index_S, , , drop = FALSE]
  pars <- lapply(seq_len(nrow(dat$pars)), function(i)
    dat$predict$transform(dat$pars[i, ]))
  rt <- lancelot_Rt_trajectories(
    dat$trajectories$step, S, pars,
    initial_step_from_parameters = TRUE,
    shared_parameters = FALSE)

  expect_error(reorder_rt_ifr(rt$beta, 1:10),
               "'x' should be an 'Rt_trajectories' or 'IFR_t_trajectories'",
               fixed = TRUE)

  expect_error(reorder_rt_ifr(rt, 1:50),
               "Unexpected length for 'rank': 50 ; should have length 10",
               fixed = TRUE)

})


test_that("reorder_rt_ifr returns expected output", {
  dat <- reference_data_lancelot_trajectories()

  index_S <- grep("^S_", names(dat$predict$index))
  S <- dat$trajectories$state[index_S, , , drop = FALSE]
  pars <- lapply(seq_len(nrow(dat$pars)), function(i)
    dat$predict$transform(dat$pars[i, ]))
  rt <- lancelot_Rt_trajectories(
    dat$trajectories$step, S, pars,
    initial_step_from_parameters = TRUE,
    shared_parameters = FALSE)

  ## maintaining the initial order returns the same as input
  expect_equal(reorder_rt_ifr(rt, 1:10), rt)

  ## ordering and then ordering back returns the same as input
  # order by increasing cumulative incidence
  rnk1 <- get_sample_rank(dat, by = "infections")
  rt2 <- reorder_rt_ifr(rt, rnk1)
  # apply revert ordering
  rt3 <- reorder_rt_ifr(rt2, match(1:10, rnk1))
  # check we are back to initial object
  expect_equal(rt3, rt)

})
