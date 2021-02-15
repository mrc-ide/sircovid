context("simulate")

test_that("Can construct an empty object", {
  expect_equal(
    sircovid_simulate_events("2020-03-01", "2020-10-10", NULL),
    structure(list(date_from = 61L, date_to = 284, data = list(NULL)),
              class = "sircovid_simulate_events"))
})


test_that("can construct simple events objects", {
  expect_equal(
    sircovid_simulate_events("2020-03-01", "2020-10-10",
                             list("2020-04-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    structure(list(date_from = c(61, 92, 122),
                   date_to = c(92, 122, 284),
                   data = list(NULL, list(a = 1), list(b = 2))),
              class = "sircovid_simulate_events"))
})


test_that("can construct events with first event before start", {
  expect_equal(
    sircovid_simulate_events("2020-04-10", "2020-10-10",
                             list("2020-04-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    structure(list(date_from = c(101, 122),
                   date_to = c(122, 284),
                   data = list(list(a = 1), list(b = 2))),
              class = "sircovid_simulate_events"))
  expect_equal(
    sircovid_simulate_events("2020-05-10", "2020-10-10",
                             list("2020-04-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    structure(list(date_from = 131,
                   date_to = 284,
                   data = list(list(b = 2))),
              class = "sircovid_simulate_events"))
})


test_that("can construct events that end before last event", {
  expect_equal(
    sircovid_simulate_events("2020-03-01", "2020-03-10",
                             list("2020-04-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    structure(list(date_from = 61,
                   date_to = 70,
                   data = list(NULL)),
              class = "sircovid_simulate_events"))

  expect_equal(
    sircovid_simulate_events("2020-03-01", "2020-04-10",
                             list("2020-04-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    structure(list(date_from = c(61, 92),
                   date_to = c(92, 101),
                   data = list(NULL, list(a = 1))),
              class = "sircovid_simulate_events"))
})


test_that("can avoid bad events", {
  expect_error(
    sircovid_simulate_events("2021-03-01", "2020-10-10",
                             list("2020-04-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    "'date_start' must be less than 'date_end'")
  expect_error(
    sircovid_simulate_events("2020-03-01", "2020-03-01",
                             list("2020-04-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    "'date_start' must be less than 'date_end'")

  expect_error(
    sircovid_simulate_events("2020-03-01", "2020-10-10",
                             list("2020-06-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    "The dates used as 'data' names must be strictly increasing")
  expect_error(
    sircovid_simulate_events("2020-03-01", "2020-10-10",
                             list("2020-05-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    "The dates used as 'data' names must be strictly increasing")
})


test_that("Basic simulation", {
  p <- basic_parameters(sircovid_date("2020-02-07"), "london")
  np <- 10
  info <- basic$new(p, 0, np)$info()
  state <- matrix(rep(basic_initial(info, np, p)$state, np), ncol = np)

  ## This is not how we'd normally model beta, but it will do here:
  events <- sircovid_simulate_events(
    "2020-02-07", "2020-06-01",
    list("2020-04-01" = list(beta_step = 0)))
  p_base <- rep(list(p), np)

  res <- sircovid_simulate(basic, state, p_base, events, seed = 1L)
  expect_equal(dim(attr(res, "state")), dim(state))
  expect_length(attr(res, "rng_state"), np * 32)

  ## Change the beta directly in the model and we should see this agree:
  p_cmp <- basic_parameters(
    sircovid_date("2020-02-07"), "london",
    beta_date = sircovid_date(c("2020-01-01", "2020-03-31", "2020-04-01")),
    beta_value = c(p$beta_step, p$beta_step, 0))
  p_cmp$beta_step[seq_len(length(p_cmp$beta_step) - 1)] <- p$beta_step
  step <- attr(res, "step")

  obj <- basic$new(rep(list(p_cmp), np), step[[1]], 1L,
                   seed = 1L, pars_multi = TRUE)
  obj$set_state(array(state, c(nrow(state), 1, ncol(state))))
  cmp <- obj$simulate(step)
  dim(cmp) <- dim(cmp)[-2]
  expect_equal(res, cmp, check.attributes = FALSE)

  expect_equal(attr(res, "state"), array(obj$state(), dim(state)))
  expect_equal(attr(res, "rng_state"), obj$rng_state())
})


test_that("validate base parameters", {
  p <- basic_parameters(sircovid_date("2020-02-07"), "london")
  np <- 10
  info <- basic$new(p, 0, np)$info()
  state <- matrix(rep(basic_initial(info, np, p)$state, np), ncol = np)
  ## This is not how we'd normally model beta, but it will do here:
  events <- sircovid_simulate_events(
    "2020-02-07", "2020-06-01",
    list("2020-04-01" = list(dt = 0.1)))
  p_base <- rep(list(p), np)

  expect_error(
    sircovid_simulate(basic, state, NULL, events),
    "Expected 'p' to be an unnamed list")
  expect_error(
    sircovid_simulate(basic, state, p, events),
    "Expected 'p' to be an unnamed list")

  expect_error(
    sircovid_simulate(basic, state, p_base, events),
    "Events must not change 'dt'")

  p_base[[5]]$dt <- 1
  expect_error(
    sircovid_simulate(basic, state, p_base, events),
    "All entries in 'p' must have the same value of 'dt'")
})
