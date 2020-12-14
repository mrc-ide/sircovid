context("simulate")

test_that("Can construct an empty object", {
  expect_equal(
    sircovid_simulate_events("2020-03-01", "2020-10-10", NULL),
    structure(list(date_from = 61L, date_to = 284, data = list(NULL)),
              class = "sircovid_events"))
})


test_that("can construct simple events objects", {
  expect_equal(
    sircovid_simulate_events("2020-03-01", "2020-10-10",
                             list("2020-04-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    structure(list(date_from = c(61, 92, 122),
                   date_to = c(92, 122, 284),
                   data = list(NULL, list(a = 1), list(b = 2))),
              class = "sircovid_events"))
})


test_that("can construct events with first event before start", {
  expect_equal(
    sircovid_simulate_events("2020-04-10", "2020-10-10",
                             list("2020-04-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    structure(list(date_from = c(101, 122),
                   date_to = c(122, 284),
                   data = list(list(a = 1), list(b = 2))),
              class = "sircovid_events"))
  expect_equal(
    sircovid_simulate_events("2020-05-10", "2020-10-10",
                             list("2020-04-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    structure(list(date_from = 131,
                   date_to = 284,
                   data = list(list(b = 2))),
              class = "sircovid_events"))
})


test_that("can construct events that end before last event", {
  expect_equal(
    sircovid_simulate_events("2020-03-01", "2020-03-10",
                             list("2020-04-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    structure(list(date_from = 61,
                   date_to = 70,
                   data = list(NULL)),
              class = "sircovid_events"))

  expect_equal(
    sircovid_simulate_events("2020-03-01", "2020-04-10",
                             list("2020-04-01" = list(a = 1),
                                  "2020-05-01" = list(b = 2))),
    structure(list(date_from = c(61, 92),
                   date_to = c(92, 101),
                   data = list(NULL, list(a = 1))),
              class = "sircovid_events"))
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
