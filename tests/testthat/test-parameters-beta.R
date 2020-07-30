context("parameters (beta)")

test_that("can set constant beta", {
  expect_equal(
    sircovid_parameters_beta(NULL, 1, 0.25), 1)
  expect_error(
    sircovid_parameters_beta(0, 1, 0.25),
    "Need at least two dates and betas for a varying beta")
})


test_that("can use a varying beta", {
  expect_equal(
    sircovid_parameters_beta(c(15, 25), c(3, 2), 1),
    c(rep(3, 15), seq(3, 2, by = -0.1)))
  expect_equal(
    sircovid_parameters_beta(c(15, 25), c(3, 2), 0.5),
    c(rep(3, 30), seq(3, 2, by = -0.05)))
})


test_that("beta date and value must have same length", {
  expect_error(
    sircovid_parameters_beta(c(4, 7), c(1, 2, 3), 0.25),
    "'date' and 'value' must have the same length")
})


test_that("beta date and value must have same length", {
  expect_error(
    sircovid_parameters_beta(c(7, 4), c(1, 3), 0.25),
    "'date' must be strictly increasing")
})


test_that("beta date and value must have same length", {
  expect_error(
    sircovid_parameters_beta(as.Date(c("2020-02-01", "2020-02-20")),
                             c(1, 3), 0.25),
    "'date' must be a sircovid date")
})
