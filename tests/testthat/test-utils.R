context("utils")

test_that("null-or-value works", {
  expect_equal(1 %||% NULL, 1)
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% NULL, NULL)
  expect_equal(NULL %||% 2, 2)
})


test_that("Can singly quote strings", {
  expect_equal(squote("a"), "'a'")
  expect_equal(squote(c("apple", "banana")),
               c("'apple'", "'banana'"))
  expect_equal(squote(character(0)), character(0))
})


test_that("name validation identifies unexpected elements", {
  x <- list(a = 1, b = 2)
  expect_silent(verify_names(x, allow_extra = TRUE))
  expect_silent(verify_names(x, required = "a", optional = "b"))
  expect_silent(verify_names(x, required = c("a", "b")))
  expect_silent(verify_names(x, optional = c("a", "b", "c")))
  expect_error(verify_names(x),
               "Extra elements in 'x': 'a', 'b'")
})


test_that("name validation prevents duplicated entries", {
  x <- list(a = 1, b = 2, a = 3, b = 4, c = 5)
  x13 <- x[1:3]
  x14 <- x[1:4]
  expect_silent(verify_names(x[1:2], allow_extra = TRUE))
  expect_error(verify_names(x13, allow_extra = TRUE),
               "Duplicate element names in 'x13': 'a'")
  expect_error(verify_names(x14, allow_extra = TRUE),
               "Duplicate element names in 'x14': 'a', 'b'")
  expect_error(verify_names(x, allow_extra = TRUE),
               "Duplicate element names in 'x': 'a', 'b'")
})


test_that("name validation identifies missing required elements", {
  x <- list(a = 1, b = 2)
  expect_silent(verify_names(x, c("a", "b")))
  expect_error(verify_names(x, c("a", "b", "c")),
               "Elements missing from 'x': 'c'")
  expect_error(verify_names(x, c("a", "c", "d"), allow_extra = TRUE),
               "Elements missing from 'x': 'c', 'd'")
})


test_that("data_frame and read_csv do not do factor conversion", {
  d <- data_frame(a = letters, b = LETTERS, c = 1:26)
  expect_type(d$a, "character")
  expect_type(d$b, "character")
  path <- tempfile()
  on.exit(unlink(path))
  write_csv(d, path)
  expect_equal(read_csv(path), d)
})


test_that("is_integer can detect integers", {
  expect_true(is_integer(1.0))
  expect_true(is_integer(1.0 + .Machine$double.eps))
  expect_false(is_integer(1.1))
})


test_that("sircovid_file throws for missing files", {
  expect_true(file.exists(sircovid_file("odin/basic.R")))
  ## NOTE: not testing error string because it comes from base R
  expect_error(sircovid_file("odin/acidic.R"))
})


test_that("rename can rename names", {
  x <- list(a = 1, b = 2, c = 3, d = 4)
  expect_equal(
    rename(x, character(), character()),
    x)
  expect_equal(
    rename(x, "a", "A"),
    list(A = 1, b = 2, c = 3, d = 4))
  expect_equal(
    rename(x, c("a", "c"), c("A", "C")),
    list(A = 1, b = 2, C = 3, d = 4))
  expect_error(
    rename(x, c("x", "y"), c("X", "Y")),
    "Elements missing from 'x': 'x', 'y'")
})


test_that("can convert and check R's dates", {
  expect_s3_class(as_date("2020-05-02"), "Date")
  expect_true(is_date(as_date("2020-05-02")))
  expect_true(is_date(Sys.Date()))
  expect_false(is_date("2020-05-02"))
  date <- as_date("2020-05-02")
  expect_identical(as_date(date), date)
  expect_error(
    as_date("02-05-2020"),
    "Expected ISO dates or R dates - please convert")
})


test_that("matrix_pow calculates powers correctly", {

  set.seed(1)
  A <- array(runif(4 ^ 2, -2, 2), dim = c(4, 4))

  expect_error(matrix_pow(A, 0.5), "n must be an integer")
  expect_error(matrix_pow(A[1:2, ], 2), "x must be a square matrix")
  expect_error(matrix_pow(A[, 1:3], 2), "x must be a square matrix")
  expect_error(matrix_pow(array(A, c(4, 4, 4)), 2), "x must have 2 dimensions")

  expect_equal(solve(A), matrix_pow(A, -1))
  expect_equal(solve(A) %*% solve(A), matrix_pow(A, -2))
  expect_equal(diag(ncol(A)), matrix_pow(A, 0))
  expect_equal(A, matrix_pow(A, 1))
  expect_equal(A %*% A, matrix_pow(A, 2))
  expect_equal(A %*% A %*% A, matrix_pow(A, 3))
  expect_equal(A %*% A %*% A %*% A, matrix_pow(A, 4))

})


test_that("spread_integer", {
  expect_equal(spread_integer(12, 5), c(3, 3, 2, 2, 2))
  expect_equal(spread_integer(15, 5), rep(3, 5))
  expect_equal(spread_integer(0, 5), rep(0, 5))
  expect_equal(spread_integer(2, 5), c(1, 1, 0, 0, 0))
})


test_that("Gamma correction works", {
  ## tests that adjusted_gamma improves match to mean of continuous distribution
  
  ## compute the mean of a discretized Gamma/Erlang distribution
  discretized_mean <- function(gamma, k, dt) {
    max_t <- qgamma(0.999, shape = k, rate = gamma)
    tt <- seq(0, max_t + 1, dt)
    cdf <- pgamma(tt, shape = k, rate = gamma)
    sum(tt[-length(tt)] * diff(cdf))
  }
  
  test_gamma_correction <- function(gamma, k, dt)
  {
    error_no_correction <- abs(discretized_mean(gamma, k, dt) - k / gamma)
    error_with_correction <- abs(discretized_mean(adjusted_gamma(gamma, k, dt), k, dt) - k / gamma)
    expect_true(error_no_correction > error_with_correction)
  }
  
  test_gamma_correction(gamma = 1 / 2.5, k = 1, dt = 1)
  test_gamma_correction(gamma = 1 / 2.5, k = 1, dt = 1/4)
  test_gamma_correction(gamma = 1 / 2.5, k = 1, dt = 1/100)
  
  test_gamma_correction(gamma = 1 / 2.5, k = 2, dt = 1)
  test_gamma_correction(gamma = 1 / 2.5, k = 2, dt = 1/4)
  test_gamma_correction(gamma = 1 / 2.5, k = 2, dt = 1/100)
  
  test_gamma_correction(gamma = 1 / 1, k = 1, dt = 1)
  test_gamma_correction(gamma = 1 / 1, k = 1, dt = 1/4)
  test_gamma_correction(gamma = 1 / 1, k = 1, dt = 1/100)
  
  test_gamma_correction(gamma = 1 / 1, k = 2, dt = 1)
  test_gamma_correction(gamma = 1 / 1, k = 2, dt = 1/4)
  test_gamma_correction(gamma = 1 / 1, k = 2, dt = 1/100)
  
  test_gamma_correction(gamma = 1 / 10, k = 1, dt = 1)
  test_gamma_correction(gamma = 1 / 10, k = 1, dt = 1/4)
  test_gamma_correction(gamma = 1 / 10, k = 1, dt = 1/100)
  
  test_gamma_correction(gamma = 1 / 10, k = 2, dt = 1)
  test_gamma_correction(gamma = 1 / 10, k = 2, dt = 1/4)
  test_gamma_correction(gamma = 1 / 10, k = 2, dt = 1/100)
})
