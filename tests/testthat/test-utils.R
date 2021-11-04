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


test_that("can interpolate an interval over a grid", {
  x <- 1:61
  f <- function(x) sin(x / 10)
  interpolate_grid <- function(x, f, every, min) {
    xx <- interpolate_grid_x(x, every, min)
    interpolate_grid_expand_y(f(xx), list(xx))
  }

  expect_equal(interpolate_grid(x, f, 1, min = 5), f(x))

  expect_equal(list(interpolate_grid_x(x, 5, 2)),
               interpolate_grid_critical_x(x, 5, NULL, 2))

  y <- interpolate_grid(x, f, 3, min = 5)
  expect_false(identical(y, f(x)))
  i <- rep(c(TRUE, FALSE, FALSE), length.out = length(x))
  expect_equal(y[i], f(x)[i], tolerance = 1e-15)
  expect_lt(
    mean((y[i] - f(x)[i])^2),
    mean((y[!i] - f(x)[!i])^2))

  x <- 1:61
  f <- function(x) sin(x / 10)
  y1 <- interpolate_grid(x, f, 3, min = 61)
  expect_equal(y1, f(x), tolerance = 1e-15)
  y2 <- interpolate_grid(x, f, 3, min = 59)
  expect_gt(sum(abs(y2 - f(x))), 1e-5)
  expect_equal(y2, f(x), tolerance = 1e-5)

  y <- interpolate_grid(x, f, 3, min = 30)
  i <- rep(c(TRUE, FALSE, FALSE), length.out = length(x))
  expect_lt(
    mean((y[i] - f(x)[i])^2),
    mean((y[!i] - f(x)[!i])^2))
  expect_equal(y[i], f(x)[i], tolerance = 1e-15)
})


test_that("can interpolate over an interval with critical points", {
  crit <- c(10, 31)
  f <- function(x) {
    if (x < crit[[1]]) {
      x
    } else if (x < crit[[2]]) {
      cos(x / 10)
    } else {
      sin(x / 10)
    }
  }
  g <- Vectorize(f)

  ## interpolate_grid_critical_x
  ## interpolate_grid_expand_y

  x <- 1:40
  y <- g(x)

  zx_1 <- interpolate_grid_critical_x(x, 2, crit, min = 3)
  zx_2 <- interpolate_grid_critical_x(x, 2, crit, min = 10)
  z1 <- interpolate_grid_expand_y(g(unlist(zx_1)), zx_1)
  z2 <- interpolate_grid_expand_y(g(unlist(zx_2)), zx_2)

  expect_equal(z1[1:10], y[1:10], tolerance = 1e-15)
  expect_equal(z2[1:10], y[1:10], tolerance = 1e-15)
  expect_equal(z2[31:40], y[31:40], tolerance = 1e-15)

  i <- c(seq(31, 40, by = 2), 40)
  j <- setdiff(31:40, i)
  expect_equal(z1[i], y[i], tolerance = 1e-15)
  expect_equal(z1[j], y[j], tolerance = 1e-4)
  expect_gt(sum(abs(y[j] - z1[j])), 1e-9)

  i <- c(seq(10, 30, by = 2), 40)
  j <- setdiff(10:30, i)
  expect_equal(z1[11:30], z2[11:30], tolerance = 1e-15)
  expect_equal(z1[i], y[i], tolerance = 1e-15)
  expect_equal(z1[j], y[j], tolerance = 1e-4)
  expect_gt(sum(abs(y[j] - z1[j])), 1e-9)
})


test_that("block expand", {
  m <- matrix(1:4, 2, 2)
  expect_identical(block_expand(m, 1), m)

  matrix(m %o% matrix(1, 2, 2), 4, 4)
  expect_identical(
    block_expand(m, 2),
    matrix(c(1:2, 1:2, 3:4, 3:4, 1:2, 1:2, 3:4, 3:4), 4, 4))
  expect_identical(
    block_expand(m, 3),
    matrix(rep(c(rep(1:2, 3), rep(3:4, 3)), 3), 6, 6))
  expect_identical(
    block_expand(m, 10),
    matrix(rep(c(rep(1:2, 10), rep(3:4, 10)), 10), 20, 20))
})

test_that("mirror_strain/unmirror_strain", {
  expect_error(mirror_strain(runif(3)), "length 1 or 2")
  expect_error(mirror_strain(matrix(runif(3), 1, 3)), "length 1 or 2")

  expect_equal(mirror_strain(pi), pi)
  expect_equal(mirror_strain(array(pi, c(4, 1, 6))), array(pi, c(4, 1, 6)))

  x <- runif(2)

  expect_equal(unmirror_strain(mirror_strain(x)), x)
  expect_equal(
    unmirror_strain(mirror_strain(matrix(x, ncol = 2))), matrix(x, ncol = 2)
  )
  expect_equal(
    unmirror_strain(mirror_strain(array(x, c(1, 2, 1)))), array(x, c(1, 2, 1))
  )

  x <- runif(3)
  expect_equal(unmirror_strain(x), x)
  expect_equal(unmirror_strain(matrix(x, ncol = 3)), matrix(x, ncol = 3))
  expect_equal(unmirror_strain(array(x, c(1, 1, 3))), array(x, c(1, 1, 3)))
})


test_that("check_sircovid_model", {
  expect_silent(check_sircovid_model("basic"))
  expect_silent(check_sircovid_model("lancelot"))
  expect_error(check_sircovid_model("camelot"),
              "Expected 'camelot' to be one of {basic, lancelot}",
                fixed = TRUE)
  expect_error(check_sircovid_model(c("basic", "lancelot")),
              "'sircovid_model' must be a scalar")
})
