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


test_that("check_rel_susceptibility rejects out of bounds errors", {
  expect_error(
    check_rel_susceptibility(NULL),
    "At least one value required for 'rel_susceptibility'")
  expect_error(
    check_rel_susceptibility(-1),
    "All values of 'rel_susceptibility' must lie in [0, 1]",
    fixed = TRUE)
  expect_error(
    check_rel_susceptibility(1.1),
    "All values of 'rel_susceptibility' must lie in [0, 1]",
    fixed = TRUE)
  expect_error(
    check_rel_susceptibility(c(1, 1.8)),
    "All values of 'rel_susceptibility' must lie in [0, 1]",
    fixed = TRUE)
  expect_error(
    check_rel_susceptibility(c(0.9, 0.8)),
    "First value of 'rel_susceptibility' must be 1")
})


test_that("check_rel_susceptibility allows sensible inputs", {
  expect_silent(
    check_rel_susceptibility(c(1, 0.5, 0.7)))
  expect_silent(
    check_rel_susceptibility(c(1, 0.7, 0.5)))
  expect_silent(
    check_rel_susceptibility(1))
  expect_silent(
    check_rel_susceptibility(c(1, 1)))
  expect_silent(
    check_rel_susceptibility(c(1, 0)))
  expect_silent(
    check_rel_susceptibility(c(1, 0, 1)))
})

test_that("build_rel_susceptibility rejects objects of the wrong dimension", {
  expect_error(
    build_rel_susceptibility(
    rel_susceptibility = matrix(c(1, 0.5, 1, 0.7), nrow = 2, byrow = TRUE),
    N_age = 3),
    "'rel_susceptibility' should have as many rows as age groups")
})

test_that("build_vaccination_rate rejects insensible inputs", {
  expect_error(
    build_vaccination_rate(vaccination_rate = c(1, 1), N_age = 3),
    "'vaccination_rate' should have as many elements as age groups")
  expect_error(
    build_vaccination_rate(vaccination_rate = c(1, 1, -1), N_age = 3),
    "'vaccination_rate' must have only non-negative values")
})

test_that("build_vaccination_rate allows sensible inputs & works as expected", {
  expect_silent(
    build_vaccination_rate(vaccination_rate = c(1, 1), N_age = 2))
  expect_silent(
    build_vaccination_rate(vaccination_rate = 1, N_age = 2))
  expect_equal(build_vaccination_rate(vaccination_rate = 1, N_age = 2), c(1, 1))
})


test_that("build_vaccine_progression_rate rejects insensible inputs", {
  expect_error(
    build_vaccine_progression_rate(vaccine_progression_rate = c(1, 1),
                                   N_age = 5, N_vacc_classes = 1),
    "'N_vacc_classes' should be at least 3")
  expect_error(
    build_vaccine_progression_rate(vaccine_progression_rate = 
                                     rbind(c(-1, 1), c(-1, 1)),
                                   N_age = 2, N_vacc_classes = 4),
  "'vaccine_progression_rate' must have only non-negative values")
  msg1 <- "'vaccine_progression_rate' must be either:"
  msg2 <- "a vector of length 'N_vacc_classes - 2'"
  msg3 <- "or a matrix with 'N_age' rows and 'N_vacc_classes - 2' columns"
  expect_error(
    build_vaccine_progression_rate(vaccine_progression_rate = c(1, 1),
                                   N_age = 5, N_vacc_classes = 3),
    paste(msg1, msg2, msg3))
  expect_error(
    build_vaccine_progression_rate(vaccine_progression_rate = c(1, 1, -1),
                                   N_age = 10, N_vacc_classes = 5),
    "'vaccine_progression_rate' must have only non-negative values")
  expect_error(
    build_vaccine_progression_rate(vaccine_progression_rate = matrix(1, 10, 5),
                                   N_age = 10, N_vacc_classes = 5),
    "'vaccine_progression_rate' must have 'N_vacc_classes - 2' columns")
  expect_error(
    build_vaccine_progression_rate(vaccine_progression_rate = matrix(1, 9, 3),
                                   N_age = 10, N_vacc_classes = 5),
    "'vaccine_progression_rate' must have as many rows as age groups")
})

test_that("build_vaccine_progression_rate allows sensible inputs and works", {
  expect_silent(
    build_vaccine_progression_rate(vaccine_progression_rate = 1,
                                   N_age = 5, N_vacc_classes = 3))
  expect_equal(
    build_vaccine_progression_rate(vaccine_progression_rate = 1,
                                   N_age = 5, N_vacc_classes = 3),
    matrix(rep(1, 5), nrow = 5))
  expect_silent(
    build_vaccine_progression_rate(vaccine_progression_rate = matrix(1, 10, 3),
                                   N_age = 10, N_vacc_classes = 5))
  expect_silent(
    build_vaccine_progression_rate(vaccine_progression_rate = NULL,
                                   N_age = 10, N_vacc_classes = 5))
  expect_equal(
    build_vaccine_progression_rate(vaccine_progression_rate = NULL,
                                   N_age = 10, N_vacc_classes = 5),
    matrix(0, 10, 3))
})
