context("utils (assert)")

test_that("assert_is", {
  expect_error(assert_is("x", "foo"), "must be a foo")
  expect_silent(assert_is(structure("x", class = "foo"), "foo"))
})


test_that("assert_increasing", {
  x <- 1:10
  y <- c(1:5, 5:10)
  expect_silent(assert_increasing(x))
  expect_silent(assert_increasing(x[1]))
  expect_silent(assert_increasing(x[0]))

  expect_error(assert_increasing(y),
               "'y' must be strictly increasing")
  expect_error(assert_increasing(rev(x)),
               "must be strictly increasing")
})


test_that("assert_date_string", {
  expect_silent(assert_date_string("2020-01-01"))
  expect_silent(assert_date_string(rep("2020-01-01", 3)))
  expect_equal(assert_date_string(as.Date("2020-01-01")),
               "2020-01-01")
  expect_error(assert_date_string("June 24, 2020"),
               "Expected ISO dates or R dates for")
})


test_that("assert_scalar", {
  expect_silent(assert_scalar(1))
  object <- NULL
  expect_error(assert_scalar(object),
               "'object' must be a scalar")
  expect_error(assert_scalar(1:5),
               "'object' must be a scalar")
})
