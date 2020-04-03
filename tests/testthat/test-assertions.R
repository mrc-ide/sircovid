context("assertions")

#------------------------------------------------
test_that("assert_numeric working correctly", {
  expect_true(assert_numeric(5))
  expect_true(assert_numeric(-5:5))
  
  expect_error(assert_numeric(NULL))
  expect_error(assert_numeric("foo"))
  expect_error(assert_numeric(c(1, "foo")))
})

#------------------------------------------------
test_that("assert_int working correctly", {
  expect_true(assert_int(5))
  expect_true(assert_int(-5))
  expect_true(assert_int(-5:5))
  expect_true(assert_int(c(a = 5)))
  
  expect_error(assert_int(NULL))
  expect_error(assert_int(0.5))
  expect_error(assert_int("foo"))
  expect_error(assert_int(c(5,"foo")))
})

#------------------------------------------------
test_that("assert_pos working correctly", {
  expect_true(assert_pos(5))
  expect_true(assert_pos(seq(1, 5, 0.5)))
  expect_true(assert_pos(seq(0, 5, 0.5), zero_allowed = TRUE))
  
  expect_error(assert_pos(NULL))
  expect_error(assert_pos(-5))
  expect_error(assert_pos(seq(-1, -5, -0.5)))
  expect_error(assert_pos(seq(0, 5, 0.5), zero_allowed = FALSE))
  expect_error(assert_pos(seq(-5, 5, 0.5), zero_allowed = TRUE))
  expect_error(assert_pos("foo"))
})

#------------------------------------------------
test_that("assert_pos_int working correctly", {
  expect_true(assert_pos_int(5))
  expect_true(assert_pos_int(0, zero_allowed = TRUE))
  expect_true(assert_pos_int(1:5))
  expect_true(assert_pos_int(0:5, zero_allowed = TRUE))
  
  expect_error(assert_pos_int(NULL))
  expect_error(assert_pos_int(-5))
  expect_error(assert_pos_int(-1:-5))
  expect_error(assert_pos_int(0:5, zero_allowed = FALSE))
  expect_error(assert_pos_int(-5:5, zero_allowed = TRUE))
  expect_error(assert_pos_int("foo"))
})

#------------------------------------------------
test_that("assert_is", {
  expect_error(assert_is("x", "foo"), "must be a foo")
  expect_silent(assert_is(structure("x", class = "foo"), "foo"))
})