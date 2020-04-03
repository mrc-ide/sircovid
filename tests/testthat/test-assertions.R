context("assertions")

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