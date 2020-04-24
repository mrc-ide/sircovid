context("util")

test_that("null-or-value works", {
  expect_equal(1 %||% NULL, 1)
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% NULL, NULL)
  expect_equal(NULL %||% 2, 2)
})

test_that("zero bounday works for matrices and arrays", {

  ## Works for matrix
  expect_true(zero_boundary(matrix(0.01, nrow = 3, ncol = 4), 0.1))
  ## Works for arrays
  expect_true(zero_boundary(array(0.01, dim = c(3, 4, 5)), 0.1))

})
