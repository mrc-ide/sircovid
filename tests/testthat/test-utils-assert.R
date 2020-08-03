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
