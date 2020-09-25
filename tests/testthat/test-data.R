context("data")

test_that("can create augmented data frame with sircovid_data", {
  start <- as_date("2020-01-15")
  from <- as_date("2020-02-01")
  to <- as_date("2020-03-01")
  d <- data_frame(date = seq(from, to, by = 1),
                  x = runif(to - from + 1))
  res <- sircovid_data(d, start_date = start, 1 / 4)
  expect_setequal(
    names(res),
    c("date_start", "date_end", "step_start", "step_end", "x"))
  expect_equal(res$date_start, c(15, 32:60))
  expect_equal(res$date_end, 32:61)
  expect_equal(res$step_start, res$date_start * 4)
  expect_equal(res$step_end, res$date_end * 4)
  expect_equal(res$x, d$x)
})
