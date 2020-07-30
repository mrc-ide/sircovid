context("basic")

## TODO: this test should come out and be replaced by something that
## does not load black-box parameters. This version just runs the model
test_that("can run the basic model", {
  d <- readRDS("basic.rds")
  mod <- basic$new(d$user, 0, 10)
  expect_equal(names(mod$index()), names(d$index)[-1])
  mod$set_index(integer(0))
  res <- mod$run(300)
  expect_equal(res, matrix(numeric(0), 0, 10))
})


test_that("can run the basic model", {
  p <- basic_parameters(sircovid_date("2020-02-07"), "england")
  mod <- basic$new(p, 0, 10)
  end <- sircovid_date("2020-07-31") / p$dt

  initial <- basic_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  mod$set_index(basic_index(mod$info()))
  res <- mod$run(end)
  expected <-
    rbind(c(1352, 2102, 1714, 1498, 878, 729, 604, 1037, 1206, 1424),
          c(275842, 273785, 275310, 275167, 276831, 276762, 277848, 276574,
            276184, 275348))
  expect_equal(res, expected)
})
