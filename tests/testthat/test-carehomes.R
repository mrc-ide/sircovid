context("carehomes")

## TODO: this test should come out and be replaced by something that
## does not load black-box parameters. This version just runs the model
test_that("can run the carehomes model", {
  d <- readRDS("carehomes.rds")
  mod <- carehomes$new(d$user, 0, 10)
  expect_equal(names(mod$index()), names(d$index)[-1])
  mod$set_state(d$state)
  mod$set_index(integer(0))
  res <- mod$run(300)
  expect_equal(res, matrix(numeric(0), 0, 10))
})
