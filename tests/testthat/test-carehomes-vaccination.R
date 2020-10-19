context("carehomes (vaccination)")

# same test as test for beta = 0 in test-carehomes-check.R
test_that("there are no infections if everyone is vaccinated with a vaccine
          preventing 100% of acquisition", {
  p <- carehomes_parameters(0, "england", rel_susc = 0)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  s <- dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$S)

  ## Susceptible population is never drawn down:
  expect_equal(s, array(s[, , 1], c(19, 1, 101)))
})
