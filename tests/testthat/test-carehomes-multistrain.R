context("carehomes (multistrain)")

test_that("carehomes_parameters_strain works as expected", {
  expect_error(
    carehomes_parameters_strain(NULL),
    "At least one value required for 'strain_transmission'")
  expect_error(
    carehomes_parameters_strain(-1),
    "'strain_transmission' must have only non-negative values",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(c(1, -1)),
    "'strain_transmission' must have only non-negative values",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(rep(0.5, 2)),
    "'strain_transmission[1]' must be 1",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(rep(0.5, 1)),
    "'strain_transmission[1]' must be 1",
    fixed = TRUE)
  expect_equal(
    carehomes_parameters_strain(1),
    list(n_strains = 1,
         strain_transmission = 1))
  expect_equal(
    carehomes_parameters_strain(c(1, 1)),
    list(n_strains = 2,
         strain_transmission = c(1, 1)))
  expect_equal(
    carehomes_parameters_strain(c(1, 2)),
    list(n_strains = 2,
         strain_transmission = c(1, 2)))
})


test_that("Adding empty strains makes no difference", {
  p1 <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  p2 <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                             strain_transmission = c(1, 0))
  p3 <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                             strain_transmission = c(1, 1, 1))
  np <- 10
  mod1 <- carehomes$new(p1, 0, np, seed = 1L)
  mod2 <- carehomes$new(p2, 0, np, seed = 1L)
  mod3 <- carehomes$new(p3, 0, np, seed = 1L)
  end <- sircovid_date("2020-03-31") / p1$dt

  mod1$set_index(carehomes_index(mod1$info())$run)
  mod2$set_index(carehomes_index(mod2$info())$run)
  mod3$set_index(carehomes_index(mod3$info())$run)

  initial1 <- carehomes_initial(mod1$info(), 1, p1)
  initial2 <- carehomes_initial(mod1$info(), 1, p2)
  initial3 <- carehomes_initial(mod1$info(), 1, p3)

  res1 <- mod1$run(end)
  res2 <- mod2$run(end)
  res3 <- mod3$run(end)

  expect_equal(res2, res1)
  expect_equal(res3, res1)
})
