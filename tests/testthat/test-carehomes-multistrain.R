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
