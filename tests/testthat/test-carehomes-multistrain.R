context("carehomes (multistrain)")

test_that("carehomes_parameters_strain works as expected", {
  expect_error(
    carehomes_parameters_strain(NULL, NULL, NULL, 1),
    "At least one value required for 'strain_transmission'")
  expect_error(
    carehomes_parameters_strain(-1, NULL, NULL, 1),
    "'strain_transmission' must have only non-negative values",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(c(1, -1), NULL, NULL, 1),
    "'strain_transmission' must have only non-negative values",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(rep(0.5, 2), NULL, NULL, 1),
    "'strain_transmission[1]' must be 1",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(rep(0.5, 1), NULL, NULL, 1),
    "'strain_transmission[1]' must be 1",
    fixed = TRUE)
  expect_equal(
    carehomes_parameters_strain(1, NULL, NULL, 1),
    list(n_strains = 1,
         strain_transmission = 1,
         strain_seed_step = 0))
  expect_equal(
    carehomes_parameters_strain(c(1, 1), NULL, NULL, 1),
    list(n_strains = 2,
         strain_transmission = c(1, 1),
         strain_seed_step = 0))
  expect_equal(
    carehomes_parameters_strain(c(1, 2), NULL, NULL, 1),
    list(n_strains = 2,
         strain_transmission = c(1, 2),
         strain_seed_step = 0))
})


test_that("Can seed with one-day window", {
  date <- c("2020-03-01", "2020-03-01")
  value <- 100
  p <- carehomes_parameters_strain(c(1, 1), sircovid_date(date), value, 1 / 4)
  expect_equal(sum(p$strain_seed_step), 100)
  expect_equal(tail(p$strain_seed_step, 6), c(0, 25, 25, 25, 25, 0))
  expect_equal(sircovid_date_as_date(length(p$strain_seed_step) / 4),
               as.Date("2020-03-02"))
})


test_that("Can seed with multiple-day window", {
  date <- c("2020-03-01", "2020-03-10")
  value <- 100
  p <- carehomes_parameters_strain(c(1, 1), sircovid_date(date), value, 1 / 4)
  expect_equal(sum(p$strain_seed_step), 100 * 10)
  expect_equal(tail(p$strain_seed_step, 6), c(25, 25, 25, 25, 25, 0))
  expect_equal(sircovid_date_as_date(length(p$strain_seed_step) / 4),
               as.Date("2020-03-11"))
})


test_that("Adding empty strains makes no difference", {
  p1 <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  p2 <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                             strain_transmission = c(1, 0))
  p3 <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                             strain_transmission = c(1, 0.5))
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


test_that("Seeding of second strain generates an epidemic", {
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1), 
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding)),
                            strain_seed_value = n_seeded_new_strain_inf)
  
  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))
  # did the seeded cases go on to infect other people?
  expect_true(sum(y$I_asympt[, , , 2, ]) > n_seeded_new_strain_inf) 
  
  # check the epidemic of the second strain starts when we expect
  steps <- seq(0, 400, by = 4)
  date <- sircovid_date_as_date(steps / 4)
  # no cases before seeding
  expect_true(all(y$E[, , , 2, date < date_seeding] == 0))
  # no cases on seeding day other than in 4th age group
  expect_true(all(y$E[-4, , , 2, date == date_seeding] == 0))
  # some cases on seeding day in 4th age group
  expect_true(y$E[4, 1, , 2, date == date_seeding] > 0)
  # some cases on all days after seeding day
  expect_true(all(colSums(y$E[, 1, , 2, date >= date_seeding]) > 0))
})
