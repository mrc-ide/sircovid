context("testing time-varying beta works correctly")

test_that("One-level beta works in odin as expected", {
  sircovid_model <- hospital_model()

  pars_model <- generate_parameters(
    sircovid_model = sircovid_model,
    transmission_model = "POLYMOD",
    beta = 0.1,
    beta_times = sircovid_date('2020-01-01'),
    hosp_transmission = 0,
    ICU_transmission = 0,
    trans_profile = 1,
    trans_increase = 1,
    dt = 1/4
  )

  mod <- sircovid_model$odin_model(user = pars_model)

  ## As in the old version:
  pars_model$beta_t <- 0
  pars_model$beta_y <- 0.1

  t_max <- max(pars_model$beta_t)+50
  t <- seq(from = 1, to = t_max)
  tmp <- mod$run(t)
  results <- mod$transform_variables(tmp)

  expect_equal(results$beta, rep(0.1, length(t)))
})

test_that("Two-level beta works in odin as expected", {
  
  sircovid_model <- hospital_model()
  
  beta = sircovid_model$generate_beta_func(beta_start = 3,
                                           beta_end = 1.5)
  pars_model = generate_parameters(sircovid_model = sircovid_model,
                               beta = beta$beta,
                               beta_times = beta$beta_times)
  mod <- sircovid_model$odin_model(user = pars_model)

  ## As in the old version:
  pars_model$beta_t <- normalise_beta(beta$beta_times, pars_model$dt)
  pars_model$beta_y <- beta$beta

  t_max <- max(pars_model$beta_t)+50
  t <- seq(from = 1, to = t_max)
  tmp <- mod$run(t)
  results <- mod$transform_variables(tmp)
  
  expect_equal(results$beta,pars_model$beta_y[sapply(t,FUN=function(t){
    max(which(pars_model$beta_t<=t))
  })])
  
  
})

test_that("Two-level beta works the same with generate_beta/generate_parameters as with update_beta", {
  
  sircovid_model <- hospital_model()
  
  beta = sircovid_model$generate_beta_func(start_date = sircovid_date("2020-02-06"),
                                           beta_start = 3,
                                           beta_end = 1.5)
  pars_model = generate_parameters(sircovid_model = sircovid_model,
                                   beta = beta$beta,
                                   beta_times = beta$beta_times,
                                   dt = 0.25)
  
  beta_update <- update_beta(sircovid_model = sircovid_model, 
                             beta_start = 3,
                             beta_end = 1.5,
                             beta_pl = NULL,
                             start_date = sircovid_date("2020-02-06"),
                             dt = 0.25)
  
  expect_equal(pars_model$beta2,beta_update$beta2)

})

test_that("Three-level beta works in odin as expected", {
  
  sircovid_model <- hospital_model()
  
  beta = sircovid_model$generate_beta_func(beta_start = 3,
                                           beta_end = 1.5,
                                           beta_pl = 1.8)
  pars_model = generate_parameters(sircovid_model = sircovid_model,
                                   beta = beta$beta,
                                   beta_times = beta$beta_times)
  mod <- sircovid_model$odin_model(user = pars_model)

  ## As in the old version:
  pars_model$beta_t <- normalise_beta(beta$beta_times, pars_model$dt)
  pars_model$beta_y <- beta$beta

  t_max <- max(pars_model$beta_t)+50
  t <- seq(from = 1, to = t_max)
  tmp <- mod$run(t)
  results <- mod$transform_variables(tmp)
  
  expect_equal(results$beta,pars_model$beta_y[sapply(t,FUN=function(t){
    max(which(pars_model$beta_t<=t))
  })])
  
  
})

test_that("Three-level beta works the same with generate_beta/generate_parameters as with update_beta", {
  
  sircovid_model <- hospital_model()
  
  beta = sircovid_model$generate_beta_func(start_date = sircovid_date("2020-02-06"),
                                           beta_start = 3,
                                           beta_end = 1.5,
                                           beta_pl = 1.8)
  pars_model = generate_parameters(sircovid_model = sircovid_model,
                                   beta = beta$beta,
                                   beta_times = beta$beta_times,
                                   dt = 0.25)
  
  beta_update <- update_beta(sircovid_model = sircovid_model, 
                             beta_start = 3,
                             beta_end = 1.5,
                             beta_pl = 1.8,
                             start_date = sircovid_date("2020-02-06"),
                             dt = 0.25)
  
  expect_equal(pars_model$beta2,beta_update$beta2)
  
})
