context("testing time-varying beta works correctly")

test_that("Two-level beta works in odin as expected", {
  
  sircovid_model <- hospital_model()
  
  beta = sircovid_model$generate_beta_func(beta_start = 3,
                                           beta_end = 1.5)
  pars_model = generate_parameters(sircovid_model = sircovid_model,
                               beta = beta$beta,
                               beta_times = beta$beta_times)
  mod <- sircovid_model$odin_model(user = pars_model)
  t_max <- max(pars_model$beta_t)+50
  t <- seq(from = 1, to = t_max)
  tmp <- mod$run(t)
  results <- mod$transform_variables(tmp)
  
  expect_equal(results$beta,pars_model$beta_y[sapply(t,FUN=function(t){
    max(which(pars_model$beta_t<=t))
  })])
  
  
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
  t_max <- max(pars_model$beta_t)+50
  t <- seq(from = 1, to = t_max)
  tmp <- mod$run(t)
  results <- mod$transform_variables(tmp)
  
  expect_equal(results$beta,pars_model$beta_y[sapply(t,FUN=function(t){
    max(which(pars_model$beta_t<=t))
  })])
  
  
})