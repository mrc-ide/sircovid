context("new hospital model for covid transmission")

## Sanity Checks

test_that("there are no infections when beta is 0", {

    pars_model <- generate_parameters_new_hospital_model(beta = 0, beta_times = "2020-02-01")
    mod <- new_hospital_model(user = pars_model)
    t_max <- 150
    t <- seq(from = 1, to = t_max)
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
     #should be TRUE
    expect_true(all(results$S[t_max,,1] == pars_model$S0))

  }
)

test_that("everyone is infected when beta is Inf", {
    pars_model <- generate_parameters_new_hospital_model(beta = 1e100, beta_times = "2020-02-01")
    mod <- new_hospital_model(user = pars_model)
    t_max <- 150
    t <- seq(from = 1, to = t_max)
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
     #should be TRUE
    expect_true(all(results$S[2:t_max,,1] == 0))
 }

)


test_that("No one is infected if I and E are 0 at t = 0", {
    pars_model <- generate_parameters_new_hospital_model(beta = 0.042, beta_times = "2020-02-01")
    pars_model$E0[,,]<- 0
    pars_model$I0_asympt[,,]<- 0
    pars_model$I0_mild[,,]<- 0
    pars_model$I0_ILI[,,] <- 0
    pars_model$I0_hosp_R[,,] <- 0
    pars_model$I0_hosp_D[,,] <- 0
    pars_model$I0_ICU_R[,,] <- 0
    pars_model$I0_ICU_D[,,] <- 0
    pars_model$I0_triage[,,] <- 0

    mod <- new_hospital_model(user = pars_model)
    t_max <- 150
    t <- seq(from = 1, to = t_max)
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    expect_true(all(results$S[t_max,,1] == pars_model$S0))

 }
)

test_that("No one is hospitalised if p_sympt_ILI is 0", {
    pars_model <- generate_parameters_new_hospital_model(beta = 0.042, beta_times = "2020-02-01")
    pars_model$p_sympt_ILI[]<-0
    mod <- new_hospital_model(user = pars_model)
    t_max <- 150
    t <- seq(from = 1, to = t_max)
    check_cases <- FALSE
    max_iter <- 10
    iter <- 0
    while (!check_cases && iter<= max_iter){
      iter <- iter + 1
      tmp <- mod$run(t, replicate = 1)
      results <- mod$transform_variables(tmp)
      if (any(results$E[,,,,1] > 0)){
        check_cases <- TRUE
      }
    }
    #check there were people infected (if yes, should be TRUE)
    expect_true(all(results$I_ILI[,,,,1] == 0))
    expect_true(all(results$I_hosp_R[,,,,1] == 0))
    expect_true(all(results$I_hosp_D[,,,,1] == 0))
    expect_true(all(results$I_ICU_R[,,,,1] == 0))
    expect_true(all(results$I_ICU_D[,,,,1] == 0))
    expect_true(all(results$I_triage[,,,,1] == 0))
    expect_true(all(results$R_stepdown[,,,,1] == 0))
    expect_true(all(results$D[,,1]==0))
 }
)


test_that("parameters can be generated", {
  t_max <- 150
  t <- seq(from = 1, to = t_max)
  max_iter <- 10

  #p_recov_ILI=1, no-one hospitalised
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$p_recov_ILI[]<-1
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$I_ILI[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(all(results$I_hosp_R[,,,,1] == 0))
  expect_true(all(results$I_hosp_D[,,,,1] == 0))
  expect_true(all(results$I_ICU_R[,,,,1] == 0))
  expect_true(all(results$I_ICU_D[,,,,1] == 0))
  expect_true(all(results$I_triage[,,,,1] == 0))
  expect_true(all(results$R_stepdown[,,,,1] == 0))
  expect_true(all(results$D[,,1]==0))

  #p_ICU_hosp=0, p_death_hosp=0 no-one goes into ICU, no deaths
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$p_ICU_hosp[]<-0
  pars_model$p_death_hosp[]<-0
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$I_hosp_R[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(all(results$I_hosp_D[,,,,1] == 0))
  expect_true(all(results$I_ICU_R[,,,,1] == 0))
  expect_true(all(results$I_ICU_D[,,,,1] == 0))
  expect_true(all(results$I_triage[,,,,1] == 0))
  expect_true(all(results$R_stepdown[,,,,1] == 0)) 
  expect_true(all(results$D[,,1]==0)) 

  #p_death_hosp=1, p_ICU_hosp=0 no-one goes into ICU, no recovery in hospital
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$p_ICU_hosp[]<-0
  pars_model$p_death_hosp[]<-1
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$I_hosp_D[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(all(results$I_hosp_R[,,,,1] == 0))
  expect_true(all(results$I_ICU_R[,,,,1] == 0))
  expect_true(all(results$I_ICU_D[,,,,1] == 0))
  expect_true(all(results$I_triage[,,,,1] == 0))
  expect_true(all(results$R_stepdown[,,,,1] == 0))
  
  #p_death_ICU=1, p_ICU_hosp=1 no-one goes in hosp_D/hosp_R, no recovery from ICU
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$p_ICU_hosp[]<-1
  pars_model$p_death_ICU[]<-1
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$I_ICU_D[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(any(results$I_hosp_R[,,,,1] == 0))
  expect_true(all(results$I_hosp_D[,,,,1] == 0))
  expect_true(all(results$I_ICU_R[,,,,1] == 0))
  expect_true(all(results$R_stepdown[,,,,1] == 0))
  
  #p_death_ICU=0, p_ICU_hosp=1 no-one goes in hosp_D/hosp_R, no deaths
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$p_ICU_hosp[]<-1
  pars_model$p_death_ICU[]<-1
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$I_ICU_R[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(any(results$I_hosp_R[,,,,1] == 0))
  expect_true(all(results$I_hosp_D[,,,,1] == 0))
  expect_true(all(results$I_ICU_D[,,,,1] == 0))
  expect_true(all(results$D[,,,,1] == 0))
  
  #function to test that setting a given gamma to Inf causes cases in corresponding compartment to
  #progress in 1 time-step
  test_gamma_inf <- function(gamma_name,compartment_name){
    pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01",
                                                         progression_groups = list(E = 2, asympt = 2, mild = 2, ILI = 2, hosp_D = 2 , hosp_R = 2, ICU_D = 2, ICU_R = 2, triage = 2, stepdown = 2))
    pars_model[gamma_name] <- Inf
    mod <- new_hospital_model(user = pars_model)
    check_cases <- FALSE
    iter <- 0
    while (!check_cases && iter<= max_iter){
      iter <- iter + 1
      tmp <- mod$run(t, replicate = 1)
      results <- mod$transform_variables(tmp)
      if (any(results[[compartment_name]][,,,,1] > 0)){
        check_cases <- TRUE
      }
    }
    expect_true(all(results[[compartment_name]][2:t_max,,2,,]==results[[compartment_name]][1:(t_max-1),,1,,]))
  }
  
  test_gamma_inf('gamma_E','E')
  test_gamma_inf('gamma_asympt','I_asympt')
  test_gamma_inf('gamma_mild','I_mild')
  test_gamma_inf('gamma_ILI','I_ILI')
  test_gamma_inf('gamma_triage','I_triage')
  test_gamma_inf('gamma_hosp_R','I_hosp_R')
  test_gamma_inf('gamma_hosp_D','I_hosp_D')
  test_gamma_inf('gamma_ICU_R','I_ICU_R')
  test_gamma_inf('gamma_ICU_D','I_ICU_D')
  test_gamma_inf('gamma_stepdown','R_stepdown')
  
  
  #function to test that setting a given gamma to 0 causes cases in corresponding compartment to
  #stay in progression stage 1
  test_gamma_zero <- function(gamma_name,compartment_name){
    pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01",
                                                         progression_groups = list(E = 2, asympt = 2, mild = 2, ILI = 2, hosp_D = 2 , hosp_R = 2, ICU_D = 2, ICU_R = 2, triage = 2, stepdown = 2))
    pars_model[gamma_name] <- 0
    mod <- new_hospital_model(user = pars_model)
    check_cases <- FALSE
    iter <- 0
    while (!check_cases && iter<= max_iter){
      iter <- iter + 1
      tmp <- mod$run(t, replicate = 1)
      results <- mod$transform_variables(tmp)
      if (any(results[[compartment_name]][,,,,1] > 0)){
        check_cases <- TRUE
      }
    }
    expect_true(all(results[[compartment_name]][,,2,,]==0))
  }
  
  test_gamma_zero('gamma_E','E')
  test_gamma_zero('gamma_asympt','I_asympt')
  test_gamma_zero('gamma_mild','I_mild')
  test_gamma_zero('gamma_ILI','I_ILI')
  test_gamma_zero('gamma_triage','I_triage')
  test_gamma_zero('gamma_hosp_R','I_hosp_R')
  test_gamma_zero('gamma_hosp_D','I_hosp_D')
  test_gamma_zero('gamma_ICU_R','I_ICU_R')
  test_gamma_zero('gamma_ICU_D','I_ICU_D')
  test_gamma_zero('gamma_stepdown','R_stepdown')
  
})
