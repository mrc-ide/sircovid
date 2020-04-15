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
  expect_true(any(results$I_hosp_R[,,,,1] == 0))
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
  
  #gamma_E = Inf, E cases must progress in 1 time-step
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_E<-Inf
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$E[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(all(results$E[2:t_max,,2,,]==results$E[1:(t_max-1),,1,,])) 

  # Check progression groups work for these parameters, even though
  # they are no longer the default
  #gamma_asympt = Inf, I_asympt cases must progress in 1 time-step
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01",
                                    progression_groups = list(E = 2, asympt = 2, mild = 1, ILI = 1, hosp_D = 2 , hosp_R = 2, ICU_D = 2, ICU_R = 2, triage = 2, stepdown = 2))
  pars_model$gamma_asympt<-Inf
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$I_asympt[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(all(results$I_asympt[2:t_max,,2,,]==results$I_asympt[1:(t_max-1),,1,,])) 

  #gamma_mild = Inf, I_mild cases must progress in 1 time-step
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01",
                                    progression_groups = list(E = 2, asympt = 1, mild = 2, ILI = 1, hosp_D = 2 , hosp_R = 2, ICU_D = 2, ICU_R = 2, triage = 2, stepdown = 2))
  pars_model$gamma_mild<-Inf
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$I_mild[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(all(results$I_mild[2:t_max,,2,,]==results$I_mild[1:(t_max-1),,1,,])) 

  #gamma_ILI = Inf, I_ILI cases must progress in 1 time-step
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01",
                                    progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 2, hosp_D = 2 , hosp_R = 2, ICU_D = 2, ICU_R = 2, triage = 2, stepdown = 2))
  pars_model$gamma_ILI<-Inf
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
  expect_true(all(results$I_ILI[2:t_max,,2,,]==results$I_ILI[1:(t_max-1),,1,,]))
  
  #gamma_triage = Inf, I_triage cases must progress in 1 time-step
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_triage<-Inf
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$I_triage[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(all(results$I_triage[2:t_max,,2,,]==results$I_triage[1:(t_max-1),,1,,]))

  #gamma_hosp_R = Inf, I_hosp_R cases must progress in 1 time-step
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_hosp_R<-Inf
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
  expect_true(all(results$I_hosp_R[2:t_max,,2,,]==results$I_hosp_R[1:(t_max-1),,1,,]))
  
  #gamma_hosp_D = Inf, I_hosp_D cases must progress in 1 time-step
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_hosp_D<-Inf
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
  expect_true(all(results$I_hosp_D[2:t_max,,2,,]==results$I_hosp_D[1:(t_max-1),,1,,]))

  #gamma_ICU_R = Inf, I_ICU_R cases must progress in 1 time-step
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_ICU_R<-Inf
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
  expect_true(all(results$I_ICU_R[2:t_max,,2,,]==results$I_ICU_R[1:(t_max-1),,1,,]))
  
  #gamma_ICU_D = Inf, I_ICU_D cases must progress in 1 time-step
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_ICU_D<-Inf
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
  expect_true(all(results$I_ICU_D[2:t_max,,2,,]==results$I_ICU_D[1:(t_max-1),,1,,]))
  
  #gamma_stepdown = Inf, R_stepdown cases must progress in 1 time-step
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_stepdown<-Inf
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$R_stepdown[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(all(results$R_stepdown[2:t_max,,2,,]==results$R_stepdown[1:(t_max-1),,1,,])) 

  #gamma_E=0, E stay in progression stage 1
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_E<-0
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$E[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(all(results$E[,,2,,]==0)) 

  #gamma_asympt=0, I_asympt stay in progression stage 1
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01",
                                    progression_groups = list(E = 2, asympt = 2, mild = 1, ILI = 1, hosp = 2, ICU = 2, rec = 2))
  pars_model$gamma_asympt<-0
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$I_asympt[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }  
  expect_true(all(results$I_asympt[,,2,,]==0)) 

  #gamma_mild=0, I_mild stay in progression stage 1
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01",
                                    progression_groups = list(E = 2, asympt = 1, mild = 2, ILI = 1, hosp = 2, ICU = 2, rec = 2))
  pars_model$gamma_mild<-0
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$I_mild[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(all(results$I_mild[,,2,,]==0)) 

  #gamma_ILI=0, I_ILI stay in progression stage 1
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01",
                                    progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 2, hosp = 2, ICU = 2, rec = 2))
  pars_model$gamma_ILI<-0
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
  expect_true(all(results$I_ILI[2:t_max,,2,,]==0)) 

  #gamma_hosp_R=0, I_hosp_R stay in progression stage 1
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_hosp_R<-0
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
  expect_true(all(results$I_hosp_R[2:t_max,,2,,]==0))
  
  #gamma_hosp_D=0, I_hosp_D stay in progression stage 1
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_hosp_D<-0
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
  expect_true(all(results$I_hosp_D[2:t_max,,2,,]==0)) 

  #gamma_ICU_R=0, I_ICU_R stay in progression stage 1
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_ICU_R<-0
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
  expect_true(all(results$I_ICU_R[,,2,,]==0)) 
  
  #gamma_ICU_D=0, I_ICU_D stay in progression stage 1
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_ICU_D<-0
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
  expect_true(all(results$I_ICU_D[2:t_max,,2,,]==0))
  
  #gamma_triage=0, I_triage stay in progression stage 1
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_triage<-0
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$I_triage[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(all(results$I_triage[2:t_max,,2,,]==0))

  #gamma_stepdown=0, R_stepdown stay in progression stage 1
  pars_model <- generate_parameters_new_hospital_model(beta = 0.126, beta_times = "2020-02-01")
  pars_model$gamma_stepdown<-0
  mod <- new_hospital_model(user = pars_model)
  check_cases <- FALSE
  iter <- 0
  while (!check_cases && iter<= max_iter){
    iter <- iter + 1
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    if (any(results$R_stepdown[,,,,1] > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(all(results$R_stepdown[,,2,,]==0)) 
})
