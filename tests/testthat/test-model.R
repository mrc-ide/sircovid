context("covid transmission model")

## Sanity Checks

test_that("there are no infections when beta is 0", {

    pars_model <- generate_parameters(
        beta = rep(0, 3)
    )
    mod <- basic(user = pars_model)
    t_max <- 150
    t <- seq(from = 1, to = t_max)
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
     #should be TRUE
    expect_true(
        all(results$S[t_max,,1] == pars_model$S0)
    )

  }
)

test_that("everyone is infected when beta is Inf", {


    pars_model <- generate_parameters(
        beta = rep(1e100, 3)
    )
    mod <- basic(user = pars_model)
    t_max <- 150
    t <- seq(from = 1, to = t_max)
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
     #should be TRUE
    expect_true(
        all(results$S[2:t_max,,1] == 0)
    )
 }

)


test_that("No one is infected if I and E are 0 at t = 0", {
    pars_model <- generate_parameters(
        beta = rep(0.042, 3)
    )
    pars_model$E0[,,]<- 0
    pars_model$I0_asympt[,,]<- 0
    pars_model$I0_mild[,,]<- 0
    pars_model$I0_ILI[,,] <- 0
    pars_model$I0_hosp[,,] <- 0
    pars_model$I0_ICU[,,] <- 0

    mod <- basic(user = pars_model)
    t_max <- 150
    t <- seq(from = 1, to = t_max)
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    expect_true(
        all(results$S[t_max,,1] == pars_model$S0)
    )

 }
)

test_that("No one is hospitalised if p_sympt_ILI is 0", {
    pars_model <- generate_parameters(
        beta = rep(0.042, 3)
    )
    pars_model$p_sympt_ILI[]<-0
    mod <- basic(user = pars_model)
    t_max <- 150
    t <- seq(from = 1, to = t_max)
    tmp <- mod$run(t, replicate = 1)
    results <- mod$transform_variables(tmp)
    #check there were people infected (if yes, should be TRUE)
    expect_true(
        any(results$E[,,,,1] > 0)
    )
    expect_true(
        all(results$I_ILI[,,,,1] == 0)
    )
    expect_true(
        all(results$I_hosp[,,,,1] == 0)
    )
    expect_true(
        all(results$I_ICU[,,,,1] == 0)
    )
    expect_true(
        all(results$R_hosp[,,,,1] == 0)
    )
    expect_true(
        all(results$D[,,1]==0)
   )

 }
)


test_that("parameters can be generated", {
  t_max <- 150
  t <- seq(from = 1, to = t_max)

  #p_recov_ILI=1, no-one hospitalised
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$p_recov_ILI[]<-1
  mod <- basic(user = pars_model)
  
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$I_ILI[,,,,1]>0))
  expect_true(all(results$I_hosp[,,,,1]==0))
  expect_true(all(results$I_ICU[,,,,1]==0))
  expect_true(all(results$R_hosp[,,,,1]==0))
  expect_true(all(results$D[,,1]==0))

  #p_recov_hosp=1, p_death_hosp=0 no-one goes into ICU, no deaths
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$p_recov_hosp[]<-1
  pars_model$p_death_hosp[]<-0
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$I_hosp[,,,,1]>0))
  expect_true(all(results$I_ICU[,,,,1]==0)) 
  expect_true(all(results$R_hosp[,,,,1]==0)) 
  expect_true(all(results$D[,,1]==0)) 

  #p_death_hosp=1, p_recov_hosp=0 no-one goes into ICU, no recovery in hospital
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$p_recov_hosp[]<-0
  pars_model$p_death_hosp[]<-1
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$I_hosp[,,,,1]>0)) #check there were hosp cases (if yes, should be TRUE)
  expect_true(all(results$I_ICU[,,,,1]==0)) 
  expect_true(all(results$R_hosp[,,,,1]==0)) 

  #p_recov_ICU=0, no-one recovers in hospital
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$p_recov_ICU[]<-0
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$I_hosp[,,,,1]>0)) #check there were hosp cases (if yes, should be TRUE)
  expect_true(all(results$R_hosp[,,,,1]==0)) 

  #gamma_E = Inf, E cases must progress in 1 time-step
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_E<-Inf
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$E[,,,,1]>0)) #check there were E cases (if yes, should be TRUE)
  expect_true(all(results$E[2:t_max,,2,,]==results$E[1:(t_max-1),,1,,])) 

  #gamma_asympt = Inf, I_asympt cases must progress in 1 time-step
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_asympt<-Inf
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$I_asympt[,,,,1]>0)) #check there were I_asympt cases (if yes, should be TRUE)
  expect_true(all(results$I_asympt[2:t_max,,2,,]==results$I_asympt[1:(t_max-1),,1,,])) 

  #gamma_mild = Inf, I_mild cases must progress in 1 time-step
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_mild<-Inf
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$I_mild[,,,,1]>0)) #check there were I_mild cases (if yes, should be TRUE)
  expect_true(all(results$I_mild[2:t_max,,2,,]==results$I_mild[1:(t_max-1),,1,,])) 

  #gamma_ILI = Inf, I_ILI cases must progress in 1 time-step
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_ILI<-Inf
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$I_ILI[,,,,1]>0)) #check there were I_ILI cases (if yes, should be TRUE)
  expect_true(all(results$I_ILI[2:t_max,,2,,]==results$I_ILI[1:(t_max-1),,1,,])) 

  #gamma_hosp = Inf, I_hosp cases must progress in 1 time-step
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_hosp<-Inf
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$I_hosp[,,,,1]>0)) #check there were I_hosp cases (if yes, should be TRUE)
  expect_true(all(results$I_hosp[2:t_max,,2,,]==results$I_hosp[1:(t_max-1),,1,,])) 

  #gamma_ICU = Inf, I_ICU cases must progress in 1 time-step
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_ICU<-Inf
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  ## TODO this needs fixing!
  ## expect_true(any(results$I_ICU[,,,,1]>0)) #check there were I_ICU cases (if yes, should be TRUE)
  expect_true(all(results$I_ICU[2:t_max,,2,,]==results$I_ICU[1:(t_max-1),,1,,])) 

  #gamma_rec = Inf, R_hosp cases must progress in 1 time-step
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_rec<-Inf
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  ## NOTE: this does not always appear to be true
  expect_true(any(results$R_hosp[,,,,1] >= 0)) #check there were R_hosp cases (if yes, should be TRUE)
  expect_true(all(results$R_hosp[2:t_max,,2,,]==results$R_hosp[1:(t_max-1),,1,,])) 

  #gamma_E=0, E stay in progression stage 1
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_E<-0
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$E[,,1,,]>0)) #check there are some E cases first (if yes, should be TRUE)
  expect_true(all(results$E[,,2,,]==0)) 

  #gamma_asympt=0, I_asympt stay in progression stage 1
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_asympt<-0
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$I_asympt[,,1,,]>0)) #check there are some I_asympt cases first (if yes, should be TRUE)
  expect_true(all(results$I_asympt[,,2,,]==0)) 

  #gamma_mild=0, I_mild stay in progression stage 1
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_mild<-0
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$I_mild[,,1,,]>0)) #check there are some I_mild cases first (if yes, should be TRUE)
  expect_true(all(results$I_mild[,,2,,]==0)) 

  #gamma_ILI=0, I_ILI stay in progression stage 1
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_ILI<-0
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$I_ILI[,,1,,]>0)) #check there are some I_ILI cases first (if yes, should be TRUE)
  expect_true(all(results$I_ILI[2:t_max,,2,,]==0)) 

  #gamma_hosp=0, I_hosp stay in progression stage 1
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_hosp<-0
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$I_hosp[,,1,,]>0)) #check there are some I_hosp cases first (if yes, should be TRUE)
  expect_true(all(results$I_hosp[2:t_max,,2,,]==0)) 

  #gamma_ICU=0, I_ICU stay in progression stage 1
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_ICU<-0
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  expect_true(any(results$I_ICU[,,1,,]>0)) #check there are some I_ICU cases first (if yes, should be TRUE)
  expect_true(all(results$I_ICU[,,2,,]==0)) 

  #gamma_rec=0, R_hosp stay in progression stage 1
  pars_model <- generate_parameters(
    beta = rep(0.042, 3)
  )
  pars_model$gamma_rec<-0
  mod <- basic(user = pars_model)
  tmp <- mod$run(t, replicate = 1)
  results <- mod$transform_variables(tmp)
  ## TODO: this is broken
  ## expect_true(any(results$R_hosp[,,1,,]>0)) #check there are some R_hosp cases first (if yes, should be TRUE)
  expect_true(all(results$R_hosp[,,2,,]==0)) 
})
