context("serology model for covid transmission")

## Sanity Checks

test_that("N_tot and N_tot2 stay constant", {
  
  sircovid_model <- serology_model()
  pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0.1)
  mod <- sircovid_model$odin_model(user = pars_model)
  t_max <- 400
  t <- seq(from = 1, to = t_max)
  tmp <- mod$run(t)
  results <- mod$transform_variables(tmp)
  #should be TRUE
  expect_true(all(t(results$N_tot) == results$N_tot[1,]))
  expect_true(all(results$N_tot2 == results$N_tot2[1]))
  expect_true(results$N_tot2[1] == sum(results$N_tot[1,]))
}
)

test_that("there are no infections when beta is 0", {
    
    sircovid_model <- serology_model()
    pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0)
    mod <- sircovid_model$odin_model(user = pars_model)
    t_max <- 150
    t <- seq(from = 1, to = t_max)
    tmp <- mod$run(t)
    results <- mod$transform_variables(tmp)
     #should be TRUE
    expect_true(all(results$S[t_max,] == pars_model$S0))

  }
)

test_that("everyone is infected when beta is Inf", {
  
    sircovid_model <- serology_model()
    pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = Inf)
    mod <- sircovid_model$odin_model(user = pars_model)
    t_max <- 150
    t <- seq(from = 1, to = t_max)
    tmp <- mod$run(t)
    results <- mod$transform_variables(tmp)
     #should be TRUE
    expect_true(all(results$S[2:t_max,] == 0))
 }

)


test_that("No one is infected if I and E are 0 at t = 0", {
  
    sircovid_model <- serology_model()
    pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0.042)
    pars_model$E0[,,]<- 0
    pars_model$I0_asympt[,,]<- 0
    pars_model$I0_mild[,,]<- 0
    pars_model$I0_ILI[,,] <- 0
    pars_model$I0_hosp_R_unconf[,,] <- 0
    pars_model$I0_hosp_R_conf[,,] <- 0
    pars_model$I0_hosp_D_unconf[,,] <- 0
    pars_model$I0_hosp_D_conf[,,] <- 0
    pars_model$I0_ICU_R_unconf[,,] <- 0
    pars_model$I0_ICU_R_conf[,,] <- 0
    pars_model$I0_ICU_D_unconf[,,] <- 0
    pars_model$I0_ICU_R_unconf[,,] <- 0
    pars_model$I0_triage_D_unconf[,,] <- 0
    pars_model$I0_triage_D_conf[,,] <- 0
    pars_model$I0_triage_R_unconf[,,] <- 0
    pars_model$I0_triage_R_conf[,,] <- 0
    pars_model$I0_comm_D[,,] <- 0

    mod <- sircovid_model$odin_model(user = pars_model)
    t_max <- 150
    t <- seq(from = 1, to = t_max)
    tmp <- mod$run(t)
    results <- mod$transform_variables(tmp)
    expect_true(all(results$S[t_max,] == pars_model$S0))

 }
)

test_that("No one is hospitalised, no-one dies if p_sympt_ILI is 0", {
    
    sircovid_model <- serology_model()
    pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0.042)
    pars_model$p_sympt_ILI[]<-0
    mod <- sircovid_model$odin_model(user = pars_model)
    t_max <- 150
    t <- seq(from = 1, to = t_max)
    check_cases <- FALSE
    max_iter <- 10
    iter <- 0
    while (!check_cases && iter<= max_iter){
      #We want to check that no-one moves from compartment E to compartment
      #I_ILI (and other subsequent compartments). It's possible that no-one
      #gets infected (so no-one is ever in compartment E), so we re-run the 
      #model until there are cases in E (or max_iter is reached)
      iter <- iter + 1
      tmp <- mod$run(t)
      results <- mod$transform_variables(tmp)
      if (any(results$E > 0)){
        check_cases <- TRUE
      }
    }
    expect_true(any(results$E > 0))
    expect_true(all(results$I_ILI == 0))
    expect_true(all(results$I_hosp_R_unconf == 0))
    expect_true(all(results$I_hosp_R_conf == 0))
    expect_true(all(results$I_hosp_D_unconf == 0))
    expect_true(all(results$I_hosp_D_conf == 0))
    expect_true(all(results$I_ICU_R_unconf == 0))
    expect_true(all(results$I_ICU_R_conf == 0))
    expect_true(all(results$I_ICU_D_unconf == 0))
    expect_true(all(results$I_ICU_D_conf == 0))
    expect_true(all(results$I_triage_R_unconf == 0))
    expect_true(all(results$I_triage_R_conf == 0))
    expect_true(all(results$I_triage_D_unconf == 0))
    expect_true(all(results$I_triage_D_conf == 0))
    expect_true(all(results$R_stepdown_unconf == 0))
    expect_true(all(results$R_stepdown_conf == 0))
    expect_true(all(results$D_hosp == 0))
    expect_true(all(results$I_comm_D == 0))
    expect_true(all(results$D_comm == 0))
 }
)


test_that("No one is hospitalised, no-one dies if p_hosp_ILI is 0, p_death_comm is 0", {
  sircovid_model <- serology_model()
  pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0.126)
  pars_model$p_hosp_ILI[] <- 0
  pars_model$p_death_comm[] <- 0
  mod <- sircovid_model$odin_model(user = pars_model)
  t_max <- 150
  t <- seq(from = 1, to = t_max)
  check_cases <- FALSE
  max_iter <- 10
  iter <- 0
  while (!check_cases && iter<= max_iter){
    #We want to check that no-one moves from compartment I_ILI to 
    #any hospital compartments. It's possible that no-one ends up in
    #I_ILI, so we re-run the model until there are cases in I_ILI
    #(or max_iter is reached)
    iter <- iter + 1
    tmp <- mod$run(t)
    results <- mod$transform_variables(tmp)
    if (any(results$I_ILI > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(any(results$I_ILI > 0))
  expect_true(all(results$I_hosp_R_unconf == 0))
  expect_true(all(results$I_hosp_R_conf == 0))
  expect_true(all(results$I_hosp_D_unconf == 0))
  expect_true(all(results$I_hosp_D_conf == 0))
  expect_true(all(results$I_ICU_R_unconf == 0))
  expect_true(all(results$I_ICU_R_conf == 0))
  expect_true(all(results$I_ICU_D_unconf == 0))
  expect_true(all(results$I_ICU_D_conf == 0))
  expect_true(all(results$I_triage_R_unconf == 0))
  expect_true(all(results$I_triage_R_conf == 0))
  expect_true(all(results$I_triage_D_unconf == 0))
  expect_true(all(results$I_triage_D_conf == 0))
  expect_true(all(results$R_stepdown_unconf == 0))
  expect_true(all(results$R_stepdown_conf == 0))
  expect_true(all(results$D_hosp == 0))
  expect_true(all(results$I_comm_D == 0))
  expect_true(all(results$D_comm == 0))
}
)

test_that("No one is hospitalised, no-one recovers if p_hosp_ILI is 0, p_death_comm is 1, p_asympt = 0, p_sympt_ILI = 1", {
  sircovid_model <- serology_model()
  pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0.126)
  #move initial infectives to ILI
  pars_model$I0_ILI <- pars_model$I0_asympt
  pars_model$I0_asympt[] <- 0
  pars_model$p_sympt_ILI[] <- 1
  pars_model$p_hosp_ILI[] <- 0
  pars_model$p_death_comm[] <- 1
  pars_model$p_asympt[] <- 0
  mod <- sircovid_model$odin_model(user = pars_model)
  t_max <- 150
  t <- seq(from = 1, to = t_max)
  check_cases <- FALSE
  max_iter <- 10
  iter <- 0
  while (!check_cases && iter<= max_iter){
    #We want to check that no-one moves from compartment I_ILI to 
    #any hospital compartments. It's possible that no-one ends up in
    #I_ILI, so we re-run the model until there are cases in I_ILI
    #(or max_iter is reached)
    iter <- iter + 1
    tmp <- mod$run(t)
    results <- mod$transform_variables(tmp)
    if (any(results$I_ILI > 0)){
      check_cases <- TRUE
    }
  }
  expect_true(any(results$I_ILI > 0))
  expect_true(all(results$I_hosp_R_unconf == 0))
  expect_true(all(results$I_hosp_R_conf == 0))
  expect_true(all(results$I_hosp_D_unconf == 0))
  expect_true(all(results$I_hosp_D_conf == 0))
  expect_true(all(results$I_ICU_R_unconf == 0))
  expect_true(all(results$I_ICU_R_conf == 0))
  expect_true(all(results$I_ICU_D_unconf == 0))
  expect_true(all(results$I_ICU_D_conf == 0))
  expect_true(all(results$I_triage_R_unconf == 0))
  expect_true(all(results$I_triage_R_conf == 0))
  expect_true(all(results$I_triage_D_unconf == 0))
  expect_true(all(results$I_triage_D_conf == 0))
  expect_true(all(results$R_stepdown_unconf == 0))
  expect_true(all(results$R_stepdown_conf == 0))
  expect_true(all(results$D_hosp == 0))
  }
)

test_that("setting hospital route probabilities to 0 or 1 result in correct path", {
  
  t_max <- 150
  t <- seq(from = 1, to = t_max)
  max_iter <- 10
  
  sircovid_model <- serology_model()
  
  check_hosp_probs <- function(
    p_ICU_hosp,
    p_death_ICU = NULL,
    p_death_hosp_D = NULL,
    cases, #name of comparment where we expect cases
    zeroes #names of compartments where we expect zeroes
    ){
    pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0.126)
    pars_model$p_ICU_hosp[] <- p_ICU_hosp
    if (!is.null(p_death_hosp_D)){
      pars_model$p_death_hosp_D[] <- p_death_hosp_D
    }
    if (!is.null(p_death_ICU)){
      pars_model$p_death_ICU[] <- p_death_ICU
    }
    mod <- sircovid_model$odin_model(user = pars_model)
    check_cases <- FALSE
    iter <- 0
    while (!check_cases && iter<=max_iter){
      #We want to check that certain probabilities result in only one possible route
      #in hospital. It's possible that no-one enters hospital, so we re-run the model
      #until we have cases in hospital (or max_iter is reached)
      iter <- iter + 1
      tmp <- mod$run(t)
      results <- mod$transform_variables(tmp)
      if (any(results$n_ILI_to_hosp > 0)){
        check_cases <- TRUE
      }
    }
    
    results$I_hosp_R <- results$I_hosp_R_unconf + results$I_hosp_R_conf
    results$I_hosp_D <- results$I_hosp_D_unconf + results$I_hosp_D_conf
    results$I_triage_R <- results$I_triage_R_unconf + results$I_triage_R_conf
    results$I_triage_D <- results$I_triage_D_unconf + results$I_triage_D_conf
    results$I_ICU_R <- results$I_ICU_R_unconf + results$I_ICU_R_conf
    results$I_ICU_D <- results$I_ICU_D_unconf + results$I_ICU_D_conf
    results$R_stepdown <- results$R_stepdown_unconf + results$R_stepdown_conf
    
    #check that there are hospital cases in the right compartment
    expect_true(any(results[[cases]] > 0))
    #check that all compartments on other hospital routes are empty
    zero_true <- rep(FALSE,length(zeroes))
    for (i in seq_len(length(zeroes))){
      zero_true[i] <- all(results[[zeroes[i]]] == 0)
    }
    expect_true(all(zero_true))
  }
  
  #p_ICU_hosp=0, p_death_hosp_D=0 no-one goes into ICU, no deaths
  check_hosp_probs(p_ICU_hosp = 0,
                   p_death_hosp_D = 0,
                   cases = "I_hosp_R",
                   zeroes = c("I_hosp_D","I_ICU_R","I_ICU_D","I_triage_R","I_triage_D","R_stepdown","D_hosp"))
  
  
  #p_death_hosp=1, p_ICU_hosp=0 no-one goes into ICU, no recovery in hospital
  check_hosp_probs(p_ICU_hosp = 0,
                   p_death_hosp_D = 1,
                   cases = "I_hosp_D",
                   zeroes = c("I_hosp_R","I_ICU_R","I_ICU_D","I_triage_R","I_triage_D","R_stepdown"))

  #p_death_ICU=1, p_ICU_hosp=1 no-one goes in hosp_D/hosp_R, no recovery from ICU
  check_hosp_probs(p_ICU_hosp = 1,
                   p_death_ICU = 1,
                   cases = "I_ICU_D",
                   zeroes = c("I_hosp_R","I_hosp_D","I_ICU_R","R_stepdown"))
  
  #p_death_ICU=0, p_ICU_hosp=1 no-one goes in hosp_D/hosp_R, no deaths
  check_hosp_probs(p_ICU_hosp = 1,
                   p_death_ICU = 0,
                   cases = "I_ICU_R",
                   zeroes = c("I_hosp_R","I_hosp_D","I_ICU_D","D_hosp"))
  
})

test_that("No one seroconverts if p_seroconversion is 0", {
    sircovid_model <- serology_model()
  pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0.126)
  pars_model$p_seroconversion[] <- 0
  mod <- sircovid_model$odin_model(user = pars_model)
  t_max <- 150
  t <- seq(from = 1, to = t_max)
  check_cases <- FALSE
  max_iter <- 10
  iter <- 0
  while (!check_cases && iter<= max_iter){
    #We want to check that no-one seroconverts. It's possible 
    #that no-one recovers, so we re-run the model until there 
    #are cases in R_pos or R_neg (or max_iter is reached)
    iter <- iter + 1
    tmp <- mod$run(t)
    results <- mod$transform_variables(tmp)
    if (any(results$R_pos + results$R_neg > 0)){
      check_cases <- TRUE
    }
  }
  
  expect_true(any(results$R_neg > 0))
  expect_true(all(results$R_pos == 0))
}
)

test_that("No one doesn't seroconvert if p_seroconversion is 1", {
  sircovid_model <- serology_model()
  pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0.126)
  pars_model$p_seroconversion[] <- 1
  mod <- sircovid_model$odin_model(user = pars_model)
  t_max <- 150
  t <- seq(from = 1, to = t_max)
  check_cases <- FALSE
  max_iter <- 10
  iter <- 0
  while (!check_cases && iter<= max_iter){
    #We want to check that no-one seroconverts. It's possible 
    #that no-one recovers, so we re-run the model until there 
    #are cases in R_pos or R_neg (or max_iter is reached)
    iter <- iter + 1
    tmp <- mod$run(t)
    results <- mod$transform_variables(tmp)
    if (any(results$R_pos + results$R_neg > 0)){
      check_cases <- TRUE
    }
  }
  
  expect_true(any(results$R_pos > 0))
  expect_true(all(results$R_neg == 0))
}
)


test_that("setting a gamma to Inf results in progress in corresponding compartment in 1 time-step", {    
  t_max <- 150
  t <- seq(from = 1, to = t_max)
  max_iter <- 10
  
  sircovid_model <- serology_model(use_fitted_parameters = FALSE ,progression_groups = list(E = 2, asympt = 2, mild = 2, ILI = 2, comm_D = 2, hosp_D = 2 , hosp_R = 2, ICU_D = 2, ICU_R = 2, triage = 2, stepdown = 2, R_pre = 2))
  
  #function to test that setting a given gamma to Inf causes cases in corresponding compartment to
  #progress in 1 time-step
  test_gamma_inf <- function(gamma_name,compartment_name){
    pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0.126)
    pars_model[gamma_name] <- Inf
    mod <- sircovid_model$odin_model(user = pars_model)
    check_cases <- FALSE
    iter <- 0
    while (!check_cases && iter<= max_iter){
      #We want to check that if a gamma is Inf, then cases progress in 1 time-step within 
      #the corresponding compartment. We re-run the model until cases are observed in that
      #compartment (or until max-iter is reached)
      iter <- iter + 1
      tmp <- mod$run(t)
      results <- mod$transform_variables(tmp)
      results$I_hosp_R <- results$I_hosp_R_unconf + results$I_hosp_R_conf
      results$I_hosp_D <- results$I_hosp_D_unconf + results$I_hosp_D_conf
      results$I_triage_R <- results$I_triage_R_unconf + results$I_triage_R_conf
      results$I_triage_D <- results$I_triage_D_unconf + results$I_triage_D_conf
      results$I_ICU_R <- results$I_ICU_R_unconf + results$I_ICU_R_conf
      results$I_ICU_D <- results$I_ICU_D_unconf + results$I_ICU_D_conf
      results$R_stepdown <- results$R_stepdown_unconf + results$R_stepdown_conf
      if (any(results[[compartment_name]] > 0)){
        check_cases <- TRUE
      }
    }
    expect_true(any(results[[compartment_name]] > 0))
    if(!compartment_name %in% c('R_stepdown','R_pre')){
      expect_true(all(results[[compartment_name]][2:t_max,,2,]==results[[compartment_name]][1:(t_max-1),,1,]))
    } else {
      expect_true(all(results[[compartment_name]][2:t_max,,2]==results[[compartment_name]][1:(t_max-1),,1]))
    }
    
  }
  
  test_gamma_inf('gamma_E','E')
  test_gamma_inf('gamma_asympt','I_asympt')
  test_gamma_inf('gamma_mild','I_mild')
  test_gamma_inf('gamma_ILI','I_ILI')
  test_gamma_inf('gamma_triage','I_triage_R')
  test_gamma_inf('gamma_triage','I_triage_D')
  test_gamma_inf('gamma_hosp_R','I_hosp_R')
  test_gamma_inf('gamma_hosp_D','I_hosp_D')
  test_gamma_inf('gamma_ICU_R','I_ICU_R')
  test_gamma_inf('gamma_ICU_D','I_ICU_D')
  test_gamma_inf('gamma_comm_D','I_comm_D')
  test_gamma_inf('gamma_stepdown','R_stepdown')
  test_gamma_inf('gamma_R_pre','R_pre')
})
  
test_that("setting a gamma to 0 results in cases in corresponding compartment to stay in progression stage 1", {  
  t_max <- 150
  t <- seq(from = 1, to = t_max)
  max_iter <- 10
  
  sircovid_model <- serology_model(use_fitted_parameters = FALSE, progression_groups = list(E = 2, asympt = 2, mild = 2, ILI = 2, comm_D = 2, hosp_D = 2 , hosp_R = 2, ICU_D = 2, ICU_R = 2, triage = 2, stepdown = 2, R_pre = 2))
  
  #function to test that setting a given gamma to 0 causes cases in corresponding compartment to
  #stay in progression stage 1
  test_gamma_zero <- function(gamma_name,compartment_name){
    pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0.126)
    pars_model[gamma_name] <- 0
    mod <- sircovid_model$odin_model(user = pars_model)
    check_cases <- FALSE
    iter <- 0
    while (!check_cases && iter<= max_iter){
      #We want to check that if a gamma is 0, then cases in the corresponding compartment stay in
      #progression stage 1. We re-run the model until cases are observed in that compartment (or 
      #until max-iter is reached)
      iter <- iter + 1
      tmp <- mod$run(t)
      results <- mod$transform_variables(tmp)
      results$I_hosp_R <- results$I_hosp_R_unconf + results$I_hosp_R_conf
      results$I_hosp_D <- results$I_hosp_D_unconf + results$I_hosp_D_conf
      results$I_triage_R <- results$I_triage_R_unconf + results$I_triage_R_conf
      results$I_triage_D <- results$I_triage_D_unconf + results$I_triage_D_conf
      results$I_ICU_R <- results$I_ICU_R_unconf + results$I_ICU_R_conf
      results$I_ICU_D <- results$I_ICU_D_unconf + results$I_ICU_D_conf
      results$R_stepdown <- results$R_stepdown_unconf + results$R_stepdown_conf
      if (any(results[[compartment_name]] > 0)){
        check_cases <- TRUE
      }
    }
    expect_true(any(results[[compartment_name]] > 0))
    if (!compartment_name %in% c('R_stepdown','R_pre')){
      expect_true(all(results[[compartment_name]][,,2,]==0))
    } else{
      expect_true(all(results[[compartment_name]][,,2]==0))
    }
  }
  
  test_gamma_zero('gamma_E','E')
  test_gamma_zero('gamma_asympt','I_asympt')
  test_gamma_zero('gamma_mild','I_mild')
  test_gamma_zero('gamma_ILI','I_ILI')
  test_gamma_zero('gamma_triage','I_triage_R')
  test_gamma_zero('gamma_triage','I_triage_D')
  test_gamma_zero('gamma_hosp_R','I_hosp_R')
  test_gamma_zero('gamma_hosp_D','I_hosp_D')
  test_gamma_zero('gamma_ICU_R','I_ICU_R')
  test_gamma_zero('gamma_ICU_D','I_ICU_D')
  test_gamma_zero('gamma_stepdown','R_stepdown')
  test_gamma_zero('gamma_R_pre','R_pre')
  test_gamma_zero('gamma_comm_D','I_comm_D')
  
})

test_that("No one is unconfirmed, if p_admit_conf = 1", {
  sircovid_model <- serology_model()
  pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0.126)
  pars_model$p_admit_conf[] <- 1
  mod <- sircovid_model$odin_model(user = pars_model)
  t_max <- 150
  t <- seq(from = 1, to = t_max)
  check_cases <- FALSE
  max_iter <- 10
  iter <- 0
  while (!check_cases && iter<= max_iter){
    #We want to check that no-one is unconfirmed. It's possible that no-one
    #is confirmed either, so we re-run the model until there are cases in 
    #confirmed compartments (or max_iter is reached)
    iter <- iter + 1
    tmp <- mod$run(t)
    results <- mod$transform_variables(tmp)
    if (all(c(any(results$I_hosp_R_conf > 0),
              any(results$I_hosp_D_conf > 0),
              any(results$I_triage_R_conf > 0),
              any(results$I_triage_D_conf > 0),
              any(results$I_ICU_R_conf > 0),
              any(results$I_ICU_D_conf > 0),
              any(results$R_stepdown_conf > 0)))){
      check_cases <- TRUE
    }
  }
  expect_true(all(results$I_hosp_R_unconf == 0))
  expect_true(any(results$I_hosp_R_conf > 0))
  expect_true(all(results$I_hosp_D_unconf == 0))
  expect_true(any(results$I_hosp_D_conf > 0))
  expect_true(all(results$I_triage_R_unconf == 0))
  expect_true(any(results$I_triage_R_conf > 0))
  expect_true(all(results$I_triage_D_unconf == 0))
  expect_true(any(results$I_triage_D_conf > 0))
  expect_true(all(results$I_ICU_R_unconf == 0))
  expect_true(any(results$I_ICU_R_conf > 0))
  expect_true(all(results$I_ICU_D_unconf == 0))
  expect_true(any(results$I_ICU_D_conf > 0))
  expect_true(all(results$R_stepdown_unconf == 0))
  expect_true(any(results$R_stepdown_conf > 0))
}
)

test_that("No one is confirmed, if p_admit_conf = 0 and gamma_test = 0", {
  sircovid_model <- serology_model()
  pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0.126)
  pars_model$p_admit_conf[] <- 0
  pars_model$gamma_test <- 0
  mod <- sircovid_model$odin_model(user = pars_model)
  t_max <- 150
  t <- seq(from = 1, to = t_max)
  check_cases <- FALSE
  max_iter <- 10
  iter <- 0
  while (!check_cases && iter<= max_iter){
    #We want to check that no-one is confirmed. It's possible that no-one
    #is unconfirmed either, so we re-run the model until there are cases in 
    #unconfirmed compartments (or max_iter is reached)
    iter <- iter + 1
    tmp <- mod$run(t)
    results <- mod$transform_variables(tmp)
    if (all(c(any(results$I_hosp_R_unconf > 0),
              any(results$I_hosp_D_unconf > 0),
              any(results$I_triage_R_unconf > 0),
              any(results$I_triage_D_unconf > 0),
              any(results$I_ICU_R_unconf > 0),
              any(results$I_ICU_D_unconf > 0),
              any(results$R_stepdown_unconf > 0)))){
      check_cases <- TRUE
    }
  }
  expect_true(any(results$I_hosp_R_unconf > 0))
  expect_true(all(results$I_hosp_R_conf == 0))
  expect_true(any(results$I_hosp_D_unconf > 0))
  expect_true(all(results$I_hosp_D_conf == 0))
  expect_true(any(results$I_triage_R_unconf > 0))
  expect_true(all(results$I_triage_R_conf == 0))
  expect_true(any(results$I_triage_D_unconf > 0))
  expect_true(all(results$I_triage_D_conf == 0))
  expect_true(any(results$I_ICU_R_unconf > 0))
  expect_true(all(results$I_ICU_R_conf == 0))
  expect_true(any(results$I_ICU_D_unconf > 0))
  expect_true(all(results$I_ICU_D_conf == 0))
  expect_true(any(results$R_stepdown_unconf > 0))
  expect_true(all(results$R_stepdown_conf == 0))
}
)

test_that("Confirmation in one time-step, if p_admit_conf = 0 and gamma_test = Inf", {
  sircovid_model <- serology_model(use_fitted_parameters = FALSE, progression_groups = list(E = 2, asympt = 2, mild = 2, ILI = 2, comm_D = 2, hosp_D = 2 , hosp_R = 2, ICU_D = 2, ICU_R = 2, triage = 2, stepdown = 2, R_pre = 2))
  pars_model <- generate_parameters(sircovid_model = sircovid_model, beta = 0.126)
  pars_model$p_admit_conf[] <- 0
  pars_model$gamma_test <- Inf
  pars_model$gamma_triage <- Inf
  pars_model$gamma_hosp_R <- Inf
  pars_model$gamma_hosp_D <- Inf
  pars_model$gamma_ICU_R <- Inf
  pars_model$gamma_ICU_D <- Inf
  pars_model$gamma_stepdown <- Inf
  pars_model$I0_ICU_R_unconf[,1,] <- 50 
  pars_model$I0_ICU_D_unconf[,1,] <- 50 
  pars_model$R0_stepdown_unconf[,1] <- 50
  mod <- sircovid_model$odin_model(user = pars_model)
  t_max <- 150
  t <- seq(from = 1, to = t_max)
  check_cases <- FALSE
  max_iter <- 10
  iter <- 0
  while (!check_cases && iter<= max_iter){
    #We want to check there are cases in compartments first (or max_iter is reached)
    iter <- iter + 1
    tmp <- mod$run(t)
    results <- mod$transform_variables(tmp)
    if (all(c(any(results$I_hosp_R_unconf > 0),
              any(results$I_hosp_D_unconf > 0),
              any(results$I_triage_R_unconf > 0),
              any(results$I_triage_D_unconf > 0)
              ))){
      check_cases <- TRUE
    }
  }
  #check hosp_R
  expect_true(all(results$I_hosp_R_conf[,,1,] == 0))
  expect_true(all(results$I_hosp_R_conf[2:t_max,,2,] == results$I_hosp_R_unconf[1:(t_max-1),,1,]))
  expect_true(all(results$I_hosp_R_unconf[,,2,] == 0))
  #check hosp_D
  expect_true(all(results$I_hosp_D_conf[,,1,] == 0))
  expect_true(all(results$I_hosp_D_conf[2:t_max,,2,] == results$I_hosp_D_unconf[1:(t_max-1),,1,]))
  expect_true(all(results$I_hosp_D_unconf[,,2,] == 0))
  #check triage_R
  expect_true(all(results$I_triage_R_conf[,,1,] == 0))
  expect_true(all(results$I_triage_R_conf[2:t_max,,2,] == results$I_triage_R_unconf[1:(t_max-1),,1,]))
  expect_true(all(results$I_triage_R_unconf[,,2,] == 0))
  #check triage_D
  expect_true(all(results$I_triage_D_conf[,,1,] == 0))
  expect_true(all(results$I_triage_D_conf[2:t_max,,2,] == results$I_triage_D_unconf[1:(t_max-1),,1,]))
  expect_true(all(results$I_triage_D_unconf[,,2,] == 0))
  #check ICU_R
  expect_true(all(results$I_ICU_R_conf[2,,2,] == results$I_ICU_R_unconf[1,,1,]))
  expect_true(all(results$I_ICU_R_conf[2:t_max,,1,] == results$I_triage_R_conf[1:(t_max-1),,2,]))
  expect_true(all(results$I_ICU_R_unconf[,,2,] == 0))
  #check ICU_D
  expect_true(all(results$I_ICU_D_conf[2,,2,] == results$I_ICU_D_unconf[1,,1,]))
  expect_true(all(results$I_ICU_D_conf[2:t_max,,1,] == results$I_triage_D_conf[1:(t_max-1),,2,]))
  expect_true(all(results$I_ICU_D_unconf[,,2,] == 0))
  #check stepdown
  expect_true(all(results$R_stepdown_conf[2,,2] == results$R_stepdown_unconf[1,,1]))
  I_ICU_R_conf <- apply(results$I_ICU_R_conf,c(1,2,3),sum)
  expect_true(all(results$R_stepdown_conf[2:t_max,,1] == I_ICU_R_conf[1:(t_max-1),,2]))
  expect_true(all(results$R_stepdown_unconf[,,2] == 0))
}
)