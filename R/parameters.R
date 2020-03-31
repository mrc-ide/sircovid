##' Create parameters for use with the model
##' @title Create parameters
##' @param transmission_model Model type
##' @param country Country name
##' @param age_limits Vector of age
##' @param progression_parameters Progression parameters
##' @param beta Beta, obvs.
##' @param dt Time-step to run the model in days
##' @param survey_pop A survey population perhaps
##' @param output_parameter_table The table of output parameters
##' @export
##' @import socialmixr
parameters <- function(
  transmission_model = "POLYMOD",
  country="United Kingdom",
  #check about 10 etc. being in first or second age group
  age_limits = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
  progression_parameters = "SPI-M-Feb-2009",
  beta = 0.1,
  dt = 0.25,
  survey_pop = NULL,
  output_parameter_table = TRUE){

  N_age <- length(age_limits)

  if(progression_parameters == "SPI-M-Feb-2009"){
    path <- sircovid_file("extdata/Final_COVID_severity.csv")
    prog_par <- utils::read.csv(file = path)
    #prog_par <- as.numeric(prog_par_table[,-1])

    #Age band from SPI-M model
    age.lim <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)

    #Proportion of asymptomatic
    p_asympt <- 1-as.numeric(prog_par[1,2:10])

    #Proportion seeking healthcare
    p_sympt_ILI <- as.numeric(prog_par[1,2:10])*as.numeric(prog_par[6,2:10])

    #Latent period mu= 4.59, k=2
    s_E <- 2
    gamma_E <- 1/(4.59/2)

    #Parameters of the I_asympt classes
    s_asympt <- 2
    gamma_asympt <- 1

    #Parameters of the I_mild classes
    s_mild <- 2
    gamma_mild <- 1

    #Parameters of the I_ILI classes
    s_ILI <- 2
    gamma_ILI <- 1

    #Parameters of the I_hosp classes
    s_hosp <- 2
    gamma_hosp <- 2/7

    #Parameters of the I_ICU classes
    s_ICU <- 2
    gamma_ICU <- 2/7
    p_recov_ICU <- 1-as.numeric(prog_par[10,2:10])

    #Parameters of the R_hosp classes
    s_rec <- 2
    gamma_rec <- 2/3

    #Proportion of ILI who recover without hospitalisation
    p_recov_ILI <- 1-as.numeric(prog_par[8,2:10])

    #Proportion of hospitalised cases who recover without needing ICU
    p_recov_hosp <- 1-as.numeric(prog_par[9,2:10])-as.numeric(prog_par[11,2:10])

    #Proportion of hospitalised cases who die without receiveing critical care
    p_death_hosp <- as.numeric(prog_par[11,2:10])

    #If survey_pop is not passed as an argument, get it from the package
    if(is.null(survey_pop)){
      survey_pop <- default_age_distribution()
      
      pop <- survey_pop$population
      
      #Get the contact matrix from socialmixr; problem no data in POLYMOD
      #for over 80+
      c_m <- contact_matrix(
        socialmixr::polymod,
        countries = country,
        age.limits = c(0, 10, 20, 30, 40, 50, 60, 70),
        symmetric = TRUE
      )
      
      #transform the matrix in (symetrical) transmission matrix rather than the contact matrix
      m <- t(t(c_m$matrix)/c_m$demography$population)
      
      #assumes that the probability of contact remains as in POLYMOD
      #and that contacts are the same in 70+ and 80+
      m <- cbind(m,m[,8])
      m <- rbind(m,m[8,])
      names_m_col <- colnames(m)
      colnames(m) <- c(names_m_col[1:7],"[70,80)","80+")
    } else
    {
      ##TODO this needs to be rewritten when a survey_pop is passed as an argument
      survey_pop_subset <-
        survey_pop[survey_pop$lower.age.limit %in%
                     c(0, 10, 20, 30, 40, 50, 60, 70), ]
      
      c_m <- contact_matrix(
        socialmixr::polymod,
        countries = country,
        age.limits = c(0, 10, 20, 30, 40, 50, 60, 70),
        symmetric = TRUE,
        survey.pop = survey_pop_subset
      )
      
      pop <- survey_pop$population
      
      m <- t(t(c_m$matrix)/pop)
      
    }
  
  }

  #Set up the heterogeneous offspring distribution #dnbinom(x=0,mu=2.2, size = 0.16)
  trans_classes <- 3
  trans_profile <- array(c(rep(.65,N_age),rep(.2,N_age),rep(.15,N_age)), c(N_age,trans_classes))
  trans_increase <- array(c(rep(0,N_age),rep(1,N_age),rep(10,N_age)), c(N_age,trans_classes))

  #TODO: flexible seeding
  #Set the initial conditions
  S0 <- pop
  E0 <- array(0, dim = c(N_age,s_E,trans_classes))
  I0_asympt <- array(0, dim = c(N_age,s_asympt,trans_classes))
  seed_SSP <- 10
  I0_asympt[4,1,3] <- seed_SSP
  S0[4] <- S0[4] - seed_SSP
  I0_mild <- array(0, dim = c(N_age,s_mild,trans_classes))
  I0_ILI <- array(0, dim = c(N_age,s_ILI,trans_classes))
  I0_hosp <- array(0, dim = c(N_age,s_hosp,trans_classes))
  I0_ICU <- array(0, dim = c(N_age,s_ICU,trans_classes))
  R0_hosp <- array(0, dim = c(N_age,s_rec,trans_classes))
  R0 <- rep(0,N_age)
  D0 <- rep(0,N_age)

  hosp_transmission <- 0.1
  ICU_transmission <- 0.05

  parameter_list <- list(N_age = N_age,
                       trans_classes = trans_classes,
                       dt = dt,
                       S0 = S0,
                       E0 = E0,
                       I0_asympt = I0_asympt,
                       I0_mild = I0_mild,
                       I0_ILI = I0_ILI,
                       I0_hosp = I0_hosp,
                       I0_ICU = I0_ICU,
                       R0_hosp = R0_hosp,
                       R0 = R0,
                       D0 = R0,
                       trans_increase = trans_increase,
                       trans_profile = trans_profile,
                       beta = beta,
                       s_E = s_E,
                       gamma_E = gamma_E,
                       s_asympt = s_asympt,
                       gamma_asympt = gamma_asympt,
                       s_mild = s_mild,
                       gamma_mild = gamma_mild,
                       s_ILI = s_ILI,
                       gamma_ILI = gamma_ILI,
                       s_hosp = s_hosp,
                       gamma_hosp = gamma_hosp,
                       s_ICU = s_ICU,
                       gamma_ICU = gamma_ICU,
                       s_rec = s_rec,
                       gamma_rec = gamma_rec,
                       m = m,
                       p_recov_hosp = p_recov_hosp,
                       p_death_hosp = p_death_hosp,
                       p_recov_ILI = p_recov_ILI,
                       p_recov_ICU = p_recov_ICU,
                       hosp_transmission = hosp_transmission,
                       ICU_transmission = ICU_transmission,
                       p_asympt = p_asympt,
                       p_sympt_ILI = p_sympt_ILI
  )


  return(parameter_list)
}


##' Parameters for the rtm
##' @title Parameters for the rtm
##' @param beta Beta
##' @param seed_SSP seed
##' @param dt dt
##' @export
generate_parameter_rtm <- function(
  beta = 0.1,
  seed_SSP = 10,
  dt = 0.25){
  N_age <- 17
  
  age_lim <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80)
  age_lim_polym <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70)
  
  ## Severity parameters from spreadsheet
  path <- sircovid_file("extdata/severity.csv")
  severity_par <- utils::read.csv(path, header = TRUE, sep = ",")
  
  pop <- as.numeric(severity_par[1,2:18])
  
  c_m <- contact_matrix(
    socialmixr::polymod,
    countries = "United Kingdom",
    age.limits = age_lim_polym,
    symmetric = TRUE
  )
  
  ## transform the matrix in (symetrical) transmission matrix rather
  ## than the contact matrix
  m <- t(t(c_m$matrix)/c_m$demography$population)
  
  ## assumes that the probability of contact remains as in POLYMOD and
  ## that contacts are the same in 70+ and 80+
  m <- cbind(m,m[,15])
  m <- rbind(m,m[15,])
  m <- cbind(m,m[,16])
  m <- rbind(m,m[16,])
  names_m_col <- colnames(m)
  colnames(m) <- c(names_m_col[1:14],"[70,75)","[75,80)","80+")
  
  ## Proportion of asymptomatic
  p_asympt <- 1-as.numeric(severity_par[2,2:18])
  
  ## Carried over from the initial NHS meeting
  prop_symp_seek_HC <- c(0.3570377550, 0.3570377550, 0.3712946230,0.3712946230,	0.420792849,0.420792849,
                         0.459552523,0.459552523,	0.488704572,0.488704572,	0.578769171,0.578769171,	0.65754772,0.65754772,	0.73278164,0.73278164,0.76501082)
  
  #Proportion seeking healthcare
  p_sympt_ILI <- as.numeric(severity_par[2,2:18])*prop_symp_seek_HC
  
  #Latent period to match parameters from Neil's IBM
  s_E <- 2
  gamma_E <- 1/2.5
  
  #Parameters of the I_asympt classes
  s_asympt <- 1
  gamma_asympt <- 1/2.09
  
  #Parameters of the I_mild classes
  s_mild <- 1
  gamma_mild <- 1/2.09
  
  #Parameters of the I_ILI classes
  s_ILI <- 1
  gamma_ILI <- 1/4
  
  #Parameters of the I_hosp classes
  s_hosp <- 2
  gamma_hosp <- 2/1
  
  #Parameters of the I_ICU classes
  s_ICU <- 2
  gamma_ICU <- 2/5
  p_recov_ICU <- 1-as.numeric(severity_par[12,2:18])
  
  #Parameters of the R_hosp classes
  s_rec <- 2
  gamma_rec <- 2/5
  
  #Proportion of ILI who recover without hospitalisation
  p_recov_ILI <- 1-as.numeric(severity_par[10,2:18])/prop_symp_seek_HC
  
  #Proportion of hospitalised cases who recover without needing ICU
  p_recov_hosp <- 1-as.numeric(severity_par[11,2:18])-as.numeric(severity_par[13,2:18])
  
  #Proportion of hospitalised cases who die without receiveing critical care
  p_death_hosp <- as.numeric(severity_par[13,2:18])
  
  #Set up the heterogeneous offspring distribution #dnbinom(x=0,mu=2.2, size = 0.16)
  #trans_classes <- 3
  #trans_profile <- array(c(rep(.65,N_age),rep(.2,N_age),rep(.15,N_age)), c(N_age,trans_classes))
  #trans_increase <- array(c(rep(0,N_age),rep(1,N_age),rep(10,N_age)), c(N_age,trans_classes))
  
  trans_classes <- 1
  trans_profile <- array(c(rep(1,N_age)), c(N_age,trans_classes))
  trans_increase <- array(c(rep(1,N_age)), c(N_age,trans_classes))
  
  #Set the initial conditions
  S0 <- pop
  E0 <- array(0, dim = c(N_age,s_E,trans_classes))
  I0_asympt <- array(0, dim = c(N_age,s_asympt,trans_classes))
  I0_asympt[4,1,1] <- seed_SSP
  S0[4] <- S0[4] - seed_SSP
  I0_mild <- array(0, dim = c(N_age,s_mild,trans_classes))
  I0_ILI <- array(0, dim = c(N_age,s_ILI,trans_classes))
  I0_hosp <- array(0, dim = c(N_age,s_hosp,trans_classes))
  I0_ICU <- array(0, dim = c(N_age,s_ICU,trans_classes))
  R0_hosp <- array(0, dim = c(N_age,s_rec,trans_classes))
  R0 <- rep(0,N_age)
  D0 <- rep(0,N_age)
  
  hosp_transmission <- 0
  ICU_transmission <- 0
  
  parameter_list <- list(N_age = N_age,
                         trans_classes = trans_classes,
                         dt = dt,
                         S0 = S0,
                         E0 = E0,
                         I0_asympt = I0_asympt,
                         I0_mild = I0_mild,
                         I0_ILI = I0_ILI,
                         I0_hosp = I0_hosp,
                         I0_ICU = I0_ICU,
                         R0_hosp = R0_hosp,
                         R0 = R0,
                         D0 = R0,
                         trans_increase = trans_increase,
                         trans_profile = trans_profile,
                         beta = beta,
                         s_E = s_E,
                         gamma_E = gamma_E,
                         s_asympt = s_asympt,
                         gamma_asympt = gamma_asympt,
                         s_mild = s_mild,
                         gamma_mild = gamma_mild,
                         s_ILI = s_ILI,
                         gamma_ILI = gamma_ILI,
                         s_hosp = s_hosp,
                         gamma_hosp = gamma_hosp,
                         s_ICU = s_ICU,
                         gamma_ICU = gamma_ICU,
                         s_rec = s_rec,
                         gamma_rec = gamma_rec,
                         m = m,
                         p_recov_hosp = p_recov_hosp,
                         p_death_hosp = p_death_hosp,
                         p_recov_ILI = p_recov_ILI,
                         p_recov_ICU = p_recov_ICU,
                         hosp_transmission = hosp_transmission,
                         ICU_transmission = ICU_transmission,
                         p_asympt = p_asympt,
                         p_sympt_ILI = p_sympt_ILI
  )
  
  return(parameter_list)
}
