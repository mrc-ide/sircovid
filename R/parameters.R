##' Create parameters for use with the model
##' @title Create parameters
##' @param transmission_model Model type
##' @param country Country name
##' @param age_limits Vector of age
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
  beta = c(0.1, 0.1, 0.1),
  beta_times = c('23-Feb-2020', '23-Mar-2020', '28-Apr-2020')
  dt = 0.25,
  survey_pop = NULL,
  output_parameter_table = TRUE){


  severity_params <- read_severity(
    severity_file = "extdata/Final_COVID_severity.csv",
    age_limits = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
    slopes = list(E = 2, asympt = 2, mild = 2, ILI = 2, hosp = 2, ICU = 2, rec = 2),
    gammas = list(E = 1/(4.59/2), asympt = 1, mild = 1, ILI = 1, hosp = 2/7, ICU = 2/7, rec = 2/3),
    survey_pop = NULL
  )

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

# Check age bins match with input file
# Warns if mismatch
# TODO - could average these when alternative bins are provided
match_age <- function(
  age_limits,
  column_headers,
  severity_file
) {
  parsed_header <- column_headers %>% stringr::str_match("(\\d+) to (\\d+)")
  bin_start <- as.numeric(parsed_header[,2])
  bin_start[-1] <- bin_start[-1] - 1 # offset of one in definitions, except for 0
  bin_end <-  as.numeric(parsed_header[,3])
  
  if (any(bin_start != age_limits)) {
    warning_message <- paste0('Passed age bins do not match those in ',
                              severity_file)
    warning(warning_message)
  }
  if (any(head(bin_end, -1)+1 != bin_start[-1]))
  {
    warning_message <- 'Passed age bins do not overlap correctly'
    warning(warning_message)
  }
  
  return(bin_start)
}

read_severity <- function(
  severity_file = "extdata/Final_COVID_severity.csv",
  age_limits = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
  slopes = list(E = 2, asympt = 2, mild = 2, ILI = 2, hosp = 2, ICU = 2, rec = 2),
  gammas = list(E = 1/(4.59/2), asympt = 1, mild = 1, ILI = 1, hosp = 2/7, ICU = 2/7, rec = 2/3),
  survey_pop = NULL
) {
  
  # Set up severity file into table
  severity_file_in <- sircovid_file(severity_file)
  severity_params <- readr::read_csv(file = severity_file)
  colnames(severity_params)[1] <- "age"
  
  # Transpose so columns are parameters, rownames are age groups
  severity_params <- t(severity_params)
  colnames(severity_params) <- severity_params[1,]
  severity_params <- severity_params[-1,]

  # Check passed bins match those in file
  rownames(severity_params) = match_age(age_limits, rownames(severity_params), severity_file)
  
  # Proportion of symptomatic
  p_sympt <- 1 - as.numeric(severity_params[,'Proportion with any symptoms'])

  #Proportion seeking healthcare
  p_sympt_ILI <- as.numeric(severity_params[,'Proportion with any symptoms']) *
    as.numeric(severity_params[,'Proportion of infections needing critical care'])
  
  # Latent period (default mu= 4.59, k=2)
  s_E <- slopes$E
  gamma_E <- gamma$E # QUESTION: divisor of s_E, latent_period_k or 2 correct?
  
  #Parameters of the I_asympt classes
  s_asympt <- slopes$asympt
  gamma_asympt <- gammas$asympt
  
  #Parameters of the I_mild classes
  s_mild <- slopes$mild
  gamma_mild <- gammas$mild
  
  #Parameters of the I_ILI classes
  s_ILI <- slopes$ILI
  gamma_ILI <- gammas$mild
  
  #Parameters of the I_hosp classes
  s_hosp <- slopes$hosp
  gamma_hosp <- gammas$hosp
  
  #Parameters of the I_ICU classes
  s_ICU <- slopes$ICU
  gamma_ICU <- gammas$ICU
  p_recov_ICU <- 1 - 
    as.numeric(severity_params[,'Proportion of hospitalised cases needing critical care'])
  
  #Parameters of the R_hosp classes
  s_rec <- slopes$rec
  gamma_rec <- gammas$rec
  
  #Proportion of ILI who recover without hospitalisation
  p_recov_ILI <- 1 - as.numeric(severity_params[,'Proportion of symptomatic cases hospitalised'])
  
  #Proportion of hospitalised cases who recover without needing ICU
  p_recov_hosp <- 1 - as.numeric(severity_params["Proportion of cases seeking healthcare who are hospitalised"])
                    - as.numeric(severity_params["Proportion of critical cases dying"])
  
  #Proportion of hospitalised cases who die without receiveing critical care
  p_death_hosp <- as.numeric(severity_params["Proportion of critical cases dying"])
  
  #If survey_pop is not passed as an argument, get it from the package
  if(is.null(survey_pop)){
    survey_pop <- default_age_distribution()
    
    pop <- survey_pop$population
    
    #Get the contact matrix from socialmixr; problem no data in POLYMOD
    #for over 80+
    c_m <- socialmixr::contact_matrix(
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
  }
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
  
  c_m <- socialmixr::contact_matrix(
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
