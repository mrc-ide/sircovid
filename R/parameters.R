##' Create parameters for use with the model
##' @title Create parameters
##' @param transmission_model Model type
##' @param country Country name
##' @param age_limits Vector of age
##' @param progression_parameters Progression parameters
##' @param beta Beta, obvs.
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
  survey_pop = NULL,
  output_parameter_table = TRUE){

  N_age <- length(age_limits)

  if(progression_parameters == "SPI-M-Feb-2009"){
    path <- sircovid_file("data/Final_COVID_severity.csv")
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


    ## proportion based on demographic data from 2019
    proportion_70_80_vs_80_plus <-
        survey_pop$population[survey_pop$lower.age.limit == 70] /
        (survey_pop$population[survey_pop$lower.age.limit == 70] +
         survey_pop$population[survey_pop$lower.age.limit == 80])

    ## if you want to use the latest UK demographic data
    if (! is.null(survey_pop)) {

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

    } else {

        c_m <- contact_matrix(
            socialmixr::polymod,
            countries = country,
            age.limits = c(0, 10, 20, 30, 40, 50, 60, 70),
            symmetric = TRUE
        )
    }
    ## This is adjusting for the fact that socialmixr doesn't have
    ## contact data on 80+ for UK.
    m <- rbind(c_m$matrix,c_m$matrix[8,])
    m <- cbind(m, m[,8])
    m[8,] <- m[8,] * proportion_70_80_vs_80_plus
    m[9,] <- m[9,] * (1 - proportion_70_80_vs_80_plus)
    m[,8] <- m[,8] * proportion_70_80_vs_80_plus
    m[,9] <- m[,9] * (1 - proportion_70_80_vs_80_plus)
    m[9,9] <- 0.6 ## this seems random. Check with Marc?

    pop <- c_m$demography$population
    pop <- c(pop, pop[8] * (1 - proportion_70_80_vs_80_plus))
    pop[8] <- pop[8] - pop[9]

    m <- t(t(m)/pop)

  }

  #Set up the heterogeneous offspring distribution #dnbinom(x=0,mu=2.2, size = 0.16)
  trans_classes <- 3
  trans_profile <- array(c(rep(.65,N_age),rep(.2,N_age),rep(.15,N_age)), c(N_age,trans_classes))
  trans_increase <- array(c(rep(0,N_age),rep(1,N_age),rep(10,N_age)), c(N_age,trans_classes))

  # if(transmission_model=="POLYMOD"){
  #   c_m <- contact_matrix(polymod, countries = country, age.limits = age_limits, symmetric = TRUE)
  #
  # }



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
