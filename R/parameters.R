##' Create parameters for use with the model
##' 
##' @title Create parameters
##' 
##' @param transmission_model Model type. Only 'POLYMOD' currently supported
##' 
##' @param country Country name
##' 
##' @param age_limits Vector of age bin (starting age of each bin)
##' 
##' @param slopes List of number of progression groups in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp', 'ICU', 'rec'
##'   
##' @param gammas List of exponential distribution rates for time in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp', 'ICU', 'rec'
##' 
##' @param infection_seeding Vector of initial cases in each age bin. These are seeded
##'   into the most infective group
##' 
##' @param beta Beta, for each time step in \code{beta_times}
##' 
##' @param beta_times Dates \code{beta} changes, in format yyyy-mm-dd
##' 
##' @param trans_profile Proportion in each infectivity group
##' 
##' @param trans_increase Relative infectivity of each group
##' 
##' @param hosp_transmission Transmission rate of hospital cases
##'
##' @param ICU_transmission Tranmissions rate of ICU cases
##' 
##' @param dt Time-step to run the model in days
##' 
##' @param survey_pop A survey population. Set to NULL to
##'   use \code{default_age_distribution()}
##' 
##' @export
##' @import socialmixr
generate_parameters <- function(
  transmission_model = "POLYMOD",
  country="United Kingdom",
  age_limits = c(0, 10, 20, 30, 40, 50, 60, 70, 80),   #check about 10 etc. being in first or second age group
  progression_groups = list(E = 2, asympt = 2, mild = 2, ILI = 2, hosp = 2, ICU = 2, rec = 2),
  gammas = list(E = 1/(4.59/2), asympt = 1, mild = 1, ILI = 1, hosp = 2/7, ICU = 2/7, rec = 2/3),
  infection_seeding = c(0, 0, 0, 10, 0, 0, 0, 0, 0),
  beta = c(0.1, 0.1, 0.1),
  beta_times = c('2020-02-01', '2020-03-01', '2020-04-01'),
  trans_profile = c(0.65, 0.2, 0.15),
  trans_increase = c(0, 1, 10),
  hosp_transmission = 0.1,
  ICU_transmission = 0.05,
  dt = 0.25,
  survey_pop = NULL
  ) 
{
  #
  # Input checks
  #
  # Currently only POLYMOD possible, throws otherwise
  if (length(trans_profile) != length(trans_increase)) {
    stop("Lengths of transmission class arguments mismatching")
  }
  
  if (length(beta) != length(beta_times))
  {
    stop("Length of beta mismatching with length of transition times")
  }
  
  if (length(age_limits) != length(infection_seeding))
  {
    stop("Length of age bins with length of infection seeds")
  }
  
  if (sum(trans_profile) != 1)
  {
    stop("trans_profile proportions must sum to 1")
  }
  
  if (transmission_model == "POLYMOD") {
    contact_survey = socialmixr::polymod
  } else {
    stop("Only POLYMOD transmission model implemented")
  }
  
  partitions = c('E', 'asympt', 'mild', 'ILI', 'hosp', 'ICU', 'rec')
  if (any(!(partitions %in% names(slopes)))) {
    stop("Slopes need to be defined for all partitions")
  }
  if (any(!(partitions %in% names(gammas)))) {
    stop("Gammas need to be defined for all partitions")
  }
  
  # 
  # This section defines proportions between partitions
  # derived from the Final_COVID_severity.csv file
  #
  severity_params <- read_severity(
    severity_file = "extdata/Final_COVID_severity.csv",
    age_limits = age_limits
  )
  
  #
  # Set up the transmission matrix
  #
  survey_pop <- get_survey_pop(survey_pop)
  transmission_matrix <- get_transmission_matrix(
    survey_pop = survey_pop,
    contact_survey = contact_survey,
    country = "United Kingdom",
    age_limits = age_limits
  )

  #
  # Set up the heterogeneous offspring distribution 
  #
  N_age_bins = length(age_limits)
  N_trans_classes <- length(trans_profile)

  # Makes an array with rows as age bins, columns as classes
  # All age bins have the same transmission changes
  trans_profile_array <- array(unlist(lapply(trans_profile, rep, N_age_bins)), 
                         c(N_age_bins, N_trans_classes))
  trans_increase_array <- array(unlist(lapply(trans_increase, rep, N_age_bins)), 
                          c(N_age_bins, N_trans_classes))

  #
  # Set up time-varying beta
  # Times are in days from first day supplied
  beta_dates = sapply(beta_times, as.Date) - as.numeric(as.Date(beta_times[1]))
  # Checks all dates are positive and ascending
  if (any(beta_dates < 0) || 
      any((beta_dates[-1] - head(beta_dates, n=-1) < 0)))
  {
    stop("Supplied dates not increasing")
  }

  #
  # Set the initial conditions for each partition
  # S0 = N, everything else zero
  #
  S0 <- survey_pop$population
  E0 <- array(0, dim = c(N_age_bins, progression_groups$E, N_trans_classes))
  I0_asympt <- array(0, dim = c(N_age_bins, progression_groups$asympt, N_trans_classes))
  I0_mild <- array(0, dim = c(N_age_bins, progression_groups$mild, N_trans_classes))
  I0_ILI <- array(0, dim = c(N_age_bins, progression_groups$ILI, N_trans_classes))
  I0_hosp <- array(0, dim = c(N_age_bins, progression_groups$hosp, N_trans_classes))
  I0_ICU <- array(0, dim = c(N_age_bins, progression_groups$ICU, N_trans_classes))
  R0_hosp <- array(0, dim = c(N_age_bins, progression_groups$rec, N_trans_classes))
  R0 <- rep(0, N_age_bins)
  D0 <- rep(0, N_age_bins)
  
  # Add some initial infections to I_asympt
  # Move from susceptible to infected
  I0_asympt[,1,N_trans_classes] <- infection_seeding
  S0 <- S0 - infection_seeding

  #
  # Returns parameters
  #
  parameter_list <- list(N_age = N_age_bins,
                         trans_classes = N_trans_classes,
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
                         trans_increase = trans_increase_array,
                         trans_profile = trans_profile_array,
                         beta_list = beta,
                         beta = beta[1], # When odin is updated, set this list
                         beta_dates = beta_dates,
                         s_E = progression_groups$E,
                         gamma_E = gammas$E,
                         s_asympt = progression_groups$asympt,
                         gamma_asympt = gammas$asympt,
                         s_mild = progression_groups$mild,
                         gamma_mild = gammas$mild,
                         s_ILI = progression_groups$ILI,
                         gamma_ILI = gammas$ILI,
                         s_hosp = progression_groups$hosp,
                         gamma_hosp = gammas$hosp,
                         s_ICU = progression_groups$ICU,
                         gamma_ICU = gammas$ICU,
                         s_rec = progression_groups$rec,
                         gamma_rec = gammas$rec,
                         m = transmission_matrix,
                         p_recov_hosp = severity_params$recov_hosp,
                         p_death_hosp = severity_params$death_hosp,
                         p_recov_ILI = severity_params$recov_ILI,
                         p_recov_ICU = severity_params$recov_ICU,
                         p_asympt = severity_params$asympt,
                         p_sympt_ILI = severity_params$sympt_ILI,
                         hosp_transmission = hosp_transmission,
                         ICU_transmission = ICU_transmission)

  parameter_list
}

#
# Internal functions
#

## Check age bins match with input file - WARNING if mismatch
## TODO - could average these when alternative bins are provided
match_age_bins <- function(
  age_limits,
  column_headers,
  severity_file
) 
{
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
  
  bin_start
}

## Sets proportion parameters using severity CSV file
read_severity <- function(
  severity_file = "extdata/Final_COVID_severity.csv",
  age_limits = c(0, 10, 20, 30, 40, 50, 60, 70, 80)
) 
{
  # Set up severity file into table
  # Use readr to avoid horrible column names
  severity_file_in <- sircovid_file(severity_file)
  severity_params <- readr::read_csv(file = severity_file_in, col_types = readr::cols())
  colnames(severity_params)[1] <- "age"
  
  # Transpose so columns are parameters, rownames are age groups
  severity_params <- t(severity_params)
  colnames(severity_params) <- severity_params[1,]
  severity_params <- severity_params[-1,]

  # Check passed bins match those in file
  age_bin_ends <- match_age_bins(age_limits, rownames(severity_params), severity_file)
  rownames(severity_params) <- age_bin_ends
  
  # Proportion of symptomatic
  p_sympt <- 1 - as.numeric(severity_params[,'Proportion with any symptoms'])

  #Proportion seeking healthcare
  p_sympt_ILI <- as.numeric(severity_params[,'Proportion with any symptoms']) *
    as.numeric(severity_params[,'Proportion of infections needing critical care'])
  
  p_recov_ICU <-
    1 - as.numeric(severity_params[,'Proportion of hospitalised cases needing critical care'])
  
  #Proportion of ILI who recover without hospitalisation
  p_recov_ILI <- 
    1 - as.numeric(severity_params[,'Proportion of symptomatic cases hospitalised'])
  
  #Proportion of hospitalised cases who recover without needing ICU
  p_recov_hosp <- 
    1 - as.numeric(severity_params[,"Proportion of cases seeking healthcare who are hospitalised"])
      - as.numeric(severity_params[,"Proportion of critical cases dying"])
  
  #Proportion of hospitalised cases who die without receiveing critical care
  p_death_hosp <- as.numeric(severity_params[,"Proportion of critical cases dying"])
  
  props = list(
    sympt = p_sympt,
    sympt_ILI = p_sympt_ILI,
    recov_ICU = p_recov_ICU,
    recov_ILI = p_recov_ILI,
    recov_hosp = p_recov_hosp
    
  )
}

## Gets the population age distribution
get_survey_pop <- function(survey_pop_in) {
  # Use the default, if no contact survey passed
  if (is.null(survey_pop_in)) {
    survey_pop <- default_age_distribution()
  } else {
    # TODO this may need looking at for other age bins
    survey_pop <-
      survey_pop_in[survey_pop_in$lower.age.limit %in%
                   c(0, 10, 20, 30, 40, 50, 60, 70), ]
  }
  survey_pop
}

## Gets a transmission matrix from requested model
get_transmission_matrix <- function(
  survey_pop,
  contact_survey = socialmixr::polymod,
  country = "United Kingdom",
  age_limits = c(0, 10, 20, 30, 40, 50, 60, 70, 80)
) 
{
  # Get the contact matrix from socialmixr; polymod only
  # goes up to age 70
  contact_matrix <- socialmixr::contact_matrix(
    contact_survey,
    countries = country,
    age.limits = age_limits[age_limits <= 70],
    survey.pop = survey_pop,
    symmetric = TRUE
  )
  
  #transform the matrix to the (symetrical) transmission matrix rather than the contact matrix
  transmission_matrix <- t(t(contact_matrix$matrix)/contact_matrix$demography$population)
  
  #assumes that the probability of contact remains as in POLYMOD
  #and that contacts are the same in 70+ and 80+
  if (any(age_limits > 70))
  {
    for (extra_age in age_limits[age_limits > 70]) 
    {
      m <- cbind(transmission_matrix, transmission_matrix[,nrow(transmission_matrix)])
      m <- rbind(transmission_matrix, transmission_matrix[nrow(transmission_matrix),])
      
      # Adds columns in the format [70,80),80+
      prev_end_age <- stringr::str_match(tail(colnames(transmission_matrix), n = 1), '(\\d+)\\+')[1,2]
      new_colnames <- c(paste0("[", prev_end_age, ",", extra_age, ")"), paste0(extra_age, "+"))
      colnames(transmission_matrix) <- c(colnames(transmission_matrix)[1:(ncol(transmission_matrix) - 2)], 
                                         new_colnames)
    }
  }
  
  transmission_matrix
}



