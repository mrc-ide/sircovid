##' Create parameters for use with the model
##' 
##' @title Create parameters
##' 
##' @param transmission_model Model type. Only 'POLYMOD' currently supported
##' 
##' @param country Country name
##' 
##' @param severity_data_file Location of file with severity data
##' 
##' @param age_limits Vector of age bin (starting age of each bin)
##' 
##' @param progression_groups List of number of progression groups in each partition
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
##' @param use_polymod_pop Set to ignore \code{survey_pop_in}
##'   and use the population from polymod when estimating the
##'   transmission matrix
##' 
##' @export
##' @import socialmixr
generate_parameters <- function(
  transmission_model = "POLYMOD",
  country="United Kingdom",
  severity_data_file=NULL,
  age_limits = seq(0, 80, by = 5),
  progression_groups = list(E = 2, asympt = 2, mild = 2, ILI = 2, hosp = 2, ICU = 2, rec = 2),
  gammas = list(E = 1/(4.59/2), asympt = 1, mild = 1, ILI = 1, hosp = 2/7, ICU = 2/7, rec = 2/3),
  infection_seeding = c(0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  beta = c(0.1, 0.1, 0.1),
  beta_times = c("2020-02-02", "2020-03-01", "2020-04-01"),
  trans_profile = c(0.65, 0.2, 0.15),
  trans_increase = c(0, 1, 10),
  hosp_transmission = 0.1,
  ICU_transmission = 0.05,
  dt = 0.25,
  use_polymod_pop = FALSE) {
  #
  # Input checks
  #
  # Currently only POLYMOD possible, throws otherwise
  if (length(trans_profile) != length(trans_increase)) {
    stop("Lengths of transmission class arguments mismatching")
  }
  
  if (length(beta) != length(beta_times)) {
    stop("Length of beta mismatching with length of transition times")
  }
  
  if (length(age_limits) != length(infection_seeding)) {
    stop("Length of age bins mismatches length of infection seeds")
  }
  
  if (any(!is.wholenumber(age_limits))){
    stop("Not yet implemented decimal age bins")
  }
  
  if (sum(trans_profile) != 1) {
    stop("trans_profile proportions must sum to 1")
  }
  
  if (transmission_model == "POLYMOD") {
    contact_survey = socialmixr::polymod
  } else {
    stop("Only POLYMOD transmission model implemented")
  }
  
  partitions <- c("E", "asympt", "mild", "ILI", "hosp", "ICU", "rec")
  if (any(!(partitions %in% names(progression_groups)))) {
    stop("progression_groups need to be defined for all partitions")
  }
  if (any(!(partitions %in% names(gammas)))) {
    stop("gammas need to be defined for all partitions")
  }
  
  #
  # Set up time-varying beta
  # Times are in days from first day supplied
  beta_dates <- as.Date(beta_times)
  beta_dates <- as.numeric(beta_dates - beta_dates[[1]])
  # Checks all dates are positive and ascending
  if (any(diff(beta_dates) < 0)) {
    stop("Supplied dates are not increasing")
  }
  
  # 
  # This section defines proportions between partitions
  # derived from the Final_COVID_severity.csv file
  #
  severity_params <- read_severity(
    severity_file_in = severity_data_file,
    age_limits = age_limits)
  
  #
  # Set up the transmission matrix
  #
  survey_pop <- get_survey_pop(severity_params$population, age_limits)
  if (use_polymod_pop) {
    transmission_matrix <- get_transmission_matrix(
      survey_pop = NULL,
      contact_survey = contact_survey,
      country = "United Kingdom",
      age_limits = age_limits)
  } else {
    transmission_matrix <- get_transmission_matrix(
      survey_pop = survey_pop,
      contact_survey = contact_survey,
      country = "United Kingdom",
      age_limits = age_limits)
  }

  #
  # Set up the heterogeneous offspring distribution 
  #
  N_age_bins <- length(age_limits)
  N_trans_classes <- length(trans_profile)

  # Makes an array with rows as age bins, columns as classes
  # All age bins have the same transmission changes
  trans_profile_array <- array(unlist(lapply(trans_profile, rep, N_age_bins)), 
                         c(N_age_bins, N_trans_classes))
  trans_increase_array <- array(unlist(lapply(trans_increase, rep, N_age_bins)), 
                          c(N_age_bins, N_trans_classes))

  #
  # Set the initial conditions for each partition
  # S0 = N, everything else zero
  #
  S0 <- severity_params$population
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
                         beta = beta[1], # When odin code is updated, set this list
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

## Check age bins match with input file - MESSAGE if mismatch
## TODO - could average these when alternative bins are provided
match_age_bins <- function(age_limits, age_headers, severity_file) {
  bin_start <- gsub("(\\d+) to (\\d+)", "\\1", age_headers)
  bin_start <- as.numeric(bin_start)
  bin_end <- gsub("(\\d+) to (\\d+)", "\\2", age_headers)
  bin_end <-  as.numeric(bin_end)

  if (any(bin_start != age_limits)) {
    warning_message <- paste0("Passed age bins do not match those in ",
                              severity_file)
    stop(warning_message)
  }
  if (any(head(bin_end, -1)+1 != bin_start[-1])) {
    warning_message <- "Passed age bins intervals do not overlap correctly"
    message(warning_message)
  }
  
  bin_start
}

## Sets proportion parameters using severity CSV file
read_severity <- function(severity_file_in = NULL, age_limits) {
  if (is.null(severity_file_in)) {
    severity_file <- sircovid_file("extdata/severity.csv")
  } else {
    stop("This is so unlikely to work that we will look at it later")
  }

  # Set up severity file into table
  severity_params <- read.csv(severity_file, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(severity_params)[1] <- "age"

  # Transpose so columns are parameters, rownames are age groups
  severity_data <- t(as.matrix(severity_params[-1]))
  colnames(severity_data) <- severity_params[[1]]
  severity_data <- cbind(age = rownames(severity_data),
             data.frame(severity_data, check.names = FALSE),
             stringsAsFactors = FALSE)
  rownames(severity_data) <- NULL

  population <- severity_data[["Size of England population"]]
  
  # Check passed bins match those in file
  age_bin_ends <- match_age_bins(age_limits, severity_data[["age"]], severity_file)
  severity_data[["age"]] <- age_bin_ends

  prop_symp_seek_HC <- severity_data[["Proportion of symptomatic cases seeking healthcare"]]
  
  # Proportion of asymptomatic
  p_asympt <- 1 - severity_data[["Proportion with symptoms"]]

  #Proportion seeking healthcare
  p_sympt_ILI <- severity_data[["Proportion with symptoms"]] *
    prop_symp_seek_HC
  
  p_recov_ICU <-
    1 - severity_data[["Proportion of critical cases dying"]]

  p_recov_ILI <- 1 - severity_data[["Proportion of symptomatic cases hospitalised"]] /
    prop_symp_seek_HC
  
  p_recov_hosp <- 
    1 - severity_data[["Proportion of hospitalised cases needing critical care"]] -
        severity_data[["Proportion of non-critical care cases dying"]]
  
  #Proportion of hospitalised cases who die without receiveing critical care
  p_death_hosp <- severity_data[["Proportion of non-critical care cases dying"]]
  
  list(
    population = population, # TODO: should be integers
    asympt = p_asympt,
    sympt_ILI = p_sympt_ILI,
    recov_ICU = p_recov_ICU,
    recov_ILI = p_recov_ILI,
    recov_hosp = p_recov_hosp,
    death_hosp = p_death_hosp)
}

## Gets the population age distribution
get_survey_pop <- function(survey_pop_in, age_bins) {
  # Use the default, if no contact survey passed
  if (is.null(survey_pop_in)) {
    survey_pop <- default_age_distribution()
  } else {
    survey_pop <- data.frame(lower.age.limit = age_bins,
                             population = survey_pop_in)
  }
  survey_pop
}

## Gets a transmission matrix from requested model
get_transmission_matrix <- function(
  survey_pop = NULL,
  contact_survey = socialmixr::polymod,
  country = "United Kingdom",
  age_limits = seq(0, 80, 5)) {
  # Get the contact matrix from socialmixr; polymod only
  # goes up to age 70
  if (is.null(survey_pop)) {
    contact_matrix <- socialmixr::contact_matrix(
      contact_survey,
      countries = country,
      age.limits = age_limits[age_limits <= 70],
      symmetric = TRUE)
  } else {
    contact_matrix <- socialmixr::contact_matrix(
      contact_survey,
      survey.pop = survey_pop,
      countries = country,
      age.limits = age_limits[age_limits <= 70],
      symmetric = TRUE)
  }
  
  #transform the matrix to the (symetrical) transmission matrix rather than the contact matrix
  transmission_matrix <- 
    contact_matrix$matrix /
    rep(contact_matrix$demography$population,
        each = ncol(contact_matrix$matrix))
  
  # POLYMOD has a max age of 70, so older bins need to be filled in
  # assumes that the probability of contact remains as in POLYMOD
  # and that contacts are the same in 70+ and 80+
  if (any(age_limits > 70)) {
    for (extra_age in age_limits[age_limits > 70]) {
      transmission_matrix <- cbind(transmission_matrix, transmission_matrix[,nrow(transmission_matrix)])
      transmission_matrix <- rbind(transmission_matrix, transmission_matrix[nrow(transmission_matrix),])
      
      # Adds columns in the format [70,80),80+
      n <- ncol(transmission_matrix)
      prev_end_age <- sub("\\+$", "", colnames(transmission_matrix)[n - 1L])
      new_colnames <- c(sprintf("[%s, %s)", prev_end_age, extra_age),
                        sprintf("%s+", extra_age))
      colnames(transmission_matrix)[n - (1:0)] <- new_colnames
    }
  }
  
  transmission_matrix
}

# from is.integer() docs
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
