##' Create a linear time-varying beta for input into 
##' \code{create_parameters()}
##' 
##' @title Generate beta
##' 
##' @param beta_start Value of beta before intervention
##' 
##' @param start_date Start date for simualtions
##' 
##' @param reduction_start Start of reduction in beta
##' 
##' @param beta_reduction Factor by which beta is reduced
##'   default is 0.238, as estimated by Jarvis et al
##'   beta_reduced = beta - beta * (1 - beta_reduction)
##' 
##' @param reduction_period Time, in days, over which the
##' reduction in beta is achieved. Default 10 days
##' 
##' @export
generate_beta <- function(beta_start, 
                          start_date = "2020-02-02",
                          reduction_start = "2020-03-16",
                          beta_reduction = 0.238, 
                          reduction_period = 10) {
  if (as.Date(start_date) > as.Date(reduction_start)) {
    stop("Start date must be earlier than intervention date")
  }
  if (beta_reduction < 0 || beta_reduction > 1) {
    stop("beta must decrease, and cannot be reduced below 0")
  }
  if (reduction_period > 100) {
    message("Reduction period over 100 days - is this correct?")
  }
  
  beta_times <- c(as.Date(start_date), 
            seq(as.Date(reduction_start), as.Date(reduction_start) + reduction_period - 1, by=1))
  
  # Corresponding change in beta
  # 
  beta_slope <- beta_start * (1 - (1 - beta_reduction ) * (seq(0, reduction_period - 1) / (reduction_period - 1)))
  
  list(beta=c(beta_start, beta_slope),
       beta_times=as.character(beta_times))
}

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
##' @param progression_groups List of number of progression groups in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp', 'ICU', 'rec'
##'   
##' @param gammas List of exponential distribution rates for time in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp', 'ICU', 'rec'
##' 
##' @param infection_seeding List of vector \code{values} of how many cases to seed 
##'  in a correspond vector of age bins in \code{bins}.
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
##' @return List of parameters for use with \code{sircovid}
##' 
##' @export
##' @import socialmixr
generate_parameters <- function(
  transmission_model = "POLYMOD",
  country="United Kingdom",
  severity_data_file=NULL,
  progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 1, hosp = 2, ICU = 2, rec = 2),
  gammas = list(E = 1/(4.59/2), asympt = 1/2.09, mild = 1/2.09, ILI = 1/4, hosp = 2, ICU = 2/5, rec = 2/5),
  infection_seeding = list(values=c(10),
                           bins=c('15 to 19')),
  beta = c(0.1, 0.1, 0.1),
  beta_times = c("2020-02-02", "2020-03-01", "2020-04-01"),
  trans_profile = c(1),
  trans_increase = c(1),
  hosp_transmission = 0.1,
  ICU_transmission = 0.05,
  dt = 0.25,
  use_polymod_pop = FALSE) {
  #
  # Input checks
  #
  # Currently only POLYMOD possible, throws otherwise
  if (length(trans_profile) != length(trans_increase)) {
    stop("Lengths of transmissibility class arguments mismatching")
  }
  
  if (length(beta) != length(beta_times)) {
    stop("Length of beta mismatching with length of transition times")
  }
  
  if (length(infection_seeding$values) != length(infection_seeding$bins)) {
    stop("Each infection seeding value must correspond to one bin")
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
  beta_t <- normalise_beta(beta_times, dt)

  # 
  # This section defines proportions between partitions
  # derived from the severity.csv file
  #
  severity_params <- read_severity(severity_data_file)
  
  #
  # Set up the transmission matrix
  #
  survey_pop <- get_survey_pop(severity_params$population, severity_params$age_bin_starts)
  if (use_polymod_pop) {
    transmission_matrix <- get_transmission_matrix(
      survey_pop = NULL,
      contact_survey = contact_survey,
      country = "United Kingdom",
      age_limits = severity_params$age_bin_starts)
  } else {
    transmission_matrix <- get_transmission_matrix(
      survey_pop = survey_pop,
      contact_survey = contact_survey,
      country = "United Kingdom",
      age_limits = severity_params$age_bin_starts)
  }

  #
  # Set up the heterogeneous offspring distribution 
  #
  N_age_bins <- length(severity_params$age_bin_starts)
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
  seed_bins <- parse_age_bins(infection_seeding$bins)
  seed_idx <- match(seed_bins$bin_start, severity_params$age_bin_starts)
  if (any(duplicated(seed_idx))) {
    stop("Seeding is into the same bin multiple times")
  }
  I0_asympt[seed_idx,1,N_trans_classes] <- infection_seeding$values
  S0[seed_idx] <- S0[seed_idx] - infection_seeding$values

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
                         D0 = D0,
                         trans_increase = trans_increase_array,
                         trans_profile = trans_profile_array,
                         beta_y = beta,
                         beta_t = beta_t,
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


##' The \code{generate_parameters_new_hospital_model} function is a
##' helper while we develop the new hospital model.
##'
##' @export
##' @rdname generate_parameters
generate_parameters_new_hospital_model <- function() {
  param_mod <- generate_parameters()
  
  ## these next two lines are because scan_beta_date will give an error if beta_t and beta_y are vectors
  param_mod$beta_t <- param_mod$beta_t[1]
  param_mod$beta_y <- param_mod$beta_y[1]

  param_mod$I0_hosp_D <- param_mod$I0_hosp
  param_mod$I0_hosp_R <- param_mod$I0_hosp
  param_mod$I0_ICU_R <- param_mod$I0_hosp
  param_mod$I0_ICU_D <- param_mod$I0_hosp
  param_mod$I0_triage <- param_mod$I0_hosp
  param_mod$R0_stepdown <- param_mod$I0_hosp
  
  param_mod$R0_hosp <- NULL
  param_mod$s_hosp <- NULL
  param_mod$gamma_hosp <- NULL
  param_mod$s_ICU <- NULL
  param_mod$gamma_ICU <- NULL
  param_mod$s_rec <- NULL
  param_mod$gamma_rec <- NULL
  
  param_mod$s_hosp_D <- 2
  param_mod$gamma_hosp_D <- 2/5
  param_mod$s_hosp_R <- 2
  param_mod$gamma_hosp_R <- 2/10
  
  param_mod$s_ICU_D <- 2
  param_mod$gamma_ICU_D <- 2/5
  param_mod$s_ICU_R <- 2
  param_mod$gamma_ICU_R <- 2/10
  
  param_mod$s_stepdown <- 2
  param_mod$gamma_stepdown <- 2/5
  
  param_mod$s_triage <- 2
  param_mod$gamma_triage <- 2
  
  param_mod$p_death_ICU <- 1-param_mod$p_recov_ICU
  param_mod$p_ICU_hosp <- 1- param_mod$p_recov_hosp - param_mod$p_death_hosp
  param_mod$p_recov_ICU <- NULL
  
  param_mod$I0_hosp <- NULL
  param_mod$I0_ICU <- NULL
  param_mod$p_recov_hosp <- NULL
  
  param_mod
}


#
# Internal functions
#

# Generates vector of times for beta changes, in the same
# terms as the odin code
normalise_beta <- function(beta_times, dt) {
  # Times are in days from first day supplied
  beta_dates <- as.Date(beta_times)
  beta_dates <- as.numeric(beta_dates - beta_dates[[1]])
  # Checks all dates are positive and ascending
  if (any(diff(beta_dates) < 0)) {
    stop("Supplied dates are not increasing")
  }
  beta_t <- beta_dates/dt
}


parse_age_bins <- function(age_bin_strings) {
  bin_start <- gsub("(\\d+) to (\\d+)", "\\1", age_bin_strings)
  bin_start <- as.numeric(bin_start)
  bin_end <- gsub("(\\d+) to (\\d+)", "\\2", age_bin_strings)
  bin_end <-  as.numeric(bin_end)
  
  if (any(!is.wholenumber(bin_start)) || any(!is.wholenumber(bin_end))){
    stop("Not yet implemented decimal age bins")
  }
  
  list(bin_start=bin_start,
       bin_end=bin_end)
}


## Check age bins match with input file - MESSAGE if mismatch
## TODO - could average these when alternative bins are provided
check_age_bins <- function(age_headers) {
  bins = parse_age_bins(age_headers)

  if (any(utils::head(bins$bin_end, - 1) + 1 != bins$bin_start[-1])) {
    warning_message <- "Passed age bins intervals do not overlap correctly"
    message(warning_message)
  }
  
  bins
}

## Sets proportion parameters using severity CSV file
read_severity <- function(severity_file_in = NULL, age_limits) {
  if (is.null(severity_file_in)) {
    severity_file <- sircovid_file("extdata/severity.csv")
  } else {
    severity_file <- sircovid_file(severity_file_in)
  }

  # Set up severity file into table
  severity_params <- utils::read.csv(severity_file, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(severity_params)[1] <- "age"

  # Transpose so columns are parameters, rownames are age groups
  severity_data <- t(as.matrix(severity_params[-1]))
  colnames(severity_data) <- severity_params[[1]]
  severity_data <- cbind(age = rownames(severity_data),
             data.frame(severity_data, check.names = FALSE),
             stringsAsFactors = FALSE)
  rownames(severity_data) <- NULL

  # Check file format
  expected_cols <- c(
    "Size of England population",
    "Proportion of symptomatic cases seeking healthcare",
    "Proportion with symptoms",
    "Unadjusted_IFR",
    "Adjusted_IFR to 1%",
    "IFR relative to 80",
    "Age specific scaling of ifr to give hospitalisation",
    "Proportion of infections hospitalised compared to age 80",
    "Proportion of infections hospitalised",
    "Proportion of infections needing critical care",
    "Proportion of symptomatic cases hospitalised",
    "Proportion of hospitalised cases needing critical care",
    "Proportion of critical cases dying",
    "Proportion of non-critical care cases dying")
  if (any(!(expected_cols %in% colnames(severity_data)))) {
    missing <- expected_cols[which(!(expected_cols %in% colnames(severity_data)))]
    error_message <- paste("Could not find the following rows in the severity file:", 
                           missing, sep=" ")
    stop(error_message)
  }
  
  # Parse the age bins. Useful to keep both start and end depending on what
  # function expects
  age_bins <- check_age_bins(severity_data[["age"]])
  
  population <- round(severity_data[["Size of England population"]])

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
    (1 - severity_data[["Proportion of hospitalised cases needing critical care"]]) *
    (1 - severity_data[["Proportion of non-critical care cases dying"]])
  
  #Proportion of hospitalised cases who die without receiveing critical care
  p_death_hosp <- (1 - severity_data[["Proportion of hospitalised cases needing critical care"]]) *
                   severity_data[["Proportion of non-critical care cases dying"]]
  
  list(
    population = population,
    age_bin_starts = age_bins$bin_start,
    age_bin_ends = age_bins$bin_end,
    asympt = p_asympt,
    sympt_ILI = p_sympt_ILI,
    recov_ICU = p_recov_ICU,
    recov_ILI = p_recov_ILI,
    recov_hosp = p_recov_hosp,
    death_hosp = p_death_hosp)
}

## Gets the population age distribution
get_survey_pop <- function(survey_pop_in, age_bins) {
  # Always use the population size from the severity file
  survey_pop <- data.frame(lower.age.limit = age_bins,
                           population = survey_pop_in)

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
