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
##' @param beta_end Overrides beta_reduction, giving the
##'   final value of beta by the end of the reduction
##'   period
##' 
##' @param reduction_period Time, in days, over which the
##' reduction in beta is achieved. Default 10 days
##' 
##' @param beta_pl Value of beta after lockdown is eased
##' 
##' @param pl_start Start of easing of lockdown
##' 
##' @param pl_transition_period Time, in days, over which
##' the change in beta after lockdown easing is achieved.
##' Default 10 days
##' 
##' @export
generate_beta <- function(beta_start, 
                          start_date = sircovid_date("2020-02-02"),
                          reduction_start = sircovid_date("2020-03-16"),
                          beta_reduction = 0.238,
                          beta_end = NULL,
                          reduction_period = 10,
                          beta_pl = NULL,
                          pl_start = sircovid_date("2020-05-11"),
                          pl_transition_period = 7) {

  if (start_date > reduction_start) {
    stop("Start date must be earlier than intervention date")
  }
  if (beta_start < 0){
    stop("beta_start must be non-negative")
  }
  if (!is.null(beta_end)) {
    if (beta_end < 0) {
      stop("beta_end must be non-negative")
    }
  } else {
    if (beta_reduction < 0) {
      stop("beta cannot be reduced below 0")
    }
  }
  if (!is.null(beta_pl)) {
    if (beta_pl < 0) {
      stop("beta_pl must be non-negative")
    }
  }
  
  if (reduction_period > 100) {
    message("Reduction period over 100 days - is this correct?")
  }
  
  beta_times <- c(start_date, seq(reduction_start, reduction_start + reduction_period - 1, by=1))

  # Corresponding change in beta
  if (!is.null(beta_end)) {
    beta_reduction <- beta_end / beta_start
  }
  beta_slope <- beta_start * (1 - (1 - beta_reduction ) * (seq(0, reduction_period - 1) / (reduction_period - 1)))
  
  beta <- c(beta_start, beta_slope)
  
  if (!(is.null(beta_pl))) {
    beta_reduction_pl <- beta_pl / beta_end
    beta_slope_pl <- beta_end * (1 - (1 - beta_reduction_pl ) * (seq(0, pl_transition_period - 1) / (pl_transition_period - 1)))
    
    beta <- c(beta,beta_slope_pl)
    beta_times <- c(beta_times, seq(pl_start, pl_start + pl_transition_period - 1, by=1))
  }
  
  list(beta=beta,
       beta_times=beta_times)
}


##' Create a piecewise linear time-varying beta for input into 
##' \code{create_parameters()}
##' 
##' @title Generate beta piecewise linear
##' 
##' @param beta_k Value of beta at changepoints
##' 
##' @param t_k Time of changepoints
##' 
##' @param start_date Start date for simulations
##' 
##' @param dt step size
##' 
##' @export
generate_beta_piecewise_linear <- function(beta_k,
                                           t_k,
                                           start_date = sircovid_date("2020-02-02"),
                                           dt) {
  
  if (any(start_date > t_k)) {
    stop("Start date must be earlier than dates in t_k")
  }
  if (any(beta_k < 0)){
    stop("beta_k must be non-negative")
  }
  if (is.unsorted(t_k,strictly = TRUE)){
    stop("t_k must be strictly increasing")
  }
  if (!length(t_k)==length(beta_k)){
    stop("t_k and beta_k must be of same length")
  }
  
  beta_times <- c(start_date, seq(t_k[1], t_k[length(t_k)], by=dt))
  
  #set beta as flat until first changepoint
  t_k <- c(start_date,t_k)
  beta_k <- c(beta_k[1],beta_k)
  
  beta_fun <- function(t){
    if (t %in% t_k){
      b <- beta_k[which(t_k==t)]
    } else {
      k <- max(which(t_k<=t))
      b <- (t - t_k[k])/(t_k[k+1] - t_k[k])*beta_k[k+1] + (t_k[k+1] - t)/(t_k[k+1] - t_k[k])*beta_k[k]
    }
    b
  }
  
  beta <- sapply(beta_times,beta_fun)
  
  list(beta=beta,
       beta_times=beta_times)
}


##' Create parameters for use with the model
##' 
##' @title Create parameters
##' 
##' @param sircovid_model Model to use, \code{basic_model()} or
##'   \code{hospital_model()}
##' 
##' @param transmission_model Model type. Only 'POLYMOD' currently supported
##' 
##' @param country Country name
##' 
##' @param severity_data_file Location of file with severity data
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
##' @param comm_D_transmission Transmission rate of cases dying in the community
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
generate_parameters <- function(
  sircovid_model = basic_model(),
  transmission_model = "POLYMOD",
  country="United Kingdom",
  severity_data_file=NULL,
  infection_seeding = list(values=c(10),
                           bins=c('15 to 19')),
  beta = 0.1,
  beta_times = sircovid_date("2020-02-02"),
  trans_profile = c(1),
  trans_increase = c(1),
  hosp_transmission = 0.1,
  ICU_transmission = 0.05,
  comm_D_transmission = 0.05,
  dt = 0.25,
  use_polymod_pop = FALSE) {
  
  if (length(infection_seeding$values) != length(infection_seeding$bins)) {
    stop("Each infection seeding value must correspond to one bin")
  }
  
  # Generate parameters used by all models
  parameter_list <- generate_parameters_base(transmission_model = transmission_model,
                                             country=country,
                                             severity_data_file=severity_data_file,
                                             beta = beta,
                                             beta_times = beta_times,
                                             trans_profile = trans_profile,
                                             trans_increase = trans_increase,
                                             dt = dt,
                                             use_polymod_pop = use_polymod_pop)
  
  #
  # Set the initial conditions for each partition
  # S0 = N (set in generate_parameters_base), everything else zero
  #
  if ("sircovid_basic" %in% class(sircovid_model)) {
    parameter_list$E0 <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$E, parameter_list$trans_classes))
    parameter_list$I0_asympt <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$asympt, parameter_list$trans_classes))
    parameter_list$I0_mild <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$mild, parameter_list$trans_classes))
    parameter_list$I0_ILI <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$ILI, parameter_list$trans_classes))
    parameter_list$I0_hosp <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$hosp, parameter_list$trans_classes))
    parameter_list$I0_ICU <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$ICU, parameter_list$trans_classes))
    parameter_list$R0_hosp <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$rec, parameter_list$trans_classes))
    parameter_list$R0 <- rep(0, parameter_list$N_age)
    parameter_list$D0 <- rep(0, parameter_list$N_age)
  } else if ("sircovid_hospital" %in% class(sircovid_model)) {
    parameter_list$E0 <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$E, parameter_list$trans_classes))
    parameter_list$I0_asympt <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$asympt, parameter_list$trans_classes))
    parameter_list$I0_mild <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$mild, parameter_list$trans_classes))
    parameter_list$I0_ILI <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$ILI, parameter_list$trans_classes))
    parameter_list$R0 <- rep(0, parameter_list$N_age)
    if ("sircovid_serology" %in% class(sircovid_model)) {
      parameter_list$R0_stepdown_unconf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$stepdown))
      parameter_list$R0_stepdown_conf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$stepdown))
      parameter_list$R0_neg <- rep(0, parameter_list$N_age)
      parameter_list$R0_pre <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$R_pre))
      parameter_list$R0_pos <- rep(0, parameter_list$N_age)
      parameter_list$I0_comm_D <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$comm_D, parameter_list$trans_classes))
      parameter_list$I0_hosp_D_unconf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$hosp_D, parameter_list$trans_classes))
      parameter_list$I0_hosp_D_conf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$hosp_D, parameter_list$trans_classes))
      parameter_list$I0_hosp_R_unconf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$hosp_R, parameter_list$trans_classes))
      parameter_list$I0_hosp_R_conf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$hosp_R, parameter_list$trans_classes))
      parameter_list$I0_triage_R_unconf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$triage, parameter_list$trans_classes))
      parameter_list$I0_triage_R_conf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$triage, parameter_list$trans_classes))
      parameter_list$I0_triage_D_unconf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$triage, parameter_list$trans_classes))
      parameter_list$I0_triage_D_conf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$triage, parameter_list$trans_classes))
      parameter_list$I0_ICU_R_unconf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$ICU_R, parameter_list$trans_classes))
      parameter_list$I0_ICU_R_conf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$ICU_R, parameter_list$trans_classes))
      parameter_list$I0_ICU_D_unconf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$ICU_D, parameter_list$trans_classes))
      parameter_list$I0_ICU_D_conf <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$ICU_D, parameter_list$trans_classes))
      parameter_list$D0_hosp <- rep(0, parameter_list$N_age)
      parameter_list$D0_comm <- rep(0, parameter_list$N_age)
      parameter_list$cum0_admit_conf <- 0
      parameter_list$cum0_new_conf <- 0
    } else {
      parameter_list$R0_stepdown <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$stepdown, parameter_list$trans_classes))
      parameter_list$I0_triage <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$triage, parameter_list$trans_classes))
      parameter_list$I0_hosp_D <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$hosp_D, parameter_list$trans_classes))
      parameter_list$I0_hosp_R <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$hosp_R, parameter_list$trans_classes))
      parameter_list$I0_ICU_R <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$ICU_R, parameter_list$trans_classes))
      parameter_list$I0_ICU_D <- array(0, dim = c(parameter_list$N_age, sircovid_model$progression_groups$ICU_D, parameter_list$trans_classes))
      parameter_list$D0 <- rep(0, parameter_list$N_age)
    }
  } else {
    stop("Model name not supported")
  }
  
  # Add some initial infections to I_asympt
  # Move from susceptible to infected
  seed_bins <- parse_age_bins(infection_seeding$bins)
  seed_idx <- match(seed_bins$bin_start, parameter_list$age_bin_starts)
  if (any(duplicated(seed_idx))) {
    stop("Seeding is into the same bin multiple times")
  }
  parameter_list$I0_asympt[seed_idx,1,parameter_list$trans_classes] <- infection_seeding$values
  parameter_list$S0[seed_idx] <- parameter_list$S0[seed_idx] - infection_seeding$values
  if ("sircovid_serology" %in% class(sircovid_model)){
   #Put the initial infectives into the sero flow
   parameter_list$R0_pre[,1] <- parameter_list$I0_asympt[,1,parameter_list$trans_classes]
  }
  
  
  #
  # Set gamma and groups
  #
  for (partition in partition_names(sircovid_model)) {
    parameter_list[paste0("s_", partition)] <- sircovid_model$progression_groups[[partition]]
    parameter_list[paste0("gamma_", partition)] <- sircovid_model$gammas[[partition]]
  }
  parameter_list$hosp_transmission <- hosp_transmission
  parameter_list$ICU_transmission <- ICU_transmission
  parameter_list$comm_D_transmission <- comm_D_transmission
  
  if ("sircovid_serology" %in% class(sircovid_model)){
   parameter_list$gamma_test <- sircovid_model$gammas[["test"]] 
  }

  # Remove parameters unused by odin
  parameter_list$age_bin_starts <- NULL
  if ("sircovid_basic" %in% class(sircovid_model)) {
    parameter_list$p_death_ICU <- NULL
    parameter_list$p_ICU_hosp <- NULL
    parameter_list$p_death_hosp_D <- NULL
    parameter_list$p_seroconversion <- NULL
    parameter_list$p_death_comm <- NULL
    parameter_list$comm_D_transmission <- NULL
    parameter_list$p_admit_conf <- NULL
    parameter_list$p_hosp_ILI <- NULL
  } else if ("sircovid_hospital" %in% class(sircovid_model)) {
    parameter_list$p_recov_hosp <- NULL
    parameter_list$p_recov_ICU <- NULL
    parameter_list$p_death_hosp <- NULL
    if (!("sircovid_serology" %in% class(sircovid_model))) {
        parameter_list$p_seroconversion <- NULL
        parameter_list$p_death_comm <- NULL
        parameter_list$comm_D_transmission <- NULL
        parameter_list$p_admit_conf <- NULL
        parameter_list$p_hosp_ILI <- NULL
    } else {
      parameter_list$p_recov_ILI <- NULL
    }
  }         

  parameter_list
}


##' Convert a date into the representation used in the sircovid package,
##' which is days since the start of 2020 (i.e. 2020-01-01 = 1)
##' 
##' Only apply this once, to data or input parameters
##' 
##' @title sircovid date
##' 
##' @param date Date as a string, or a format otherwise understood by lubridate
##' 
##' @return days that \code{date} is after the beginning of 2020
##' 
##' @import lubridate
##'
##' @export 
sircovid_date <- function(date) {
  days_into_2020 <- as.numeric(lubridate::as_date(date) - lubridate::as_date('2019-12-31'))
  if (any(days_into_2020 < 0)) {
    stop("Negative dates, sircovid_date likely applied twice")
  }
  days_into_2020
}

#
# Internal functions
#

##' General parameter creation for use with any model
##' 
##' @title Create parameters base
##' 
##' @param transmission_model Model type. Only 'POLYMOD' currently supported
##' 
##' @param country Country name
##' 
##' @param severity_data_file Location of file with severity data
##' 
##' @param beta Beta, for each time step in \code{beta_times}
##' 
##' @param beta_times Dates \code{beta} changes, in format yyyy-mm-dd
##' 
##' @param trans_profile Proportion in each infectivity group
##' 
##' @param trans_increase Relative infectivity of each group
##' 
##' @param dt Time-step to run the model in days
##' 
##' @param hospital_fitted If TRUE the parameters describing the progression
##' through the hospital compartments of the model are pulled from a set fitted
##' to a hospital model
##'   
##' @param use_polymod_pop Set to ignore \code{survey_pop_in}
##'   and use the population from polymod when estimating the
##'   transmission matrix
##' 
##' @return List of parameters for use with \code{sircovid}
##' 
##' @import socialmixr
##' 
generate_parameters_base <- function(
  transmission_model,
  country,
  severity_data_file,
  beta,
  beta_times,
  trans_profile,
  trans_increase,
  dt,
  use_polymod_pop,
  hospital_fitted) {
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
  
  if (sum(trans_profile) != 1) {
    stop("trans_profile proportions must sum to 1")
  }
  
  if (transmission_model == "POLYMOD") {
    contact_survey = socialmixr::polymod
  } else {
    stop("Only POLYMOD transmission model implemented")
  }
  
  #
  # Set up time-varying beta
  # Times are in days from first day supplied
  beta_t <- normalise_beta(beta_times, dt)

  # 
  # This section defines proportions between partitions
  # derived from the severity.csv file
   
  severity_params <- read_severity(severity_file_in=severity_data_file)
  
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
  # Returns parameters
  #
  parameter_list <- list(N_age = N_age_bins,
                         age_bin_starts = severity_params$age_bin_starts,
                         trans_classes = N_trans_classes,
                         dt = dt,
                         S0 = severity_params$population,
                         trans_increase = trans_increase_array,
                         trans_profile = trans_profile_array,
                         beta_y = beta,
                         beta_t = beta_t,
                         m = transmission_matrix,
                         p_recov_hosp = severity_params$recov_hosp,
                         p_death_hosp = severity_params$death_hosp,
                         p_death_hosp_D = severity_params$death_hosp_D,
                         p_ICU_hosp = severity_params$ICU_hosp,
                         p_hosp_ILI = severity_params$hosp_ILI,
                         p_recov_ILI = severity_params$recov_ILI,
                         p_recov_ICU = severity_params$recov_ICU,
                         p_death_ICU = severity_params$death_ICU,
                         p_asympt = severity_params$asympt,
                         p_sympt_ILI = severity_params$sympt_ILI,
                         p_seroconversion = severity_params$seroconversion,
                         p_death_comm = severity_params$death_comm,
                         p_admit_conf = severity_params$admit_conf)

  parameter_list
}


##' Convert a date into the representation used in the sircovid package,
##' which is days since the start of 2020 (i.e. 2020-01-01 = 1)
##' 
##' Only apply this once, to data or input parameters
##' 
##' @title sircovid date
##' 
##' @param date Date as a string, or a format otherwise understood by lubridate
##' 
##' @return days that \code{date} is after the beginning of 2020
##' 
##' @import lubridate
##' @export
sircovid_date <- function(date) {
  days_into_2020 <- as.numeric(lubridate::as_date(date) - lubridate::as_date('2019-12-31'))
  if (any(days_into_2020 < 0)) {
    stop("Negative dates, sircovid_date likely applied twice")
  }
  days_into_2020
}



##' Convert a date from the representation used in the sircovid package,
##' which is days since the start of 2020 (i.e. 2020-01-01 = 1) into 
##' a lubridate 
##' 
##' @title sircovid_date_as_Date
##' 
##' @param sircovid_date numeric days since the start of 2020 (i.e. 2020-01-01 = 1)
##' 
##' @return \code{date} 
##' 
##' @import lubridate
##' @export

sircovid_date_as_Date <- function(sircovid_date)
{
  if (class(sircovid_date) != "numeric") {
    stop("sircovid_date must be numeric")
  }

  lubridate::as_date('2019-12-31') + sircovid_date
}

#
# Internal functions
#
# Generates vector of times for beta changes, in the same
# terms as the odin code
normalise_beta <- function(beta_times, dt) {
  # Times are in days from first day supplied
  beta_dates <- beta_times - beta_times[[1]]
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
    severity_file <- sircovid_file("extdata/severity_first.csv")
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
    "Adjusted_IFR to 1%",
    "IFR relative to 80",
    "Age specific scaling of ifr to give hospitalisation",
    "Proportion of infections hospitalised compared to age 80",
    "Proportion of infections hospitalised",
    "Proportion of infections needing critical care",
    "Proportion of symptomatic cases hospitalised",
    "Proportion of hospitalised cases getting critical care",
    "Proportion of critical cases dying",
    "Proportion of non-critical care cases dying",
    "Proportion of cases that seroconvert",
    "Proportion of severe cases dying in the community")
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
  
  p_death_hosp_D <- severity_data[["Proportion of non-critical care cases dying"]]
  
  p_ICU_hosp <- severity_data[["Proportion of hospitalised cases getting critical care"]]
  
  p_death_ICU <- severity_data[["Proportion of critical cases dying"]]

  p_recov_ILI <- 1 - severity_data[["Proportion of symptomatic cases hospitalised"]] /
    prop_symp_seek_HC
  
  p_hosp_ILI <- 1 - p_recov_ILI
  
  p_recov_hosp <- 
    (1 - severity_data[["Proportion of hospitalised cases getting critical care"]]) *
    (1 - severity_data[["Proportion of non-critical care cases dying"]])
  
  p_recov_ICU <- 1-p_death_ICU
  
  p_death_hosp <- (1 - severity_data[["Proportion of hospitalised cases getting critical care"]]) *
    severity_data[["Proportion of non-critical care cases dying"]]
  
  p_seroconversion <- severity_data[["Proportion of cases that seroconvert"]]
  
  p_death_comm <- severity_data[["Proportion of severe cases dying in the community"]]
  
  p_admit_conf <- severity_data[["Proportion of hospitalised cases admitted as confirmed"]]
  
  list(
    population = population,
    age_bin_starts = age_bins$bin_start,
    age_bin_ends = age_bins$bin_end,
    asympt = p_asympt,
    sympt_ILI = p_sympt_ILI,
    recov_ICU = p_recov_ICU,
    recov_ILI = p_recov_ILI,
    hosp_ILI = p_hosp_ILI,
    recov_hosp = p_recov_hosp,
    death_hosp = p_death_hosp,
    death_hosp_D = p_death_hosp_D,
    death_ICU = p_death_ICU,
    ICU_hosp = p_ICU_hosp,
    seroconversion = p_seroconversion,
    death_comm = p_death_comm,
    admit_conf = p_admit_conf)
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
