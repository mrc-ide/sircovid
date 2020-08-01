## These could be moved to be defaults within the models
sircovid_parameters_shared <- function(start_date, region,
                                       beta_date, beta_value) {
  dt <- 0.25
  assert_sircovid_date(start_date)
  beta_step <- sircovid_parameters_beta(beta_date, beta_value %||% 0.08, dt)
  list(hosp_transmission = 0.1,
       ICU_transmission = 0.05,
       comm_D_transmission = 0.05,
       dt = dt,
       initial_step = start_date / dt,
       N_age = length(sircovid_age_bins()$start),
       beta_step = beta_step,
       population = sircovid_population(region))
}


##' Construct the `beta` (contact rate) over time array for use within
##' sircovid models.
##'
##' @title Construct beta array
##'
##' @param date Either `NULL`, if one value of beta will be used for
##'   all time steps, or a vector of times that will be used as change
##'   points. Must be provided as a [sircovid_date()], i.e., days into
##'   2020.
##'
##' @param value A vector of values to use for beta - either a scalar
##'   (if `date` is `NULL`) or a vector the same length as `date`.
##'
##' @param dt The timestep that will be used in the simulation. This
##'   must be of the form `1 / n` where `n` is an integer representing
##'   the number of steps per day. Ordinarily this is set by sircovid
##'   internally to be `0.25` but this will become tuneable in a
##'   future version.
##'
##' @return Returns a vector of beta values, one per timestep, until
##'   the beta values stabilise.  After this point beta is assumed to
##'   be constant.
##'
##' @export
##' # If "date" is NULL, then beta is constant and this function is
##' # trivial:
##' sircovid_parameters_beta(NULL, 0.1, 0.25)
##'
##' # If a vector of dates is provided then, it's more complex. We'll
##' # use dt of 1 here as it's easier to visualise
##' t <- seq(0, 100, by = 1)
##'
##' date <- sircovid_date(c("2020-02-01", "2020-02-14", "2020-03-15"))
##' value <- c(3, 1, 2)
##' beta <- sircovid_parameters_beta(date, value, 1)
##'
##' # The implied time series looks like this:
##' t <- seq(0, date[[3]])
##' plot(t, beta, type = "o")
##' points(date, value, pch = 19, col = "red")
##'
##' # After 2020-03-15, the beta value will be fixed at 2, the value
##' # that it reached at that date.
##'
##' # If dt is less than 1, this is scaled, but the pattern of beta
##' # change is the same
##' beta <- sircovid_parameters_beta(date, value, 0.5)
##' t <- seq(0, date[[3]], by = 0.5)
##' plot(t, beta, type = "o", cex = 0.25)
##' points(date, value, pch = 19, col = "red")
sircovid_parameters_beta <- function(date, value, dt) {
  if (is.null(date)) {
    if (length(value) != 1L) {
      stop("As 'date' is NULL, expected single value")
    }
    return(value)
  }
  if (length(date) != length(value)) {
    stop("'date' and 'value' must have the same length")
  }
  if (length(date) < 2) {
    stop("Need at least two dates and betas for a varying beta")
  }
  assert_sircovid_date(date)
  assert_increasing(date)
  stats::approx(c(0, date), c(value[[1]], value),
                seq(0, date[[length(date)]], by = dt))$y
}


sircovid_parameters_severity <- function(params) {
  ## Set up severity file into table
  if (is.null(params)) {
    params <- severity_default()
  }

  ## Transpose so columns are parameters, rownames are age groups
  data <- t(as.matrix(params[-1L]))
  colnames(data) <- params[[1]]
  data <- cbind(age = rownames(data),
                data.frame(data, check.names = FALSE),
                stringsAsFactors = FALSE)
  rownames(data) <- NULL

  required <- c(
    population = "Size of England population",
    p_sympt_seek_hc = "Proportion of symptomatic cases seeking healthcare",
    p_sympt = "Proportion with symptoms",
    p_sympt_hosp = "Proportion of symptomatic cases hospitalised",
    p_ICU_hosp = "Proportion of hospitalised cases getting critical care",
    p_death_ICU = "Proportion of critical cases dying",
    p_death_hosp_D = "Proportion of non-critical care cases dying",
    p_seroconversion = "Proportion of cases that seroconvert",
    p_death_comm = "Proportion of severe cases dying in the community",
    p_admit_conf = "Proportion of hospitalised cases admitted as confirmed")
  data <- rename(data, required, names(required))

  # Parse the age bins. Useful to keep both start and end depending on what
  # function expects
  age_bins <- check_age_bins(data[["age"]])

  p_recov_ILI <- 1 - data[["p_sympt_hosp"]] / data[["p_sympt_seek_hc"]]

  list(
    p_asympt = 1 - data[["p_sympt"]],
    p_sympt_ILI = data[["p_sympt"]] * data[["p_sympt_seek_hc"]],
    p_recov_ICU = 1 - data[["p_death_ICU"]],
    p_recov_ILI = p_recov_ILI,
    p_hosp_ILI = 1 - p_recov_ILI,
    p_recov_hosp = (1 - data[["p_ICU_hosp"]]) * (1 - data[["p_death_hosp_D"]]),
    p_death_hosp = (1 - data[["p_ICU_hosp"]]) * data[["p_death_hosp_D"]],
    p_death_hosp_D = data[["p_death_hosp_D"]],
    p_death_ICU = data[["p_death_ICU"]],
    p_ICU_hosp = data[["p_ICU_hosp"]],
    p_seroconversion = data[["p_seroconversion"]],
    p_death_comm = data[["p_death_comm"]],
    p_admit_conf = data[["p_admit_conf"]])
}
