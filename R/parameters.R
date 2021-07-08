## These could be moved to be defaults within the models
sircovid_parameters_shared <- function(start_date, region,
                                       beta_date, beta_value,
                                       beta_type = 'piecewise-linear',
                                       population = NULL) {
  steps_per_day <- 4
  dt <- 1 / steps_per_day
  assert_sircovid_date(start_date)
  if (beta_type == "piecewise-linear") {
    beta_step <- sircovid_parameters_piecewise_linear(beta_date,
                                                      beta_value %||% 0.08, dt)
  } else if (beta_type == "piecewise-constant") {
    beta_step <- sircovid_parameters_piecewise_constant(beta_date,
                                                        beta_value %||% 0.08,
                                                        dt)
  } else {
    stop("'beta_type' must be 'piecewise-linear' or 'piecewise-constant'")
  }


  population <- population %||% sircovid_population(region)

  list(hosp_transmission = 0,
       ICU_transmission = 0,
       G_D_transmission = 0,
       dt = dt,
       steps_per_day = steps_per_day,
       initial_step = start_date / dt,
       n_age_groups = length(sircovid_age_bins()$start),
       beta_step = beta_step,
       population = population)
}


##' Construct a piecewise linear quantity over time array for use within
##' sircovid models.
##'
##' @title Construct piecewise linear array
##'
##' @param date Either `NULL`, if one value of the quantity will be used for
##'   all time steps, or a vector of times that will be used as change
##'   points. Must be provided as a [sircovid_date()], i.e., days into
##'   2020. The first date must be 0.
##'
##' @param value A vector of values to use for the quantity - either a scalar
##'   (if `date` is `NULL`) or a vector the same length as `date`.
##'
##' @param dt The timestep that will be used in the simulation. This
##'   must be of the form `1 / n` where `n` is an integer representing
##'   the number of steps per day. Ordinarily this is set by sircovid
##'   internally to be `0.25` but this will become tuneable in a
##'   future version.
##'
##' @return Returns a vector of piecewise linear values, one per timestep,
##'   until the values stabilise.  After this point the quantity is assumed to
##'   be constant.
##'
##' @seealso [sircovid_parameters_expand_step()] - see examples below
##'
##' @export
##' @examples
##' # If "date" is NULL, then the quantity is constant and this function is
##' # trivial:
##' sircovid::sircovid_parameters_piecewise_linear(NULL, 0.1, 0.25)
##'
##' date <- sircovid::sircovid_date(
##'    c("2020-02-01", "2020-02-14", "2020-03-15"))
##' value <- c(3, 1, 2)
##' y <- sircovid::sircovid_parameters_piecewise_linear(date, value, 1)
##'
##' # The implied time series looks like this:
##' t <- seq(0, date[[3]])
##' plot(t, y, type = "o")
##' points(date, value, pch = 19, col = "red")
##'
##' # After 2020-03-15, the quantity value will be fixed at 2, the value
##' # that it reached at that date.
##'
##' # You can see this using sircovid_parameters_expand_step
##' # If a vector of dates is provided then, it's more complex. We'll
##' # use dt of 1 here as it's easier to visualise
##' t <- seq(0, 100, by = 1)
##' sircovid::sircovid_parameters_expand_step(t, y)
##' plot(t, sircovid::sircovid_parameters_expand_step(t, y), type = "o")
##' points(date, value, pch = 19, col = "red")
##'
##' # If dt is less than 1, this is scaled, but the pattern of
##' # change is the same
##' y <- sircovid::sircovid_parameters_piecewise_linear(date, value, 0.5)
##' t <- seq(0, date[[3]], by = 0.5)
##' plot(t, y, type = "o", cex = 0.25)
##' points(date, value, pch = 19, col = "red")
sircovid_parameters_piecewise_linear <- function(date, value, dt) {
  ## for one strain coerce numeric to matrix, will coerce back later
  if (!inherits(value, "matrix")) {
    value <- matrix(value, ncol = 1)
  }
  if (is.null(date)) {
    if (nrow(value) != 1L) {
      stop("As 'date' is NULL, expected single value")
    }
    ## coerce back
    if (ncol(value) == 1) {
      value <- as.numeric(value)
    }
    return(value)
  }
  if (length(date) != nrow(value)) {
    stop("'date' and 'value' must have the same length")
  }
  if (length(date) < 2) {
    stop("Need at least two dates and values for a varying piecewise linear")
  }
  assert_sircovid_date(date)
  assert_increasing(date)

  if (date[[1]] != 0) {
    date <- c(0, date)
    value <- rbind(value[1, ], value)
  }

  value <- apply(value, 2, function(x) {
    stats::approx(date, x, seq(0, date[[length(date)]], by = dt))$y
  })

  ## coerce back
  if (ncol(value) == 1) {
    value <- as.numeric(value)
  }

  value
}


##' Construct a piecewise constant quantity over time array for use within
##' sircovid models.
##'
##' @title Construct piecewise constant array
##'
##' @param date Either `NULL`, if one value of the quantity will be used for
##'   all time steps, or a vector of times that will be used as change
##'   points. Must be provided as a [sircovid_date()], i.e., days into
##'   2020. The first date must be 0.
##'
##' @param value A vector of values to use for the quantity - either a scalar
##'   (if `date` is `NULL`) or a vector the same length as `date`.
##'
##' @param dt The timestep that will be used in the simulation. This
##'   must be of the form `1 / n` where `n` is an integer representing
##'   the number of steps per day. Ordinarily this is set by sircovid
##'   internally to be `0.25` but this will become tuneable in a
##'   future version.
##'
##' @return Returns a vector of piecewise constant values, one per timestep,
##'   until the values stabilise.  After this point the quantity is assumed to
##'   be constant.
##'
##' @export
##' @examples
##' # If "date" is NULL, then the quantity is constant and this function is
##' # trivial:
##' sircovid::sircovid_parameters_piecewise_constant(NULL, 0.1, 0.25)
##'
##' date <- sircovid::sircovid_date(
##'    c("2019-12-31", "2020-02-01", "2020-02-14", "2020-03-15"))
##' value <- c(0, 3, 1, 2)
##' y <- sircovid::sircovid_parameters_piecewise_constant(date, value, 1)
##'
##' # The implied time series looks like this:
##' t <- seq(0, date[[4]])
##' plot(t, y, type = "o")
##' points(date, value, pch = 19, col = "red")
##'
##' # After 2020-03-15, the quantity value will be fixed at 2, the value
##' # that it reached at that date.
##'
##' # You can see this using sircovid_parameters_expand_step
##' # If a vector of dates is provided then, it's more complex. We'll
##' # use dt of 1 here as it's easier to visualise
##' t <- seq(0, 100, by = 1)
##' sircovid::sircovid_parameters_expand_step(t, y)
##' plot(t, sircovid::sircovid_parameters_expand_step(t, y), type = "o")
##' points(date, value, pch = 19, col = "red")
##'
##' # If dt is less than 1, this is scaled, but the pattern of
##' # change is the same
##' y <- sircovid::sircovid_parameters_piecewise_constant(date, value, 0.5)
##' t <- seq(0, date[[4]], by = 0.5)
##' plot(t, y, type = "o", cex = 0.25)
##' points(date, value, pch = 19, col = "red")
sircovid_parameters_piecewise_constant <- function(date, value, dt) {
  if (is.null(date)) {
    if (length(value) != 1L) {
      stop("As 'date' is NULL, expected single value")
    }
    return(value)
  }
  if (length(date) != length(value)) {
    stop("'date' and 'value' must have the same length")
  }
  assert_sircovid_date(date)
  assert_increasing(date)
  if (!is.null(date)) {
    if (date[1L] != 0) {
      stop("As 'date' is not NULL, first date should be 0")
    }
  }

  stats::approx(date, value, method = "constant",
                xout = seq(0, date[[length(date)]], by = dt))$y

}


##' Expand `value_step` based on a series of `step`s.  Use this to
##' convert between the values passed to
##' [sircovid_parameters_piecewise_linear()] and the actual values
##' for a given set of steps.
##'
##' @title Expand beta steps
##'
##' @param step A vector of steps
##'
##' @param value_step A vector of values
##'
##' @return A numeric vector the same length as `step`
##'
##' @export
sircovid_parameters_expand_step <- function(step, value_step) {
  value_step[pmin(step, length(value_step) - 1L) + 1L]
}


##' Process severity data
##' @title Process severity data
##'
##' @param params Severity data, via Bob Verity's `markovid`
##'   package. This needs to be `NULL` (use the default bundled data
##'   version in the package), a [data.frame] object (for raw severity
##'   data) or a list (for data that has already been processed by
##'   `sircovid` for use).  New severity data comes from Bob Verity
##'   via the markovid package, and needs to be carefully calibrated
##'   with the progression parameters.
##'
##' @return A list of age-structured probabilities
##' @export
sircovid_parameters_severity <- function(params) {
  if (is.null(params)) {
    params <- severity_default()
  } else if (!is.data.frame(params)) {
    expected <- c("p_star", "p_C", "p_G_D",
                  "p_H_D", "p_ICU_D", "p_W_D",
                  "p_ICU", "p_H",
                  "p_sero_pos_1", "p_sero_pos_2")
    verify_names(params, expected)
    return(params)
  }

  ## Transpose so columns are parameters, rownames are age groups
  data <- t(as.matrix(params[-1L]))
  colnames(data) <- params[[1]]
  data <- cbind(age = rownames(data),
                data.frame(data, check.names = FALSE),
                stringsAsFactors = FALSE)
  rownames(data) <- NULL

  required <- c(
    p_C = "p_C",
    p_H = "p_H",
    p_ICU = "p_ICU",
    p_ICU_D = "p_ICU_D",
    p_H_D = "p_H_D",
    p_W_D = "p_W_D",
    p_sero_pos_1 = "p_sero_pos_1",
    p_sero_pos_2 = "p_sero_pos_2",
    p_G_D = "p_G_D",
    p_star = "p_star")
  data <- rename(data, required, names(required))

  list(
    p_star = data[["p_star"]],
    p_C = data[["p_C"]],
    p_G_D = data[["p_G_D"]],
    p_H_D = data[["p_H_D"]],
    p_ICU_D = data[["p_ICU_D"]],
    p_W_D = data[["p_W_D"]],
    p_ICU = data[["p_ICU"]],
    p_sero_pos_1 = data[["p_sero_pos_1"]],
    p_sero_pos_2 = data[["p_sero_pos_2"]],
    p_H = data[["p_H"]])
}
