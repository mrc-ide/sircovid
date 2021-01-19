##' Create parameters that include beta values correspond to different
##' Rt scenarios in the future.
##'
##' This function is called `add_future_betas` because updates beta
##' values, however it does this by applying Rt changes in the future
##' via `future_Rt`.
##'
##' There are two basic ways that new beta values in the future might
##' be added. First, as relative to some previous day in the past. For
##' example to use a value of beta that corresponds to 1.5x the Rt
##' recorded on the 28th of October, we might write
##'
##' ```
##' sircovid::future_Rt(1.5, "2020-10-28")
##' ```
##'
##' Second, as an absolute Rt value. For example to use a value of
##' beta corresponding to an Rt of 0.9 we might write
##'
##' ```
##' sircovid::future_Rt(0.9)
##' ```
##'
##' These new levels need to be associated with dates at which they
##' occur, and will result in a piece-wise constant pattern with 1
##' day's transition between each level.
##'
##' @title Create future betas
##'
##' @param sample A `mcstate_pmcmc` object
##'
##' @param rt The results of [sircovid::carehomes_Rt_trajectories]
##'   with matrix elements `Rt_general` and `date` which will be used to
##'   calculate relative Rt values.
##'
##' @param future A named list of [sircovid::future_Rt] values. Each
##'   element name corresponds to the date that the change takes
##'   effect (in ISO-8601 "YYYY-MM-DD" format) and must be strictly
##'   increasing with at least two days separating changes.
##'
##'@param rt_type A string giving the entry of `rt` to modify,
##' defaults to `Rt_general`
##' 
##' @export
##' @author Richard Fitzjohn
add_future_betas <- function(sample, rt, future, rt_type = "Rt_general") {
  assert_is(sample, "mcstate_pmcmc")
  assert_is(rt, "Rt_trajectories")
  sample <- drop_trajectory_predicted(sample)

  ## We never want the future projections here, so exclude them
  n_pars <- nrow(sample$pars)
  n_time <- length(sample$trajectories$step)

  ## Check that our trajectories are consistent with the sample:
  i <- match(sample$trajectories$step, rt$step[, 1])
  rt_date <- rt$date[i, 1]

  rt_value <- rt[[rt_type]][i, , drop = FALSE]

  beta_scaling <- future_relative_beta(future, rt_date, rt_value, "fbs_")
  sample$pars <- cbind(sample$pars, beta_scaling$value)

  ## Capture the transformation of base parameters here as we will
  ## overwrite this, but the new transformation function will pick
  ## this up due to lexical scope.
  base_transform <- sample$predict$transform

  ## The new transform will do the old transform but also extend the
  ## beta trace with new values that correspond to our new scenario
  sample$predict$transform <- function(p) {
    ret <- base_transform(p)
    ret$beta_step <- add_future_beta(ret,
                                     beta_scaling$date,
                                     p[colnames(beta_scaling$value)])
    ret
  }

  sample
}


##' @rdname add_future_betas
##' @export
##'
##' @param value A value to add in the future. If `relative_to` is
##'   `NULL`, then this is an *absolute* Rt value, otherwise it is a
##'   relative value.
##'
##' @param relative_to Optionally an ISO 8601 (YYYY-MM-DD) format date
##'   string or [Date] object indicating the date that the value
##'   should be taken relative to.
future_Rt <- function(value, relative_to = NULL) {
  assert_scalar(value)
  if (is.null(relative_to)) {
    ret <- list(value = value,
                relative_value = NA_real_, relative_to = NA_character_)
  } else {
    ret <- list(value = NA_real_,
                relative_value = value,
                relative_to = assert_date_string(relative_to))
  }
  class(ret) <- "future_Rt"
  ret
}


## This can be split into validate (turning 'future' into a data.frame
## or similar) and 'calculate'. There's still a *lot* of logic in
## here, but the final product is fairly easily understood.
future_relative_beta <- function(future, rt_date, rt_value, prefix = NULL) {
  if (length(future) == 0) {
    stop("Expected at least one element in 'future'")
  }
  if (!all(vlapply(future, inherits, "future_Rt"))) {
    stop("Expected all elements of 'future' to be 'future_Rt' objects")
  }
  if (is.null(names(future))) {
    stop("Expected 'future' to be named")
  }
  future_date <- as_date(names(future))

  if (!all(diff(future_date) > 0)) {
    stop("Future change dates must be increasing")
  }
  if (!all(diff(future_date) > 1)) {
    stop("At least one full date required between all future change dates")
  }

  date <- last(rt_date)
  if (sircovid_date(future_date[[1]]) - date <= 1) {
    stop(sprintf("The first future date must be at least %s (but was %s)",
                 sircovid_date_as_date(date) + 2, future_date[[1]]))
  }

  current_rt <- rt_value[nrow(rt_value), ]
  future_value <- vnapply(future, "[[", "value", USE.NAMES = FALSE)
  is_relative <- is.na(future_value)

  relative_to <- vcapply(future[is_relative], "[[", "relative_to",
                         USE.NAMES = FALSE)
  relative_value <- vnapply(future[is_relative], "[[", "relative_value",
                            USE.NAMES = FALSE)
  relative_to_index <- match(sircovid_date(relative_to), rt_date)
  if (any(is.na(relative_to_index))) {
    stop("Relative date not found in rt set: ",
         paste(relative_to[is.na(relative_to_index)], collapse = ", "))
  }

  value <- matrix(NA_real_, length(future), ncol(rt_value))
  value[is_relative, ] <- relative_value * rt_value[relative_to_index, ] /
    rep(current_rt, each = length(relative_to_index))
  value[!is_relative, ] <- outer(future_value[!is_relative], current_rt, "/")
  value <- t(value)

  date <- sircovid_date(future_date)

  ## We need to do a little transformation here so that we end up with
  ## a piecewise-constant-with-one-day-ease pattern
  date_use <- c(rbind(date - 1, date))
  n <- length(future)
  value_use <- cbind(1, value[, c(rep(seq_len(n - 1), each = 2), n)])

  if (!is.null(prefix)) {
    colnames(value_use) <- paste0(prefix, seq_len(ncol(value_use)))
  }

  list(date = date_use, value = value_use)
}


add_future_beta <- function(p, beta_date, beta_scaling) {
  dt <- p$dt
  beta_step <- p$beta_step
  beta_last <- last(beta_step)

  beta_value <- beta_last * beta_scaling
  beta_extra <- sircovid_parameters_beta(beta_date, beta_value, dt)

  beta_index <- seq(beta_date[[1]] / dt + 1L, length(beta_extra))
  beta_step[beta_index] <- beta_extra[beta_index]

  ## Values between the final change and first forecast change will
  ## end up as an NA here, so replace with a constant
  beta_step[is.na(beta_step)] <- beta_last
  beta_step
}
