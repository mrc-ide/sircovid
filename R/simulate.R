##' Simulate from a sircovid model ([basic] or [carehomes]) given a
##' set of starting positions and parameters. This forms the core of
##' our projections support, and is typically used after fitting a
##' model.
##'
##' @title Simulate a sircovid model
##'
##' @param mod The model object to simulate from
##'
##' @param state A matrix of state; rows correspond to state
##'   variables, columns correspond to different parameter sets or
##'   realisations.
##'
##' @param p An unnamed list of parameter objects; must be the same
##'   length as `state` has columns.
##'
##' @param events A [sircovid_simulate_events()] object with "events"
##'   data, representing changes in the parameters, and the overall
##'   timescale of the simulation.
##'
##' @param index A (possibly named) integer vector representing an
##'   index into the model state that should be returned (see
##'   [dust::dust_simulate])
##'
##' @param n_threads The number of threads to use in the
##'   simulation. This only has an effect if your sircovid was built
##'   with OpenMP support (possibly not the case on macOS).  Each
##'   simulation will be run on potentially a different thread.
##'
##' @param seed Optional seed (see [dust::dust_simulate])
##'
##' @param export
##'
##' @return A 3-d array with dimensions representing state (along
##'   `index`), parameter set (along `p`) and time. The array has
##'   attributes `step` and `date` corresponding to the step (relative
##'   to 2020-01-01) and date.
sircovid_simulate <- function(mod, state, p, events,
                              index = NULL, n_threads = 1L, seed = NULL) {
  assert_is(events, "sircovid_simulate_events")

  if (!is.list(p) || !is.null(names(p))) {
    stop("Expected 'p' to be an unnamed list")
  }

  dt <- vnapply(p, "[[", "dt")
  if (length(unique(dt)) > 1L) {
    stop("All entries in 'p' must have the same value of 'dt'")
  }
  dt <- dt[[1L]]

  for (e in events$data) {
    if ("dt" %in% names(e) && e$dt != dt) {
      stop("Events must not change 'dt'")
    }
  }

  n_epoch <- length(events$data)
  step_from <- events$date_from / dt
  step_to <- events$date_to / dt
  steps <- seq(step_from[[1]], step_to[[n_epoch]], by = 1 / dt)

  n_state <- if (is.null(index)) nrow(state) else length(index)
  res <- array(NA_real_, c(n_state, length(p), length(steps)))

  for (i in seq_len(n_epoch)) {
    p_i <- p
    p_new <- events$data[[i]]
    for (j in seq_along(p_i)) {
      p_i[[j]][names(p_new)] <- p_new
    }
    i_step <- steps >= step_from[[i]] & steps <= step_to[[i]]
    y <- dust::dust_simulate(mod, steps[i_step], p_i, state,
                             index = index, n_threads = n_threads,
                             seed = seed, return_state = TRUE)
    res[, , i_step] <- y
    state <- attr(y, "state")
    seed <- attr(y, "rng_state")
  }

  rownames(res) <- names(index)

  ## We're basically mimicing the interface here of
  ## dust::dust_simulate which probably should return a list. So to
  ## make downstream use a bit easier we'll use attributes to return
  ## time domain information.
  attr(res, "step") <- steps
  attr(res, "date") <- sircovid_date_as_date(steps * dt)

  res
}


##' Events for use with [sircovid_simulate()]. This captures the
##' overall timescape of the simulation and a set of "events" which
##' correspond to changes in parameters.
##'
##' @title Simulation events
##'
##' @param date_start The start date for the simulation, as an ISO
##'   date or [sircovid_date]
##'
##' @param date_end The end date for the simulation, as an ISO
##'   date or [sircovid_date]
##'
##' @param data Events data; this must be a named list, with each name
##'   being in ISO 8601 (YYYY-MM-DD) format, and value being a named
##'   set of parameters to apply to our base set. So a value of
##'   `list("2021-03-01" = list(vaccine_daily_doses = 10000))` would
##'   represent a change to the `vaccine_daily_doses` parameter from
##'   the 1st of March 2021.  Each entry can change as many parameters
##'   as wanted. These changes are *not cumulative*, so if you want a
##'   change to persist across multiple events it must exist in each
##'   of these events.
##'
##' @return A `sircovid_simulate_events` object; consider its contents
##'   opaque for now as these may change.
##'
##' @export
sircovid_simulate_events <- function(date_start, date_end, data) {
  date_start <- sircovid_date(date_start)
  date_end <- sircovid_date(date_end)
  date_vary <- sircovid_date(names(data))

  if (date_start >= date_end) {
    stop("'date_start' must be less than 'date_end'")
  }

  if (!all(diff(date_vary) > 0)) {
    stop("The dates used as 'data' names must be strictly increasing")
  }

  ## Events we will not hit are relatively easy to fix; filter them off:
  if (any(date_vary >= date_end)) {
    i <- date_vary < date_end
    data <- data[i]
    date_vary <- date_vary[i]
  }

  if (any(date_vary <= date_start)) {
    i <- seq(max(which(date_vary < date_start)), length(data))
    data <- unname(data[i])
    date_from <- c(date_start, date_vary[i][-1])
    date_to <- c(date_vary[i][-1], date_end)
  } else {
    data <- c(list(NULL), unname(data))
    date_from <- c(date_start, date_vary)
    date_to <- c(date_vary, date_end)
  }

  ret <- list(date_from = date_from, date_to = date_to, data = data)
  class(ret) <- "sircovid_simulate_events"
  ret
}
