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
##'   index into the model state that should be returned; will be set
##'   into the model with `$set_index`
##'
##' @param n_threads The number of threads to use in the
##'   simulation. This only has an effect if your sircovid was built
##'   with OpenMP support (possibly not the case on macOS).  Each
##'   simulation will be run on potentially a different thread.
##'
##' @param seed Optional seed to create the model with
##'
##' @export
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

  ## NOTE: this is basically the same as mcstate::pmcmc_predict, which
  ## suggests that's where this belongs (see mcstate issue #83)
  obj <- mod$new(p, steps[[1]], 1L, n_threads = n_threads,
                 seed = seed, pars_multi = TRUE)
  ## TODO: update to mcstate::array_reshape(state, 2, c(1, length(p)))
  dim(state) <- c(nrow(state), 1, length(p))
  obj$update_state(state = state)
  if (!is.null(index)) {
    obj$set_index(index)
  }

  for (i in seq_len(n_epoch)) {
    p_i <- p
    p_new <- events$data[[i]]
    for (j in seq_along(p_i)) {
      p_i[[j]][names(p_new)] <- p_new
    }
    obj$update_state(pars = p_i, set_initial_state = FALSE)
    i_step <- steps >= step_from[[i]] & steps <= step_to[[i]]
    res[, , i_step] <- obj$simulate(steps[i_step])
  }

  rownames(res) <- names(index)

  ## TODO: use mcstate::array_drop(state, 2)
  state <- obj$state()
  dim(state) <- dim(state)[c(1, 3)]

  ## TODO: This mimicked the interface of dust::dust_simulate, but
  ## that function is deprecated. If we retain this function in the
  ## package we should check what is really wanted here.
  attr(res, "step") <- steps
  attr(res, "date") <- sircovid_date_as_date(steps * dt)
  attr(res, "rng_state") <- obj$rng_state()
  attr(res, "state") <- state

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
