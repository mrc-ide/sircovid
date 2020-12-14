sircovid_simulate <- function(mod, state, p, events,
                              index = NULL, n_threads = 1L, seed = NULL) {
  assert_is(events, "sircovid_events")

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
  class(ret) <- "sircovid_events"
  ret
}
