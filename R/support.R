## General things that we need but that aren't that interesting.

##' Generate lists of regions that we use for various tasks.
##'
##' @title Regions
##'
##' @param type The name of a region type; must be one of "all",
##'   "england", or "nations"
##'
##' @return A character vector
##' @export
##' @examples
##' sircovid::regions("england")
regions <- function(type) {
  regions_england <- c("east_of_england",
                       "london",
                       "midlands",
                       "north_east_and_yorkshire",
                       "north_west",
                       "south_east",
                       "south_west")
  celtic_nations <- c("scotland", "wales", "northern_ireland")

  switch(type,
         "all" = c(regions_england, celtic_nations),
         "england" = regions_england,
         "nations" = c("england", celtic_nations),
         stop(sprintf("Unknown region type '%s'", type)))
}


## We always use these age bands, so rather than detect them, we will
## check that things conform to them.
sircovid_age_bins <- function() {
  end <- c(seq(4, 79, by = 5), 100)
  start <- c(0, end[-length(end)] + 1L)
  list(start = start, end = end)
}


check_age_bins <- function(age_headers) {
  bins <- sircovid_age_bins()
  expected <- sprintf("%d to %d", bins$start, bins$end)
  if (!identical(age_headers, expected)) {
    stop(sprintf("Incorrect age bands:\nexpected: %s\ngiven: %s",
                 paste(squote(expected), collapse = ", "),
                 paste(squote(age_headers), collapse = ", ")))
  }
  bins
}


sircovid_population <- function(region) {
  if (is.null(region)) {
    stop("'region' must not be NULL")
  }

  if (is.null(cache$population)) {
    cache$population <- read_csv(sircovid_file("extdata/population.csv"))
  }

  population <- cache$population[[tolower(region)]]
  if (is.null(population)) {
    valid <- setdiff(names(cache$population), "age")
    stop(sprintf("Population not found for '%s': must be one of %s",
                 region, paste(squote(valid), collapse = ", ")))
  }

  population
}


sircovid_carehome_beds <- function(region) {
  if (is.null(region)) {
    stop("'region' must not be NULL")
  }

  if (is.null(cache$carehomes)) {
    cache$carehomes <- read_csv(sircovid_file("extdata/carehomes.csv"))
  }

  i <- match(tolower(region), cache$carehomes$region)
  if (is.na(i)) {
    valid <- paste(squote(cache$carehomes$region), collapse = ", ")
    stop(sprintf("Carehome beds not found for '%s': must be one of %s",
                 region, valid))
  }

  cache$carehomes$carehome_beds[[i]]
}


##' @importFrom stats dnbinom rexp
ll_nbinom <- function(data, model, kappa, exp_noise) {
  if (is.na(data)) {
    return(numeric(length(model)))
  }
  mu <- model + rexp(length(model), rate = exp_noise)
  dnbinom(data, kappa, mu = mu, log = TRUE)
}


##' @importFrom stats dbinom
ll_binom <- function(data_x, data_size, model_prob) {
  if (is.na(data_x) || is.na(data_size)) {
    return(numeric(length(model_prob)))
  }
  dbinom(data_x, data_size, model_prob, log = TRUE)
}


ll_betabinom <- function(data_x, data_size, model_prob, rho) {
  if (is.na(data_x) || is.na(data_size)) {
    return(numeric(length(model_prob)))
  }
  dbetabinom(data_x, data_size, model_prob, rho, log = TRUE)
}


dbetabinom <- function(x, size, prob, rho, log = FALSE) {

  a <- prob * (1 / rho - 1)
  b <- (1 - prob) * (1 / rho - 1)

  out <- lchoose(size, x) + lbeta(x + a, size - x + b) - lbeta(a, b)

  if (!log) {
    out <- exp(out)
  }

  out
}


##' @importFrom stats rexp
test_prob_pos <- function(pos, neg, sensitivity, specificity, exp_noise) {

  ## We add some exponential noise to the number of positives and negatives
  ## to help ensure prob_pos is not 0 or 1. If e.g. prob_pos were 0 and there
  ## were individuals who tested positive, this would result in a weight of 0
  ## for a particle. If all particles have weights of 0, the particle filter
  ## breaks. The exponential noise produces small non-zero weights in these
  ## circumstances to prevent the particle filter from breaking.

  pos <- pos + rexp(length(pos), exp_noise)
  neg <- neg + rexp(length(neg), exp_noise)

  prob_pos <- (sensitivity * pos + (1 - specificity) * neg) / (pos + neg)
  prob_pos
}


sircovid_transmission_matrix <- function(region, population = NULL) {
  if (is.null(cache$transmission_matrix[[region]])) {
    if (is.null(cache$transmission_matrix)) {
      cache$transmission_matrix <- list()
    }
    age_bins <- sircovid_age_bins()

    max_polymod_age <- 70

    population <- population %||% sircovid_population(region)

    ## Survey population in socialmixr format
    survey_pop <- data_frame(
      "lower.age.limit" = sircovid_age_bins()$start,
      population = population)

    ## Get the contact matrix from socialmixr; polymod only
    ## goes up to age 70
    contact <- suppressMessages(socialmixr::contact_matrix(
      survey.pop = survey_pop,
      socialmixr::polymod,
      countries = "United Kingdom",
      age.limits = age_bins$start[age_bins$start <= max_polymod_age],
      symmetric = TRUE))

    ## Transform the matrix to the (symetrical) transmission matrix
    ## rather than the contact matrix
    transmission <- contact$matrix /
      rep(contact$demography$population, each = ncol(contact$matrix))

    ## POLYMOD has a max age of 70, so older bins need to be filled in
    ## assumes that the probability of contact remains as in POLYMOD
    ## and that contacts are the same in 70+ and 80+
    extra <- age_bins$start[age_bins$start > max_polymod_age]
    i <- c(seq_len(nrow(transmission)), rep(nrow(transmission), length(extra)))
    transmission <- transmission[i, i]

    ## Rename for cosmetic reasons only
    nms <- sprintf("[%d,%d)", age_bins$start, age_bins$end)
    dimnames(transmission) <- list(nms, nms)

    cache$transmission_matrix[[region]] <- transmission
  }

  cache$transmission_matrix[[region]]
}


severity_default <- function() {
  if (is.null(cache$severity_default)) {
    cache$severity_default <-
      read_csv(sircovid_file("extdata/severity_default.csv"))
  }
  cache$severity_default
}


##' Augment (or remove) trajectories from a `mcstate_trajectories` or
##' `mcstate_pmcmc` object.
##'
##' @title Add incidence to trajectories
##'
##' @param obj A `mcstate_trajectories` or `mcstate_pmcmc` object to
##'   update
##'
##' @param states A vector of cumulative states to compute incidence
##'   for
##'
##' @param suffix A string to append to the input states to create the
##'   incidence variables
##'
##' @export
add_trajectory_incidence <- function(obj, states, suffix = "_inc") {
  if (inherits(obj, "mcstate_pmcmc")) {
    obj$trajectories <- add_trajectory_incidence(obj$trajectories,
                                                 states, suffix)
    return(obj)
  }
  assert_is(obj, "mcstate_trajectories")
  if (length(states) == 0) {
    return(obj)
  }

  if (length(dim(obj$state)) == 3) {
    add_trajectory_incidence_single(obj, states, suffix)
  } else {
    add_trajectory_incidence_nested(obj, states, suffix)
  }
}

## In order to compute incidence we have to add two NA values; one
## is the usual one dropped in a rolling difference, the other is
## dropped as the first time interval is potentially much longer
## than a day.
trajectory_incidence <- function(x) {
  c(NA, NA, diff(x[-1L]))
}

add_trajectory_incidence_single <- function(obj, states, suffix) {
  ## This is less complicated than it looks, but takes a diff over the
  ## time dimension and converts back into the correct array dimension
  ## order.
  traj_inc <- aperm(
    apply(obj$state[states, , , drop = FALSE], c(1, 2), trajectory_incidence),
    c(2, 3, 1))
  rownames(traj_inc) <- paste0(states, suffix)
  obj$state <- abind1(obj$state, traj_inc)

  obj
}


add_trajectory_incidence_nested <- function(obj, states, suffix) {
  add_incidence <- function(a, b) {

    traj_inc <- aperm(apply(a[b, , , drop = FALSE],
                            c(1, 2), trajectory_incidence), c(2, 3, 1))
    rownames(traj_inc) <- paste0(b, suffix)
    abind1(a, traj_inc)
  }

  new_state <- array(NA,
                     c(length(states) + nrow(obj$state), dim(obj$state)[2:4]),
                     c(list(c(rownames(obj$state), paste0(states, suffix))),
                       dimnames(obj$state)[2:4]))

  for (i in seq_len(dim(obj$state)[[3L]])) {
    new_state[, , i, ] <- add_incidence(obj$state[, , i, ], states)
  }

  obj$state <- new_state

  obj
}

##' @rdname add_trajectory_incidence
##' @export
drop_trajectory_incidence <- function(obj) {
  if (inherits(obj, "mcstate_pmcmc")) {
    obj$trajectories <- drop_trajectory_incidence(obj$trajectories)
    return(obj)
  }

  assert_is(obj, "mcstate_trajectories")
  k <- grep("_inc$", rownames(obj$state))

  if (length(dim(obj$state)) == 3) {
    obj$state <- obj$state[-k, , ]
  } else {
    obj$state <- obj$state[-k, , , ]
  }

  obj
}


##' Combine trajectories across multiple runs
##'
##' @title Combine trajectories
##'
##' @param samples A list of samples from [mcstate::pmcmc]
##'
##' @param rank Logical, indicating if trajectories should be ranked
##'   before combination.
##'
##' @return A set of combined trajectories (not a combined samples object).
##'
##' @export
combine_trajectories <- function(samples, rank = TRUE) {
  states <- lapply(samples, function(x) x$trajectories$state)

  ## Ensure all trajectories are the same length
  dates <- lapply(samples, function(x) x$trajectories$date)
  dates_keep <- dates[[1]]
  for (i in seq_along(dates)[-1]) {
    dates_keep <- intersect(dates_keep, dates[[i]])
  }

  ## This really only needs doing if length(unique(dates)) > 1 but is
  ## relatively harmless otherwise.
  states <- lapply(seq_along(dates), function(i)
    states[[i]][, , dates[[i]] %in% dates_keep, drop = FALSE])

  if (rank) {
    states <- lapply(states, rank_trajectories)
  }

  ## This is a sum over regions, somewhat equivalent to a do.call over
  ## this list
  state <- Reduce(`+`, states)

  ## We count a date as "predicted" if any of the entries have a
  ## prediction on that day.
  predicted <- vapply(seq_along(dates), function(i)
    samples[[i]]$trajectories$predicted[dates[[i]] %in% dates_keep],
    logical(length(dates_keep)))
  predicted <- apply(predicted, 1, any)

  ## Create a mcstate_trajectories object based on the first element,
  ## filtered by the the dates. The state is updated but everything
  ## else is valid.
  trajectories <- samples[[1]]$trajectories
  trajectories$predicted <- predicted
  trajectories$date <- trajectories$date[dates[[1]] %in% dates_keep]
  trajectories$state <- state

  trajectories
}


rank_trajectories <- function(state) {
  ## Reorder the trajectories by the area under each curve
  for (i in seq_len(nrow(state))) {
    rank_state <- order(rowSums(state[i, , ], na.rm = TRUE))
    state[i, , ] <- state[i, rank_state, ]
  }

  state
}

##' Get the ranking of sample
##'
##' @title Get the ranking of sample by a given variable, e.g. infections
##'
##' @param sample An `mcstate_pmcmc` object
##'
##' @param by The name of the variable used for ranking. The ranking is
##' computed based on the value of this variable at the last time step.
##' Default is set to "infections", which will give the rank based on the
##' cumulative number of infections.
##'
##' @return A vector of integers giving the ranks
##'
##' @export
get_sample_rank <- function(sample, by = "infections") {
  if (!inherits(sample, "mcstate_pmcmc")) {
    stop("'sample' should be an 'mcstate_pmcmc' object")
  }
  if (!(by %in% rownames(sample$trajectories$state))) {
    stop(paste("Unkwnown 'by' argument. Should be one of: \n  -",
               paste(rownames(sample$trajectories$state), collapse = "\n  - ")))
  }
  last_time <- dim(sample$trajectories$state)[3]
  cum_inc <- sample$trajectories$state[by, , last_time]
  order(cum_inc, decreasing = FALSE)
}


##' Reorder samples according to a predefined ranking
##'
##' @title Reorder samples
##'
##' @param sample An `mcstate_pmcmc` object
##'
##' @param rank A vector of ranks to reorder by
##'
##' @return An `mcstate_pmcmc` object with appropriately reordered elements
##'
##' @export
reorder_sample <- function(sample, rank) {
  if (!inherits(sample, "mcstate_pmcmc")) {
    stop("'sample' should be an 'mcstate_pmcmc' object")
  }
  if (length(sample$iteration) != length(rank)) {
    stop(paste("Unexpected length for 'rank':", length(rank),
               "; should have length", length(sample$iteration)))
  }
  if (!is.null(sample$chain)) {
    sample$chain <- sample$chain[rank]
  }
  ## This looks like a no-op, but acts to drop any "full" parameters
  ## (using mcstate 0.8.2 or later).  Without this, the
  ## pars/probability subsetting will behave unexpectedly.
  sample <- mcstate::pmcmc_thin(sample)
  sample$iteration <- sample$iteration[rank]
  sample$pars <- sample$pars[rank, ]
  sample$probabilities <- sample$probabilities[rank, ]
  sample$state <- sample$state[, rank]
  sample$trajectories$state <- sample$trajectories$state[, rank, ]
  sample
}


##' Reorder Rt or IFR trajectories
##'
##' @title Reorder Rt or IFR trajectories
##'
##' @param x An `Rt_trajectories` object, as returned
##'   by [lancelot_Rt_trajectories()]
##'
##' @param rank A vector of ranks to reorder by
##'
##' @return An `Rt_trajectories` object with
##'   appropriately reordered elements
##'
##' @export
reorder_rt_ifr <- function(x, rank) {
  if (!(inherits(x, "Rt_trajectories") || inherits(x, "IFR_t_trajectories"))) {
    stop("'x' should be an 'Rt_trajectories' or 'IFR_t_trajectories' object")
  }
  what <- setdiff(names(x), c("time", "date"))
  if (ncol(x[[what[[1]]]]) != length(rank)) {
    stop(paste("Unexpected length for 'rank':", length(rank),
               "; should have length", ncol(x[[what[[1]]]])))
  }
  for (i in what) {
    x[[i]] <- x[[i]][, rank]
  }
  x
}


##' Combine Rt across multiple runs.
##'
##' @title Combine Rt estimates
##'
##' @param rt A list of Rt calculations from
##'   [lancelot_Rt_trajectories()] (though any Rt calculation that
##'   confirms to this will work)
##'
##' @param samples A list of samples from [mcstate::pmcmc]
##'
##' @param rank A boolean deciding whether to rank trajectories by increasing
##' incidence or not before combining Rt estimates
##'
##' @param weight A string specifying what trajectory to use for weighting,
##'   defaults to "infections_inc" (used for Rt calculation)
##'
##' @return A list of Rt output in the same structure as the first
##'   element of `rt`. All Rt estimates will be aggregated across
##'   regions (or whatever else you are aggregating on) based on the
##'   parameters in `samples`.
##'
##' @export
combine_rt <- function(rt, samples, rank = TRUE, weight = "infections_inc") {
  ## Ensure all trajectories are the same length
  ret <- rt[[1L]]
  what <- setdiff(names(ret), c("date", "time"))
  ret[what] <- lapply(what, combine_rt1, rt, samples, rank, weight)
  ret
}

##' Combine Rt estimates from EpiEstim across multiple runs.
##'
##' @title Combine Rt estimates from EpiEstim
##'
##' @param rt A list of Rt calculations from
##'   [lancelot_rt_trajectories_epiestim()] (though any Rt calculation that
##'   confirms to this will work)
##'
##' @param samples A list of samples from [mcstate::pmcmc]
##'
##' @param q A vector of quantiles to return values for
##'
##' @param rank A boolean deciding whether to rank trajectories by increasing
##' incidence or not before combining Rt estimates
##'
##' @return A list of Rt output in the same structure as the first
##'   element of `rt`. Rt estimates will be aggregated across
##'   regions (or whatever else you are aggregating on) based on the
##'   parameters in `samples`.
##'
##' @export
combine_rt_epiestim <- function(rt, samples, q = NULL, rank = TRUE) {
  q <- q %||% c(0.025, 0.5, 0.975)
  ## Ensure all trajectories are the same length
  ret <- rt[[1L]]
  if (!("Rt" %in% names(ret))) {
    stop(paste("rt$Rt missing. Did you forget 'save_all_Rt_sample = TRUE'",
               "in 'lancelot_EpiEstim_Rt_trajectories'?"))
  }
  ret$Rt <- combine_rt1_epiestim("Rt", rt, samples, rank)
  summary_R <- apply(ret$Rt, 2,
                     stats::quantile, q, na.rm = TRUE)
  mean_R <- apply(ret$Rt, 2, mean, na.rm = TRUE)
  summary_R <- rbind(summary_R, mean_R)
  ret$Rt_summary <- summary_R
  ret
}


combine_rt1 <- function(what, rt, samples, rank, weight) {
  dates <- lapply(samples, function(x) x$trajectories$date)
  dates_keep <- dates[[1]]
  for (i in seq_along(dates)[-1]) {
    dates_keep <- intersect(dates_keep, dates[[i]])
  }

  idx <- lapply(seq_along(samples), function(i) dates[[i]] %in% dates_keep)

  incidence <- Map(function(s, i)
    t(s$trajectories$state[weight, , i]), samples, idx)
    rt_what <- Map(function(r, i) r[[what]][i, ], rt, idx)

  if (rank) {
    ## Calculate rank of particles by area under Rt curve
    rank_x <- lapply(rt_what, function(x) order(colSums(x)))

    ## Rank based on the transpose (vs when this is done in the trajectories).
    reorder_by_rank <- function(x, rank_x) {
      Map(function(x, rank) x[, rank], x, rank_x)
    }

    ## Rank Rt and incidence
    x <- reorder_by_rank(rt_what, rank_x)
    w <- reorder_by_rank(incidence, rank_x)
  } else {
    x <- rt_what
    w <- incidence
  }
  sum_w <- Reduce(`+`, w)

  ## Weight Rt by incidence to combine regions
  ret <- Map(function(x, w) x * w / sum_w, x, w)

  Reduce(`+`, ret)
}


combine_rt1_epiestim <- function(what, rt, samples, rank) {

  dates <- lapply(samples, function(x) x$trajectories$date)
  dates_keep <- dates[[1]]
  for (i in seq_along(dates)[-1]) {
    dates_keep <- intersect(dates_keep, dates[[i]])
    dates_keep <- intersect(dates_keep, rt[[i]]$t_end)
  }
  idx <- lapply(seq_along(rt), function(i) which(rt[[i]]$t_end %in% dates_keep))

  incidence <- Map(function(s, i)
    s$trajectories$state["infections_inc", , i], samples, idx)

  rt_what <- Map(function(r, i) r[[what]][, i], rt, idx)

  n_R_per_traj <- round(nrow(rt_what[[1]]) / nrow(incidence[[1]]))
  incidence_rep <- lapply(incidence, function(m)
    m[rep(seq_len(nrow(m)), each = n_R_per_traj), ])

  if (rank) {
    ## Calculate rank of particles by area under Rt curve
    ## average across sets of n_R_per_traj Rs
    rank_x <- lapply(rt_what, function(x) order(rowSums(x, na.rm = TRUE)))

    ## Rank based on the transpose (vs when this is done in the trajectories).
    reorder_by_rank <- function(x, rank_x) {
      Map(function(x, rank) x[rank, ], x, rank_x)
    }

    ## Rank Rt and incidence
    x <- reorder_by_rank(rt_what, rank_x)
    w <- reorder_by_rank(incidence_rep, rank_x)
  } else {
    x <- rt_what
    w <- incidence_rep
  }
  sum_w <- Reduce(`+`, w)

  ## Weight Rt by incidence to combine regions
  ret <- Map(function(x, w) x * w / sum_w, x, w)

  Reduce(`+`, ret)
}


combine_steps_groups <- function(step, n_groups, n_time_steps, n_strains,
                                 n_vacc_classes, p_step, rel_p,
                                 strain_rel_p) {

  ret <- vapply(
    seq_len(n_groups),
    function(i) {
      outer(
        sircovid_parameters_expand_step(step, p_step[, i]),
        rel_p[i, , ] * strain_rel_p
      )
    },
    array(0, c(n_time_steps, n_strains, n_vacc_classes))
  )

  ## If these dimensions are all 1 then the above reduces to a vector
  ## rather than an array, so here we just reshape into an array
  if (all(c(n_time_steps, n_strains, n_vacc_classes) == 1)) {
    ret <- array(ret, c(1, 1, 1, length(ret)))
  }

  ret <- pmin(ret, 1)
  ret <- aperm(ret, c(4, 2, 3, 1))

  ret
}


unmirror_pars <- function(p) {
  which <-
    vapply(p,
           function(x)
             (is.matrix(x) && length(dim(x)) == 2 && ncol(x) == 4 &&
                identical(x[, 1:2], x[, 4:3])) ||
             (length(x) == 4 && identical(x[1:2], x[4:3])) ||
             (inherits(x, "array") && length(dim(x)) == 3 &&
                ncol(x) == 4 &&
                identical(x[, 1:2, ], x[, 4:3, ])),
           logical(1))
  p[which] <- lapply(p[which], unmirror_strain)
  p
}
