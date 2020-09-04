## General things that we need but that aren't that interesting.

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


##' @importFrom stats dnbinom rexp
ll_nbinom <- function(data, model, k, exp_noise) {
  if (is.na(data)) {
    return(numeric(length(model)))
  }
  mu <- model + rexp(length(model), rate = exp_noise)
  dnbinom(data, k, mu = mu, log = TRUE)
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


sircovid_transmission_matrix <- function(region) {
  if (is.null(cache$transmission_matrix[[region]])) {
    if (is.null(cache$transmission_matrix)) {
      cache$transmission_matrix <- list()
    }
    age_bins <- sircovid_age_bins()

    max_polymod_age <- 70

    ## Survey population in socialmixr format
    survey_pop <- data_frame(
      "lower.age.limit" = sircovid_age_bins()$start,
      population = sircovid_population(region))

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
