sircovid_infection_seeding <- function(values, bins) {
  if (length(values) != length(bins)) {
    stop("Each infection seeding value must correspond to one bin")
  }
  ret <- list(values = values, bins)
  class(ret) <- "sircovid_infection_seeding"
  ret
}


## We always use these age bands, so rather than detect them, we will
## check that things conform to them.
sircovid_age_bins <- function() {
  end <- c(seq(4, 79, by = 5), 100)
  start <- c(0, end[-length(end)] + 1L)
  list(start = start, end = end)
}


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


sircovid_population <- function(region) {
  if (is.null(region)) {
    stop("'region' must not be NULL")
  }

  ## TODO: cache this file read as it's constant within a session
  data <- read_csv(sircovid_file("extdata/population.csv"))
  population <- data[[tolower(region)]]
  if (is.null(population)) {
    valid <- paste(squote(setdiff(names(data), "age")), collapse = ", ")
    stop(sprintf("Population not found for '%s': must be one of %s",
                 region, valid))
  }

  population
}
