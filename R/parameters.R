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
sircovid_parameters_shared <- function() {
  list(psi = 0.1,
       hosp_transmission = 0.1,
       ICU_transmission = 0.05,
       dt = 0.25,
       N_age = length(sircovid_age_bins()$start))
}
