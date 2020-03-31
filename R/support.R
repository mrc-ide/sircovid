default_age_distribution <- function() {
  if (is.null(cache$age_distribution)) {
    cache$age_distribution <-
      read.csv(sircovid_file("extdata/age_distribution.csv"),
               stringsAsFactors = FALSE)
  }
  cache$age_distribution
}
