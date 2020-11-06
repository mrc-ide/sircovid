build_waning_rate <- function(waning_rate) {
  if (length(waning_rate) == 0) {
    stop("At least one value required for 'rel_susceptibility'")
  }
  if (any(waning_rate < 0)) {
    stop("'waning_rate' must have only non-negative values")
  }
  if (length(waning_rate) == 1) {
    waning_rate <- rep(waning_rate, carehomes_n_groups())
  } else if (length(waning_rate) != carehomes_n_groups()) {
    stop(
      "'waning_rate' should have as many elements as age groups or be a scalar")
  }
  waning_rate
}