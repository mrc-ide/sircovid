sircovid_parameters_beta <- function(date, value, dt) {
  if (is.null(date)) {
    if (length(value) != 1L) {
      stop("As 'date' is NULL, expected single value")
    }
    return(value)
  }
  if (length(date) != length(value)) {
    stop("'date' and 'value' must have the same length")
  }
  if (length(date) < 2) {
    stop("Need at least two dates and betas for a varying beta")
  }
  if (!is.numeric(date) || any(date < 0)) {
    stop("'date' must be a sircovid_date")
  }
  assert_increasing(date)
  approx(c(0, date), c(value[[1]], value),
         seq(0, date[[length(date)]], by = dt))$y
}
