sircovid_date <- function(date) {
  days_into_2020 <- as.numeric(as_date(date) - as_date("2019-12-31"))
  if (any(days_into_2020 < 0)) {
    stop("Negative dates, sircovid_date likely applied twice")
  }
  days_into_2020
}


sircovid_date_as_date <- function(date) {
  assert_sircovid_date(date)
  as_date("2019-12-31") + date
}


assert_sircovid_date <- function(date) {
  if (!is.numeric(date)) {
    stop("'date' must be numeric - did you forget sircovid_date()?")
  }
  date
}


as_sircovid_date <- function(date) {
  if (is.character(date)) {
    sircovid_date(as_date(date))
  } else if (is_date(date)) {
    sircovid_date(date)
  } else {
    assert_sircovid_date(date)
  }
}
