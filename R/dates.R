##' We need to map "dates" onto [`dust::dust`]'s concept of model
##' "step" and we do this by mapping a date such as `2020-03-02` into
##' the number of days into 2020 (62 here, with the 1st of January
##' being day 1). We call this integer number a "sircovid date".
##'
##' There are several related functions here
##'
##' * `sircovid_date` converts its argument into an R `Date` object,
##'   then applies this tranformation. If the argument is not a `Date`
##'   object or a string representing one, an error will be thrown.
##'
##' * `sircovid_date_to_date` does the reverse conversion to
##'   `sircovid_date`, converting an integer sircovid date into an R
##'   `Date`
##'
##' * `as_sircovid_date` does the same conversion as `sircovid_date`
##'   but will assume that an integer *already* represents a sircovid
##'   date and will return it unmodified rather than erroring.
##'
##' * `as_date` does a string to date conversion, using [as.Date()]
##'   but requiring the dates are in ISO 8601 (YYYY-MM-DD) format (it
##'   is a helper that avoids conversion to `NA`, instead throwing an
##'   error)
##'
##' @title Date handling for sircovid2
##'
##' @param date A Date object, or something that can be converted to
##'   one, or a "sircovid date"; see Details
##'
##' @return An integer, being the number of days into 2020
##' @export
##' @examples
##' # Convert dates into sircovid dates:
##' sircovid2::sircovid_date("2020-01-01")
##' sircovid2::sircovid_date(c("2020-03-01", "2020-10-01"))
##'
##' # Reverse the conversion:
##' sircovid2::sircovid_date_as_date(1)
##' sircovid2::sircovid_date_as_date(c(61, 275))
##'
##' # Double conversion not possible with sircovid_date...
##' try(sircovid2::sircovid_date(61))
##' # ...but allowed with as_sircovid_date
##' sircovid2::as_sircovid_date(61)
##'
##' # Strict date conversion with as_date
##' sircovid2::as_date("2020-03-01")
##' try(sircovid2::as_date("03-01-2020"))
sircovid_date <- function(date) {
  days_into_2020 <- as.numeric(as_date(date) - as_date("2019-12-31"))
  if (any(days_into_2020 < 0)) {
    stop("Negative dates, sircovid_date likely applied twice")
  }
  days_into_2020
}


##' @export
##' @rdname sircovid_date
sircovid_date_as_date <- function(date) {
  assert_sircovid_date(date)
  as_date("2019-12-31") + date
}


##' @export
##' @rdname sircovid_date
as_sircovid_date <- function(date) {
  if (is.character(date) || is_date(date)) {
    sircovid_date(as_date(date))
  } else {
    assert_sircovid_date(date)
  }
}

##' @export
##' @rdname sircovid_date
as_date <- function(date) {
  if (is_date(date)) {
    return(date)
  }
  if (!all(grepl("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", date))) {
    stop("Expected ISO dates or R dates - please convert")
  }
  as.Date(date)
}


assert_sircovid_date <- function(date) {
  if (!is.numeric(date)) {
    stop("'date' must be numeric - did you forget sircovid_date()?")
  }
  if (any(date < 0)) {
    stop("Negative dates, sircovid_date likely applied twice")
  }
  date
}

is_date <- function(x) {
  inherits(x, "Date")
}
