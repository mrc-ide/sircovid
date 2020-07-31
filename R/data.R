##' Convert data into mcstate format; this is a thin wrapper around
##' [mcstate::particle_filter_data()] which adds a dummy step in front
##' of the first data point so that we can use the previous state and
##' the current states to convert cumulative measures into net daily
##' changes.
##'
##' @title Prepare data for mcstate
##'
##' @param data A `data.frame` object suitable for
##'   [mcstate::particle_filter_data()]
##'
##' @param start_date The start date, as a [sircovid_date()], R "Date"
##'   object or a string in ISO 8601 format (YYYY-MM-DD)
##'
##' @param dt The time step (fraction of a day that each step
##'   represents) as used to create the model object
##'
##' @return A data.frame suitable for use with `mcstate` functions
##'   such as [mcstate::particle_filter()] and [mcstate::pmcmc()]
##'
##' @export
##' @examples
##' # A data sert that has data from the first of February to the first of
##' # March (one column of data called 'x')
##' from <- as.Date("2020-02-01")
##' to <- as.Date("2020-03-01")
##' d <- data.frame(date = seq(from, to, by = 1),
##'                 x = runif(to - from + 1),
##'                 stringsAsFactors = FALSE)
##'
##' # Get this ready for sircovid/mcstate assuming the seeding starts on
##' # the 15th of January and we take 4 steps per day.
##' sircovid_data(d, start_date = "2020-01-15", 1 / 4)
sircovid_data <- function(data, start_date, dt) {
  start_date <- as_sircovid_date(start_date)
  ## Some horrid off-by-one unpleasantness lurking here. See this commit:
  ##   https://github.com/mrc-ide/mcstate/commit/97e68ad
  ## for for more details, and the accompanying PR.
  ##
  ## To make this work, we've manually inserted a fake reporting
  ## period at the first row of the file so that our compare works
  ## correctly; this should be something that mcstate can do for us.
  data$date <- sircovid_date(data$date)
  rate <- 1 / dt
  data <- mcstate::particle_filter_data(data, "date", rate, start_date)
  data
}
