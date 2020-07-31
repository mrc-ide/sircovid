sircovid_data <- function(data, start_date, dt) {
  ## Some horrid off-by-one still lurking here. See here for more
  ## details, and the accompanying PR
  ## https://github.com/mrc-ide/mcstate/commit/97e68ade560c9028204e691bc7b57ef2ef2ef557
  ## To make this work, we've manually inserted a fake reporting
  ## period at the first row of the file so that our compare works
  ## correctly; this should be something that mcstate can do for us.
  data$date <- sircovid_date(data$date)
  rate <- 1 / pars$dt
  data <- mcstate::particle_filter_data(data, "date", rate, start_date)
  data
}
