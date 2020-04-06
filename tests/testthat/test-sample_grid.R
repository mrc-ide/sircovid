context("sample_grid_scan")

# Only tests that a grid search can be run
test_that("sample_grid_scan works", {
  set.seed(1)
  
  # grab the data
  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)
  
  
  # filter and tidy it
  # Don't start the data until it is after thid date
  data <- data[data$date > as.Date("2020-02-29"),]
  
  # which dates should be set to NA
  na_dates <- as.Date(c("2020-03-01","2020-03-02","2020-03-03","2020-03-04",
                        "2020-03-05","2020-03-06","2020-03-07","2020-03-08",
                        "2020-03-15","2020-03-16","2020-03-19"))
  data$itu[data$date %in% na_dates] <- NA
  
  # bring the deaths back by 2 days as ECDC is out of sync
  data$deaths[head(seq_len(nrow(data)), -2)] <- tail(data$deaths, -2)
  data$deaths[tail(seq_len(nrow(data)), 2)] <- NA
  
  
  # Parameters for run
  min_beta <- 0.10
  max_beta <- 0.18
  beta_step <- 0.08
  first_start_date <- "2020-01-29"
  last_start_date <- "2020-02-14"
  day_step <- 12
  
  scan_results <- scan_beta_date(
    min_beta = min_beta,
    max_beta = max_beta,
    beta_step = beta_step,
    first_start_date = first_start_date, 
    last_start_date = last_start_date, 
    day_step = day_step,
    data = data)
  
  n_sample_pairs <- 4 
  res <- sample_grid_scan(scan_results = scan_results,
                                       n_sample_pairs = n_sample_pairs, 
                                       n_particles = 10)
  
  # check length based on model and dates
  days_between <- length( min(as.Date(res$param_grid$start_date)) : as.Date(tail(rownames(res$trajectories[,,1]),1)))
  expect_equal(dim(res$trajectories), c(days_between, length(res$inputs$model$initial()), n_sample_pairs))
  
  ## Testing plotting is always a nightmare
  if (TRUE) {
    plot(res, what = "ICU")
    plot(res, what = "Deaths")
  }
  
})
