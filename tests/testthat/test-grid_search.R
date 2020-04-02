context("grid_search")

# Only tests that a grid search can be run
test_that("Small grid search works", {
  set.seed(1)

  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)
 
  # Parameters for run
  min_beta <- 0.1
  max_beta <- 0.2
  beta_step <- 0.05
  first_start_date <- "2020-01-21"
  last_start_date <- "2020-01-22"
  day_step <- 1
  
  scan_results = scan_beta_date(
    min_beta = min_beta,
    max_beta = max_beta,
    beta_step = beta_step,
    first_start_date = first_start_date, 
    last_start_date = last_start_date, 
    day_step = day_step,
    data = data)
   
  expect_is(scan_results, "sircovid_scan")
  
  beta_grid = seq(min_beta, max_beta, beta_step)
  date_grid = seq(as.Date(first_start_date), as.Date(last_start_date), day_step)
  expect_equal(scan_results$x, beta_grid)
  expect_equal(scan_results$y, date_grid)
  expect_equal(dim(scan_results$renorm_mat_LL), dim(scan_results$mat_log_ll))
  expect_equal(dim(scan_results$renorm_mat_LL), c(length(beta_grid), length(date_grid)))
  expect_true(all(scan_results$renorm_mat_LL < 1 & scan_results$renorm_mat_LL > 0))
  
  # Plots run, but not checked
  plot(scan_results)
})

test_that("Transmission is more likely", {
  set.seed(1)
  
  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)
  
  # Parameters for run
  min_beta <- 0
  max_beta <- 0.1
  beta_step <- 0.1
  first_start_date <- "2020-01-21"
  last_start_date <- "2020-01-21"
  day_step <- 1
  
  scan_results = scan_beta_date(
    min_beta = min_beta,
    max_beta = max_beta,
    beta_step = beta_step,
    first_start_date = first_start_date, 
    last_start_date = last_start_date, 
    day_step = day_step,
    data = data)
  
  # No transmission b = 0 much less likely than some transmission b = 0.1
  expect_lt(scan_results$renorm_mat_LL[[1]], scan_results$renorm_mat_LL[[2]])
  
})

test_that("Unreasonable start dates are less likely", {
  set.seed(1)
  
  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)
  
  # Parameters for run
  min_beta <- 0.1
  max_beta <- 0.1
  beta_step <- 0.1
  first_start_date <- "2020-01-01"
  last_start_date <- "2020-02-29"
  day_step <- 20
  
  scan_results = scan_beta_date(
    min_beta = min_beta,
    max_beta = max_beta,
    beta_step = beta_step,
    first_start_date = first_start_date, 
    last_start_date = last_start_date, 
    day_step = day_step,
    data = data)
  
  # Mid Jan start most likely
  expect_gt(scan_results$renorm_mat_LL[[2]], scan_results$renorm_mat_LL[[1]])
  expect_gt(scan_results$renorm_mat_LL[[2]], scan_results$renorm_mat_LL[[3]])
  
})
