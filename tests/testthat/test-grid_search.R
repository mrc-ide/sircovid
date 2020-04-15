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
    data = data,
    sircovid_model = )
   
  expect_is(scan_results, "sircovid_scan")
  expect_true("inputs" %in% names(scan_results))
  expect_setequal(names(scan_results$inputs),
                  c("model", "model_params", "pars_obs", "data"))
  
  beta_grid <- seq(min_beta, max_beta, beta_step)
  date_grid <- seq(as.Date(first_start_date), as.Date(last_start_date), day_step)
  expect_equal(scan_results$x, beta_grid)
  expect_equal(scan_results$y, date_grid)
  expect_equal(dim(scan_results$renorm_mat_LL), dim(scan_results$mat_log_ll))
  expect_equal(dim(scan_results$renorm_mat_LL), c(length(beta_grid), length(date_grid)))
  expect_true(all(scan_results$renorm_mat_LL < 1 & scan_results$renorm_mat_LL > 0))
  
  # Plots run, but not checked
  plot(scan_results)
})

test_that("Small grid search works with new model", {
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
  
  skip("TODO: Add patched parameters here")
  scan_results = scan_beta_date(
    min_beta = min_beta,
    max_beta = max_beta,
    beta_step = beta_step,
    first_start_date = first_start_date, 
    last_start_date = last_start_date, 
    day_step = day_step,
    data = data,
    model = "new_hospital_model")

  expect_is(scan_results, "sircovid_scan")
  expect_true("inputs" %in% names(scan_results))
  expect_setequal(names(scan_results$inputs),
                  c("model", "model_params", "pars_obs", "data"))
  
  beta_grid <- seq(min_beta, max_beta, beta_step)
  date_grid <- seq(as.Date(first_start_date), as.Date(last_start_date), day_step)
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

test_that("Varying beta is set in the right place", {
  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)
  
  # Parameters for run
  min_beta <- 0.1
  max_beta <- 0.1
  beta_step <- 0.1
  first_start_date <- "2020-01-01"
  last_start_date <- "2020-02-29"
  day_step <- 20
  
  model_params <- generate_parameters(
    transmission_model = "POLYMOD",
    progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 1, hosp = 2, ICU = 2, rec = 2),
    gammas = list(E = 1/2.5, asympt = 1/2.09, mild = 1/2.09, ILI = 1/4, hosp = 2/1, ICU = 2/5, rec = 2/5),
    hosp_transmission = 0,
    ICU_transmission = 0,
    trans_profile = 1,
    trans_increase = 1,
    dt = 0.25)
  
  expect_error(scan_beta_date(
    min_beta = min_beta,
    max_beta = max_beta,
    model_params = model_params,
    beta_step = beta_step,
    first_start_date = first_start_date, 
    last_start_date = last_start_date, 
    day_step = day_step,
    data = data),
    "Set beta variation through generate_beta, not model_params")
})
