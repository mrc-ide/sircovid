context("sample_grid_scan")

# Only tests that a grid search can be run
test_that("sample_grid_scan works", {
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
  
  trajectories <- sample_grid_scan(scan_results = scan_results,
                                       n_sample_pairs = 10, 
                                       n_particles = 100)
  
  expect_equal(dim(sample_grid_scan), c(69, 238, 10))
  
  ## Testing plotting is always a nightmare
  if (FALSE) {
    mod <- sircovid(params = scan_results$inputs$model_params)
    index <- c(odin_index(mod)$I_ICU) - 1L
    particles <- apply(trajectories[, index, ], c(1, 3), sum)
    plot_particles(particles, ylab = "ICU")
    points(as.Date(data$date), data$itu / pars_obs$phi_ICU, pch = 19)
  }
  
})
