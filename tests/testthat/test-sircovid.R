context("sircovid")

test_that("Can be run on real data", {
  set.seed(1)
  time_steps_per_day <- 4
  
  data <- generate_data(death_data_file = "covid_cases_2020_4_3.csv",
                        admissions_data_file = "combin_time_series.csv")
  
  vary_beta <- generate_beta(0.1)
  model_params <- generate_parameters(beta = vary_beta$beta,
                                      beta_times = vary_beta$beta_times,
                                      dt = 1/time_steps_per_day)
  
  results <- run_particle_filter(data = data,
                                 model_params = model_params,
                                 obs_params = list(phi_ICU = 0.95,
                                                   k_ICU = 2,
                                                   phi_death = 1789/1651,
                                                   k_death = 2,
                                                   exp_noise = 1e6),
                                 n_particles = 1000)
  # No check of correctness  
  expect_equal(results$log_likelihood, -269.9478, tolerance=1e-3)
})

test_that("Poor formatting of real data errors", {
  # Try some bad files and dates
  expect_error(generate_data(death_data_file = "covid_cases_2020_4_3.csv",
                             admissions_data_file = "combin_time_series.csv",
                             itu_colname = "itu"),
               "itu_colname not found in time_series_file")
  expect_error(generate_data(death_data_file = "covid_cases_2020_4_3.csv",
                             admissions_data_file = "combin_time_series.csv",
                             first_data_date = "2022-01-01"),
               "No case data matching country and date filtering criteria")
  
  data <- generate_data(death_data_file = "covid_cases_2020_4_3.csv",
                        admissions_data_file = "combin_time_series.csv",
                        first_data_date = "2020-01-01")
  expect_error(run_particle_filter(data = data,
                                   model_start_date = "2020-03-02",
                                   obs_params = list(phi_ICU = 0.95,
                                                     k_ICU = 2,
                                                     phi_death = 1789/1651,
                                                     k_death = 2,
                                                     exp_noise = 1e6),
                                   n_particles = 1000),
               "Model start date is later than data start date")
})

