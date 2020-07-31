## These could be moved to be defaults within the models
sircovid_parameters_shared <- function(start_date, region,
                                       beta_date, beta_value) {
  dt <- 0.25
  assert_sircovid_date(start_date)
  beta_step <- sircovid_parameters_beta(beta_date, beta_value %||% 0.08, dt)
  list(hosp_transmission = 0.1,
       ICU_transmission = 0.05,
       comm_D_transmission = 0.05,
       dt = dt,
       initial_step = start_date / dt,
       N_age = length(sircovid_age_bins()$start),
       beta_step = beta_step,
       population = sircovid_population(region))
}
