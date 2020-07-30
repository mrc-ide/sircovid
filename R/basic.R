basic_parameters <- function(beta_date = NULL, beta_value = NULL,
                             severity_path = NULL) {
  ret <- sircovid_parameters_shared()
  ret$beta_step <-
    sircovid_parameters_beta(beta_date, beta_value %||% 0.08, ret$dt)
  ret$m <- sircovid_transmission_matrix()
  c(ret,
    sircovid_parameters_severity(severity_path),
    basic_parameters_progression())
}


basic_parameters_progression <- function() {
  ## These need to be aligned with Bob's severity outputs, and we will
  ## come up with a better way of correlating the two.

  ## The s_ parameters are the scaling parameters for the Erlang
  ## distibution (a.k.a 'k'), while the gamma parameters are the gamma
  ## parameters of that distribution.
  list(s_E = 2,
       s_asympt = 1,
       s_mild = 1,
       s_ILI = 1,
       s_hosp = 2,
       s_ICU = 2,
       s_rec = 2,

       gamma_E = 1 / (4.59 / 2),
       gamma_asympt = 1 / 2.09,
       gamma_mild = 1 / 2.09,
       gamma_ILI = 1 / 4,
       gamma_hosp = 2,
       gamma_ICU = 2 / 5,
       gamma_rec = 2 / 5)
}
