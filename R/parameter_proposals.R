# Generates beta_step correctly for use in the odin model
update_beta <- function(sircovid_model, 
                        beta_start,
                        beta_end,
                        beta_pl,
                        start_date,
                        dt) {
  new_beta <- sircovid_model$generate_beta_func(beta_start = beta_start,
                                                beta_end = beta_end, 
                                                beta_pl = beta_pl,
                                                start_date = start_date) 

  beta_t <- normalise_beta(new_beta$beta_times, dt)

  beta_fun <- cinterpolate::interpolation_function(beta_t, new_beta$beta,
                                                   "constant")
  beta_step <- beta_fun(seq.int(min(beta_t), max(beta_t)))
  
  list(beta_step = beta_step)
}

# Generates beta_step correctly for use in the odin model
update_beta_piecewise_linear <- function(sircovid_model, 
                                         beta_k,
                                         t_k,
                                         start_date,
                                         dt) {
  new_beta <- sircovid_model$generate_beta_func(beta_k = beta_k,
                                                t_k = t_k, 
                                                start_date = start_date,
                                                dt = dt) 
  
  beta_t <- normalise_beta(new_beta$beta_times, dt)

  beta_fun <- cinterpolate::interpolation_function(beta_t, new_beta$beta,
                                                   "constant")
  beta_step <- beta_fun(seq.int(min(beta_t), max(beta_t)))
  
  list(beta_step = beta_step)
}
