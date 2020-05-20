# Generates beta_y and beta_t correctly for use in the odin model
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
  
  list(beta_y = new_beta$beta,
       beta_t = beta_t)
}


