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

# Converts dates from data into a numeric offset as used in the MCMC
# Automatically converts type
start_date_to_offset <- function(first_data_date, start_date)
{
  if (class(first_data_date) == "Date" || class(start_date) == "Date") {
    stop("Do not use Date objects with offsets")
  }
  
  # Convert any strings to days since start of 2020
  if (class(first_data_date) == "character" || class(first_data_date) == "factor") {
    first_data_date = sircovid_date(first_data_date)
  }
  if (class(start_date) == "character" || class(start_date) == "factor") {
    start_date = sircovid_date(start_date)
  }

  first_data_date - start_date
}

# Converts dates from from numeric offset as used in the MCMC to a Date
# Automatically converts type
offset_to_start_date <- function(first_data_date, offset)
{
  if (class(start_date) != "numeric") {
    stop("Offset start date must be numeric")
  }

  # Convert any strings to Dates
  if (class(first_data_date) == "character" || class(first_data_date) == "factor") {
    sircovid_date(first_data_date)
  } else if (class(first_data_date) != "numeric") {
    stop("Start date must be numeric")
  }

  first_data_date - offset
}