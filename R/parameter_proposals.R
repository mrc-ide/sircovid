# Generates beta_y and beta_t correctly for use in the odin model
update_beta <- function(sircovid_model, 
                        beta_start,
                        beta_end,
                        start_date,
                        dt) {
  if(!is.null(beta_end)) {
    
    new_beta <- sircovid_model$generate_beta_func(beta_start = beta_start,
                                                  beta_end = beta_end, 
                                                  start_date = start_date) 
  } else {
    
    new_beta <- sircovid_model$generate_beta_func(beta_start = beta_start, 
                                                  start_date = start_date)

  }
  
  beta_t <- normalise_beta(new_beta$beta_times, dt)
  
  list(beta_y = new_beta$beta,
       beta_t = beta_t)
}

# Converts dates from data into a numeric offset as used in the MCMC
# Automatically converts type
start_date_to_offset <- function(first_data_date, start_date)
{
  # Format conversion cascades as required
  # string -> Date -> numeric

  # Convert any strings to Dates
  if (class(first_data_date) == "character") {
    first_data_date = as.Date(first_data_date)
  }
  if (class(start_date) == "character") {
    start_date = as.Date(start_date)
  }

  # Convert any Dates to numerics
  if (class(first_data_date) == "Date") {
    first_data_date = as.numeric(first_data_date)
  }
  if (class(start_date) == "Date") {
    start_date = as.numeric(start_date)
  }

  first_data_date - start_date
}

# Converts dates from data into a numeric offset as used in the MCMC
# Automatically converts type
offset_to_start_date <- function(first_data_date, start_date)
{
  # Format conversion cascades as required
  # string -> Date -> numeric

  # Convert any strings to Dates
  if (class(first_data_date) == "character") {
    first_data_date = as.Date(first_data_date)
  }
  if (class(start_date) == "character") {
    start_date = as.Date(start_date)
  }

  # Convert any Dates to numerics
  if (class(start_date) == "Date") {
    start_date = as.numeric(start_date)
  }

  as.Date(-start_date, origin=first_data_date)
}