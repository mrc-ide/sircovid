##' Compare the model to ICU data for use with the particle filter
##' 
##' @title Compare model to ICU data
##' 
##' @param model An \code{odin_model} object
##' 
##' @param pars_obs Parameters for the observations
##' 
##' @param data The data to be compared against
##' 
##' @export
compare_output <- function(model, pars_obs, data, type="sircovid_basic") {
  index <- odin_index(model)

  ## Unpack things that we will use repeatedly
  phi_ICU <- pars_obs$phi_ICU
  k_ICU <- pars_obs$k_ICU
  phi_death <- pars_obs$phi_death
  k_death <- pars_obs$k_death
  exp_noise <- pars_obs$exp_noise
  if (type == "sircovid_basic") {
    index_ICU <- c(index$I_ICU) - 1L
    index_D <- c(index$D) - 1L
  } else if (type == "sircovid_hospital") {
    phi_general <- pars_obs$phi_general
    k_general <- pars_obs$k_general
    index_general <- c(c(index$I_triage),c(index$I_hosp_R),c(index$I_hosp_D),c(index$R_stepdown)) - 1L
    index_ICU <- c(c(index$I_ICU_R),c(index$I_ICU_D)) - 1L
    index_D <- c(index$D) - 1L
  }

  force(data)

  ## This returns a closure, with the above variables bound, the
  ## sampler will provide the arguments below.
  function(t, state, prev_state) {
    if (is.na(data$itu[t] && is.na(data$deaths[t]))) {
      return(NULL)
    }

    log_weights <- rep(0, ncol(state))

    if (!is.na(data$itu[t])) {
      ## sum model ITU cases output across ages/infectivities
      model_icu <- colSums(state[index_ICU, ])
      log_weights <- log_weights +
        ll_nbinom(data$itu[t], model_icu, phi_ICU, k_ICU, exp_noise)
    }
    
    if (type == "sircovid_hospital" && !is.na(data$general[t])) {
      ## sum model output across ages/infectivities
      model_general <- colSums(state[index_general, ])
      log_weights <- log_weights +
        ll_nbinom(data$general[t], model_general, phi_general, k_general, exp_noise)
    }

    if (!is.na(data$deaths[t])) {
      ## new deaths summed across ages/infectivities
      model_deaths <- colSums(state[index_D, ]) -
        colSums(prev_state[index_D, ])
      log_weights <- log_weights +
        ll_nbinom(data$deaths[t], model_deaths, phi_death, k_death, exp_noise)
    }

    log_weights
  }
}
