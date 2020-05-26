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
##' @param type The class of the model, either
##'   \code{"sircovid_basic"},  \code{"sircovid_hospital"} or \code{"sircovid_serology"}
##' 
##' @export
##' 
##' @importFrom stats dbinom
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
  } else if (type == "sircovid_serology") {
    phi_general <- pars_obs$phi_general
    k_general <- pars_obs$k_general
    phi_admitted <- pars_obs$phi_admitted
    k_admitted <- pars_obs$k_admitted
    phi_new <- pars_obs$phi_new
    k_new <- pars_obs$k_new
    p_specificity <- pars_obs$p_specificity
    index_general <- c(c(index$I_triage_R_conf),c(index$I_triage_D_conf),c(index$I_hosp_R_conf),c(index$I_hosp_D_conf),c(index$R_stepdown_conf)) - 1L
    index_admit <- c(index$cum_admit_conf) - 1L
    index_new <- c(index$cum_new_conf) - 1L
    index_ICU <- c(c(index$I_ICU_R_conf),c(index$I_ICU_D_conf)) - 1L
    index_D <- c(index$D_hosp,c(index$D_comm)) - 1L
    index_R_pos <- c(index$R_pos) - 1L
    index_R_neg <- c(index$R_neg) - 1L
    index_R_pre <- index$R_pre - 1L #need to retain array structure for this
    
    tmp <- model$run(1)
    res <- model$transform_variables(tmp)
    N_tot <- res$N_tot[1,]
    
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
    
    if (type %in% c("sircovid_hospital","sircovid_serology")  && !is.na(data$general[t])) {
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
    
    if (type %in% c("sircovid_serology")  && !is.na(data$admitted[t])) {
      ## sum model output across ages/infectivities
      model_admitted <- colSums(state[index_admit, ,drop = FALSE]) -
        colSums(prev_state[index_admit, ,drop = FALSE])
      log_weights <- log_weights +
        ll_nbinom(data$admitted[t], model_admitted, phi_admitted, k_admitted, exp_noise)
    }
    
    if (type %in% c("sircovid_serology")  && !is.na(data$new[t])) {
      ## sum model output across ages/infectivities
      model_new <- colSums(state[index_new, ,drop = FALSE]) -
        colSums(prev_state[index_new, ,drop = FALSE])
      log_weights <- log_weights +
        ll_nbinom(data$new[t], model_new, phi_new, k_new, exp_noise)
    }
    
    if (type %in% c("sircovid_serology")  && !is.na(data$ntot_0_14[t]) && !is.na(data$npos_0_14[t])) {
      agegroups <- seq.int(1,3)
      prob_true_pos <- colSums(state[index_R_pos[agegroups], ,drop = FALSE]) / (sum(N_tot[agegroups]) - colSums(state[index_D[agegroups], ,drop = FALSE]))
      prob_false_pos <- (1 - p_specificity) * (1 - colSums(state[c(index_R_pos[agegroups], index_R_neg[agegroups], c(index_R_pre[agegroups,])), ,drop = FALSE]) / (sum(N_tot[agegroups]) - colSums(state[index_D[agegroups], ,drop = FALSE])))
        
      log_weights <- log_weights +
        dbinom(data$npos_0_14[t], size = data$ntot_0_14[t], prob = prob_true_pos + prob_false_pos, log = TRUE)
    }
    
    if (type %in% c("sircovid_serology")  && !is.na(data$ntot_0_14[t]) && !is.na(data$npos_0_14[t])) {
      agegroups <- seq.int(4,13)
      prob_true_pos <- colSums(state[index_R_pos[agegroups], ,drop = FALSE]) / (sum(N_tot[agegroups]) - colSums(state[index_D[agegroups], ,drop = FALSE]))
      prob_false_pos <- (1 - p_specificity) * (1 - colSums(state[c(index_R_pos[agegroups], index_R_neg[agegroups], c(index_R_pre[agegroups,])), ,drop = FALSE]) / (sum(N_tot[agegroups]) - colSums(state[index_D[agegroups], ,drop = FALSE])))
      
      log_weights <- log_weights +
        dbinom(data$npos_15_64[t], size = data$ntot_15_64[t], prob = prob_true_pos + prob_false_pos, log = TRUE)
    }
    
    if (type %in% c("sircovid_serology")  && !is.na(data$ntot_0_14[t]) && !is.na(data$npos_0_14[t])) {
      agegroups <- seq.int(14,17)
      prob_true_pos <- colSums(state[index_R_pos[agegroups], ,drop = FALSE]) / (sum(N_tot[agegroups]) - colSums(state[index_D[agegroups], ,drop = FALSE]))
      prob_false_pos <- (1 - p_specificity) * (1 - colSums(state[c(index_R_pos[agegroups], index_R_neg[agegroups], c(index_R_pre[agegroups,])), ,drop = FALSE]) / (sum(N_tot[agegroups]) - colSums(state[index_D[agegroups], ,drop = FALSE])))
      
      log_weights <- log_weights +
        dbinom(data$npos_65plus[t], size = data$ntot_65plus[t], prob = prob_true_pos + prob_false_pos, log = TRUE)
    }
    
    log_weights
  }
}