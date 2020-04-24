##' Seeding function for use with the particle filter
##' 
##' @title Compare model to ICU data
##' 
##' @param model An \code{odin_model} object
##' 
##' @param data The data, must be constructed with
##'   \code{particle_filter_data}
##'
##' @param pars_seeding A list of parameters used for seeding,
##' e.g. pars_seeding = list(lambda = 20) 
##' 
##' @importFrom extraDistr rtpois
##' 
##' @export
seeding_function <- function(model,
                             data,
                             pars_seeding = NULL
){
  
  steps_per_day <- attr(data, "steps_per_day")
  min_seeding_step <- data$step_start[1L] #earliest possible seeding time
  max_seeding_step <- data$step_end[1L] #latest possible seeding time (should be 1 day before first_data_step)
  first_data_step <- data$step_end[2L] #first data point time
  
  function(n_particles){
    if (is.null(pars_seeding)){
    
      #if pars_seeding is NULL, then we seed all particles at start_date with same state  
      step <- rep(min_seeding_step,n_particles)
      initial_state <- model$initial()
      state <- array(initial_state,dim=c(length(initial_state),n_particles))
    
    } else {
    
      #variable seeding time
      #sample a poisson number of time steps before the first data time point. Truncated such that
      #min_seeding_step <= step <= max_seeding_step
      step <- first_data_step - rtpois(n = n_particles, lambda = pars_seeding$lambda*steps_per_day,
                                       a = first_data_step - max_seeding_step - 1, b = first_data_step - min_seeding_step)
    
      #this does the same as when pars_seeding is NULL, but the seeding state could also be varied between particles
      initial_state <- model$initial()
      state <- array(initial_state,dim=c(length(initial_state),n_particles))
    }
  
    list(step = step, state = state)
  }

}