##' Seeding function for use with the particle filter
##' 
##' @title Compare model to ICU data
##' 
##' @param model An \code{odin_model} object
##' 
##' @param data The data, must be constructed with
##'   \code{particle_filter_data}
##'
##' @param pars_seeding Parameters for seeding
##' 
##' @importFrom stats rpois
##' 
##' @export
seeding_function <- function(model,
                             data,
                             pars_seeding = NULL
){
  
  steps_per_day <- attr(data, "steps_per_day")
  
  function(n_particles){
    if (is.null(pars_seeding)){
    
      #if pars_seeding is NULL, then we seed all particles at start_date with same state  
      step <- rep(data$step_start[1L],n_particles)
      initial_state <- model$initial()
      state <- array(initial_state,dim=c(length(initial_state),n_particles))
    
    } else {
    
      #variable seeding time  
      step <- data$step_end[1L] - rpois(n_particles,pars_seeding$lambda*steps_per_day)
      while (any(step < data$step_start[1L])){
        #make sure no seeding time is before start_date
        k <- which(step < data$step_start[1L])
        step[k] <- data$step_end[1L] - rpois(length(k),pars_seeding$lambda*steps_per_day)
      }
    
      #this does the same as when pars_seeding is NULL, but the seeding state could also be varied between particles
      initial_state <- model$initial()
      state <- array(initial_state,dim=c(length(initial_state),n_particles))
    }
  
    list(step = step, state = state)
  }

}