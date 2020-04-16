##' Create a basic model 
##' 
##' @title Basic model
##' 
##' @param progression_groups List of number of progression groups in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp', 'ICU', 'rec'
##'   
##' @param gammas List of exponential distribution rates for time in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp', 'ICU', 'rec'
##'   
##' @export
basic_model <- function(progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 1, hosp = 2, ICU = 2, rec = 2),
                        gammas = list(E = 1/(4.59/2), asympt = 1/2.09, mild = 1/2.09, ILI = 1/4, hosp = 2, ICU = 2/5, rec = 2/5)) {
  model_class <- "sircovid_basic"
  odin_model <- load_odin_model("basic")
  generate_beta_func <- generate_beta
  compare_model <- function(model, pars_obs, data) {compare_output(model, pars_obs, data, type=model_class)}
  
  model_partitions <- c("E", "asympt", "mild", "ILI", "hosp", "ICU", "rec")
  if (any(!(model_partitions %in% names(progression_groups)))) {
    stop("progression_groups need to be defined for all partitions")
  }
  if (any(!(model_partitions %in% names(gammas)))) {
    stop("gammas need to be defined for all partitions")
  }
  
  basic <- list(odin_model = odin_model,
                generate_beta_func = generate_beta_func,
                progression_groups = progression_groups,
                gammas = gammas,
                compare_model = compare_model)
  class(basic) <- (model_class)
  basic
}

##' Create a hosptial model
##' 
##' @title Hosptial model
##' 
##' @param progression_groups List of number of progression groups in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp_D', 'hosp_R', 'ICU_D', 'ICU_R', 'triage', 'stepdown', 'rec'
##'   
##' @param gammas List of exponential distribution rates for time in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp_D', 'hosp_R', 'ICU_D', 'ICU_R', 'triage', 'stepdown', 'rec'
##'   
##' @export
hospital_model <- function(use_fitted_parameters = TRUE,
                           progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 1, hosp_D = 2 , hosp_R = 2, ICU_D = 2, ICU_R = 2, triage = 2, stepdown = 2),
                           gammas = list(E = 1/(4.59/2), asympt = 1/2.09, mild = 1/2.09, ILI = 1/4, hosp_D = 2/5, hosp_R = 2/10, ICU_D = 2/5, ICU_R = 2/10, triage = 2, stepdown = 2/5)) {
  model_class <- "sircovid_hospital"
  odin_model <- load_odin_model("new_hospital_model")
  generate_beta_func <- generate_beta
  compare_model <- function(model, pars_obs, data) {compare_output(model, pars_obs, data, type=model_class)}
  
  if (use_fitted_parameters) {
    fitted_parameters <- read_fitted_parameters()
    progression_groups <- fitted_parameters$progression_groups 
    gammas <- fitted_parameters$gammas
  }
  
  model_partitions <- partition_names(model_class)
  if (any(!(model_partitions %in% names(progression_groups)))) {
    stop("progression_groups need to be defined for all partitions")
  }
  if (any(!(model_partitions %in% names(gammas)))) {
    stop("gammas need to be defined for all partitions")
  }
  
  hospital <- list(odin_model = odin_model,
                generate_beta_func = generate_beta_func,
                progression_groups = progression_groups,
                gammas = gammas,
                compare_model = compare_model)
  class(hospital) <- c(model_class, "sircovid_basic")
  hospital
}

#
# Internal functions
#

# Definitions of partition names
partition_names <- function(model_name) {
  if (inherits(model_name, "sircovid_basic")) {
    model_partitions <- names(model_name$progression_groups)
  } else {
    if (model_name == "sircovid_basic") {
      model_partitions <- c("E", "asympt", "mild", "ILI", "hosp", "ICU", "rec")
    } else if (model_name == "sircovid_hospital") {
      model_partitions <- c("E", "asympt", "mild", "ILI", "hosp_D", "hosp_R", "ICU_D", "ICU_R", "triage", "stepdown")
    } else {
      stop("Unknown model name")
    }
  }
  model_partitions
}

# pull the fitted value from fitted_parameters.csv 
# this is calculated externally at the moment
read_fitted_parameters <- function(parameter_file = "extdata/fitted_parameters.csv") {
    fitted_parameter_file <- sircovid_file(parameter_file)
    
    # Read file with fitted parameters (from Bob Verity's hospital model)
    fitted_parameters <- read.csv(file = fitted_parameter_file)
    progression_groups <- list()

    s_hosp_D <- fitted_parameters[fitted_parameters$parameter=="s_hosp_D","value"]
    gamma_hosp_D <- fitted_parameters[fitted_parameters$parameter=="gamma_hosp_D","value"]
    s_hosp_R <- fitted_parameters[fitted_parameters$parameter=="s_hosp_R","value"]
    gamma_hosp_R <- fitted_parameters[fitted_parameters$parameter=="gamma_hosp_R","value"]
    s_ICU_D <- fitted_parameters[fitted_parameters$parameter=="s_ICU_D","value"]
    gamma_ICU_D <- fitted_parameters[fitted_parameters$parameter=="gamma_ICU_D","value"]
    s_ICU_R <- fitted_parameters[fitted_parameters$parameter=="s_ICU_R","value"]
    gamma_ICU_R <- fitted_parameters[fitted_parameters$parameter=="gamma_ICU_R","value"]
    s_triage <- fitted_parameters[fitted_parameters$parameter=="s_triage","value"]
    gamma_triage <- fitted_parameters[fitted_parameters$parameter=="gamma_triage","value"]
    s_stepdown <- fitted_parameters[fitted_parameters$parameter=="s_stepdown","value"]
    gamma_stepdown <- fitted_parameters[fitted_parameters$parameter=="gamma_stepdown","value"]

    parameters <- (list(progression_groups = list(hosp_D = s_hosp_D,
                                          hosp_R = s_hosp_R,
                                          ICU_D = s_ICU_D,
                                          ICU_R = s_ICU_R,
                                          triage = s_triage,
                                          stepdown = s_stepdown),
                gammas = list(hosp_D = gamma_hosp_D,
                              hosp_R = gamma_hosp_R,
                              ICU_D = gamma_ICU_D,
                              ICU_R = gamma_ICU_R,
                              triage = gamma_triage,
                              stepdown = gamma_stepdown)))
              
    parameters
  }

# Loads a model by its name
# Must be in inst/odin
load_odin_model <- function(x) {
  if (inherits(x, "sircovid_basic")) {
    return(x)
  }
  if (is.null(x) || x == "basic") {
    model <- basic
  } else {
    path <- system.file("odin", package = "sircovid", mustWork = TRUE)
    possible <- sub("\\.json$", "", dir(path, pattern = "\\.json$"))
    if (x %in% possible) {
      env <- asNamespace("sircovid")
      model <- get(x, envir = env, mode = "function", inherits = FALSE)
    } else {
      stop("Unknown model: ", x)
    }
  }
  model
}
