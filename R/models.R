##' Create a basic model 
##' 
##' @title Basic model
##'
##' @param use_fitted_parameters Override progression_groups and gammas with fitted
##'   parameters loaded by \code{read_fitted_parameters()}
##'
##' @param progression_groups List of number of progression groups in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp', 'ICU', 'rec'
##'   
##' @param gammas List of exponential distribution rates for time in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp', 'ICU', 'rec'
##'   
##' @export
basic_model <- function(use_fitted_parameters = FALSE,
                        progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 1, hosp = 2, ICU = 2, rec = 2),
                        gammas = list(E = 1/(4.59/2), asympt = 1/2.09, mild = 1/2.09, ILI = 1/4, hosp = 2, ICU = 2/5, rec = 2/5)) {
  basic_model <- model_constructor("sircovid_basic", "basic", 
                                   use_fitted_parameters, progression_groups, gammas)
  basic_model
}

##' Create a hospital model
##' 
##' @title Hospital model
##' 
##' @param use_fitted_parameters Override progression_groups and gammas with fitted
##'   parameters loaded by \code{read_fitted_parameters()}
##' 
##' @param progression_groups List of number of progression groups in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp_D', 'hosp_R', 'ICU_D', 'ICU_R', 'triage', 'stepdown'
##'   
##' @param gammas List of exponential distribution rates for time in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp_D', 'hosp_R', 'ICU_D', 'ICU_R', 'triage', 'stepdown'
##'   
##' @export
hospital_model <- function(use_fitted_parameters = TRUE,
                           progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 1, hosp_D = 2 , hosp_R = 2, ICU_D = 2, ICU_R = 2, triage = 2, stepdown = 2),
                           gammas = list(E = 1/(4.59/2), asympt = 1/2.09, mild = 1/2.09, ILI = 1/4, hosp_D = 2/5, hosp_R = 2/10, ICU_D = 2/5, ICU_R = 2/10, triage = 2, stepdown = 2/5)) {
  hospital_model <- model_constructor("sircovid_hospital", "new_hospital_model", 
                                      use_fitted_parameters, progression_groups, gammas)
  
  hospital_model
}

##' Create a serology model
##' 
##' @title Serology model
##' 
##' @param use_fitted_parameters Override progression_groups and gammas with fitted
##'   parameters loaded by \code{read_fitted_parameters()}
##' 
##' @param progression_groups List of number of progression groups in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp_D', 'hosp_R', 'ICU_D', 'ICU_R', 'triage', 'stepdown', 'R_pre'
##'   
##' @param gammas List of exponential distribution rates for time in each partition
##'   needs 'E', 'asympt', 'mild', 'ILI', 'hosp_D', 'hosp_R', 'ICU_D', 'ICU_R', 'triage', 'stepdown', 'R_pre', 'test'
##'   
##' @export
serology_model <- function(use_fitted_parameters = TRUE,
                           progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 1, comm_D =2, hosp_D = 2 , hosp_R = 2, ICU_D = 2, ICU_R = 2, triage = 2, stepdown = 2, R_pre = 2),
                           gammas = list(E = 1/(4.59/2), asympt = 1/2.09, mild = 1/2.09, ILI = 1/4, comm_D = 2/5, hosp_D = 2/5, hosp_R = 2/10, ICU_D = 2/5, ICU_R = 2/10, triage = 2, stepdown = 2/5, R_pre = 1/5, test = 3/10)) {
  model_class <- "sircovid_serology" 
  serology_model <- model_constructor(model_class, "hospital_with_serology_testing", 
                                      use_fitted_parameters, progression_groups, gammas)
  
  # This inherits from the sircovid_hospital model, it only adds new partitions/parameters
  class(serology_model) <- c(model_class, "sircovid_hospital")
  serology_model
}


#
# Internal functions
#

# Base constructor that can be overridden by child classes
model_constructor <- function(model_class,
                              odin_model_name,
                              use_fitted_parameters,
                              progression_groups,
                              gammas) {
  odin_model <- load_odin_model(odin_model_name)
  generate_beta_func <- generate_beta
  seeding_model <- seeding_function
  compare_model <- function(model, pars_obs, data) {compare_output(model, pars_obs, data, type=model_class)}

  # Overwrite with fitted parameters if loaded
  if (use_fitted_parameters) {
    fitted_parameters <- read_fitted_parameters()
    progression_groups[names(fitted_parameters$progression_groups)] <- fitted_parameters$progression_groups 
    gammas[names(fitted_parameters$gammas)] <- fitted_parameters$gammas 
  }

  model_partitions <- partition_names(model_class)
  if (any(!(model_partitions %in% names(progression_groups)))) {
    stop("progression_groups need to be defined for all partitions")
  } else {
    progression_groups <- progression_groups[model_partitions] 
  }
  if (any(!(model_partitions %in% names(gammas)))) {
    stop("gammas need to be defined for all partitions")
  } else {
    if(model_class == "sircovid_serology"){
      if (!("test" %in% names(gammas))){
        stop("gamma needs to be defined for test")
      } else {
        gammas <- gammas[c(model_partitions,"test")] 
      }
    } else {
      gammas <- gammas[model_partitions] 
    }
  }
  
  model_object <- list(odin_model = odin_model,
                generate_beta_func = generate_beta_func,
                progression_groups = progression_groups,
                gammas = gammas,
                seeding_model = seeding_model,
                compare_model = compare_model)
  class(model_object) <- (model_class)
  model_object
}

# Definitions of partition names
partition_names <- function(model_name) {
  if (inherits(model_name, "sircovid_basic") || inherits(model_name, "sircovid_hospital")) {
    model_partitions <- names(model_name$progression_groups)
  } else {
    if (model_name == "sircovid_basic") {
      model_partitions <- c("E", "asympt", "mild", "ILI", "hosp", "ICU", "rec")
    } else if (model_name == "sircovid_hospital") {
      model_partitions <- c("E", "asympt", "mild", "ILI", "hosp_D", "hosp_R", "ICU_D", "ICU_R", "triage", "stepdown")
    } else if (model_name == "sircovid_serology") {
      model_partitions <- c("E", "asympt", "mild", "ILI", "comm_D", "hosp_D", "hosp_R", "ICU_D", "ICU_R", "triage", "stepdown", "R_pre")
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
    s_R_pre <- fitted_parameters[fitted_parameters$parameter=="s_R_pre","value"]
    gamma_R_pre <- fitted_parameters[fitted_parameters$parameter=="gamma_R_pre","value"]
    s_comm_D <- fitted_parameters[fitted_parameters$parameter=="s_comm_D","value"]
    gamma_comm_D <- fitted_parameters[fitted_parameters$parameter=="gamma_comm_D","value"]
    gamma_test <- fitted_parameters[fitted_parameters$parameter=="gamma_test","value"]

    parameters <- (list(progression_groups = list(hosp_D = s_hosp_D,
                                          hosp_R = s_hosp_R,
                                          ICU_D = s_ICU_D,
                                          ICU_R = s_ICU_R,
                                          triage = s_triage,
                                          stepdown = s_stepdown,
                                          R_pre = s_R_pre,
                                          comm_D = s_comm_D),
                        gammas = list(hosp_D = gamma_hosp_D,
                                      hosp_R = gamma_hosp_R,
                                      ICU_D = gamma_ICU_D,
                                      ICU_R = gamma_ICU_R,
                                      triage = gamma_triage,
                                      stepdown = gamma_stepdown,
                                      R_pre = gamma_R_pre,
                                      comm_D = gamma_comm_D,
                                      test = gamma_test)))
              
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
