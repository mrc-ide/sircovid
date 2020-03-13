##' Create a model
##' @title Create a model
##' @export
sircovid <- function() {
  mod <- basic(user = pars_model)
  survey_pop <- age_distribution()
  pars_model <- generate_parameters_nCov_severity_model(
    beta = 0.042,
    survey_pop = survey_pop)
  discrete_structured_SIR_with_severity(user = pars_model)
}
