##' Create a model
##' @title Create a model
##' @export
##' @param params List of parameters of the model
sircovid <- function(params) {
  survey_pop <- default_age_distribution()
  ## pars_model <- generate_parameters(
  ##   beta = rep(0.042,3))
  ##basic(user = pars_model)
  basic(user = params)
}
