##' Create a model
##' @title Create a model
##' @export
sircovid <- function() {
  survey_pop <- default_age_distribution()
  pars_model <- parameters(
    beta = 0.042,
    survey_pop = survey_pop)
  basic(user = pars_model)
}
