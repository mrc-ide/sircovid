library(socialmixr)
source("uk_latest_demography.R")
source("generate_parameters_nCov_severity_model.R")
discrete_structured_SIR_with_severity <- odin::odin("discrete_structured_SIR_with_severity.R")
survey.pop <- get_age_distribution()
pars_model <- generate_parameters_nCov_severity_model(
    beta = .042,
    survey.pop = survey.pop
)
mod <- discrete_structured_SIR_with_severity(user = pars_model)
