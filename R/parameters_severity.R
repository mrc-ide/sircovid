sircovid_parameters_severity <- function(path) {
  ## Set up severity file into table
  path <- path %||% sircovid_file("extdata/severity_default.csv")
  params <- read_csv(path)

  ## Transpose so columns are parameters, rownames are age groups
  data <- t(as.matrix(params[-1L]))
  colnames(data) <- params[[1]]
  data <- cbind(age = rownames(data),
                data.frame(data, check.names = FALSE),
                stringsAsFactors = FALSE)
  rownames(data) <- NULL

  required <- c(
    population = "Size of England population",
    p_sympt_seek_hc = "Proportion of symptomatic cases seeking healthcare",
    p_sympt = "Proportion with symptoms",
    p_sympt_hosp = "Proportion of symptomatic cases hospitalised",
    p_ICU_hosp = "Proportion of hospitalised cases getting critical care",
    p_death_ICU = "Proportion of critical cases dying",
    p_death_hosp_D = "Proportion of non-critical care cases dying",
    p_seroconversion = "Proportion of cases that seroconvert",
    p_death_comm = "Proportion of severe cases dying in the community",
    p_admit_conf = "Proportion of hospitalised cases admitted as confirmed")
  data <- rename(data, required, names(required))

  # Parse the age bins. Useful to keep both start and end depending on what
  # function expects
  age_bins <- check_age_bins(data[["age"]])

  p_recov_ILI <- 1 - data[["p_sympt_hosp"]] / data[["p_sympt_seek_hc"]]

  list(
    p_asympt = 1 - data[["p_sympt"]],
    p_sympt_ILI = data[["p_sympt"]] * data[["p_sympt_seek_hc"]],
    p_recov_ICU = 1 - data[["p_death_ICU"]],
    p_recov_ILI = p_recov_ILI,
    p_hosp_ILI = 1 - p_recov_ILI,
    p_recov_hosp = (1 - data[["p_ICU_hosp"]]) * (1 - data[["p_death_hosp_D"]]),
    p_death_hosp = (1 - data[["p_ICU_hosp"]]) * data[["p_death_hosp_D"]],
    p_death_hosp_D = data[["p_death_hosp_D"]],
    p_death_ICU = data[["p_death_ICU"]],
    p_ICU_hosp = data[["p_ICU_hosp"]],
    p_seroconversion = data[["p_seroconversion"]],
    p_death_comm = data[["p_death_comm"]],
    p_admit_conf = data[["p_admit_conf"]])
}


check_age_bins <- function(age_headers) {
  bins <- sircovid_age_bins()
  expected <- sprintf("%d to %d", bins$start, bins$end)
  if (!identical(age_headers, expected)) {
    stop(sprintf("Incorrect age bands:\nexpected: %s\ngiven: %s",
                 expected, age_headers))
  }
  bins
}
