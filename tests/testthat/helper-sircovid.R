if (!file.exists("survey_pop.csv")) {
  survey_pop <- default_age_distribution()
  write.csv(survey_pop, "survey_pop.csv", row.names = FALSE)
}

survey_pop <- read.csv("survey_pop.csv")
