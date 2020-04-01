context("grid_search")

## This is not a real test but simply tries to run the model
test_that("Small grid search works") {
  set.seed(1)

  data <- read.csv("combin_time_series.csv")
  mortality_data <- readRDS("mortaility_data.Rds")
  
  names(data)[names(data) == "number_of_confirmed_covid_19_patients_in_hdu_itu_at_0800_total"] <- "itu"
  data$date <- as.Date(data$date,format="%Y-%m-%d")
  data <- subset(data,date >= as.Date("2020-03-01"), select=c(date, itu))
  data <- merge(data, mortality_data, by.x="date", by.y="dateRep", all.x=TRUE, all.y=TRUE)

  #ignoring ITU on 15/03
  data$itu[which(data$date == as.Date("2020-03-15"))] <- NA
 
  scan_results = scan_beta_date(
    min_beta = 0.1,
    max_beta = 0.2,
    beta_step = 0.05,
    start_date = "2020-01-21", 
    end_date = "2020-01-22", 
    day_step = 1,
    data = data)
   
}