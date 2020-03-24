default_age_distribution <- function() {
  if (is.null(cache$age_distribution)) {
    message("Computing age distribution")
    cache$age_distribution <- age_distribution()
  }
  cache$age_distribution
}

## source of age distribution:
## https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/EXCEL_FILES/5_Interpolated/WPP2019_INT_F03_1_POPULATION_BY_AGE_ANNUAL_BOTH_SEXES.xlsx
age_distribution <- function(infile = NULL,
                             country = "United Kingdom",
                             year = 2019,
                             age.limits = c(0, 10, 20, 30, 40, 50, 60, 70, 80)
                             ) {
    if (is.null(infile)) {
      infile <- sircovid_file("extdata/WPP2019_INT_F03_1_POPULATION_BY_AGE_ANNUAL_BOTH_SEXES.xlsx")
    }

    estimates <- readxl::read_xlsx(infile, sheet = "ESTIMATES", skip = 16)
    estimates <- estimates[estimates$`Region, subregion, country or area *` == country, ]
    estimates <- estimates[estimates[["Reference date (as of 1 July)"]] == year, ]
    pop_by_age <- dplyr::select(estimates,`0`:`98`)
    pop_by_age <-tidyr::gather(pop_by_age, key = "age", value = "pop", convert = TRUE)
    pop_by_age$pop <- as.numeric(pop_by_age$pop)

    age_bands_ll <- age.limits
    age_bands_ul <- age.limits[-1] - 1
    age_bands_ul <- append(age_bands_ul, max(pop_by_age$age))

    ## In the structure required by socialmixr
    out <- data.frame(
        lower.age.limit = age_bands_ll
    )

    out$population <- unlist(
        purrr::map2(
            age_bands_ll,
            age_bands_ul,
            function(ll, ul) sum(pop_by_age$pop[pop_by_age$age >= ll & pop_by_age$age <= ul]))
    ) * 1000

    out
}
