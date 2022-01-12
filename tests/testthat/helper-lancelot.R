helper_lancelot_data <- function(data, start_date, dt) {
  expected <- c(deaths_hosp = NA_real_, deaths_comm = NA_real_,
                deaths_carehomes = NA_real_, deaths_non_hosp = NA_real_,
                deaths_hosp_0_49 = NA_real_, deaths_hosp_50_54 = NA_real_,
                deaths_hosp_55_59 = NA_real_, deaths_hosp_60_64 = NA_real_,
                deaths_hosp_65_69 = NA_real_, deaths_hosp_70_74 = NA_real_,
                deaths_hosp_75_79 = NA_real_, deaths_hosp_80_plus = NA_real_,
                icu = NA_real_, general = NA_real_, hosp = NA_real_,
                deaths = NA_real_, admitted = NA_real_, diagnoses = NA_real_,
                all_admission = NA_real_, sero_pos_15_64_1 = NA_real_,
                sero_tot_15_64_1 = NA_real_, sero_pos_15_64_2 = NA_real_,
                sero_tot_15_64_2 = NA_real_, pillar2_tot = NA_real_,
                pillar2_pos = NA_real_, pillar2_cases = NA_real_,
                pillar2_over25_tot = NA_real_, pillar2_under15_tot = NA_real_,
                pillar2_15_24_tot = NA_real_, pillar2_25_49_tot = NA_real_,
                pillar2_50_64_tot = NA_real_, pillar2_65_79_tot = NA_real_,
                pillar2_80_plus_tot = NA_real_, pillar2_over25_pos = NA_real_,
                pillar2_under15_pos = NA_real_, pillar2_15_24_pos = NA_real_,
                pillar2_25_49_pos = NA_real_, pillar2_50_64_pos = NA_real_,
                pillar2_65_79_pos = NA_real_, pillar2_80_plus_pos = NA_real_,
                pillar2_over25_cases = NA_real_,
                pillar2_under15_cases = NA_real_,
                pillar2_15_24_cases = NA_real_, pillar2_25_49_cases = NA_real_,
                pillar2_50_64_cases = NA_real_, pillar2_65_79_cases = NA_real_,
                pillar2_80_plus_cases = NA_real_, react_pos = NA_real_,
                react_tot = NA_real_, strain_non_variant = NA_real_,
                strain_tot = NA_real_, strain_over25_non_variant = NA_real_,
                strain_over25_tot = NA_real_)
  data <- helper_sircovid_data(data, start_date, dt, expected)
  lancelot_check_data(data)
  data
}


helper_lancelot_particle_filter <- function(data, n_particles,
                                            n_threads = 1L, seed = NULL,
                                            compiled_compare = FALSE) {

  fit_index <- function(info) {
    lancelot_index(info, rt = TRUE)
  }

  lancelot_check_data(data)
  mcstate::particle_filter$new(
    data,
    lancelot,
    n_particles,
    if (compiled_compare) NULL else lancelot_compare,
    fit_index,
    lancelot_initial,
    n_threads,
    seed)
}


## similar to above but really simple; load a data set and add a bunch
## of empty data where we require it.
lancelot_simple_data <- function(data) {
  extra <- c("deaths_hosp",
             "deaths_hosp_0_49", "deaths_hosp_50_54",
             "deaths_hosp_55_59", "deaths_hosp_60_64",
             "deaths_hosp_65_69", "deaths_hosp_70_74",
             "deaths_hosp_75_79", "deaths_hosp_80_plus",
             "deaths_carehomes", "deaths_comm", "deaths_non_hosp",
             "general", "hosp", "admitted", "diagnoses",
             "all_admission", "sero_pos_15_64_1", "sero_tot_15_64_1",
             "sero_pos_15_64_2", "sero_tot_15_64_2", "pillar2_pos",
             "pillar2_tot", "pillar2_cases", "pillar2_over25_pos",
             "pillar2_over25_tot", "pillar2_over25_cases",
             "react_pos", "react_tot", "strain_non_variant",
             "strain_tot", "strain_over25_non_variant",
             "strain_over25_tot", "pillar2_under15_cases",
             "pillar2_15_24_cases", "pillar2_25_49_cases",
             "pillar2_50_64_cases", "pillar2_65_79_cases",
             "pillar2_80_plus_cases", "pillar2_under15_tot",
             "pillar2_15_24_tot", "pillar2_25_49_tot",
             "pillar2_50_64_tot", "pillar2_65_79_tot",
             "pillar2_80_plus_tot", "pillar2_under15_pos",
             "pillar2_15_24_pos", "pillar2_25_49_pos",
             "pillar2_50_64_pos", "pillar2_65_79_pos",
             "pillar2_80_plus_pos")
  for (v in setdiff(extra, names(data))) {
    data[[v]] <- NA_real_
  }
  data
}
