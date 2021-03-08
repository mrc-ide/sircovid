dput_str <- function(x, name, prefix, width) {
  txt <- trimws(deparse(x, width.cutoff = width))
  if (last(txt) == ")") {
    txt[length(txt) - 1] <- paste0(txt[length(txt) - 1], ")")
    txt <- txt[-length(txt)]
  }
  pad1 <- strrep(" ", prefix - nchar(name))
  txt[[1]] <- sprintf("%s%s = %s", name, pad1, txt[[1]])
  if (length(txt) > 1) {
    txt[-1] <- paste0(strrep(" ", prefix + 5), txt[-1])
  }
  paste0(strrep(" ", 10), txt, collapse = "\n")
}


dput_named_matrix <- function(m, width = 55) {
  prefix1 <- "   rbind("
  prefix2 <- "         "
  nms <- rownames(m)
  width <- width - max(nchar(nms))
  contents <- paste(vcapply(seq_len(nrow(m)), function(i)
    dput_str(m[i, ], nms[i], max(nchar(nms)), width)),
    collapse = ",\n")
  res <- sprintf("    rbind(%s)", trimws(contents, "left"))
  writeLines(res)
  invisible(res)
}


on_ci <- function() {
  isTRUE(as.logical(Sys.getenv("CI")))
}


on_windows <- function() {
  tolower(Sys.info()[["sysname"]]) == "windows"
}


on_mac <- function() {
  tolower(Sys.info()[["sysname"]]) == "darwin"
}


skip_on_windows_gha <- function() {
  if (on_ci() && on_windows()) {
    testthat::skip("On Windows Github Actions")
  }
}


skip_on_mac_gha <- function() {
  if (on_ci() && on_mac()) {
    testthat::skip("On macOS Github Actions")
  }
}


test_example_uptake <- function() {
  c(rep(0, 3), # no vaccination in <15
  (2 / 5) * 0.75, # no vaccination in 15-17yo
  rep(0.75, 6),
  rep(0.85, 6),
  0.95,
  0.85,
  0.95)
}


test_vaccine_schedule <- function(daily_doses = 20000, region = "london",
                                  n_days = 365,
                                  uptake = test_example_uptake(),
                                  mean_days_between_doses = 12 * 7) {
  daily_doses <- rep(daily_doses, n_days)
  n <- vaccine_priority_population(region, uptake)
  vaccine_schedule_future(
    0, daily_doses, mean_days_between_doses, n)
}


test_vaccine_data <- function() {
  age_start <- seq(15, 95, by = 5)
  data <- data_frame(
    date = rep(seq(as.Date("2021-03-05"), length.out = 25, by = 1),
               each = length(age_start)),
    age_band_min = age_start)
  data$dose1 <- rpois(nrow(data), 200)
  data$dose2 <- rpois(nrow(data), 100)
  data
}


expect_approx_equal <- function(x1, x2, rel_tol = 0.05) {
  x1_zeros <- x1 == 0
  x2_zeros <- x2 == 0
  expect_true(all(abs(x1[!x1_zeros] - x2[!x1_zeros]) / x1[!x1_zeros] < rel_tol))
  expect_true(all(abs(x1[x1_zeros & !x2_zeros] - x2[x1_zeros & !x2_zeros]) /
                    x2[x1_zeros & !x2_zeros] < rel_tol))
}
