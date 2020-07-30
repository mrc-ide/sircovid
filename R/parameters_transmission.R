sircovid_transmission_matrix <- function() {
  age_bins <- sircovid_age_bins()
  max_polymod_age <- 70

  # Get the contact matrix from socialmixr; polymod only
  # goes up to age 70
  contact <- suppressMessages(socialmixr::contact_matrix(
    socialmixr::polymod,
    countries = "United Kingdom",
    age.limits = age_bins$start[age_bins$start <= max_polymod_age],
    symmetric = TRUE))

  # transform the matrix to the (symetrical) transmission matrix
  # rather than the contact matrix
  transmission <- contact$matrix /
    rep(contact$demography$population, each = ncol(contact$matrix))

  # POLYMOD has a max age of 70, so older bins need to be filled in
  # assumes that the probability of contact remains as in POLYMOD
  # and that contacts are the same in 70+ and 80+
  extra <- age_bins$start[age_bins$start > max_polymod_age]
  i <- c(seq_len(nrow(transmission)), rep(nrow(transmission), length(extra)))
  transmission <- transmission[i, i]

  # Rename for cosmetic reasons only
  nms <- sprintf("[%d,%d)", age_bins$start, age_bins$end)
  dimnames(transmission) <- list(nms, nms)

  transmission
}
