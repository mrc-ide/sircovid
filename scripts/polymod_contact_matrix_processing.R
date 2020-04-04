# Loading in Raw Per-Capita Contact Matrices Generated in polymod_contact_matrix_extraction.R
overall_matrix <- read_rds("scripts/overall_contact_matrix.rds")
home_matrix <- read_rds("scripts/home_contact_matrix.rds")
work_matrix <- read_rds("scripts/work_contact_matrix.rds")
school_matrix <- read_rds("scripts/school_contact_matrix.rds")
other_matrix <- read_rds("scripts/other_contact_matrix.rds")

# Get Age Distribution of the United Kingdom 
age_dist <- age_distribution(infile = NULL, country = "United Kingdom", 
                             year = 2019, age.limits = seq(0, 75, 5))
population <- age_dist$population

# Generating the Symmetric Version of these Raw Matrices
# This is equivalent to calling "contact_matrix" from socialmixr with argument symmetric = TRUE
balance_matrices <- function(contact_matrix, population) {
  if (ncol(contact_matrix) != length(population)) {
    stop("Number of age groups doesn't match number of columns in contact matrix")
  }
  processed_matrix <- matrix(nrow = nrow(contact_matrix), ncol = ncol(contact_matrix))
  for (i in 1:nrow(contact_matrix)) {
    for (j in 1:ncol(contact_matrix)) {
      pop_age_i <- population[i]
      pop_age_j <- population[j]
      processed_matrix[i, j] <- (1/(2 * pop_age_i)) * (contact_matrix[i, j] * pop_age_i + 
                                                         contact_matrix[j, i] * pop_age_j)
    }
  }
  return(processed_matrix)
}

# Transforming the contact matrices to the (symetrical) transmission matrix 
balanced_overall_matrix <- balance_matrices(overall_matrix, population)
trans_overall_matrix <- balanced_overall_matrix / rep(population, each = ncol(balanced_overall_matrix))
saveRDS(trans_overall_matrix, file = "inst/extdata/overall_trans_matrix.rds")

balanced_home_matrix <- balance_matrices(home_matrix, population)
trans_home_matrix <- balanced_home_matrix / rep(population, each = ncol(balanced_home_matrix))
saveRDS(trans_overall_matrix, file = "inst/extdata/overall_home_matrix.rds")

balanced_work_matrix <- balance_matrices(work_matrix, population)
trans_work_matrix <- balanced_work_matrix / rep(population, each = ncol(balanced_work_matrix))
saveRDS(trans_overall_matrix, file = "inst/extdata/overall_work_matrix.rds")

balanced_school_matrix <- balance_matrices(school_matrix, population)
trans_school_matrix <- balanced_school_matrix / rep(population, each = ncol(balanced_school_matrix))
saveRDS(trans_overall_matrix, file = "inst/extdata/overall_school_matrix.rds")

balanced_other_matrix <- balance_matrices(other_matrix, population)
trans_other_matrix <- balanced_other_matrix / rep(population, each = ncol(balanced_other_matrix))
saveRDS(trans_overall_matrix, file = "inst/extdata/overall_othermatrix.rds")

