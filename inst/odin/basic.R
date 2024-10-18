## E and R stage indexed by i, j, k with
## i for the age group
## j for the progression (not exponential latent and infectious period)
## k for the infectivity group

## Number of age classes & number of transmissibility classes
n_age_groups <- parameter(type = "integer", constant = TRUE)
n_trans_classes <- parameter(1, type = "integer", constant = TRUE)

## Seeding of first wave; this will happen on the S->E flow
seed <- interpolate(seed_time, seed_value, "constant")
seed_age_band <- as.integer(4) # 15-19y band

seed_time <- parameter()
seed_value <- parameter()
dim(seed_time) <- parameter(rank = 1)
dim(seed_value) <- parameter(rank = 1)

## Core equations for transitions between compartments:
update(S[]) <- S[i] - n_SE[i]

update(E[, , ]) <- E[i, j, k] + delta_E[i, j, k]
update(I_A[, , ]) <- I_A[i, j, k] + delta_I_A[i, j, k]
update(I_C[, , ]) <- I_C[i, j, k] + delta_I_C[i, j, k]
update(R[]) <- R[i] + delta_R[i]
update(I_hosp[, , ]) <- I_hosp[i, j, k] + delta_I_hosp[i, j, k]
update(I_ICU[, , ]) <- new_I_ICU[i, j, k]
update(R_hosp[, , ]) <- R_hosp[i, j, k] + delta_R_hosp[i, j, k]

update(D[]) <- D[i] + new_D[i]

## Stuff we want to track in addition to number of individuals
## in each compartment.

## Individual probabilities of transition:
p_SE[] <- 1 - exp(-lambda[i] * dt) # S to I - age dependent
p_EE <- 1 - exp(-gamma_E * dt) # progression of latent period
p_II_A <- 1 - exp(-gamma_A * dt) # progression of infectious period
p_II_C <- 1 - exp(-gamma_C * dt)
p_II_hosp <- 1 - exp(-gamma_hosp * dt)
p_II_ICU <- 1 - exp(-gamma_ICU * dt)
p_R_hosp <- 1 - exp(-gamma_rec * dt)

## Draws from binomial distributions for numbers changing between
## compartments:
n_SE[] <- Binomial(S[i], p_SE[i])
n_SE[seed_age_band] <- min(S[i], n_SE[i] + Poisson(seed))

n_EE[, , ] <- Binomial(E[i, j, k], p_EE)
n_II_A[, , ] <- Binomial(I_A[i, j, k], p_II_A)
n_II_C[, , ] <- Binomial(I_C[i, j, k], p_II_C)
n_II_hosp[, , ] <- Binomial(I_hosp[i, j, k], p_II_hosp)
n_II_ICU[, , ] <- Binomial(I_ICU[i, j, k], p_II_ICU)
n_R_hosp[, , ] <- Binomial(R_hosp[i, j, k], p_R_hosp)

## Computes the number of asymptomatic
n_EI_A[, ] <- Binomial(n_EE[i, k_E, j], 1 - p_C[i])


## Computes the number of symptomatic cases
n_EI_C[, ] <- n_EE[i, k_E, j] - n_EI_A[i, j]

## Compute the aux_p_bin matrix of binom nested coeff
aux_p_bin[, 1] <- trans_profile[i, 1]
aux_p_bin[, 2:(n_trans_classes - 1)] <-
  trans_profile[i, j] / sum(trans_profile[i, j:n_trans_classes])

## Implementation of multinom via nested binomial
aux_EE[, 1, 1] <- Binomial(n_SE[i], aux_p_bin[i, 1])
aux_EE[, 1, 2:(n_trans_classes - 1)] <- Binomial(n_SE[i] - sum(aux_EE[i, 1, 1:(k - 1)]), aux_p_bin[i, k])
aux_EE[, 1, n_trans_classes] <-
  n_SE[i] - sum(aux_EE[i, 1, 1:(n_trans_classes - 1)])

## Work out the E->E transitions
aux_EE[, 2:k_E, ] <- n_EE[i, j - 1, k]
aux_EE[, 1:k_E, ] <- aux_EE[i, j, k] - n_EE[i, j, k]
delta_E[, , ] <- aux_EE[i, j, k]

## Work out the I_A->I_A transitions
aux_II_A[, 1, ] <- n_EI_A[i, k]
aux_II_A[, 2:k_A, ] <- n_II_A[i, j - 1, k]
aux_II_A[, 1:k_A, ] <- aux_II_A[i, j, k] - n_II_A[i, j, k]
delta_I_A[, , ] <- aux_II_A[i, j, k]

## Work out the I_C->I_C transitions
aux_II_C[, 1, ] <- n_EI_C[i, k]
aux_II_C[, 2:k_C, ] <- n_II_C[i, j - 1, k]
aux_II_C[, 1:k_C, ] <- aux_II_C[i, j, k] - n_II_C[i, j, k]
delta_I_C[, , ] <- aux_II_C[i, j, k]

## Work out the I_hosp->I_hosp transitions
n_sympt_to_hosp[, ] <- Binomial(n_II_C[i, k_C, j], 1 - p_recov_sympt[i])
aux_II_hosp[, 1, ] <- n_sympt_to_hosp[i, k]
aux_II_hosp[, 2:k_hosp, ] <- n_II_hosp[i, j - 1, k]
aux_II_hosp[, 1:k_hosp, ] <- aux_II_hosp[i, j, k] - n_II_hosp[i, j, k]
delta_I_hosp[, , ] <- aux_II_hosp[i, j, k]

## Work out the death in hospital without getting critical care
n_death_hosp[, ] <- Binomial(n_II_hosp[i, k_hosp, j], p_death_hosp[i])

## Work out the I_ICU -> I_ICU transitions
n_hosp_to_ICU[, ] <- Binomial(n_II_hosp[i, k_hosp, j] - n_death_hosp[i, j], 1 - p_recov_hosp[i] - p_death_hosp[i])
aux_II_ICU[, 1, ] <- n_hosp_to_ICU[i, k]
aux_II_ICU[, 2:k_ICU, ] <- n_II_ICU[i, j - 1, k]
aux_II_ICU[, 1:k_ICU, ] <- aux_II_ICU[i, j, k] - n_II_ICU[i, j, k]
delta_I_ICU[, , ] <- aux_II_ICU[i, j, k]

new_I_ICU[, , ] <- I_ICU[i, j, k] + delta_I_ICU[i, j, k]

## Work out the R_hosp -> R_hosp transitions
n_ICU_to_R_hosp[, ] <- Binomial(n_II_ICU[i, k_ICU, j], p_recov_ICU[i])
aux_R_hosp[, 1, ] <- n_ICU_to_R_hosp[i, k]
aux_R_hosp[, 2:k_rec, ] <- n_R_hosp[i, j - 1, k]
aux_R_hosp[, 1:k_rec, ] <- aux_R_hosp[i, j, k] - n_R_hosp[i, j, k]
delta_R_hosp[, , ] <- aux_R_hosp[i, j, k]

## Work out the number of deaths
delta_D[] <- sum(n_II_ICU[i, k_ICU, ]) - sum(n_ICU_to_R_hosp[i, ]) +
  sum(n_death_hosp[i, ])

new_D[] <- delta_D[i]

## Work out the number of recovery
delta_R[] <-
  sum(n_II_A[i, k_A, ]) +
  sum(n_II_C[i, k_C, ]) -
  sum(n_sympt_to_hosp[i, ]) +
  sum(n_II_hosp[i, k_hosp, ]) -
  sum(n_hosp_to_ICU[i, ]) -
  sum(n_death_hosp[i, ]) +
  sum(n_R_hosp[i, k_rec, ])

## Compute the force of infection
I_with_diff_trans[, ] <- trans_increase[i, j] * (
  sum(I_A[i, , j]) +
  sum(I_C[i, , j]) +
  hosp_transmission * sum(I_hosp[i, , j]) +
  ICU_transmission * sum(I_ICU[i, , j]))

s_ij[, ] <- m[i, j] * sum(I_with_diff_trans[j, ])
lambda[] <- beta * sum(s_ij[i, ])

## Initial states are all zerod as we will provide a state vector
## setting S and I based on the seeding model.
initial(S[]) <- 0
initial(E[, , ]) <- 0
initial(I_A[, , ]) <- 0
initial(I_C[, , ]) <- 0
initial(I_hosp[, , ]) <- 0
initial(I_ICU[, , ]) <- 0
initial(R_hosp[, , ]) <- 0
initial(R[]) <- 0
initial(D[]) <- 0
initial(N_tot) <- 0

## aggregate our reporting statistics
initial(I_ICU_tot) <- 0
update(I_ICU_tot) <- sum(new_I_ICU)

tot_new_D <- sum(new_D)
initial(D_tot) <- 0
update(D_tot) <- D_tot + tot_new_D
initial(D_inc, zero_every = 1) <- 0
update(D_inc) <- D_inc + tot_new_D


## User defined parameters - default in parentheses:

## Parameters of the E classes
k_E <- parameter(type = "integer", constant = TRUE)
gamma_E <- parameter(0.1)

## Probability of transitioning from the E to the asymptomatic
## class, the rest goes in the symptomatic class
p_C <- parameter()

## Parameters of the I_A classes
k_A <- parameter(type = "integer", constant = TRUE)
gamma_A <- parameter(0.1)

## Parameters of the I_C classes
k_C <- parameter(type = "integer", constant = TRUE)
gamma_C <- parameter(0.1)
p_recov_sympt <- parameter()

## Parameters of the I_hosp classes
k_hosp <- parameter(type = "integer", constant = TRUE)
gamma_hosp <- parameter(0.1)
p_recov_hosp <- parameter()
p_death_hosp <- parameter()

## Parameters of the I_ICU classes
k_ICU <- parameter(type = "integer", constant = TRUE)
gamma_ICU <- parameter(0.1)
p_recov_ICU <- parameter()

## Parameters of the R_hosp classes
k_rec <- parameter(type = "integer", constant = TRUE)
gamma_rec <- parameter(0.1)

## Parameters of the age stratified transmission
beta_time <- parameter()
beta_value <- parameter()
dim(beta_time) <- parameter(rank = 1)
dim(beta_value) <- parameter(rank = 1)
beta <- interpolate(beta_time, beta_value, "constant")

## Useful for debugging
initial(beta_out) <- beta_value[1]
update(beta_out) <- beta

m <- parameter()
## TODO: trans_profile and trans_increase can be removed as not used;
## this will required removing one layer off many variables so be
## careful.
trans_profile[, ] <- 1
trans_increase[, ] <- 1
hosp_transmission <- parameter()
ICU_transmission <- parameter()

## Dimensions of the different "vectors" here vectors stand for
## multi-dimensional arrays

## Vectors handling the S class
dim(S) <- n_age_groups

## Vectors handling the E class
dim(E) <- c(n_age_groups, k_E, n_trans_classes)
dim(aux_EE) <- c(n_age_groups, k_E, n_trans_classes)
dim(delta_E) <- c(n_age_groups, k_E, n_trans_classes)
dim(n_EE) <- c(n_age_groups, k_E, n_trans_classes)

## Vectors handling the I_A class
dim(I_A) <- c(n_age_groups, k_A, n_trans_classes)
dim(aux_II_A) <- c(n_age_groups, k_A, n_trans_classes)
dim(delta_I_A) <- c(n_age_groups, k_A, n_trans_classes)
dim(n_II_A) <- c(n_age_groups, k_A, n_trans_classes)

## Vectors handling the I_C class
dim(I_C) <- c(n_age_groups, k_C, n_trans_classes)
dim(aux_II_C) <- c(n_age_groups, k_C, n_trans_classes)
dim(delta_I_C) <- c(n_age_groups, k_C, n_trans_classes)
dim(n_II_C) <- c(n_age_groups, k_C, n_trans_classes)
dim(p_recov_sympt) <- c(n_age_groups)

## Vectors handling the I_hosp class
dim(I_hosp) <- c(n_age_groups, k_hosp, n_trans_classes)
dim(aux_II_hosp) <- c(n_age_groups, k_hosp, n_trans_classes)
dim(delta_I_hosp) <- c(n_age_groups, k_hosp, n_trans_classes)
dim(n_II_hosp) <- c(n_age_groups, k_hosp, n_trans_classes)
dim(p_recov_hosp) <- c(n_age_groups)
dim(n_death_hosp) <- c(n_age_groups, n_trans_classes)

## Vectors handling the I_ICU class
dim(I_ICU) <- c(n_age_groups, k_ICU, n_trans_classes)
dim(aux_II_ICU) <- c(n_age_groups, k_ICU, n_trans_classes)
dim(delta_I_ICU) <- c(n_age_groups, k_ICU, n_trans_classes)
dim(n_II_ICU) <- c(n_age_groups, k_ICU, n_trans_classes)
dim(p_recov_ICU) <- c(n_age_groups)
dim(new_I_ICU) <- c(n_age_groups, k_ICU, n_trans_classes)

## Vectors handling the R_hosp class
dim(R_hosp) <- c(n_age_groups, k_rec, n_trans_classes)
dim(aux_R_hosp) <- c(n_age_groups, k_rec, n_trans_classes)
dim(delta_R_hosp) <- c(n_age_groups, k_rec, n_trans_classes)
dim(n_R_hosp) <- c(n_age_groups, k_rec, n_trans_classes)

## Vectors handling the R class
dim(R) <- c(n_age_groups)
dim(delta_R) <- c(n_age_groups)

## Vectors handling the D class
dim(D) <- c(n_age_groups)
dim(delta_D) <- c(n_age_groups)
dim(new_D) <- c(n_age_groups)

## Vectors handling the S->E transition where infected are split
## between level of infectivity
dim(p_SE) <- n_age_groups
dim(n_SE) <- n_age_groups
dim(aux_p_bin) <- c(n_age_groups, n_trans_classes)

## Vectors handling the E->I transition where newly infectious cases
## are split between level of severity
dim(n_EI_A) <- c(n_age_groups, n_trans_classes)
dim(n_EI_C) <- c(n_age_groups, n_trans_classes)

## Vectors handling number of new hospitalisations, ICU admissions and
## recoveries in hospital
dim(n_hosp_to_ICU) <- c(n_age_groups, n_trans_classes)
dim(n_sympt_to_hosp) <- c(n_age_groups, n_trans_classes)
dim(n_ICU_to_R_hosp) <- c(n_age_groups, n_trans_classes)

## Vectors handling the severity profile
dim(p_C) <- c(n_age_groups)

## Vectors handling the potential death in hospital before critical care
dim(p_death_hosp) <- c(n_age_groups)

## Vectors handling the age specific heterogeneous transmission process
dim(lambda) <- n_age_groups
dim(s_ij) <- c(n_age_groups, n_age_groups)
dim(m) <- c(n_age_groups, n_age_groups)
dim(trans_profile) <- c(n_age_groups, n_trans_classes)
dim(trans_increase) <- c(n_age_groups, n_trans_classes)
dim(I_with_diff_trans) <- c(n_age_groups, n_trans_classes)

## Used for error checking - population should be constant
update(N_tot) <- sum(S) + sum(R) + sum(D) + sum(E) + sum(I_A) +
  sum(I_C) + sum(I_hosp) + sum(I_ICU) + sum(R_hosp)


## compare
exp_noise <- parameter()
phi_ICU <- parameter()
kappa_ICU <- parameter()
phi_death <- parameter()
kappa_death <- parameter()

icu <- data()
deaths <- data()

icu_with_noise <- phi_ICU * I_ICU_tot + Exponential(exp_noise)
icu ~ NegativeBinomial(kappa_ICU, mu = icu_with_noise)

deaths_with_noise <- phi_death * D_inc + Exponential(exp_noise)
deaths ~ NegativeBinomial(kappa_death, mu = deaths_with_noise)
