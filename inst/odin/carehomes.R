## E and R stage indexed by i, j, k, l with
## i for the age group
## j for the strain
## k for the progression (not exponential latent and infectious period)
## l for the vacc. group

## Parameter to turn on or off the serology flows

model_pcr_and_serology_user <- user(1)
model_pcr_and_serology <- if (model_pcr_and_serology_user == 1) 1 else 0

## Number of classes (age & vaccination)

## Number of "groups", being the age classes, Carehome workers and
## Carehome residents. This will be 19 in all but experimental uses.
n_age_groups <- user()
n_groups <- user()

## Definition of the time-step and output as "time"
steps_per_day <- user(integer = TRUE)
dt <- 1 / steps_per_day
initial(time) <- 0
update(time) <- (step + 1) * dt

## output number of individuals vaccinated by age and vaccine stage
## For example, for E, we sum over n_E_next_vacc_class (those moving vaccine
## stage without progressing disease stages) and also n_EE_next_vacc_class
## (those moving vaccine stage and also progressing disease stages)
## vaccinated S
n_S_vaccinated[, ] <-
  n_S_next_vacc_class[i, j] + sum(n_SE_next_vacc_class[i, , j])
dim(n_S_vaccinated) <- c(n_groups, n_vacc_classes)
n_E_vaccinated[, ] <-
  sum(n_E_next_vacc_class[i, , , j]) + sum(n_EE_next_vacc_class[i, , , j])
dim(n_E_vaccinated) <- c(n_groups, n_vacc_classes)
n_I_A_vaccinated[, ] <-
  sum(n_I_A_next_vacc_class[i, , , j]) +
  sum(n_II_A_next_vacc_class[i, , , j])
dim(n_I_A_vaccinated) <- c(n_groups, n_vacc_classes)
n_I_P_vaccinated[, ] <-
  sum(n_I_P_next_vacc_class[i, , , j]) +
  sum(n_II_P_next_vacc_class[i, , , j])
dim(n_I_P_vaccinated) <- c(n_groups, n_vacc_classes)
n_R_vaccinated[, ] <-
  sum(n_R_next_vacc_class[i, , j]) + sum(n_RS_next_vacc_class[i, , j])
dim(n_R_vaccinated) <- c(n_groups, n_vacc_classes)

initial(cum_n_S_vaccinated[, ]) <- 0
update(cum_n_S_vaccinated[, ]) <-
  cum_n_S_vaccinated[i, j] + n_S_vaccinated[i, j]
dim(cum_n_S_vaccinated) <- c(n_groups, n_vacc_classes)
initial(cum_n_E_vaccinated[, ]) <- 0
update(cum_n_E_vaccinated[, ]) <-
  cum_n_E_vaccinated[i, j] + n_E_vaccinated[i, j]
dim(cum_n_E_vaccinated) <- c(n_groups, n_vacc_classes)
initial(cum_n_I_A_vaccinated[, ]) <- 0
update(cum_n_I_A_vaccinated[, ]) <-
  cum_n_I_A_vaccinated[i, j] + n_I_A_vaccinated[i, j]
dim(cum_n_I_A_vaccinated) <- c(n_groups, n_vacc_classes)
initial(cum_n_I_P_vaccinated[, ]) <- 0
update(cum_n_I_P_vaccinated[, ]) <-
  cum_n_I_P_vaccinated[i, j] + n_I_P_vaccinated[i, j]
dim(cum_n_I_P_vaccinated) <- c(n_groups, n_vacc_classes)
initial(cum_n_R_vaccinated[, ]) <- 0
update(cum_n_R_vaccinated[, ]) <-
  cum_n_R_vaccinated[i, j] + n_R_vaccinated[i, j]
dim(cum_n_R_vaccinated) <- c(n_groups, n_vacc_classes)

n_vaccinated[, ] <-
  n_S_vaccinated[i, j] +
  n_E_vaccinated[i, j] +
  n_I_A_vaccinated[i, j] +
  n_I_P_vaccinated[i, j] +
  n_R_vaccinated[i, j]
dim(n_vaccinated) <- c(n_groups, n_vacc_classes)

## Total number of vaccinations over S, E, I_asypmt, R for convenience
initial(cum_n_vaccinated[, ]) <- 0
update(cum_n_vaccinated[, ]) <- cum_n_vaccinated[i, j] + n_vaccinated[i, j]
dim(cum_n_vaccinated) <- c(n_groups, n_vacc_classes)

## Core equations for transitions between compartments:
update(S[, ]) <- new_S[i, j]
update(E[, , , ]) <- new_E[i, j, k, l]
update(I_A[, , , ]) <- new_I_A[i, j, k, l]
update(I_P[, , , ]) <- new_I_P[i, j, k, l]
update(I_C_1[, , , ]) <- new_I_C_1[i, j, k, l]
update(I_C_2[, , , ]) <- new_I_C_2[i, j, k, l]
update(G_D[, , , ]) <- new_G_D[i, j, k, l]
update(ICU_pre_unconf[, , , ]) <- new_ICU_pre_unconf[i, j, k, l]
update(ICU_pre_conf[, , , ]) <- new_ICU_pre_conf[i, j, k, l]
update(H_R_unconf[, , , ]) <- new_H_R_unconf[i, j, k, l]
update(H_R_conf[, , , ]) <- new_H_R_conf[i, j, k, l]
update(H_D_unconf[, , , ]) <- new_H_D_unconf[i, j, k, l]
update(H_D_conf[, , , ]) <- new_H_D_conf[i, j, k, l]
update(ICU_W_R_unconf[, , , ]) <- new_ICU_W_R_unconf[i, j, k, l]
update(ICU_W_R_conf[, , , ]) <- new_ICU_W_R_conf[i, j, k, l]
update(ICU_W_D_unconf[, , , ]) <- new_ICU_W_D_unconf[i, j, k, l]
update(ICU_W_D_conf[, , , ]) <- new_ICU_W_D_conf[i, j, k, l]
update(ICU_D_unconf[, , , ]) <- new_ICU_D_unconf[i, j, k, l]
update(ICU_D_conf[, , , ]) <- new_ICU_D_conf[i, j, k, l]
update(W_R_unconf[, , , ]) <- new_W_R_unconf[i, j, k, l]
update(W_R_conf[, , , ]) <- new_W_R_conf[i, j, k, l]
update(W_D_unconf[, , , ]) <- new_W_D_unconf[i, j, k, l]
update(W_D_conf[, , , ]) <- new_W_D_conf[i, j, k, l]
update(T_sero_pre[, , , ]) <- new_T_sero_pre[i, j, k, l]
update(T_sero_pos[, , , ]) <- new_T_sero_pos[i, j, k, l]
update(T_sero_neg[, , ]) <- new_T_sero_neg[i, j, k]
update(R[, , ]) <- new_R[i, j, k]
update(D_hosp[, ]) <- D_hosp[i, j] + delta_D_hosp[i, j]
update(D_non_hosp[]) <- D_non_hosp[i] + delta_D_non_hosp[i]
update(T_PCR_pre[, , , ]) <- new_T_PCR_pre[i, j, k, l]
update(T_PCR_pos[, , , ]) <- new_T_PCR_pos[i, j, k, l]
update(T_PCR_neg[, , ]) <- new_T_PCR_neg[i, j, k]

delta_admit_conf <-
  sum(n_I_C_2_to_H_D_conf) +
  sum(n_I_C_2_to_H_R_conf) +
  sum(n_I_C_2_to_ICU_pre_conf)
update(cum_admit_conf) <- cum_admit_conf + delta_admit_conf

delta_new_conf <-
  sum(n_H_D_unconf_to_conf) +
  sum(n_H_R_unconf_to_conf) +
  sum(n_ICU_pre_unconf_to_conf) +
  sum(n_ICU_D_unconf_to_conf) +
  sum(n_ICU_W_R_unconf_to_conf) +
  sum(n_ICU_W_D_unconf_to_conf) +
  sum(n_W_R_unconf_to_conf) +
  sum(n_W_D_unconf_to_conf)
update(cum_new_conf) <- cum_new_conf + delta_new_conf

initial(diagnoses_admitted[, ]) <- 0
update(diagnoses_admitted[, ]) <- diagnoses_admitted[i, j] +
  sum(n_I_C_2_to_H_D_conf[i, , j]) +
  sum(n_I_C_2_to_H_R_conf[i, , j]) +
  sum(n_I_C_2_to_ICU_pre_conf[i, , j]) +
  sum(n_H_D_unconf_to_conf[i, , , j]) +
  sum(n_H_R_unconf_to_conf[i, , , j]) +
  sum(n_ICU_pre_unconf_to_conf[i, , , j]) +
  sum(n_ICU_D_unconf_to_conf[i, , , j]) +
  sum(n_ICU_W_R_unconf_to_conf[i, , , j]) +
  sum(n_ICU_W_D_unconf_to_conf[i, , , j]) +
  sum(n_W_R_unconf_to_conf[i, , , j]) +
  sum(n_W_D_unconf_to_conf[i, , , j])
dim(diagnoses_admitted) <- c(n_groups, n_vacc_classes)

initial(admit_conf_inc) <- 0
update(admit_conf_inc) <- if (step %% steps_per_day == 0)
  delta_admit_conf else
    admit_conf_inc + delta_admit_conf

initial(new_conf_inc) <- 0
update(new_conf_inc) <- if (step %% steps_per_day == 0)
  delta_new_conf else new_conf_inc + delta_new_conf

update(cum_admit_by_age[]) <- cum_admit_by_age[i] + sum(n_I_C_2_to_hosp[i, , ])

## Individual probabilities of transition:

## vaccination
p_S_next_vacc_class[, ] <- vaccine_probability[i, j]
p_E_next_vacc_class[, , , ] <- vaccine_probability[i, l]
p_I_A_next_vacc_class[, , , ] <- vaccine_probability[i, l]
p_I_P_next_vacc_class[, , , ] <- vaccine_probability[i, l]
p_R_next_vacc_class[, , ] <- vaccine_probability[i, k]

## clinical progression
p_SE[, ] <- 1 - exp(-sum(lambda[i, ]) *
                      rel_susceptibility[i, j] * dt) # S to I age/vacc dependent
p_E_progress <- 1 - exp(-gamma_E * dt) # progression of latent period
p_I_A_progress[] <- 1 - exp(-gamma_A[i] * dt) # progression of infectious period
p_I_P_progress[] <- 1 - exp(-gamma_P[i] * dt)
p_I_C_1_progress[] <- 1 - exp(-gamma_C_1[i] * dt)
p_I_C_2_progress[] <- 1 - exp(-gamma_C_2[i] * dt)
p_G_D_progress <- 1 - exp(-gamma_G_D * dt)
p_ICU_pre_progress <- 1 - exp(-gamma_ICU_pre * dt)
p_H_R_progress <- 1 - exp(-gamma_H_R * dt)
p_H_D_progress <- 1 - exp(-gamma_H_D * dt)
p_ICU_W_R_progress <- 1 - exp(-gamma_ICU_W_R * dt)
p_ICU_W_D_progress <- 1 - exp(-gamma_ICU_W_D * dt)
p_ICU_D_progress <- 1 - exp(-gamma_ICU_D * dt)
p_W_R_progress <- 1 - exp(-gamma_W_R * dt)
p_W_D_progress <- 1 - exp(-gamma_W_D * dt)
p_T_sero_pre_progress[, , , ] <- 1 - exp(-gamma_sero_pre[k] * dt)
p_T_sero_pos_progress <- 1 - exp(-gamma_sero_pos * dt)
p_test <- 1 - exp(-gamma_U * dt)
p_T_PCR_pre_progress <- 1 - exp(-gamma_PCR_pre * dt)
p_T_PCR_pos_progress <- 1 - exp(-gamma_PCR_pos * dt)
p_RS[] <- 1 - exp(-waning_rate[i] * dt) # R to S age dependent

dim(p_I_A_progress) <- n_strains
dim(p_I_P_progress) <- n_strains
dim(p_I_C_1_progress) <- n_strains
dim(p_I_C_2_progress) <- n_strains

## Work out time-varying probabilities
p_ICU <- if (as.integer(step) >= length(p_ICU_step))
  p_ICU_step[length(p_ICU_step)] else p_ICU_step[step + 1]
p_ICU_by_age[] <- p_ICU * psi_ICU[i]

p_H <- if (as.integer(step) >= length(p_H_step))
  p_H_step[length(p_H_step)] else p_H_step[step + 1]
p_H_by_age[] <- p_H * psi_H[i]

p_ICU_D[] <- if (as.integer(step) >= dim(p_ICU_D_step, 1))
  p_ICU_D_step[dim(p_ICU_D_step, 1), i] else p_ICU_D_step[step + 1, i]
p_ICU_D_by_age[, ] <- p_ICU_D[j] * psi_ICU_D[i, j]

p_H_D[] <- if (as.integer(step) >= dim(p_H_D_step, 1))
  p_H_D_step[dim(p_H_D_step, 1), i] else p_H_D_step[step + 1, i]
p_H_D_by_age[, ] <- p_H_D[j] * psi_H_D[i, j]

p_W_D[] <- if (as.integer(step) >= dim(p_W_D_step, 1))
  p_W_D_step[dim(p_W_D_step, 1), i] else p_W_D_step[step + 1, i]
p_W_D_by_age[, ] <- p_W_D[j] * psi_W_D[i, j]

p_G_D[] <- if (as.integer(step) >= dim(p_G_D_step, 1))
  p_G_D_step[dim(p_G_D_step, 1), i] else p_G_D_step[step + 1, i]
p_G_D_by_age[, ] <- p_G_D[j] * psi_G_D[i, j]

p_star <- if (as.integer(step) >= length(p_star_step))
  p_star_step[length(p_star_step)] else p_star_step[step + 1]
p_star_by_age[] <- p_star * psi_star[i]

## Work out time-varying gammas for hospital durations
gamma_H_R <- if (as.integer(step) >= length(gamma_H_R_step))
  gamma_H_R_step[length(gamma_H_R_step)] else gamma_H_R_step[step + 1]

gamma_W_R <- if (as.integer(step) >= length(gamma_W_R_step))
  gamma_W_R_step[length(gamma_W_R_step)] else gamma_W_R_step[step + 1]

gamma_ICU_W_R <- if (as.integer(step) >= length(gamma_ICU_W_R_step))
  gamma_ICU_W_R_step[length(gamma_ICU_W_R_step)] else
    gamma_ICU_W_R_step[step + 1]

gamma_H_D <- if (as.integer(step) >= length(gamma_H_D_step))
  gamma_H_D_step[length(gamma_H_D_step)] else gamma_H_D_step[step + 1]

gamma_W_D <- if (as.integer(step) >= length(gamma_W_D_step))
  gamma_W_D_step[length(gamma_W_D_step)] else gamma_W_D_step[step + 1]

gamma_ICU_W_D <- if (as.integer(step) >= length(gamma_ICU_W_D_step))
  gamma_ICU_W_D_step[length(gamma_ICU_W_D_step)] else
    gamma_ICU_W_D_step[step + 1]

gamma_ICU_D <- if (as.integer(step) >= length(gamma_ICU_D_step))
  gamma_ICU_D_step[length(gamma_ICU_D_step)] else gamma_ICU_D_step[step + 1]

gamma_ICU_pre <- if (as.integer(step) >= length(gamma_ICU_pre_step))
  gamma_ICU_pre_step[length(gamma_ICU_pre_step)] else
    gamma_ICU_pre_step[step + 1]

## Draws from binomial distributions for numbers changing between
## compartments:

## modelling infections and vaccine progression, which can happen simultaneously

#### flow out of S ####

## new infections

## Compute the new infections with multiple strains using nested binomials
n_S_progress_tot[, ] <- rbinom(S[i, j], p_SE[i, j])
n_S_progress[, , ] <- if (j == 1)
  rbinom(n_S_progress_tot[i, k], lambda[i, 1] / sum(lambda[i, ])) else
    rbinom(n_S_progress_tot[i, k] - sum(n_S_progress[i, 1:(j - 1), k]),
           lambda[i, j] / sum(lambda[i, j:n_strains]))

## Introduction of new strains. n_S_progress is arranged as:
##
## [age, vaccine stage, strain infected with]
##
## As in the model initialisation we will use the teenager category,
## and only infect *unvaccinated* people. For now we will model only
## movement into the second compartment as that represents our "new"
## strain.
strain_seed_step[] <- user()
dim(strain_seed_step) <- user()
strain_rate <- (if (as.integer(step) >= length(strain_seed_step))
  strain_seed_step[length(strain_seed_step)]
  else strain_seed_step[step + 1])
strain_seed <- rpois(strain_rate)
## We must never try to move more individuals from this S category
## than are available, so need to do this with a min()
##
## NOTE: We *must* use the range 2:n_strains here even though only one
## strain variant is allowed exist, otherwise the generated code leads
## us to write out-of-bounds when running with a single strain.
n_S_progress[4, 2:n_strains, 1] <-
  min(n_S_progress[i, j, k] + strain_seed,
      n_S_progress[i, j, k] + S[i, k] - sum(n_S_progress[i, , k]))

## of those some can also be vaccinated or progress through vaccination classes
## --> number transitioning from S[k] to E[k+1] (k vaccination class)
n_SE_next_vacc_class[, , ] <-
  rbinom(n_S_progress[i, j, k], p_S_next_vacc_class[i, k])
## resulting transitions from S[k] to E[k]
## (k vaccination class)
n_SE[, , 1:n_vacc_classes] <-
  n_S_progress[i, j, k] - n_SE_next_vacc_class[i, j, k]

## vaccine progression
n_S_next_vacc_class[, ] <- rbinom(S[i, j] - sum(n_S_progress[i, , j]),
                                  p_S_next_vacc_class[i, j])

#### flow out of E ####

n_E_progress[, , , ] <- rbinom(E[i, j, k, l], p_E_progress)
## of those some can also be vaccinated or progress through vaccination classes
## --> number transitioning from E[j, k] to E[j+1, k+1] (k vaccination class)
n_EE_next_vacc_class[, , , ] <-
  rbinom(n_E_progress[i, j, k, l], p_E_next_vacc_class[i, j, k, l])
## resulting transitions from E[j, k] to E[j+1, k]
## (k vaccination class)
n_EE[, , , ] <- n_E_progress[i, j, k, l] - n_EE_next_vacc_class[i, j, k, l]

## vaccine progression
n_E_next_vacc_class[, , , ] <- rbinom(E[i, j, k, l] - n_E_progress[i, j, k, l],
                                      p_E_next_vacc_class[i, j, k, l])

#### flow out of I_A ####

n_I_A_progress[, , , ] <- rbinom(I_A[i, j, k, l], p_I_A_progress[j])
## of those some can also be vaccinated or progress through vaccination classes
## --> number transitioning from I_A[j, k] to I_A[j+1, k+1]
## (k vaccination class)
n_II_A_next_vacc_class[, , , ] <-
  rbinom(n_I_A_progress[i, j, k, l],
         p_I_A_next_vacc_class[i, j, k, l])
## resulting transitions from I_A[j, l] to I_A[j + 1, l]
## (l vaccination class)
n_II_A[, , , ] <- n_I_A_progress[i, j, k, l] -
  n_II_A_next_vacc_class[i, j, k, l]

## vaccine progression
n_I_A_next_vacc_class[, , , ] <- rbinom(
  I_A[i, j, k, l] - n_I_A_progress[i, j, k, l],
  p_I_A_next_vacc_class[i, j, k, l])


#### flow out of I_P ####

n_I_P_progress[, , , ] <- rbinom(I_P[i, j, k, l], p_I_P_progress[j])
## of those some can also be vaccinated or progress through vaccination classes
## --> number transitioning from I_P[j, k] to I_P[j+1, k+1]
## (k vaccination class)
n_II_P_next_vacc_class[, , , ] <-
  rbinom(n_I_P_progress[i, j, k, l],
         p_I_P_next_vacc_class[i, j, k, l])
## resulting transitions from I_P[j, l] to I_P[j + 1, l]
## (l vaccination class)
n_II_P[, , , ] <- n_I_P_progress[i, j, k, l] -
  n_II_P_next_vacc_class[i, j, k, l]

## vaccine progression
n_I_P_next_vacc_class[, , , ] <- rbinom(
  I_P[i, j, k, l] - n_I_P_progress[i, j, k, l],
  p_I_P_next_vacc_class[i, j, k, l])


#### flow out of R ####

n_R_progress_tmp[, , ] <- rbinom(R[i, j, k], p_RS[i])
## cap on people who can move out of R based on numbers in T_sero_neg and
## T_PCR_neg
n_R_progress_capped[, , ] <-
  min(n_R_progress_tmp[i, j, k], T_sero_neg[i, j, k], T_PCR_neg[i, j, k])
## use cap or not depending on model_pcr_and_serology value
n_R_progress[, , ] <- if (model_pcr_and_serology == 1)
  n_R_progress_capped[i, j, k] else n_R_progress_tmp[i, j, k]
## of those some can also be vaccinated or progress through vaccination classes
## --> number transitioning from R[j] to S[j+1]
## (j vaccination class)
n_RS_next_vacc_class[, , ] <-
  rbinom(n_R_progress[i, j, k], p_R_next_vacc_class[i, j, k])
## resulting transitions from R[j] to S[j]
## (j vaccination class)
n_RS[, , ] <- n_R_progress[i, j, k] - n_RS_next_vacc_class[i, j, k]

## vaccine progression
n_R_next_vacc_class_tmp[, , ] <- rbinom(
  R[i, j, k] - n_R_progress[i, j, k], p_R_next_vacc_class[i, j, k])
n_R_next_vacc_class_capped[, , ] <-
  min(n_R_next_vacc_class_tmp[i, j, k],
      T_sero_neg[i, j, k] - n_R_progress[i, j, k],
      T_PCR_neg[i, j, k] - n_R_progress[i, j, k])
n_R_next_vacc_class[, , ] <- if (model_pcr_and_serology == 1)
  n_R_next_vacc_class_capped[i, j, k] else n_R_next_vacc_class_tmp[i, j, k]

#### other transitions ####

n_I_C_1_progress[, , , ] <- rbinom(I_C_1[i, j, k, l], p_I_C_1_progress[j])
n_I_C_2_progress[, , , ] <- rbinom(I_C_2[i, j, k, l], p_I_C_2_progress[j])
n_G_D_progress[, , , ] <- rbinom(G_D[i, j, k, l], p_G_D_progress)
n_ICU_pre_unconf_progress[, , , ] <-
  rbinom(ICU_pre_unconf[i, j, k, l], p_ICU_pre_progress)
n_ICU_pre_conf_progress[, , , ] <-
  rbinom(ICU_pre_conf[i, j, k, l], p_ICU_pre_progress)
n_H_R_unconf_progress[, , , ] <- rbinom(H_R_unconf[i, j, k, l], p_H_R_progress)
n_H_R_conf_progress[, , , ] <- rbinom(H_R_conf[i, j, k, l], p_H_R_progress)
n_H_D_unconf_progress[, , , ] <- rbinom(H_D_unconf[i, j, k, l], p_H_D_progress)
n_H_D_conf_progress[, , , ] <- rbinom(H_D_conf[i, j, k, l], p_H_D_progress)
n_ICU_W_R_unconf_progress[, , , ] <-
  rbinom(ICU_W_R_unconf[i, j, k, l], p_ICU_W_R_progress)
n_ICU_W_R_conf_progress[, , , ] <-
  rbinom(ICU_W_R_conf[i, j, k, l], p_ICU_W_R_progress)
n_ICU_W_D_unconf_progress[, , , ] <-
  rbinom(ICU_W_D_unconf[i, j, k, l], p_ICU_W_D_progress)
n_ICU_W_D_conf_progress[, , , ] <-
  rbinom(ICU_W_D_conf[i, j, k, l], p_ICU_W_D_progress)
n_ICU_D_unconf_progress[, , , ] <-
  rbinom(ICU_D_unconf[i, j, k, l], p_ICU_D_progress)
n_ICU_D_conf_progress[, , , ] <-
  rbinom(ICU_D_conf[i, j, k, l], p_ICU_D_progress)
n_W_R_unconf_progress[, , , ] <- rbinom(W_R_unconf[i, j, k, l], p_W_R_progress)
n_W_R_conf_progress[, , , ] <- rbinom(W_R_conf[i, j, k, l], p_W_R_progress)
n_W_D_unconf_progress[, , , ] <- rbinom(W_D_unconf[i, j, k, l], p_W_D_progress)
n_W_D_conf_progress[, , , ] <- rbinom(W_D_conf[i, j, k, l], p_W_D_progress)
n_T_sero_pre_progress[, , , ] <-
  rbinom(T_sero_pre[i, j, k, l], p_T_sero_pre_progress[i, j, k, l])
n_T_sero_pos_progress[, , , ] <-
  rbinom(T_sero_pos[i, j, k, l], p_T_sero_pos_progress)
n_T_PCR_pre_progress[, , , ] <-
  rbinom(T_PCR_pre[i, j, k, l], p_T_PCR_pre_progress)
n_T_PCR_pos_progress[, , , ] <-
  rbinom(T_PCR_pos[i, j, k, l], p_T_PCR_pos_progress)

## Cumulative infections, summed over all age groups
initial(cum_infections) <- 0
update(cum_infections) <- cum_infections + sum(n_S_progress)

initial(cum_infections_per_strain[]) <- 0
update(cum_infections_per_strain[]) <-
  cum_infections_per_strain[i] + sum(n_S_progress[, i, ])
dim(cum_infections_per_strain) <- n_strains

## Work out the new S (i for age, j for vaccination status)
new_S[, ] <- S[i, j] + sum(n_RS[i, , j]) - sum(n_S_progress[i, , j]) -
  n_S_next_vacc_class[i, j]
new_S[, ] <- new_S[i, j] +
  (if (j == 1) n_S_next_vacc_class[i, n_vacc_classes] +
     sum(n_RS_next_vacc_class[i, , n_vacc_classes]) else
       n_S_next_vacc_class[i, j - 1] + sum(n_RS_next_vacc_class[i, , j - 1]))

## Computes the number of asymptomatic
n_EI_A[, , ] <- rbinom(n_EE[i, j, k_E, k],
                       1 - p_C[i] * rel_p_sympt[i, k])
n_EI_A_next_vacc_class[, , ] <-
  rbinom(n_EE_next_vacc_class[i, j, k_E, k],
         1 - p_C[i] * rel_p_sympt[i, k])

## Computes the number of symptomatic cases
n_EI_P[, , ] <- n_EE[i, j, k_E, k] - n_EI_A[i, j, k]
n_EI_P_next_vacc_class[, , ] <- n_EE_next_vacc_class[i, j, k_E, k] -
  n_EI_A_next_vacc_class[i, j, k]

## Work out the S->E and E->E transitions
aux_E[, , , ] <- (if (k == 1) n_SE[i, j, l] else n_EE[i, j, k - 1, l]) -
  n_EE[i, j, k, l] -
  n_EE_next_vacc_class[i, j, k, l] -
  n_E_next_vacc_class[i, j, k, l] +
  (if (l == 1) n_E_next_vacc_class[i, j, k, n_vacc_classes] else
    n_E_next_vacc_class[i, j, k, l - 1]) +
  (if (k == 1) (if (l == 1) n_SE_next_vacc_class[i, j, n_vacc_classes] else
    n_SE_next_vacc_class[i, j, l - 1]) else
      (if (l == 1) n_EE_next_vacc_class[i, j, k - 1, n_vacc_classes] else
        n_EE_next_vacc_class[i, j, k - 1, l - 1]))

new_E[, , , ] <- E[i, j, k, l] + aux_E[i, j, k, l]

## Work out the I_A->I_A transitions
aux_I_A[, , , ] <- (if (k == 1) n_EI_A[i, j, l] else n_II_A[i, j, k - 1, l]) -
  n_II_A[i, j, k, l] -
  n_II_A_next_vacc_class[i, j, k, l] -
  n_I_A_next_vacc_class[i, j, k, l] +
  (if (l == 1) n_I_A_next_vacc_class[i, j, k, n_vacc_classes] else
    n_I_A_next_vacc_class[i, j, k, l - 1]) +
  (if (k == 1) (if (l == 1) n_EI_A_next_vacc_class[i, j, n_vacc_classes] else
    n_EI_A_next_vacc_class[i, j, l - 1]) else
      (if (l == 1) n_II_A_next_vacc_class[i, j, k - 1, n_vacc_classes] else
        n_II_A_next_vacc_class[i, j, k - 1, l - 1]))

new_I_A[, , , ] <- I_A[i, j, k, l] + aux_I_A[i, j, k, l]

## Work out the I_P->I_P transitions
aux_I_P[, , , ] <- (if (k == 1) n_EI_P[i, j, l] else n_II_P[i, j, k - 1, l]) -
  n_II_P[i, j, k, l] -
  n_II_P_next_vacc_class[i, j, k, l] -
  n_I_P_next_vacc_class[i, j, k, l] +
  (if (l == 1) n_I_P_next_vacc_class[i, j, k, n_vacc_classes] else
    n_I_P_next_vacc_class[i, j, k, l - 1]) +
  (if (k == 1) (if (l == 1) n_EI_P_next_vacc_class[i, j, n_vacc_classes] else
    n_EI_P_next_vacc_class[i, j, l - 1]) else
      (if (l == 1) n_II_P_next_vacc_class[i, j, k - 1, n_vacc_classes] else
        n_II_P_next_vacc_class[i, j, k - 1, l - 1]))

new_I_P[, , , ] <- I_P[i, j, k, l] + aux_I_P[i, j, k, l]

## Work out the I_C_1->I_C_1 transitions
aux_I_C_1[, , , ] <- (if (k == 1)
  (if (l == 1) n_II_P[i, j, k_P, 1] +
     n_II_P_next_vacc_class[i, j, k_P, n_vacc_classes] else
       n_II_P[i, j, k_P, l] + n_II_P_next_vacc_class[i, j, k_P, l - 1]) else
         n_I_C_1_progress[i, j, k - 1, l]) - n_I_C_1_progress[i, j, k, l]

new_I_C_1[, , , ] <- I_C_1[i, j, k, l] + aux_I_C_1[i, j, k, l]

## Work out the I_C_2->I_C_2 transitions
aux_I_C_2[, , , ] <- (if (k == 1) n_I_C_1_progress[i, j, k_C_1, l] else
  n_I_C_2_progress[i, j, k - 1, l]) - n_I_C_2_progress[i, j, k, l]

new_I_C_2[, , , ] <- I_C_2[i, j, k, l] + aux_I_C_2[i, j, k, l]

## Work out the flow from I_C_2 -> R, G_D, hosp
n_I_C_2_to_R[, , ] <- rbinom(n_I_C_2_progress[i, j, k_C_2, k],
                             1 - p_H_by_age[i] * rel_p_hosp_if_sympt[i, k])
n_I_C_2_to_G_D[, , ] <-
  rbinom(n_I_C_2_progress[i, j, k_C_2, k] - n_I_C_2_to_R[i, j, k],
         p_G_D_by_age[i, j])
n_I_C_2_to_hosp[, , ] <- n_I_C_2_progress[i, j, k_C_2, k] -
  n_I_C_2_to_R[i, j, k] - n_I_C_2_to_G_D[i, j, k]

## Work out the G_D -> G_D transitions
aux_G_D[, , , ] <- (if (k == 1) n_I_C_2_to_G_D[i, j, l] else
  n_G_D_progress[i, j, k - 1, l]) - n_G_D_progress[i, j, k, l]

new_G_D[, , , ] <- G_D[i, j, k, l] + aux_G_D[i, j, k, l]

## Work out the split in hospitals between H_D, H_R and ICU_pre
n_I_C_2_to_ICU_pre[, , ] <- rbinom(n_I_C_2_to_hosp[i, j, k], p_ICU_by_age[i])
n_I_C_2_to_ICU_pre_conf[, , ] <- rbinom(n_I_C_2_to_ICU_pre[i, j, k],
                                        p_star_by_age[i])
n_hosp_non_ICU[, , ] <- n_I_C_2_to_hosp[i, j, k] - n_I_C_2_to_ICU_pre[i, j, k]
n_I_C_2_to_H_D[, , ] <- rbinom(n_hosp_non_ICU[i, j, k], p_H_D_by_age[i, j])
n_I_C_2_to_H_D_conf[, , ] <- rbinom(n_I_C_2_to_H_D[i, j, k],
                                    p_star_by_age[i])
n_I_C_2_to_H_R[, , ] <- n_hosp_non_ICU[i, j, k] - n_I_C_2_to_H_D[i, j, k]
n_I_C_2_to_H_R_conf[, , ] <- rbinom(n_I_C_2_to_H_R[i, j, k],
                                    p_star_by_age[i])

## Work out the ICU_pre -> ICU_pre transitions
aux_ICU_pre_unconf[, , , ] <- ICU_pre_unconf[i, j, k, l] +
  (if (k > 1) n_ICU_pre_unconf_progress[i, j, k - 1, l] else 0) -
  n_ICU_pre_unconf_progress[i, j, k, l]
aux_ICU_pre_conf[, , , ] <- ICU_pre_conf[i, j, k, l] +
  (if (k > 1) n_ICU_pre_conf_progress[i, j, k - 1, l] else 0) -
  n_ICU_pre_conf_progress[i, j, k, l]

n_ICU_pre_unconf_to_conf[, , , ] <-
  rbinom(aux_ICU_pre_unconf[i, j, k, l], p_test)

new_ICU_pre_unconf[, , , ] <-
  aux_ICU_pre_unconf[i, j, k, l] - n_ICU_pre_unconf_to_conf[i, j, k, l] +
  (if (k == 1) n_I_C_2_to_ICU_pre[i, j, l] - n_I_C_2_to_ICU_pre_conf[i, j, l]
   else 0)
new_ICU_pre_conf[, , , ] <-
  aux_ICU_pre_conf[i, j, k, l] + n_ICU_pre_unconf_to_conf[i, j, k, l] +
  (if (k == 1) n_I_C_2_to_ICU_pre_conf[i, j, l] else 0)

## Work out the H_R->H_R transitions
aux_H_R_unconf[, , , ] <- H_R_unconf[i, j, k, l] +
  (if (k > 1) n_H_R_unconf_progress[i, j, k - 1, l] else 0) -
  n_H_R_unconf_progress[i, j, k, l]
aux_H_R_conf[, , , ] <- H_R_conf[i, j, k, l] +
  (if (k > 1) n_H_R_conf_progress[i, j, k - 1, l] else 0) -
  n_H_R_conf_progress[i, j, k, l]

n_H_R_unconf_to_conf[, , , ] <-
  rbinom(aux_H_R_unconf[i, j, k, l], p_test)

new_H_R_unconf[, , , ] <-
  aux_H_R_unconf[i, j, k, l] - n_H_R_unconf_to_conf[i, j, k, l] +
  (if (k == 1) n_I_C_2_to_H_R[i, j, l] - n_I_C_2_to_H_R_conf[i, j, l] else 0)
new_H_R_conf[, , , ] <-
  aux_H_R_conf[i, j, k, l] + n_H_R_unconf_to_conf[i, j, k, l] +
  (if (k == 1) n_I_C_2_to_H_R_conf[i, j, l] else 0)

## Work out the H_D->H_D transitions
aux_H_D_unconf[, , , ] <- H_D_unconf[i, j, k, l] +
  (if (k > 1) n_H_D_unconf_progress[i, j, k - 1, l] else 0) -
  n_H_D_unconf_progress[i, j, k, l]
aux_H_D_conf[, , , ] <- H_D_conf[i, j, k, l] +
  (if (k > 1) n_H_D_conf_progress[i, j, k - 1, l] else 0) -
  n_H_D_conf_progress[i, j, k, l]

n_H_D_unconf_to_conf[, , , ] <-
  rbinom(aux_H_D_unconf[i, j, k, l], p_test)

new_H_D_unconf[, , , ] <-
  aux_H_D_unconf[i, j, k, l] - n_H_D_unconf_to_conf[i, j, k, l] +
  (if (k == 1) n_I_C_2_to_H_D[i, j, l] - n_I_C_2_to_H_D_conf[i, j, l] else 0)
new_H_D_conf[, , , ] <-
  aux_H_D_conf[i, j, k, l] + n_H_D_unconf_to_conf[i, j, k, l] +
  (if (k == 1) n_I_C_2_to_H_D_conf[i, j, l] else 0)

## Work out the ICU_pre to ICU_D, ICU_W_R and ICU_W_D splits
n_ICU_pre_unconf_to_ICU_D_unconf[, , ] <-
  rbinom(n_ICU_pre_unconf_progress[i, j, k_ICU_pre, k], p_ICU_D_by_age[i, j])
n_ICU_pre_conf_to_ICU_D_conf[, , ] <-
  rbinom(n_ICU_pre_conf_progress[i, j, k_ICU_pre, k], p_ICU_D_by_age[i, j])
n_ICU_pre_unconf_to_ICU_W_D_unconf[, , ] <-
  rbinom(n_ICU_pre_unconf_progress[i, j, k_ICU_pre, k] -
           n_ICU_pre_unconf_to_ICU_D_unconf[i, j, k],
         p_W_D_by_age[i, j])
n_ICU_pre_unconf_to_ICU_W_R_unconf[, , ] <-
  n_ICU_pre_unconf_progress[i, j, k_ICU_pre, k] -
  n_ICU_pre_unconf_to_ICU_D_unconf[i, j, k] -
  n_ICU_pre_unconf_to_ICU_W_D_unconf[i, j, k]
n_ICU_pre_conf_to_ICU_W_D_conf[, , ] <-
  rbinom(n_ICU_pre_conf_progress[i, j, k_ICU_pre, k] -
           n_ICU_pre_conf_to_ICU_D_conf[i, j, k], p_W_D_by_age[i, j])
n_ICU_pre_conf_to_ICU_W_R_conf[, , ] <-
  n_ICU_pre_conf_progress[i, j, k_ICU_pre, k] -
  n_ICU_pre_conf_to_ICU_D_conf[i, j, k] -
  n_ICU_pre_conf_to_ICU_W_D_conf[i, j, k]


## Work out the ICU_W_R->ICU_W_R transitions
aux_ICU_W_R_unconf[, , , ] <- ICU_W_R_unconf[i, j, k, l] +
  (if (k == 1) n_ICU_pre_unconf_to_ICU_W_R_unconf[i, j, l] else
    n_ICU_W_R_unconf_progress[i, j, k - 1, l]) -
  n_ICU_W_R_unconf_progress[i, j, k, l]
aux_ICU_W_R_conf[, , , ] <- ICU_W_R_conf[i, j, k, l] +
  (if (k == 1) n_ICU_pre_conf_to_ICU_W_R_conf[i, j, l] else
    n_ICU_W_R_conf_progress[i, j, k - 1, l]) -
  n_ICU_W_R_conf_progress[i, j, k, l]

n_ICU_W_R_unconf_to_conf[, , , ] <-
  rbinom(aux_ICU_W_R_unconf[i, j, k, l], p_test)
new_ICU_W_R_unconf[, , , ] <-
  aux_ICU_W_R_unconf[i, j, k, l] - n_ICU_W_R_unconf_to_conf[i, j, k, l]
new_ICU_W_R_conf[, , , ] <-
  aux_ICU_W_R_conf[i, j, k, l] + n_ICU_W_R_unconf_to_conf[i, j, k, l]

## Work out the ICU_W_D->ICU_W_D transitions
aux_ICU_W_D_unconf[, , , ] <- ICU_W_D_unconf[i, j, k, l] +
  (if (k == 1) n_ICU_pre_unconf_to_ICU_W_D_unconf[i, j, l] else
    n_ICU_W_D_unconf_progress[i, j, k - 1, l]) -
  n_ICU_W_D_unconf_progress[i, j, k, l]
aux_ICU_W_D_conf[, , , ] <- ICU_W_D_conf[i, j, k, l] +
  (if (k == 1) n_ICU_pre_conf_to_ICU_W_D_conf[i, j, l] else
    n_ICU_W_D_conf_progress[i, j, k - 1, l]) -
  n_ICU_W_D_conf_progress[i, j, k, l]

n_ICU_W_D_unconf_to_conf[, , , ] <-
  rbinom(aux_ICU_W_D_unconf[i, j, k, l], p_test)
new_ICU_W_D_unconf[, , , ] <-
  aux_ICU_W_D_unconf[i, j, k, l] - n_ICU_W_D_unconf_to_conf[i, j, k, l]
new_ICU_W_D_conf[, , , ] <-
  aux_ICU_W_D_conf[i, j, k, l] + n_ICU_W_D_unconf_to_conf[i, j, k, l]

## Work out the ICU_D->ICU_D transitions
aux_ICU_D_unconf[, , , ] <- ICU_D_unconf[i, j, k, l] +
  (if (k == 1) n_ICU_pre_unconf_to_ICU_D_unconf[i, j, l] else
    n_ICU_D_unconf_progress[i, j, k - 1, l]) -
  n_ICU_D_unconf_progress[i, j, k, l]
aux_ICU_D_conf[, , , ] <- ICU_D_conf[i, j, k, l] +
  (if (k == 1) n_ICU_pre_conf_to_ICU_D_conf[i, j, l] else
    n_ICU_D_conf_progress[i, j, k - 1, l]) -
  n_ICU_D_conf_progress[i, j, k, l]

n_ICU_D_unconf_to_conf[, , , ] <-
  rbinom(aux_ICU_D_unconf[i, j, k, l], p_test)
new_ICU_D_unconf[, , , ] <-
  aux_ICU_D_unconf[i, j, k, l] - n_ICU_D_unconf_to_conf[i, j, k, l]
new_ICU_D_conf[, , , ] <-
  aux_ICU_D_conf[i, j, k, l] + n_ICU_D_unconf_to_conf[i, j, k, l]

## Work out the W_R->W_R transitions
aux_W_R_unconf[, , , ] <- W_R_unconf[i, j, k, l] +
  (if (k == 1) n_ICU_W_R_unconf_progress[i, j, k_ICU_W_R, l] else
    n_W_R_unconf_progress[i, j, k - 1, l]) -
  n_W_R_unconf_progress[i, j, k, l]
aux_W_R_conf[, , , ] <- W_R_conf[i, j, k, l] +
  (if (k == 1) n_ICU_W_R_conf_progress[i, j, k_ICU_W_R, l] else
    n_W_R_conf_progress[i, j, k - 1, l]) -
  n_W_R_conf_progress[i, j, k, l]

n_W_R_unconf_to_conf[, , , ] <-
  rbinom(aux_W_R_unconf[i, j, k, l], p_test)
new_W_R_unconf[, , , ] <-
  aux_W_R_unconf[i, j, k, l] -
  n_W_R_unconf_to_conf[i, j, k, l]
new_W_R_conf[, , , ] <-
  aux_W_R_conf[i, j, k, l] + n_W_R_unconf_to_conf[i, j, k, l]

## Work out the W_D->W_D transitions
aux_W_D_unconf[, , , ] <- W_D_unconf[i, j, k, l] +
  (if (k == 1) n_ICU_W_D_unconf_progress[i, j, k_ICU_W_D, l] else
    n_W_D_unconf_progress[i, j, k - 1, l]) -
  n_W_D_unconf_progress[i, j, k, l]
aux_W_D_conf[, , , ] <- W_D_conf[i, j, k, l] +
  (if (k == 1) n_ICU_W_D_conf_progress[i, j, k_ICU_W_D, l] else
    n_W_D_conf_progress[i, j, k - 1, l]) -
  n_W_D_conf_progress[i, j, k, l]

n_W_D_unconf_to_conf[, , , ] <-
  rbinom(aux_W_D_unconf[i, j, k, l], p_test)
new_W_D_unconf[, , , ] <-
  aux_W_D_unconf[i, j, k, l] -
  n_W_D_unconf_to_conf[i, j, k, l]
new_W_D_conf[, , , ] <-
  aux_W_D_conf[i, j, k, l] + n_W_D_unconf_to_conf[i, j, k, l]

## Work out the number of deaths in hospital

delta_D_hosp_disag[, ] <-
  sum(n_H_D_unconf_progress[i, , k_H_D, j]) +
  sum(n_H_D_conf_progress[i, , k_H_D, j]) +
  sum(n_ICU_D_unconf_progress[i, , k_ICU_D, j]) +
  sum(n_ICU_D_conf_progress[i, , k_ICU_D, j]) +
  sum(n_W_D_unconf_progress[i, , k_W_D, j]) +
  sum(n_W_D_conf_progress[i, , k_W_D, j])
delta_D_non_hosp_disag[, ] <- sum(n_G_D_progress[i, , k_G_D, j])
dim(delta_D_hosp_disag) <- c(n_groups, n_vacc_classes)
dim(delta_D_non_hosp_disag) <- c(n_groups, n_vacc_classes)

initial(D[, ]) <- 0
update(D[, ]) <- D[i, j] +
  delta_D_hosp_disag[i, j] + delta_D_non_hosp_disag[i, j]
dim(D) <- c(n_groups, n_vacc_classes)

delta_D_hosp[, ] <- sum(delta_D_hosp_disag[i, j])

## Work out the number of deaths in the community
delta_D_non_hosp[] <- sum(delta_D_non_hosp_disag[i, ])

## Work out the number of people entering the seroconversion flow
n_com_to_T_sero_pre[, , 1, ] <- rbinom(
  n_EE[i, j, k_E, l] +
    (if (l == 1) n_EE_next_vacc_class[i, j, k_E, n_vacc_classes] else
      n_EE_next_vacc_class[i, j, k_E, l - 1]),
  p_sero_pre_1)
n_com_to_T_sero_pre[, , 2, ] <- n_EE[i, j, k_E, l] +
  (if (l == 1) n_EE_next_vacc_class[i, j, k_E, n_vacc_classes] else
    n_EE_next_vacc_class[i, j, k_E, l - 1]) -
  n_com_to_T_sero_pre[i, j, 1, l]

new_T_sero_pre[, , , ] <- T_sero_pre[i, j, k, l] +
  n_com_to_T_sero_pre[i, j, k, l] - n_T_sero_pre_progress[i, j, k, l]

## Split the seroconversion flow between people who are going to
## seroconvert and people who are not
n_T_sero_pre_to_T_sero_pos[, , ] <-
  rbinom(sum(n_T_sero_pre_progress[i, j, , k]), p_sero_pos[i])

new_T_sero_pos[, , , ] <- T_sero_pos[i, j, k, l] -
  n_T_sero_pos_progress[i, j, k, l] +
  (if (k == 1)  n_T_sero_pre_to_T_sero_pos[i, j, l] else
    n_T_sero_pos_progress[i, j, k - 1, l])

new_T_sero_neg[, , ] <- T_sero_neg[i, j, k] +
  sum(n_T_sero_pre_progress[i, j, , k]) - n_T_sero_pre_to_T_sero_pos[i, j, k] +
  n_T_sero_pos_progress[i, j, k_sero_pos, k] -
  model_pcr_and_serology * n_R_progress[i, j, k] -
  model_pcr_and_serology * n_R_next_vacc_class[i, j, k] +
  model_pcr_and_serology *
  (if (k == 1) n_R_next_vacc_class[i, j, n_vacc_classes] else
    n_R_next_vacc_class[i, j, k - 1])

## Work out the total number of recovery
new_R[, , ] <- R[i, j, k] +
  n_II_A[i, j, k_A, k] +
  n_I_C_2_to_R[i, j, k] +
  n_H_R_conf_progress[i, j, k_H_R, k] +
  n_H_R_unconf_progress[i, j, k_H_R, k] +
  n_W_R_conf_progress[i, j, k_W_R, k] +
  n_W_R_unconf_progress[i, j, k_W_R, k] -
  n_R_progress[i, j, k] -
  n_R_next_vacc_class[i, j, k] +
  (if (k == 1) n_II_A_next_vacc_class[i, j, k_A, n_vacc_classes] +
     n_R_next_vacc_class[i, j, n_vacc_classes] else
       n_II_A_next_vacc_class[i, j, k_A, k - 1] +
     n_R_next_vacc_class[i, j, k - 1])

## Work out the PCR positivity
new_T_PCR_pre[, , , ] <- T_PCR_pre[i, j, k, l] -
  n_T_PCR_pre_progress[i, j, k, l] +
  (if (k == 1) n_S_progress[i, j, l] else n_T_PCR_pre_progress[i, j, k - 1, l])

new_T_PCR_pos[, , , ] <- T_PCR_pos[i, j, k, l] -
  n_T_PCR_pos_progress[i, j, k, l] +
  (if (k == 1) n_T_PCR_pre_progress[i, j, k_PCR_pre, l] else
    n_T_PCR_pos_progress[i, j, k - 1, l])

new_T_PCR_neg[, , ] <- T_PCR_neg[i, j, k] +
  n_T_PCR_pos_progress[i, j, k_PCR_pos, k] -
  model_pcr_and_serology * n_R_progress[i, j, k] -
  model_pcr_and_serology * n_R_next_vacc_class[i, j, k] +
  model_pcr_and_serology *
  (if (k == 1) n_R_next_vacc_class[i, j, n_vacc_classes] else
    n_R_next_vacc_class[i, j, k - 1])

## Compute the force of infection

I_with_diff_trans[, , ] <-
  rel_infectivity[i, k] * strain_transmission[j] * (
    I_A_transmission * sum(I_A[i, j, , k]) +
      I_P_transmission * sum(I_P[i, j, , k]) +
      I_C_1_transmission * sum(I_C_1[i, j, , k]) +
      I_C_2_transmission * sum(I_C_2[i, j, , k]) +
      hosp_transmission * (
        sum(ICU_pre_unconf[i, j, , k]) +
          sum(ICU_pre_conf[i, j, , k]) +
          sum(H_R_unconf[i, j, , k]) +
          sum(H_R_conf[i, j, , k]) +
          sum(H_D_unconf[i, j, , k]) +
          sum(H_D_conf[i, j, , k])) +
      ICU_transmission * (
        sum(ICU_W_R_unconf[i, j, , k]) +
          sum(ICU_W_R_conf[i, j, , k]) +
          sum(ICU_W_D_unconf[i, j, , k]) +
          sum(ICU_W_D_conf[i, j, , k]) +
          sum(ICU_D_unconf[i, j, , k]) +
          sum(ICU_D_conf[i, j, , k])) +
      G_D_transmission * sum(G_D[i, j, , k]))

## NOTE: "age groups" 1-17 are age groups, 18 are CHW and 19 CHR. Here we apply
## beta to all contacts *except* within care home contacts
s_ij[, , ] <- m[i, j] * sum(I_with_diff_trans[j, k, ])
s_ij[1:n_age_groups, 1:n_groups, ] <- beta * s_ij[i, j, k]
s_ij[(n_age_groups + 1):n_groups, 1:n_age_groups, ] <- beta * s_ij[i, j, k]
lambda[, ] <- sum(s_ij[i, , j])

## Initial states are all zerod as we will provide a state vector
## setting S and I based on the seeding model.
initial(S[, ]) <- 0
initial(E[, , , ]) <- 0
initial(I_A[, , , ]) <- 0
initial(I_P[, , , ]) <- 0
initial(I_C_1[, , , ]) <- 0
initial(I_C_2[, , , ]) <- 0
initial(G_D[, , , ]) <- 0
initial(ICU_pre_unconf[, , , ]) <- 0
initial(ICU_pre_conf[, , , ]) <- 0
initial(H_R_unconf[, , , ]) <- 0
initial(H_R_conf[, , , ]) <- 0
initial(H_D_unconf[, , , ]) <- 0
initial(H_D_conf[, , , ]) <- 0
initial(ICU_W_R_unconf[, , , ]) <- 0
initial(ICU_W_R_conf[, , , ]) <- 0
initial(ICU_W_D_unconf[, , , ]) <- 0
initial(ICU_W_D_conf[, , , ]) <- 0
initial(ICU_D_unconf[, , , ]) <- 0
initial(ICU_D_conf[, , , ]) <- 0
initial(W_R_unconf[, , , ]) <- 0
initial(W_R_conf[, , , ]) <- 0
initial(W_D_unconf[, , , ]) <- 0
initial(W_D_conf[, , , ]) <- 0
initial(T_sero_pre[, , , ]) <- 0
initial(T_sero_pos[, , , ]) <- 0
initial(T_sero_neg[, , ]) <- 0
initial(R[, , ]) <- 0
initial(D_hosp[, ]) <- 0
initial(D_non_hosp[]) <- 0
initial(T_PCR_pre[, , , ]) <- 0
initial(T_PCR_pos[, , , ]) <- 0
initial(T_PCR_neg[, , ]) <- 0
initial(cum_admit_conf) <- 0
initial(cum_new_conf) <- 0
initial(cum_admit_by_age[]) <- 0

## User defined parameters - default in parentheses:

## Vaccination parameters
rel_susceptibility[, ] <- user()
dim(rel_susceptibility) <- c(n_groups, n_vacc_classes)
n_vacc_classes <- user()
rel_p_sympt[, ] <- user()
dim(rel_p_sympt) <- c(n_groups, n_vacc_classes)
rel_p_hosp_if_sympt[, ] <- user()
dim(rel_p_hosp_if_sympt) <- c(n_groups, n_vacc_classes)
rel_infectivity[, ] <- user()
dim(rel_infectivity) <- c(n_groups, n_vacc_classes)

vaccine_progression_rate_base[, ] <- user()
dim(vaccine_progression_rate_base) <- c(n_groups, n_vacc_classes)

## Parameters of the E classes
k_E <- user()
gamma_E <- user(0.1)

## Probability of transitioning from the E to the asymptomatic class,
## the rest go into the symptomatic class
p_C[] <- user()

## Parameters of the I_A classes
k_A <- user()
gamma_A[] <- user()
dim(gamma_A) <- n_strains

## Parameters of the I_P classes
k_P <- user()
gamma_P[] <- user()
dim(gamma_P) <- n_strains

## Parameters of the I_C_1 classes
k_C_1 <- user()
gamma_C_1[] <- user()
dim(gamma_C_1) <- n_strains

## Parameters of the I_C_2 classes
k_C_2 <- user()
gamma_C_2[] <- user()
dim(gamma_C_2) <- n_strains
dim(p_H_step) <- user()
p_H_step[] <- user()
psi_H[] <- user()

## Parameters of the G_D class
k_G_D <- user()
gamma_G_D <- user(0.1)
dim(p_G_D) <- n_strains
p_G_D_step[, ] <- user()
dim(p_G_D_step) <- user()
psi_G_D[, ] <- user()

## Parameters of the ICU_pre classes
k_ICU_pre <- user()
dim(gamma_ICU_pre_step) <- user()
gamma_ICU_pre_step[] <- user()

## Proportion of hospital cases progressing to ICU
dim(p_ICU_step) <- user()
p_ICU_step[] <- user()
psi_ICU[] <- user()

## Proportion of stepdown cases dying
dim(p_W_D_step) <- user()
dim(p_W_D) <- n_strains
p_W_D_step[, ] <- user()
psi_W_D[, ] <- user()

## Parameters of the H_R classes
k_H_R <- user()
dim(gamma_H_R_step) <- user()
gamma_H_R_step[] <- user()

## Parameters of the H_D classes
k_H_D <- user()
dim(gamma_H_D_step) <- user()
gamma_H_D_step[] <- user()
dim(p_H_D) <- n_strains
dim(p_H_D_step) <- user()
p_H_D_step[, ] <- user()
psi_H_D[, ] <- user()

## Parameters of the ICU_W_R classes
k_ICU_W_R <- user()
dim(gamma_ICU_W_R_step) <- user()
gamma_ICU_W_R_step[] <- user()

## Parameters of the ICU_W_D classes
k_ICU_W_D <- user()
dim(gamma_ICU_W_D_step) <- user()
gamma_ICU_W_D_step[] <- user()

## Parameters of the ICU_D classes
k_ICU_D <- user()
dim(p_ICU_D) <- n_strains
dim(gamma_ICU_D_step) <- user()
gamma_ICU_D_step[] <- user()
dim(p_ICU_D_step) <- user()
p_ICU_D_step[, ] <- user()
psi_ICU_D[, ] <- user()

## Waning of immunity
waning_rate[] <- user()
dim(waning_rate) <- n_groups

## Parameters of the W_R classes
k_W_R <- user()
dim(gamma_W_R_step) <- user()
gamma_W_R_step[] <- user()

## Parameters of the W_D classes
k_W_D <- user()
dim(gamma_W_D_step) <- user()
gamma_W_D_step[] <- user()

## Parameters of the T_sero_pre classes
gamma_sero_pre_1 <- user(0.1)
gamma_sero_pre_2 <- user(0.1)
gamma_sero_pre[1] <- gamma_sero_pre_1
gamma_sero_pre[2] <- gamma_sero_pre_2
## Governs the mixing - pretty much only makes sense at 0.5
p_sero_pre_1 <- user(0.5)
p_sero_pos[] <- user()

## Parameters of the T_sero_pos classes
k_sero_pos <- user()
gamma_sero_pos <- user(0.1)

## Parameters relating to testing
gamma_U <- user(0.1)
dim(p_star_step) <- user()
p_star_step[] <- user()
psi_star[] <- user()

## Parameters relating to PCR positivity
k_PCR_pre <- user()
gamma_PCR_pre <- user(0.1)
k_PCR_pos <- user()
gamma_PCR_pos <- user(0.1)

## Parameters of the age stratified transmission
beta_step[] <- user()
dim(beta_step) <- user()
## What we really want is min(step + 1, length(beta_step)) but that's not
## supported by odin (it could be made to support this).
beta <- if (as.integer(step) >= length(beta_step))
  beta_step[length(beta_step)] else beta_step[step + 1]

## Useful for debugging
initial(beta_out) <- beta_step[1]
update(beta_out) <- beta

m[, ] <- user()
I_A_transmission <- user()
I_P_transmission <- user()
I_C_1_transmission <- user()
I_C_2_transmission <- user()
hosp_transmission <- user()
ICU_transmission <- user()
G_D_transmission <- user()
strain_transmission[] <- user()
dim(strain_transmission) <- n_strains
n_strains <- user()

## Dimensions of the different "vectors" here vectors stand for
## multi-dimensional arrays

## Vectors handling the S class
dim(S) <- c(n_groups, n_vacc_classes)
dim(new_S) <- c(n_groups, n_vacc_classes)

## Vectors handling the E class
dim(E) <- c(n_groups, n_strains, k_E, n_vacc_classes)
dim(aux_E) <- c(n_groups, n_strains, k_E, n_vacc_classes)
dim(new_E) <- c(n_groups, n_strains, k_E, n_vacc_classes)
dim(n_EE) <- c(n_groups, n_strains, k_E, n_vacc_classes)

## Vectors handling the I_A class
dim(I_A) <- c(n_groups, n_strains, k_A, n_vacc_classes)
dim(aux_I_A) <- c(n_groups, n_strains, k_A, n_vacc_classes)
dim(new_I_A) <- c(n_groups, n_strains, k_A, n_vacc_classes)
dim(n_II_A) <- c(n_groups, n_strains, k_A, n_vacc_classes)

## Vectors handling the I_P class
dim(I_P) <- c(n_groups, n_strains, k_P, n_vacc_classes)
dim(aux_I_P) <- c(n_groups, n_strains, k_P, n_vacc_classes)
dim(new_I_P) <- c(n_groups, n_strains, k_P, n_vacc_classes)
dim(n_II_P) <- c(n_groups, n_strains, k_P, n_vacc_classes)

## Vectors handling the I_C_2 class
dim(I_C_1) <- c(n_groups, n_strains, k_C_1, n_vacc_classes)
dim(aux_I_C_1) <- c(n_groups, n_strains, k_C_1, n_vacc_classes)
dim(new_I_C_1) <- c(n_groups, n_strains, k_C_1, n_vacc_classes)
dim(n_I_C_1_progress) <- c(n_groups, n_strains, k_C_1, n_vacc_classes)

## Vectors handling the I_C_2 class
dim(I_C_2) <- c(n_groups, n_strains, k_C_2, n_vacc_classes)
dim(aux_I_C_2) <- c(n_groups, n_strains, k_C_2, n_vacc_classes)
dim(new_I_C_2) <- c(n_groups, n_strains, k_C_2, n_vacc_classes)
dim(n_I_C_2_progress) <- c(n_groups, n_strains, k_C_2, n_vacc_classes)
dim(p_H_by_age) <- n_groups
dim(psi_H) <- n_groups

## Vectors handling the G_D class
dim(G_D) <- c(n_groups, n_strains, k_G_D, n_vacc_classes)
dim(aux_G_D) <- c(n_groups, n_strains, k_G_D, n_vacc_classes)
dim(new_G_D) <- c(n_groups, n_strains, k_G_D, n_vacc_classes)
dim(n_G_D_progress) <- c(n_groups, n_strains, k_G_D, n_vacc_classes)
dim(p_G_D_by_age) <- c(n_groups, n_strains)
dim(psi_G_D) <- c(n_groups, n_strains)

## Vectors handling the ICU_pre class
dim(ICU_pre_unconf) <- c(n_groups, n_strains, k_ICU_pre, n_vacc_classes)
dim(aux_ICU_pre_unconf) <- c(n_groups, n_strains, k_ICU_pre, n_vacc_classes)
dim(new_ICU_pre_unconf) <- c(n_groups, n_strains, k_ICU_pre, n_vacc_classes)
dim(n_ICU_pre_unconf_progress) <-
  c(n_groups, n_strains, k_ICU_pre, n_vacc_classes)
dim(ICU_pre_conf) <- c(n_groups, n_strains, k_ICU_pre, n_vacc_classes)
dim(aux_ICU_pre_conf) <- c(n_groups, n_strains, k_ICU_pre, n_vacc_classes)
dim(new_ICU_pre_conf) <- c(n_groups, n_strains, k_ICU_pre, n_vacc_classes)
dim(n_ICU_pre_conf_progress) <-
  c(n_groups, n_strains, k_ICU_pre, n_vacc_classes)
dim(n_ICU_pre_unconf_to_conf) <-
  c(n_groups, n_strains, k_ICU_pre, n_vacc_classes)

## Vector handling who progress to ICU
dim(p_ICU_by_age) <- n_groups
dim(psi_ICU) <- n_groups

## Vectors handling the H_R class
dim(H_R_unconf) <- c(n_groups, n_strains, k_H_R, n_vacc_classes)
dim(aux_H_R_unconf) <- c(n_groups, n_strains, k_H_R, n_vacc_classes)
dim(new_H_R_unconf) <- c(n_groups, n_strains, k_H_R, n_vacc_classes)
dim(n_H_R_unconf_progress) <- c(n_groups, n_strains, k_H_R, n_vacc_classes)
dim(H_R_conf) <- c(n_groups, n_strains, k_H_R, n_vacc_classes)
dim(aux_H_R_conf) <- c(n_groups, n_strains, k_H_R, n_vacc_classes)
dim(new_H_R_conf) <- c(n_groups, n_strains, k_H_R, n_vacc_classes)
dim(n_H_R_conf_progress) <- c(n_groups, n_strains, k_H_R, n_vacc_classes)
dim(n_H_R_unconf_to_conf) <-
  c(n_groups, n_strains, k_H_R, n_vacc_classes)

## Vectors handling the H_D class
dim(H_D_unconf) <- c(n_groups, n_strains, k_H_D, n_vacc_classes)
dim(aux_H_D_unconf) <- c(n_groups, n_strains, k_H_D, n_vacc_classes)
dim(new_H_D_unconf) <- c(n_groups, n_strains, k_H_D, n_vacc_classes)
dim(n_H_D_unconf_progress) <- c(n_groups, n_strains, k_H_D, n_vacc_classes)
dim(H_D_conf) <- c(n_groups, n_strains, k_H_D, n_vacc_classes)
dim(aux_H_D_conf) <- c(n_groups, n_strains, k_H_D, n_vacc_classes)
dim(new_H_D_conf) <- c(n_groups, n_strains, k_H_D, n_vacc_classes)
dim(n_H_D_conf_progress) <- c(n_groups, n_strains, k_H_D, n_vacc_classes)
dim(n_H_D_unconf_to_conf) <-
  c(n_groups, n_strains, k_H_D, n_vacc_classes)

## Vectors handling the ICU_W_R class
dim(ICU_W_R_unconf) <- c(n_groups, n_strains, k_ICU_W_R, n_vacc_classes)
dim(aux_ICU_W_R_unconf) <- c(n_groups, n_strains, k_ICU_W_R, n_vacc_classes)
dim(new_ICU_W_R_unconf) <- c(n_groups, n_strains, k_ICU_W_R, n_vacc_classes)
dim(n_ICU_W_R_unconf_progress) <-
  c(n_groups, n_strains, k_ICU_W_R, n_vacc_classes)
dim(ICU_W_R_conf) <- c(n_groups, n_strains, k_ICU_W_R, n_vacc_classes)
dim(aux_ICU_W_R_conf) <- c(n_groups, n_strains, k_ICU_W_R, n_vacc_classes)
dim(new_ICU_W_R_conf) <- c(n_groups, n_strains, k_ICU_W_R, n_vacc_classes)
dim(n_ICU_W_R_conf_progress) <-
  c(n_groups, n_strains, k_ICU_W_R, n_vacc_classes)
dim(n_ICU_W_R_unconf_to_conf) <-
  c(n_groups, n_strains, k_ICU_W_R, n_vacc_classes)

## Vectors handling the ICU_W_D class
dim(ICU_W_D_unconf) <- c(n_groups, n_strains, k_ICU_W_D, n_vacc_classes)
dim(aux_ICU_W_D_unconf) <- c(n_groups, n_strains, k_ICU_W_D, n_vacc_classes)
dim(new_ICU_W_D_unconf) <- c(n_groups, n_strains, k_ICU_W_D, n_vacc_classes)
dim(n_ICU_W_D_unconf_progress) <-
  c(n_groups, n_strains, k_ICU_W_D, n_vacc_classes)
dim(ICU_W_D_conf) <- c(n_groups, n_strains, k_ICU_W_D, n_vacc_classes)
dim(aux_ICU_W_D_conf) <- c(n_groups, n_strains, k_ICU_W_D, n_vacc_classes)
dim(new_ICU_W_D_conf) <- c(n_groups, n_strains, k_ICU_W_D, n_vacc_classes)
dim(n_ICU_W_D_conf_progress) <-
  c(n_groups, n_strains, k_ICU_W_D, n_vacc_classes)
dim(n_ICU_W_D_unconf_to_conf) <-
  c(n_groups, n_strains, k_ICU_W_D, n_vacc_classes)

## Vectors handling the ICU_D class
dim(ICU_D_unconf) <- c(n_groups, n_strains, k_ICU_D, n_vacc_classes)
dim(aux_ICU_D_unconf) <- c(n_groups, n_strains, k_ICU_D, n_vacc_classes)
dim(new_ICU_D_unconf) <- c(n_groups, n_strains, k_ICU_D, n_vacc_classes)
dim(n_ICU_D_unconf_progress) <- c(n_groups, n_strains, k_ICU_D, n_vacc_classes)
dim(ICU_D_conf) <- c(n_groups, n_strains, k_ICU_D, n_vacc_classes)
dim(aux_ICU_D_conf) <- c(n_groups, n_strains, k_ICU_D, n_vacc_classes)
dim(new_ICU_D_conf) <- c(n_groups, n_strains, k_ICU_D, n_vacc_classes)
dim(n_ICU_D_conf_progress) <- c(n_groups, n_strains, k_ICU_D, n_vacc_classes)
dim(n_ICU_D_unconf_to_conf) <- c(n_groups, n_strains, k_ICU_D, n_vacc_classes)

## Vectors handling the W_R class
dim(W_R_unconf) <- c(n_groups, n_strains, k_W_R, n_vacc_classes)
dim(aux_W_R_unconf) <-
  c(n_groups, n_strains, k_W_R, n_vacc_classes)
dim(new_W_R_unconf) <-
  c(n_groups, n_strains, k_W_R, n_vacc_classes)
dim(n_W_R_unconf_progress) <-
  c(n_groups, n_strains, k_W_R, n_vacc_classes)
dim(W_R_conf) <- c(n_groups, n_strains, k_W_R, n_vacc_classes)
dim(aux_W_R_conf) <-
  c(n_groups, n_strains, k_W_R, n_vacc_classes)
dim(new_W_R_conf) <-
  c(n_groups, n_strains, k_W_R, n_vacc_classes)
dim(n_W_R_conf_progress) <- c(n_groups, n_strains, k_W_R, n_vacc_classes)
dim(n_W_R_unconf_to_conf) <-
  c(n_groups, n_strains, k_W_R, n_vacc_classes)

## Vectors handling the W_D class
dim(W_D_unconf) <- c(n_groups, n_strains, k_W_D, n_vacc_classes)
dim(aux_W_D_unconf) <-
  c(n_groups, n_strains, k_W_D, n_vacc_classes)
dim(new_W_D_unconf) <-
  c(n_groups, n_strains, k_W_D, n_vacc_classes)
dim(n_W_D_unconf_progress) <-
  c(n_groups, n_strains, k_W_D, n_vacc_classes)
dim(W_D_conf) <-
  c(n_groups, n_strains, k_W_D, n_vacc_classes)
dim(aux_W_D_conf) <-
  c(n_groups, n_strains, k_W_D, n_vacc_classes)
dim(new_W_D_conf) <-
  c(n_groups, n_strains, k_W_D, n_vacc_classes)
dim(n_W_D_conf_progress) <- c(n_groups, n_strains, k_W_D, n_vacc_classes)
dim(n_W_D_unconf_to_conf) <-
  c(n_groups, n_strains, k_W_D, n_vacc_classes)

## Vectors handling the T_sero_pos class
dim(R) <- c(n_groups, n_strains, n_vacc_classes)
dim(new_R) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the T_sero_pre class and seroconversion
dim(T_sero_pre) <- c(n_groups, n_strains, 2, n_vacc_classes)
dim(new_T_sero_pre) <- c(n_groups, n_strains, 2, n_vacc_classes)
dim(n_T_sero_pre_progress) <- c(n_groups, n_strains, 2, n_vacc_classes)
dim(gamma_sero_pre) <- 2
dim(p_T_sero_pre_progress) <- c(n_groups, n_strains, 2, n_vacc_classes)
dim(p_sero_pos) <- n_groups

## Vectors handling the T_sero_pos class
dim(T_sero_pos) <- c(n_groups, n_strains, k_sero_pos, n_vacc_classes)
dim(n_T_sero_pos_progress) <- c(n_groups, n_strains, k_sero_pos, n_vacc_classes)
dim(new_T_sero_pos) <- c(n_groups, n_strains, k_sero_pos, n_vacc_classes)
dim(n_T_sero_pre_to_T_sero_pos) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the T_sero_neg class
dim(T_sero_neg) <- c(n_groups, n_strains, n_vacc_classes)
dim(new_T_sero_neg) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the D_hosp class
dim(D_hosp) <- c(n_groups, n_vacc_classes)
dim(delta_D_hosp) <- c(n_groups, n_vacc_classes)

## Vectors handling the D_non_hosp class
dim(D_non_hosp) <- n_groups
dim(delta_D_non_hosp) <- n_groups

## Vectors handling the PCR classes
dim(T_PCR_pre) <- c(n_groups, n_strains, k_PCR_pre, n_vacc_classes)
dim(n_T_PCR_pre_progress) <- c(n_groups, n_strains, k_PCR_pre, n_vacc_classes)
dim(new_T_PCR_pre) <- c(n_groups, n_strains, k_PCR_pre, n_vacc_classes)
dim(T_PCR_pos) <- c(n_groups, n_strains, k_PCR_pos, n_vacc_classes)
dim(n_T_PCR_pos_progress) <- c(n_groups, n_strains, k_PCR_pos, n_vacc_classes)
dim(new_T_PCR_pos) <- c(n_groups, n_strains, k_PCR_pos, n_vacc_classes)
dim(T_PCR_neg) <- c(n_groups, n_strains, n_vacc_classes)
dim(new_T_PCR_neg) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the S->S transitions i.e. moving between vaccination classes
dim(p_S_next_vacc_class) <- c(n_groups, n_vacc_classes)
dim(n_S_next_vacc_class) <- c(n_groups, n_vacc_classes)

dim(p_E_next_vacc_class) <- c(n_groups, n_strains, k_E, n_vacc_classes)
dim(n_E_next_vacc_class) <- c(n_groups, n_strains, k_E, n_vacc_classes)
dim(n_E_progress) <- c(n_groups, n_strains, k_E, n_vacc_classes)
dim(n_EE_next_vacc_class) <- c(n_groups, n_strains, k_E, n_vacc_classes)

dim(p_I_A_next_vacc_class) <-
  c(n_groups, n_strains, k_A, n_vacc_classes)
dim(n_I_A_next_vacc_class) <-
  c(n_groups, n_strains, k_A, n_vacc_classes)
dim(n_I_A_progress) <- c(n_groups, n_strains, k_A, n_vacc_classes)
dim(n_II_A_next_vacc_class) <-
  c(n_groups, n_strains, k_A, n_vacc_classes)

dim(p_I_P_next_vacc_class) <-
  c(n_groups, n_strains, k_P, n_vacc_classes)
dim(n_I_P_next_vacc_class) <-
  c(n_groups, n_strains, k_P, n_vacc_classes)
dim(n_I_P_progress) <- c(n_groups, n_strains, k_P, n_vacc_classes)
dim(n_II_P_next_vacc_class) <-
  c(n_groups, n_strains, k_P, n_vacc_classes)

dim(p_R_next_vacc_class) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_R_next_vacc_class) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_R_next_vacc_class_capped) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_R_next_vacc_class_tmp) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_R_progress) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_R_progress_capped) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_R_progress_tmp) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_RS_next_vacc_class) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the S->E transition where infected are split
## between level of infectivity
dim(p_SE) <- c(n_groups, n_vacc_classes)
dim(n_SE) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_S_progress) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_S_progress_tot) <- c(n_groups, n_vacc_classes)
dim(n_SE_next_vacc_class) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the E->I transition where newly infectious cases
## are split between level of severity
dim(n_EI_A) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_EI_P) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_EI_A_next_vacc_class) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_EI_P_next_vacc_class) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling I_C_2 to R, G_D transition
dim(n_I_C_2_to_G_D) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_2_to_R) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling number of new hospitalisations, ICU admissions and
## recoveries in hospital
dim(n_I_C_2_to_hosp) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_2_to_ICU_pre) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_2_to_ICU_pre_conf) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_hosp_non_ICU) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_2_to_H_D) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_2_to_H_D_conf) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_2_to_H_R) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_2_to_H_R_conf) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_ICU_pre_unconf_to_ICU_D_unconf) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_ICU_pre_conf_to_ICU_D_conf) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_ICU_pre_unconf_to_ICU_W_R_unconf) <-
  c(n_groups, n_strains, n_vacc_classes)
dim(n_ICU_pre_conf_to_ICU_W_R_conf) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_ICU_pre_unconf_to_ICU_W_D_unconf) <-
  c(n_groups, n_strains, n_vacc_classes)
dim(n_ICU_pre_conf_to_ICU_W_D_conf) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the serology flow
dim(n_com_to_T_sero_pre) <- c(n_groups, n_strains, 2, n_vacc_classes)

## Vectors handling the severity profile
dim(p_C) <- n_groups

## Vectors handling the potential death in hospital (general beds and ICU)
dim(p_H_D_by_age) <- c(n_groups, n_strains)
dim(psi_H_D) <- c(n_groups, n_strains)
dim(p_ICU_D_by_age) <- c(n_groups, n_strains)
dim(psi_ICU_D) <- c(n_groups, n_strains)
dim(p_W_D_by_age) <- c(n_groups, n_strains)
dim(psi_W_D) <- c(n_groups, n_strains)

## Vector handling the probability of being admitted as confirmed
dim(p_star_by_age) <- n_groups
dim(psi_star) <- n_groups

dim(cum_admit_by_age) <- n_groups

## Vectors handling the age specific heterogeneous transmission process
dim(lambda) <- c(n_groups, n_strains)
dim(s_ij) <- c(n_groups, n_groups, n_strains)
dim(m) <- c(n_groups, n_groups)
dim(I_with_diff_trans) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the loss of immunity
dim(n_RS) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_RS) <- n_groups

## Total population
initial(N_tot[]) <- 0
update(N_tot[]) <- sum(S[i, ]) + sum(R[i, , ]) + sum(D_hosp[i, ]) + sum(E[i, , , ]) +
  sum(I_A[i, , , ]) + sum(I_P[i, , , ]) +
  sum(I_C_1[i, , , ]) + sum(I_C_2[i, , , ]) +
  sum(ICU_pre_conf[i, , , ]) + sum(ICU_pre_unconf[i, , , ])  +
  sum(H_R_conf[i, , , ]) + sum(H_R_unconf[i, , , ]) +
  sum(H_D_conf[i, , , ]) + sum(H_D_unconf[i, , , ]) +
  sum(ICU_W_R_conf[i, , , ]) + sum(ICU_W_R_unconf[i, , , ]) +
  sum(ICU_W_D_conf[i, , , ]) + sum(ICU_W_D_unconf[i, , , ]) +
  sum(ICU_D_conf[i, , , ]) + sum(ICU_D_unconf[i, , , ]) +
  sum(W_R_conf[i, , , ]) + sum(W_R_unconf[i, , , ]) +
  sum(W_D_conf[i, , , ]) + sum(W_D_unconf[i, , , ]) +
  sum(G_D[i, , , ]) + D_non_hosp[i]
dim(N_tot) <- n_groups

## Total population calculated with seroconversion flow
initial(N_tot2) <- 0
update(N_tot2) <- sum(S) + sum(T_sero_pre) +
  sum(T_sero_pos) + sum(T_sero_neg) + sum(E)

## Total population calculated with PCR flow
initial(N_tot3) <- 0
update(N_tot3) <- sum(S) + sum(T_PCR_pre) + sum(T_PCR_pos) + sum(T_PCR_neg)

## Aggregate our reporting statistics by summing across age (simple
## for everything except for seropositivity data, done last)
initial(ICU_tot) <- 0
new_ICU_tot <- sum(new_ICU_W_R_conf) + sum(new_ICU_W_D_conf) +
  sum(new_ICU_D_conf)
update(ICU_tot) <- new_ICU_tot

initial(general_tot) <- 0
new_general_tot <- sum(new_ICU_pre_conf) + sum(new_H_R_conf) +
  sum(new_H_D_conf) + sum(new_W_R_conf) +
  sum(new_W_D_conf)
update(general_tot) <- new_general_tot

initial(hosp_tot) <- 0
update(hosp_tot) <- new_ICU_tot + new_general_tot

initial(D_hosp_tot) <- 0
delta_D_hosp_tot <- sum(delta_D_hosp)
update(D_hosp_tot) <- D_hosp_tot + delta_D_hosp_tot

## community deaths are non-hospital deaths in groups 1 to 18
initial(D_comm_tot) <- 0
delta_D_comm_tot <- sum(delta_D_non_hosp[1:18])
update(D_comm_tot) <- D_comm_tot + delta_D_comm_tot

initial(D_comm_inc) <- 0
update(D_comm_inc) <- if (step %% steps_per_day == 0)
  delta_D_comm_tot else D_comm_inc + delta_D_comm_tot

## carehome deaths are non-hospital deaths in group 19
initial(D_carehomes_tot) <- 0
delta_D_carehomes_tot <- delta_D_non_hosp[19]
update(D_carehomes_tot) <- D_carehomes_tot + delta_D_carehomes_tot

initial(D_carehomes_inc) <- 0
update(D_carehomes_inc) <- if (step %% steps_per_day == 0)
  delta_D_carehomes_tot else D_carehomes_inc + delta_D_carehomes_tot

initial(D_hosp_inc) <- 0
update(D_hosp_inc) <- if (step %% steps_per_day == 0)
  delta_D_hosp_tot else D_hosp_inc + delta_D_hosp_tot

initial(D_tot) <- 0
update(D_tot) <- D_tot + delta_D_hosp_tot + delta_D_comm_tot +
  delta_D_carehomes_tot

## Our age groups for serology are fixed: we break them down into the
##
## * 0-14 (1, 2, 3)
## * 15-64 (4, 5, ..., 13)
## * 65-100 (14, 15, ..., 17)
##
## NOTE: this excludes CHW (18) and CHR (19) but we probably should
## sum in CHW into the figures here.
##
## To fit with the data currently available, we currently only
## consider the middle group, though this could be expanded easily by
## more statements like the ones below.
##
initial(sero_pos) <- 0
update(sero_pos) <- sum(new_T_sero_pos[4:13, , , ])

initial(cum_sympt_cases) <- 0
new_sympt_cases <- sum(n_EI_P) + sum(n_EI_P_next_vacc_class)
update(cum_sympt_cases) <- cum_sympt_cases + new_sympt_cases

## only over 25s (exclude groups 1 to 5)
initial(cum_sympt_cases_over25) <- 0
new_sympt_cases_over25 <- sum(n_EI_P[6:n_groups, , ]) +
  sum(n_EI_P_next_vacc_class[6:n_groups, , ])
update(cum_sympt_cases_over25) <- cum_sympt_cases_over25 +
  new_sympt_cases_over25

initial(cum_sympt_cases_non_variant_over25) <- 0
new_sympt_cases_non_variant_over25 <- sum(n_EI_P[6:n_groups, 1, ]) +
  sum(n_EI_P_next_vacc_class[6:n_groups, 1, ])
update(cum_sympt_cases_non_variant_over25) <-
  cum_sympt_cases_non_variant_over25 + new_sympt_cases_non_variant_over25

## And incidence:
initial(sympt_cases_inc) <- 0
update(sympt_cases_inc) <- (
  if (step %% steps_per_day == 0) new_sympt_cases
  else sympt_cases_inc + new_sympt_cases)

initial(sympt_cases_over25_inc) <- 0
update(sympt_cases_over25_inc) <- (
  if (step %% steps_per_day == 0) new_sympt_cases_over25
  else sympt_cases_over25_inc + new_sympt_cases_over25)

initial(sympt_cases_non_variant_over25_inc) <- 0
update(sympt_cases_non_variant_over25_inc) <- (
  if (step %% steps_per_day == 0) new_sympt_cases_non_variant_over25
  else sympt_cases_non_variant_over25_inc + new_sympt_cases_non_variant_over25)

## For REACT we exclude the 0-4 (1) and CHR (19) groups
initial(react_pos) <- 0
update(react_pos) <- sum(new_T_PCR_pos[2:18, , , ])


## prob_strain gives probability of an infection in group i being of strain j
initial(prob_strain[, ]) <- 0
dim(prob_strain) <- c(n_groups, n_strains)
update(prob_strain[, ]) <- lambda[i, j] / sum(lambda[i, ])


## I_weighted used in IFR calculation
initial(I_weighted[, ]) <- 0
dim(I_weighted) <- c(n_groups, n_vacc_classes)
dim(I_weighted_strain) <- c(n_groups, n_strains, n_vacc_classes)
I_weighted_strain[, , ] <-
  strain_transmission[j] * (
    I_A_transmission * sum(new_I_A[i, j, , k]) +
      I_P_transmission * sum(new_I_P[i, j, , k]) +
      I_C_1_transmission * sum(new_I_C_1[i, j, , k]) +
      I_C_2_transmission * sum(new_I_C_2[i, j, , k]) +
      hosp_transmission * (
        sum(new_ICU_pre_unconf[i, j, , k]) +
          sum(new_ICU_pre_conf[i, j, , k]) +
          sum(new_H_R_unconf[i, j, , k]) +
          sum(new_H_R_conf[i, j, , k]) +
          sum(new_H_D_unconf[i, j, , k]) +
          sum(new_H_D_conf[i, j, , k])) +
      ICU_transmission * (
        sum(new_ICU_W_R_unconf[i, j, , k]) +
          sum(new_ICU_W_R_conf[i, j, , k]) +
          sum(new_ICU_W_D_unconf[i, j, , k]) +
          sum(new_ICU_W_D_conf[i, j, , k]) +
          sum(new_ICU_D_unconf[i, j, , k]) +
          sum(new_ICU_D_conf[i, j, , k])) +
      G_D_transmission * sum(new_G_D[i, j, , k]))
update(I_weighted[, ]) <- sum(I_weighted_strain[i, , j])

## Vaccination engine
n_doses <- 2
index_dose[] <- user(integer = TRUE)
dim(index_dose) <- n_doses

vaccine_dose_step[, , ] <- user() # n_groups, n_doses, n_time
dim(vaccine_dose_step) <- user()

## First, the number of candidates
vaccine_n_candidates[, ] <-
  S[i, index_dose[j]] +
  sum(E[i, , , index_dose[j]]) +
  sum(I_A[i, , , index_dose[j]]) +
  sum(I_P[i, , , index_dose[j]]) +
  sum(R[i, , index_dose[j]])
dim(vaccine_n_candidates) <- c(n_groups, n_doses)

## Work out the vaccination probability via doses, driven by the
## schedule
vaccine_probability_doses[, ] <- min(
  vaccine_attempted_doses[i, j] / vaccine_n_candidates[i, j],
  as.numeric(1))
dim(vaccine_probability_doses) <- c(n_groups, n_doses)

vaccine_attempted_doses[, ] <- vaccine_missed_doses[i, j] + (
  if (as.integer(step) >= dim(vaccine_dose_step, 3) ||
      vaccine_n_candidates[i, j] == 0) 0
  else vaccine_dose_step[i, j, step + 1])
dim(vaccine_attempted_doses) <- c(n_groups, n_doses)

initial(vaccine_missed_doses[, ]) <- 0
update(vaccine_missed_doses[, ]) <-
  vaccine_catchup_fraction *
  max(vaccine_attempted_doses[i, j] - n_vaccinated[i, index_dose[j]],
      as.numeric(0))
dim(vaccine_missed_doses) <- c(n_groups, n_doses)

vaccine_catchup_fraction <- user(0)


## TODO: There's a divide-by-zero here causing NaNs in
## vaccine_probability. Worth fixing

## Then fix everything based on progression at a constant rate (will
## be zero for the cases that have probabilities above)
vaccine_probability[, ] <-
  1 - exp(-vaccine_progression_rate_base[i, j] * dt)
dim(vaccine_probability) <- c(n_groups, n_vacc_classes)

## This can't be automatically driven from the number of doses, so we
## have to unroll it here and write both out manually. This is the
## reason why n_doses is fixed as 2 rather than being user-supplied.
vaccine_probability[, index_dose[1]] <- vaccine_probability_doses[i, 1]
vaccine_probability[, index_dose[2]] <- vaccine_probability_doses[i, 2]

initial(tmp_vaccine_n_candidates[, ]) <- 0
update(tmp_vaccine_n_candidates[, ]) <- vaccine_n_candidates[i, j]
dim(tmp_vaccine_n_candidates) <- c(n_groups, n_doses)

initial(tmp_vaccine_probability[, ]) <- 0
update(tmp_vaccine_probability[, ]) <- vaccine_probability[i, j]
dim(tmp_vaccine_probability) <- c(n_groups, n_vacc_classes)

config(compare) <- "compare_carehomes.cpp"
## Parameters and code to support the compare function. Because these
## do not appear in any odin equation we mark them as "ignore.unused"
## so that odin doesn't complain that they appear redundant. This
## means that unwanted things can accumulate here so keep on top of
## it.
N_tot_all <- user() # ignore.unused
N_tot_over25 <- user() # ignore.unused
N_tot_react <- user() # ignore.unused
N_tot_15_64 <- user() # ignore.unused

p_NC <- user() # ignore.unused
pillar2_sensitivity <- user() # ignore.unused
pillar2_specificity <- user() # ignore.unused
react_sensitivity <- user() # ignore.unused
react_specificity <- user() # ignore.unused
sero_sensitivity <- user() # ignore.unused
sero_specificity <- user() # ignore.unused
exp_noise <- user() # ignore.unused
phi_ICU <- user() # ignore.unused
kappa_ICU <- user() # ignore.unused
phi_general <- user() # ignore.unused
kappa_general <- user() # ignore.unused
phi_hosp <- user() # ignore.unused
kappa_hosp <- user() # ignore.unused
phi_death_hosp <- user() # ignore.unused
kappa_death_hosp <- user() # ignore.unused
phi_death_comm <- user() # ignore.unused
kappa_death_comm <- user() # ignore.unused
phi_death_carehomes <- user() # ignore.unused
kappa_death_carehomes <- user() # ignore.unused
kappa_death <- user() # ignore.unused
kappa_death_non_hosp <- user() # ignore.unused
phi_admitted <- user() # ignore.unused
kappa_admitted <- user() # ignore.unused
phi_diagnoses <- user() # ignore.unused
kappa_diagnoses <- user() # ignore.unused
phi_all_admission <- user() # ignore.unused
kappa_all_admission <- user() # ignore.unused
rho_pillar2_tests <- user() # ignore.unused
phi_pillar2_cases <- user() # ignore.unused
kappa_pillar2_cases <- user() # ignore.unused
