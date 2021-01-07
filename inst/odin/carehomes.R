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
dt <- user()
initial(time) <- 0
update(time) <- (step + 1) * dt

## output number of individuals vaccinated by age and vaccine stage
## For example, for E, we sum over n_E_next_vacc_class (those moving vaccine
## stage without progressing disease stages) and also n_EE_next_vacc_class
## (those moving vaccine stage and also progressing disease stages)
## vaccinated S
initial(cum_n_S_vaccinated[, ]) <- 0
update(cum_n_S_vaccinated[, ]) <- cum_n_S_vaccinated[i, j] +
  n_S_next_vacc_class[i, j] + sum(n_SE_next_vacc_class[i, , j])
dim(cum_n_S_vaccinated) <- c(n_groups, n_vacc_classes)
## vaccinated E
initial(cum_n_E_vaccinated[, ]) <- 0
update(cum_n_E_vaccinated[, ]) <- cum_n_E_vaccinated[i, j] +
  sum(n_E_next_vacc_class[i, , , j]) + sum(n_EE_next_vacc_class[i, , , j])
dim(cum_n_E_vaccinated) <- c(n_groups, n_vacc_classes)
## vaccinated I_A
initial(cum_n_I_A_vaccinated[, ]) <- 0
update(cum_n_I_A_vaccinated[, ]) <- cum_n_I_A_vaccinated[i, j] +
  sum(n_I_A_next_vacc_class[i, , , j]) +
  sum(n_II_A_next_vacc_class[i, , , j])
dim(cum_n_I_A_vaccinated) <- c(n_groups, n_vacc_classes)
## vaccinated R
initial(cum_n_R_vaccinated[, ]) <- 0
update(cum_n_R_vaccinated[, ]) <- cum_n_R_vaccinated[i, j] +
  sum(n_R_next_vacc_class[i, , j]) + sum(n_RS_next_vacc_class[i, , j])
dim(cum_n_R_vaccinated) <- c(n_groups, n_vacc_classes)

## Total number of vaccinations over S, E, I_asypmt, R for convenience
initial(cum_n_vaccinated[, ]) <- 0
update(cum_n_vaccinated[, ]) <-
  cum_n_S_vaccinated[i, j] +
  cum_n_E_vaccinated[i, j] +
  cum_n_I_A_vaccinated[i, j] +
  cum_n_R_vaccinated[i, j]
dim(cum_n_vaccinated) <- c(n_groups, n_vacc_classes)

## Core equations for transitions between compartments:
update(S[, ]) <- new_S[i, j]
update(E[, , , ]) <- new_E[i, j, k, l]
update(I_A[, , , ]) <- new_I_A[i, j, k, l]
update(I_C[, , , ]) <- new_I_C[i, j, k, l]
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
update(R_pre[, , , ]) <- new_R_pre[i, j, k, l]
update(R_pos[, , , ]) <- new_R_pos[i, j, k, l]
update(T_sero_neg[, , ]) <- new_T_sero_neg[i, j, k]
update(R[, , ]) <- new_R[i, j, k]
update(D_hosp[]) <- new_D_hosp[i]
update(D_comm[]) <- new_D_comm[i]
update(T_PCR_pre[, , , ]) <- new_T_PCR_pre[i, j, k, l]
update(T_PCR_pos[, , , ]) <- new_T_PCR_pos[i, j, k, l]
update(T_PCR_neg[, , ]) <- new_T_PCR_neg[i, j, k]
update(cum_admit_conf) <-
  cum_admit_conf +
  sum(n_I_C_to_H_D_conf) +
  sum(n_I_C_to_H_R_conf) +
  sum(n_I_C_to_ICU_pre_conf)
update(cum_new_conf) <-
  cum_new_conf +
  sum(n_H_D_unconf_to_conf) +
  sum(n_H_R_unconf_to_conf) +
  sum(n_ICU_pre_unconf_to_conf) +
  sum(n_ICU_D_unconf_to_conf) +
  sum(n_ICU_W_R_unconf_to_conf) +
  sum(n_ICU_W_D_unconf_to_conf) +
  sum(n_W_R_unconf_to_conf) +
  sum(n_W_D_unconf_to_conf)
update(cum_admit_by_age[]) <- cum_admit_by_age[i] + sum(n_I_C_to_hosp[i, , ])

## Individual probabilities of transition:

## vaccination
p_S_next_vacc_class[, ] <- 1 - exp(-vaccine_progression_rate[i, j] * dt)
p_E_next_vacc_class[, , , ] <- 1 - exp(-vaccine_progression_rate[i, l] * dt)
p_I_A_next_vacc_class[, , , ] <-
  1 - exp(-vaccine_progression_rate[i, l] * dt)
p_R_next_vacc_class[, , ] <-
  1 - exp(-vaccine_progression_rate[i, k] * dt)
## clinical progression
p_SE[, ] <- 1 - exp(-sum(lambda[i, ]) *
                      rel_susceptibility[i, j] * dt) # S to I age/vacc dependent
p_EE <- 1 - exp(-gamma_E * dt) # progression of latent period
p_II_A <- 1 - exp(-gamma_A * dt) # progression of infectious period
p_II_C <- 1 - exp(-gamma_C * dt)
p_GG_D <- 1 - exp(-gamma_G_D * dt)
p_IICU_pre <- 1 - exp(-gamma_ICU_pre * dt)
p_HH_R <- 1 - exp(-gamma_H_R * dt)
p_HH_D <- 1 - exp(-gamma_H_D * dt)
p_IICU_W_R <- 1 - exp(-gamma_ICU_W_R * dt)
p_IICU_W_D <- 1 - exp(-gamma_ICU_W_D * dt)
p_IICU_D <- 1 - exp(-gamma_ICU_D * dt)
p_WW_R <- 1 - exp(-gamma_W_R * dt)
p_WW_D <- 1 - exp(-gamma_W_D * dt)
p_R_pre[, , , ] <- 1 - exp(-gamma_sero_pre[k] * dt)
p_R_pos <- 1 - exp(-gamma_sero_pos * dt)
p_test <- 1 - exp(-gamma_U * dt)
p_T_PCR_pre <- 1 - exp(-gamma_PCR_pre * dt)
p_T_PCR_pos <- 1 - exp(-gamma_PCR_pos * dt)
p_RS[] <- 1 - exp(-waning_rate[i] * dt) # R to S age dependent

## Work out time-varying probabilities
p_ICU <- if (as.integer(step) >= length(p_ICU_step))
  p_ICU_step[length(p_ICU_step)] else p_ICU_step[step + 1]
prob_ICU[] <- p_ICU * psi_ICU[i]

p_H <- if (as.integer(step) >= length(p_H_step))
  p_H_step[length(p_H_step)] else p_H_step[step + 1]
prob_H[] <- p_H * psi_H[i]

p_ICU_D <- if (as.integer(step) >= length(p_ICU_D_step))
  p_ICU_D_step[length(p_ICU_D_step)] else p_ICU_D_step[step + 1]
prob_ICU_D[] <- p_ICU_D * psi_ICU_D[i]

p_H_D <- if (as.integer(step) >= length(p_H_D_step))
  p_H_D_step[length(p_H_D_step)] else
    p_H_D_step[step + 1]
prob_H_D[] <- p_H_D * psi_H_D[i]

p_W_D <- if (as.integer(step) >= length(p_W_D_step))
  p_W_D_step[length(p_W_D_step)] else
    p_W_D_step[step + 1]
prob_W_D[] <- p_W_D * psi_W_D[i]

p_G_D <- if (as.integer(step) >= length(p_G_D_step))
  p_G_D_step[length(p_G_D_step)] else p_G_D_step[step + 1]
prob_G_D[] <- p_G_D * psi_G_D[i]

p_star <- if (as.integer(step) >= length(p_star_step))
  p_star_step[length(p_star_step)] else p_star_step[step + 1]
prob_admit_conf[] <- p_star * psi_star[i]

## Draws from binomial distributions for numbers changing between
## compartments:

## modelling infections and vaccine progression, which can happen simultaneously

#### flow out of S ####

## new infections

## Compute the new infections with multiple strains using nested binomials
n_S_progress_tot[, ] <- rbinom(S[i, j], p_SE[i, j])
n_S_progress[, 1, ] <-
  rbinom(n_S_progress_tot[i, k], lambda[i, 1] / sum(lambda[i, ]))
n_S_progress[, 2:n_strains, ] <-
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
strain_seed <- (if (as.integer(step) >= length(strain_seed_step))
                  strain_seed_step[length(strain_seed_step)]
                else strain_seed_step[step + 1])
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

n_E_progress[, , , ] <- rbinom(E[i, j, k, l], p_EE)
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

n_I_A_progress[, , , ] <- rbinom(I_A[i, j, k, l], p_II_A)
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


#### flow out of R ####

n_R_progress_tmp[, , ] <- rbinom(R[i, j, k], p_RS[i])
## cap on people who can move out of R based on numbers in T_sero_neg and T_PCR_neg
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
n_R_next_vacc_class_capped[, , ] <- min(n_R_next_vacc_class_tmp[i, j, k],
  T_sero_neg[i, j, k] - n_R_progress[i, j, k],
  T_PCR_neg[i, j, k] - n_R_progress[i, j, k])
n_R_next_vacc_class[, , ] <- if (model_pcr_and_serology == 1)
  n_R_next_vacc_class_capped[i, j, k] else n_R_next_vacc_class_tmp[i, j, k]
  
#### other transitions ####

n_II_C[, , , ] <- rbinom(I_C[i, j, k, l], p_II_C)
n_GG_D[, , , ] <- rbinom(G_D[i, j, k, l], p_GG_D)
n_ICU_pre_unconf[, , , ] <- rbinom(ICU_pre_unconf[i, j, k, l], p_IICU_pre)
n_ICU_pre_conf[, , , ] <- rbinom(ICU_pre_conf[i, j, k, l], p_IICU_pre)
n_HH_R_unconf[, , , ] <- rbinom(H_R_unconf[i, j, k, l], p_HH_R)
n_HH_R_conf[, , , ] <- rbinom(H_R_conf[i, j, k, l], p_HH_R)
n_HH_D_unconf[, , , ] <- rbinom(H_D_unconf[i, j, k, l], p_HH_D)
n_HH_D_conf[, , , ] <- rbinom(H_D_conf[i, j, k, l], p_HH_D)
n_IICU_W_R_unconf[, , , ] <-
  rbinom(ICU_W_R_unconf[i, j, k, l], p_IICU_W_R)
n_IICU_W_R_conf[, , , ] <- rbinom(ICU_W_R_conf[i, j, k, l], p_IICU_W_R)
n_IICU_W_D_unconf[, , , ] <-
  rbinom(ICU_W_D_unconf[i, j, k, l], p_IICU_W_D)
n_IICU_W_D_conf[, , , ] <- rbinom(ICU_W_D_conf[i, j, k, l], p_IICU_W_D)
n_IICU_D_unconf[, , , ] <- rbinom(ICU_D_unconf[i, j, k, l], p_IICU_D)
n_IICU_D_conf[, , , ] <- rbinom(ICU_D_conf[i, j, k, l], p_IICU_D)
n_WW_R_unconf[, , , ] <- rbinom(W_R_unconf[i, j, k, l], p_WW_R)
n_WW_R_conf[, , , ] <- rbinom(W_R_conf[i, j, k, l], p_WW_R)
n_WW_D_unconf[, , , ] <- rbinom(W_D_unconf[i, j, k, l], p_WW_D)
n_WW_D_conf[, , , ] <- rbinom(W_D_conf[i, j, k, l], p_WW_D)
n_R_pre[, , , ] <- rbinom(R_pre[i, j, k, l], p_R_pre[i, j, k, l])
n_R_pos[, , , ] <- rbinom(R_pos[i, j, k, l], p_R_pos)
n_T_PCR_pre[, , , ] <- rbinom(T_PCR_pre[i, j, k, l], p_T_PCR_pre)
n_T_PCR_pos[, , , ] <- rbinom(T_PCR_pos[i, j, k, l], p_T_PCR_pos)

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
new_S[, 1] <- new_S[i, 1] + n_S_next_vacc_class[i, n_vacc_classes] +
  sum(n_RS_next_vacc_class[i, , n_vacc_classes])
new_S[, 2:n_vacc_classes] <- new_S[i, j] + n_S_next_vacc_class[i, j - 1] +
  sum(n_RS_next_vacc_class[i, , j - 1])


## Computes the number of asymptomatic
n_EI_A[, , ] <- rbinom(n_EE[i, j, s_E, k],
                          1 - p_C[i] * rel_p_C[i, k])
n_EI_A_next_vacc_class[, , ] <-
  rbinom(n_EE_next_vacc_class[i, j, s_E, k],
         1 - p_C[i] * rel_p_C[i, k])

## Computes the number of symptomatic cases
n_EI_C[, , ] <- n_EE[i, j, s_E, k] - n_EI_A[i, j, k]
n_EI_C_next_vacc_class[, , ] <- n_EE_next_vacc_class[i, j, s_E, k] -
  n_EI_A_next_vacc_class[i, j, k]

## Work out the S->E and E->E transitions
aux_EE[, , 1, ] <- n_SE[i, j, l]
aux_EE[, , 2:s_E, ] <- n_EE[i, j, k - 1, l]
aux_EE[, , , ] <- aux_EE[i, j, k, l] - n_EE[i, j, k, l] -
  n_EE_next_vacc_class[i, j, k, l] -
  n_E_next_vacc_class[i, j, k, l]
aux_EE[, , , 1] <- aux_EE[i, j, k, 1] +
  n_E_next_vacc_class[i, j, k, n_vacc_classes]
aux_EE[, , , 2:n_vacc_classes] <- aux_EE[i, j, k, l] +
  n_E_next_vacc_class[i, j, k, l - 1]
aux_EE[, , 1, 1] <- aux_EE[i, j, 1, 1] +
  n_SE_next_vacc_class[i, j, n_vacc_classes]
aux_EE[, , 1, 2:n_vacc_classes] <- aux_EE[i, j, k, l] +
  n_SE_next_vacc_class[i, j, l - 1]
aux_EE[, , 2:s_E, 1] <- aux_EE[i, j, k, l] +
  n_EE_next_vacc_class[i, j, k - 1, n_vacc_classes]
aux_EE[, , 2:s_E, 2:n_vacc_classes] <- aux_EE[i, j, k, l] +
  n_EE_next_vacc_class[i, j, k - 1, l - 1]
new_E[, , , ] <- E[i, j, k, l] + aux_EE[i, j, k, l]

## Work out the I_A->I_A transitions
aux_II_A[, , 1, ] <- n_EI_A[i, j, l]
aux_II_A[, , 2:s_A, ] <- n_II_A[i, j, k - 1, l]
aux_II_A[, , , ] <- aux_II_A[i, j, k, l] - n_II_A[i, j, k, l] -
  n_II_A_next_vacc_class[i, j, k, l] -
  n_I_A_next_vacc_class[i, j, k, l]
aux_II_A[, , , 1] <- aux_II_A[i, j, k, l]  +
  n_I_A_next_vacc_class[i, j, 1, n_vacc_classes]
aux_II_A[, , , 2:n_vacc_classes] <- aux_II_A[i, j, k, l] +
  n_I_A_next_vacc_class[i, j, k, l - 1]
aux_II_A[, , 1, 1] <- aux_II_A[i, j, k, l] +
  n_EI_A_next_vacc_class[i, j, n_vacc_classes]
aux_II_A[, , 1, 2:n_vacc_classes] <- aux_II_A[i, j, k, l] +
  n_EI_A_next_vacc_class[i, j, l - 1]
aux_II_A[, , 2:s_A, 1] <- aux_II_A[i, j, k, l] +
  n_II_A_next_vacc_class[i, j, k - 1, n_vacc_classes]
aux_II_A[, , 2:s_A, 2:n_vacc_classes] <- aux_II_A[i, j, k, l] +
  n_II_A_next_vacc_class[i, j, k - 1, l - 1]
new_I_A[, , , ] <- I_A[i, j, k, l] + aux_II_A[i, j, k, l]

## Work out the I_C->I_C transitions
aux_II_C[, , 1, 1] <- n_EI_C[i, j, 1] +
  n_EI_C_next_vacc_class[i, j, n_vacc_classes]
aux_II_C[, , 1, 2:n_vacc_classes] <-
  n_EI_C[i, j, l] + n_EI_C_next_vacc_class[i, j, l - 1]

aux_II_C[, , 2:s_C, ] <- n_II_C[i, j, k - 1, l]
aux_II_C[, , 1:s_C, ] <-
  aux_II_C[i, j, k, l] - n_II_C[i, j, k, l]
new_I_C[, , , ] <- I_C[i, j, k, l] + aux_II_C[i, j, k, l]

## Work out the flow from I_C -> R, G_D, hosp
n_I_C_to_R[, , ] <- rbinom(n_II_C[i, j, s_C, k],
                         1 - prob_H[i] * rel_p_H[i, k])
n_I_C_to_G_D[, , ] <-
  rbinom(n_II_C[i, j, s_C, k] - n_I_C_to_R[i, j, k],
         prob_G_D[i])
n_I_C_to_hosp[, , ] <- n_II_C[i, j, s_C, k] - n_I_C_to_R[i, j, k] -
  n_I_C_to_G_D[i, j, k]

## Work out the G_D -> G_D transitions
aux_GG_D[, , 1, ] <- n_I_C_to_G_D[i, j, l]
aux_GG_D[, , 2:s_G_D, ] <- n_GG_D[i, j, k - 1, l]
aux_GG_D[, , 1:s_G_D, ] <-
  aux_GG_D[i, j, k, l] - n_GG_D[i, j, k, l]
new_G_D[, , , ] <- G_D[i, j, k, l] + aux_GG_D[i, j, k, l]

## Work out the split in hospitals between H_D, H_R and ICU_pre
n_I_C_to_ICU_pre[, , ] <- rbinom(n_I_C_to_hosp[i, j, k], prob_ICU[i])
n_I_C_to_ICU_pre_conf[, , ] <- rbinom(n_I_C_to_ICU_pre[i, j, k],
                                   prob_admit_conf[i])
n_hosp_non_ICU[, , ] <- n_I_C_to_hosp[i, j, k] - n_I_C_to_ICU_pre[i, j, k]
n_I_C_to_H_D[, , ] <- rbinom(n_hosp_non_ICU[i, j, k], prob_H_D[i])
n_I_C_to_H_D_conf[, , ] <- rbinom(n_I_C_to_H_D[i, j, k],
                                     prob_admit_conf[i])
n_I_C_to_H_R[, , ] <- n_hosp_non_ICU[i, j, k] - n_I_C_to_H_D[i, j, k]
n_I_C_to_H_R_conf[, , ] <- rbinom(n_I_C_to_H_R[i, j, k],
                                     prob_admit_conf[i])

## Work out the ICU_pre -> ICU_pre transitions
aux_ICU_pre_unconf[, , , ] <- ICU_pre_unconf[i, j, k, l]
aux_ICU_pre_unconf[, , 2:s_ICU_pre, ] <-
  aux_ICU_pre_unconf[i, j, k, l] + n_ICU_pre_unconf[i, j, k - 1, l]
aux_ICU_pre_unconf[, , 1:s_ICU_pre, ] <-
  aux_ICU_pre_unconf[i, j, k, l] - n_ICU_pre_unconf[i, j, k, l]
aux_ICU_pre_conf[, , , ] <-
  ICU_pre_conf[i, j, k, l]
aux_ICU_pre_conf[, , 2:s_ICU_pre, ] <-
  aux_ICU_pre_conf[i, j, k, l] + n_ICU_pre_conf[i, j, k - 1, l]
aux_ICU_pre_conf[, , 1:s_ICU_pre, ] <-
  aux_ICU_pre_conf[i, j, k, l] - n_ICU_pre_conf[i, j, k, l]
n_ICU_pre_unconf_to_conf[, , , ] <-
  rbinom(aux_ICU_pre_unconf[i, j, k, l], p_test)
new_ICU_pre_unconf[, , , ] <-
  aux_ICU_pre_unconf[i, j, k, l] - n_ICU_pre_unconf_to_conf[i, j, k, l]
new_ICU_pre_unconf[, , 1, ] <-
  new_ICU_pre_unconf[i, j, 1, l] + n_I_C_to_ICU_pre[i, j, l] -
  n_I_C_to_ICU_pre_conf[i, j, l]
new_ICU_pre_conf[, , , ] <-
  aux_ICU_pre_conf[i, j, k, l] + n_ICU_pre_unconf_to_conf[i, j, k, l]
new_ICU_pre_conf[, , 1, ] <-
  new_ICU_pre_conf[i, j, 1, l] + n_I_C_to_ICU_pre_conf[i, j, l]

## Work out the H_R->H_R transitions
aux_H_R_unconf[, , , ] <- H_R_unconf[i, j, k, l]
aux_H_R_unconf[, , 2:s_H_R, ] <-
  aux_H_R_unconf[i, j, k, l] + n_HH_R_unconf[i, j, k - 1, l]
aux_H_R_unconf[, , 1:s_H_R, ] <-
  aux_H_R_unconf[i, j, k, l] - n_HH_R_unconf[i, j, k, l]
aux_H_R_conf[, , , ] <- H_R_conf[i, j, k, l]
aux_H_R_conf[, , 2:s_H_R, ] <-
  aux_H_R_conf[i, j, k, l] + n_HH_R_conf[i, j, k - 1, l]
aux_H_R_conf[, , 1:s_H_R, ] <-
  aux_H_R_conf[i, j, k, l] - n_HH_R_conf[i, j, k, l]
n_H_R_unconf_to_conf[, , , ] <-
  rbinom(aux_H_R_unconf[i, j, k, l], p_test)
new_H_R_unconf[, , , ] <-
  aux_H_R_unconf[i, j, k, l] - n_H_R_unconf_to_conf[i, j, k, l]
new_H_R_unconf[, , 1, ] <-
  new_H_R_unconf[i, j, 1, l] + n_I_C_to_H_R[i, j, l] -
  n_I_C_to_H_R_conf[i, j, l]
new_H_R_conf[, , , ] <-
  aux_H_R_conf[i, j, k, l] + n_H_R_unconf_to_conf[i, j, k, l]
new_H_R_conf[, , 1, ] <-
  new_H_R_conf[i, j, 1, l] + n_I_C_to_H_R_conf[i, j, l]

## Work out the H_D->H_D transitions
aux_H_D_unconf[, , , ] <- H_D_unconf[i, j, k, l]
aux_H_D_unconf[, , 2:s_H_D, ] <-
  aux_H_D_unconf[i, j, k, l] + n_HH_D_unconf[i, j, k - 1, l]
aux_H_D_unconf[, , 1:s_H_D, ] <-
  aux_H_D_unconf[i, j, k, l] - n_HH_D_unconf[i, j, k, l]
aux_H_D_conf[, , , ] <- H_D_conf[i, j, k, l]
aux_H_D_conf[, , 2:s_H_D, ] <-
  aux_H_D_conf[i, j, k, l] + n_HH_D_conf[i, j, k - 1, l]
aux_H_D_conf[, , 1:s_H_D, ] <-
  aux_H_D_conf[i, j, k, l] - n_HH_D_conf[i, j, k, l]
n_H_D_unconf_to_conf[, , , ] <-
  rbinom(aux_H_D_unconf[i, j, k, l], p_test)
new_H_D_unconf[, , , ] <-
  aux_H_D_unconf[i, j, k, l] - n_H_D_unconf_to_conf[i, j, k, l]
new_H_D_unconf[, , 1, ] <-
  new_H_D_unconf[i, j, 1, l] + n_I_C_to_H_D[i, j, l] -
  n_I_C_to_H_D_conf[i, j, l]
new_H_D_conf[, , , ] <-
  aux_H_D_conf[i, j, k, l] + n_H_D_unconf_to_conf[i, j, k, l]
new_H_D_conf[, , 1, ] <-
  new_H_D_conf[i, j, 1, l] + n_I_C_to_H_D_conf[i, j, l]

## Work out the ICU_pre to ICU_D, ICU_W_R and ICU_W_D splits
n_ICU_pre_unconf_to_ICU_D_unconf[, , ] <-
  rbinom(n_ICU_pre_unconf[i, j, s_ICU_pre, k], prob_ICU_D[i])
n_ICU_pre_conf_to_ICU_D_conf[, , ] <-
  rbinom(n_ICU_pre_conf[i, j, s_ICU_pre, k], prob_ICU_D[i])
n_ICU_pre_unconf_to_ICU_W_D_unconf[, , ] <-
  rbinom(n_ICU_pre_unconf[i, j, s_ICU_pre, k] -
           n_ICU_pre_unconf_to_ICU_D_unconf[i, j, k],
         prob_W_D[i])
n_ICU_pre_unconf_to_ICU_W_R_unconf[, , ] <-
  n_ICU_pre_unconf[i, j, s_ICU_pre, k] -
  n_ICU_pre_unconf_to_ICU_D_unconf[i, j, k] -
  n_ICU_pre_unconf_to_ICU_W_D_unconf[i, j, k]
n_ICU_pre_conf_to_ICU_W_D_conf[, , ] <-
  rbinom(n_ICU_pre_conf[i, j, s_ICU_pre, k] -
           n_ICU_pre_conf_to_ICU_D_conf[i, j, k], prob_W_D[i])
n_ICU_pre_conf_to_ICU_W_R_conf[, , ] <- n_ICU_pre_conf[i, j, s_ICU_pre, k] -
  n_ICU_pre_conf_to_ICU_D_conf[i, j, k] -
  n_ICU_pre_conf_to_ICU_W_D_conf[i, j, k]


## Work out the ICU_W_R->ICU_W_R transitions
aux_ICU_W_R_unconf[, , , ] <- ICU_W_R_unconf[i, j, k, l]
aux_ICU_W_R_unconf[, , 1, ] <-
  aux_ICU_W_R_unconf[i, j, k, l] +
  n_ICU_pre_unconf_to_ICU_W_R_unconf[i, j, l]
aux_ICU_W_R_unconf[, , 2:s_ICU_W_R, ] <-
  aux_ICU_W_R_unconf[i, j, k, l] + n_IICU_W_R_unconf[i, j, k - 1, l]
aux_ICU_W_R_unconf[, , 1:s_ICU_W_R, ] <-
  aux_ICU_W_R_unconf[i, j, k, l] - n_IICU_W_R_unconf[i, j, k, l]
aux_ICU_W_R_conf[, , , ] <- ICU_W_R_conf[i, j, k, l]
aux_ICU_W_R_conf[, , 1, ] <-
  aux_ICU_W_R_conf[i, j, k, l] + n_ICU_pre_conf_to_ICU_W_R_conf[i, j, l]
aux_ICU_W_R_conf[, , 2:s_ICU_W_R, ] <-
  aux_ICU_W_R_conf[i, j, k, l] + n_IICU_W_R_conf[i, j, k - 1, l]
aux_ICU_W_R_conf[, , 1:s_ICU_W_R, ] <-
  aux_ICU_W_R_conf[i, j, k, l] - n_IICU_W_R_conf[i, j, k, l]
n_ICU_W_R_unconf_to_conf[, , , ] <-
  rbinom(aux_ICU_W_R_unconf[i, j, k, l], p_test)
new_ICU_W_R_unconf[, , , ] <-
  aux_ICU_W_R_unconf[i, j, k, l] - n_ICU_W_R_unconf_to_conf[i, j, k, l]
new_ICU_W_R_conf[, , , ] <-
  aux_ICU_W_R_conf[i, j, k, l] + n_ICU_W_R_unconf_to_conf[i, j, k, l]

## Work out the ICU_W_D->ICU_W_D transitions
aux_ICU_W_D_unconf[, , , ] <- ICU_W_D_unconf[i, j, k, l]
aux_ICU_W_D_unconf[, , 1, ] <-
  aux_ICU_W_D_unconf[i, j, k, l] +
  n_ICU_pre_unconf_to_ICU_W_D_unconf[i, j, l]
aux_ICU_W_D_unconf[, , 2:s_ICU_W_D, ] <-
  aux_ICU_W_D_unconf[i, j, k, l] + n_IICU_W_D_unconf[i, j, k - 1, l]
aux_ICU_W_D_unconf[, , 1:s_ICU_W_D, ] <-
  aux_ICU_W_D_unconf[i, j, k, l] - n_IICU_W_D_unconf[i, j, k, l]
aux_ICU_W_D_conf[, , , ] <- ICU_W_D_conf[i, j, k, l]
aux_ICU_W_D_conf[, , 1, ] <-
  aux_ICU_W_D_conf[i, j, k, l] + n_ICU_pre_conf_to_ICU_W_D_conf[i, j, l]
aux_ICU_W_D_conf[, , 2:s_ICU_W_D, ] <-
  aux_ICU_W_D_conf[i, j, k, l] + n_IICU_W_D_conf[i, j, k - 1, l]
aux_ICU_W_D_conf[, , 1:s_ICU_W_D, ] <-
  aux_ICU_W_D_conf[i, j, k, l] - n_IICU_W_D_conf[i, j, k, l]
n_ICU_W_D_unconf_to_conf[, , , ] <-
  rbinom(aux_ICU_W_D_unconf[i, j, k, l], p_test)
new_ICU_W_D_unconf[, , , ] <-
  aux_ICU_W_D_unconf[i, j, k, l] - n_ICU_W_D_unconf_to_conf[i, j, k, l]
new_ICU_W_D_conf[, , , ] <-
  aux_ICU_W_D_conf[i, j, k, l] + n_ICU_W_D_unconf_to_conf[i, j, k, l]

## Work out the ICU_D->ICU_D transitions
aux_ICU_D_unconf[, , , ] <- ICU_D_unconf[i, j, k, l]
aux_ICU_D_unconf[, , 1, ] <-
  aux_ICU_D_unconf[i, j, k, l] + n_ICU_pre_unconf_to_ICU_D_unconf[i, j, l]
aux_ICU_D_unconf[, , 2:s_ICU_D, ] <-
  aux_ICU_D_unconf[i, j, k, l] + n_IICU_D_unconf[i, j, k - 1, l]
aux_ICU_D_unconf[, , 1:s_ICU_D, ] <-
  aux_ICU_D_unconf[i, j, k, l] - n_IICU_D_unconf[i, j, k, l]
aux_ICU_D_conf[, , , ] <- ICU_D_conf[i, j, k, l]
aux_ICU_D_conf[, , 1, ] <-
  aux_ICU_D_conf[i, j, k, l] + n_ICU_pre_conf_to_ICU_D_conf[i, j, l]
aux_ICU_D_conf[, , 2:s_ICU_D, ] <-
  aux_ICU_D_conf[i, j, k, l] + n_IICU_D_conf[i, j, k - 1, l]
aux_ICU_D_conf[, , 1:s_ICU_D, ] <-
  aux_ICU_D_conf[i, j, k, l] - n_IICU_D_conf[i, j, k, l]
n_ICU_D_unconf_to_conf[, , , ] <-
  rbinom(aux_ICU_D_unconf[i, j, k, l], p_test)
new_ICU_D_unconf[, , , ] <-
  aux_ICU_D_unconf[i, j, k, l] - n_ICU_D_unconf_to_conf[i, j, k, l]
new_ICU_D_conf[, , , ] <-
  aux_ICU_D_conf[i, j, k, l] + n_ICU_D_unconf_to_conf[i, j, k, l]

## Work out the W_R->W_R transitions
aux_W_R_unconf[, , , ] <- W_R_unconf[i, j, k, l]
aux_W_R_unconf[, , 1, ] <-
  aux_W_R_unconf[i, j, k, l] + n_IICU_W_R_unconf[i, j, s_ICU_W_R, l]
aux_W_R_unconf[, , 2:s_W_R, ] <-
  aux_W_R_unconf[i, j, k, l] + n_WW_R_unconf[i, j, k - 1, l]
aux_W_R_unconf[, , 1:s_W_R, ] <-
  aux_W_R_unconf[i, j, k, l] - n_WW_R_unconf[i, j, k, l]
aux_W_R_conf[, , , ] <- W_R_conf[i, j, k, l]
aux_W_R_conf[, , 1, ] <-
  aux_W_R_conf[i, j, k, l] + n_IICU_W_R_conf[i, j, s_ICU_W_R, l]
aux_W_R_conf[, , 2:s_W_R, ] <-
  aux_W_R_conf[i, j, k, l] + n_WW_R_conf[i, j, k - 1, l]
aux_W_R_conf[, , 1:s_W_R, ] <-
  aux_W_R_conf[i, j, k, l] - n_WW_R_conf[i, j, k, l]
n_W_R_unconf_to_conf[, , , ] <-
  rbinom(aux_W_R_unconf[i, j, k, l], p_test)
new_W_R_unconf[, , , ] <-
  aux_W_R_unconf[i, j, k, l] -
  n_W_R_unconf_to_conf[i, j, k, l]
new_W_R_conf[, , , ] <-
  aux_W_R_conf[i, j, k, l] + n_W_R_unconf_to_conf[i, j, k, l]

## Work out the W_D->W_D transitions
aux_W_D_unconf[, , , ] <- W_D_unconf[i, j, k, l]
aux_W_D_unconf[, , 1, ] <-
  aux_W_D_unconf[i, j, k, l] + n_IICU_W_D_unconf[i, j, s_ICU_W_D, l]
aux_W_D_unconf[, , 2:s_W_D, ] <-
  aux_W_D_unconf[i, j, k, l] + n_WW_D_unconf[i, j, k - 1, l]
aux_W_D_unconf[, , 1:s_W_D, ] <-
  aux_W_D_unconf[i, j, k, l] - n_WW_D_unconf[i, j, k, l]
aux_W_D_conf[, , , ] <- W_D_conf[i, j, k, l]
aux_W_D_conf[, , 1, ] <-
  aux_W_D_conf[i, j, k, l] + n_IICU_W_D_conf[i, j, s_ICU_W_D, l]
aux_W_D_conf[, , 2:s_W_D, ] <-
  aux_W_D_conf[i, j, k, l] + n_WW_D_conf[i, j, k - 1, l]
aux_W_D_conf[, , 1:s_W_D, ] <-
  aux_W_D_conf[i, j, k, l] - n_WW_D_conf[i, j, k, l]
n_W_D_unconf_to_conf[, , , ] <-
  rbinom(aux_W_D_unconf[i, j, k, l], p_test)
new_W_D_unconf[, , , ] <-
  aux_W_D_unconf[i, j, k, l] -
  n_W_D_unconf_to_conf[i, j, k, l]
new_W_D_conf[, , , ] <-
  aux_W_D_conf[i, j, k, l] + n_W_D_unconf_to_conf[i, j, k, l]

## Work out the number of deaths in hospital
new_D_hosp[] <- D_hosp[i] +
  sum(n_HH_D_unconf[i, , s_H_D, ]) +
  sum(n_HH_D_conf[i, , s_H_D, ]) +
  sum(n_IICU_D_unconf[i, , s_ICU_D, ]) +
  sum(n_IICU_D_conf[i, , s_ICU_D, ]) +
  sum(n_WW_D_unconf[i, , s_W_D, ]) +
  sum(n_WW_D_conf[i, , s_W_D, ])

## Work out the number of deaths in the community
new_D_comm[] <- D_comm[i] + sum(n_GG_D[i, , s_G_D, ])

## Work out the number of people entering the seroconversion flow
n_com_to_R_pre[, , 1, 1] <- rbinom(
  n_EE[i, j, s_E, 1] + n_EE_next_vacc_class[i, j, s_E, n_vacc_classes],
  p_sero_pre_1)
n_com_to_R_pre[, , 1, 2:n_vacc_classes] <- rbinom(
  n_EE[i, j, s_E, l] + n_EE_next_vacc_class[i, j, s_E, l - 1], p_sero_pre_1)
n_com_to_R_pre[, , 2, 1] <- n_EE[i, j, s_E, 1] +
  n_EE_next_vacc_class[i, j, s_E, n_vacc_classes] - n_com_to_R_pre[i, j, 1, 1]
n_com_to_R_pre[, , 2, 2:n_vacc_classes] <- n_EE[i, j, s_E, l] +
  n_EE_next_vacc_class[i, j, s_E, l - 1] - n_com_to_R_pre[i, j, 1, l]
new_R_pre[, , , ] <-
  R_pre[i, j, k, l] + n_com_to_R_pre[i, j, k, l] - n_R_pre[i, j, k, l]


## Split the seroconversion flow between people who are going to
## seroconvert and people who are not
n_R_pre_to_R_pos[, , ] <- rbinom(sum(n_R_pre[i, j, , k]), p_sero_pos[i])

new_R_pos[, , , ] <- R_pos[i, j, k, l] - n_R_pos[i, j, k, l]
new_R_pos[, , 1, ] <- new_R_pos[i, j, 1, l] + n_R_pre_to_R_pos[i, j, l]
new_R_pos[, , 2:s_sero_pos, ] <- new_R_pos[i, j, k, l] + n_R_pos[i, j, k - 1, l]

new_T_sero_neg[, , ] <- T_sero_neg[i, j, k] + sum(n_R_pre[i, j, , k]) -
  n_R_pre_to_R_pos[i, j, k] + n_R_pos[i, j, s_sero_pos, k] -
  model_pcr_and_serology * n_R_progress[i, j, k] -
  model_pcr_and_serology * n_R_next_vacc_class[i, j, k]
new_T_sero_neg[, , 1] <- new_T_sero_neg[i, j, 1] +
  model_pcr_and_serology * n_R_next_vacc_class[i, j, n_vacc_classes]
new_T_sero_neg[, , 2:n_vacc_classes] <- new_T_sero_neg[i, j, k] +
  model_pcr_and_serology * n_R_next_vacc_class[i, j, k - 1]


## Work out the total number of recovery
new_R[, , ] <- R[i, j, k] +
  n_II_A[i, j, s_A, k] +
  n_I_C_to_R[i, j, k] +
  n_HH_R_conf[i, j, s_H_R, k] +
  n_HH_R_unconf[i, j, s_H_R, k] +
  n_WW_R_conf[i, j, s_W_R, k] +
  n_WW_R_unconf[i, j, s_W_R, k] -
  n_R_progress[i, j, k] -
  n_R_next_vacc_class[i, j, k]
new_R[, , 1] <- new_R[i, j, 1] +
  n_II_A_next_vacc_class[i, j, s_A, n_vacc_classes] +
  n_R_next_vacc_class[i, j, n_vacc_classes]
new_R[, , 2:n_vacc_classes] <- new_R[i, j, k] +
  n_II_A_next_vacc_class[i, j, s_A, k - 1] +
  n_R_next_vacc_class[i, j, k - 1]


## Work out the PCR positivity
new_T_PCR_pre[, , , ] <- T_PCR_pre[i, j, k, l] - n_T_PCR_pre[i, j, k, l]
new_T_PCR_pre[, , 1, ] <- new_T_PCR_pre[i, j, 1, l] + n_S_progress[i, j, l]
new_T_PCR_pre[, , 2:s_PCR_pre, ] <-
  new_T_PCR_pre[i, j, k, l] + n_T_PCR_pre[i, j, k - 1, l]

new_T_PCR_pos[, , , ] <- T_PCR_pos[i, j, k, l] - n_T_PCR_pos[i, j, k, l]
new_T_PCR_pos[, , 1, ] <- new_T_PCR_pos[i, j, 1, l] +
  n_T_PCR_pre[i, j, s_PCR_pre, l]
new_T_PCR_pos[, , 2:s_PCR_pos, ] <-
  new_T_PCR_pos[i, j, k, l] + n_T_PCR_pos[i, j, k - 1, l]

new_T_PCR_neg[, , ] <- T_PCR_neg[i, j, k] + n_T_PCR_pos[i, j, s_PCR_pos, k] -
  model_pcr_and_serology * n_R_progress[i, j, k] -
  model_pcr_and_serology * n_R_next_vacc_class[i, j, k]
new_T_PCR_neg[, , 1] <- new_T_PCR_neg[i, j, 1] +
  model_pcr_and_serology * n_R_next_vacc_class[i, j, n_vacc_classes]
new_T_PCR_neg[, , 2:n_vacc_classes] <- new_T_PCR_neg[i, j, k] +
  model_pcr_and_serology * n_R_next_vacc_class[i, j, k - 1]


## Compute the force of infection

I_with_diff_trans[, , ] <-
    strain_transmission[j] * (
      sum(I_A[i, j, , k]) + sum(I_C[i, j, , k]) +
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
initial(I_C[, , , ]) <- 0
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
initial(R_pre[, , , ]) <- 0
initial(R_pos[, , , ]) <- 0
initial(T_sero_neg[, , ]) <- 0
initial(R[, , ]) <- 0
initial(D_hosp[]) <- 0
initial(D_comm[]) <- 0
initial(T_PCR_pre[, , , ]) <- 0
initial(T_PCR_pos[, , , ]) <- 0
initial(T_PCR_neg[, , ]) <- 0
initial(cum_admit_conf) <- 0
initial(cum_new_conf) <- 0
initial(cum_admit_by_age[]) <- 0

## User defined parameters - default in parentheses:

## Vaccination parameters
rel_susceptibility[, ] <- user()
dim(rel_susceptibility) <- user() # use length as provided by the user
n_vacc_classes <- dim(rel_susceptibility, 2)
rel_p_C[, ] <- user()
dim(rel_p_C) <- c(n_groups, n_vacc_classes)
rel_p_H[, ] <- user()
dim(rel_p_H) <- c(n_groups, n_vacc_classes)

vaccine_progression_rate_base[, ] <- user()
dim(vaccine_progression_rate_base) <- c(n_groups, n_vacc_classes)

## Parameters of the E classes
s_E <- user()
gamma_E <- user(0.1)

## Probability of transitioning from the E to the asymptomatic class,
## the rest go into the symptomatic class
p_C[] <- user()

## Parameters of the I_A classes
s_A <- user()
gamma_A <- user(0.1)

## Parameters of the I_C classes
s_C <- user()
gamma_C <- user(0.1)
dim(p_H_step) <- user()
p_H_step[] <- user()
psi_H[] <- user()

## Parameters of the G_D class
s_G_D <- user()
gamma_G_D <- user(0.1)
dim(p_G_D_step) <- user()
p_G_D_step[] <- user()
psi_G_D[] <- user()

## Parameters of the ICU_pre classes
s_ICU_pre <- user()
gamma_ICU_pre <- user(0.1)

## Proportion of hospital cases progressing to ICU
dim(p_ICU_step) <- user()
p_ICU_step[] <- user()
psi_ICU[] <- user()

## Proportion of stepdown cases dying
dim(p_W_D_step) <- user()
p_W_D_step[] <- user()
psi_W_D[] <- user()

## Parameters of the H_R classes
s_H_R <- user()
gamma_H_R <- user(0.1)

## Parameters of the H_D classes
s_H_D <- user()
gamma_H_D <- user(0.1)
dim(p_H_D_step) <- user()
p_H_D_step[] <- user()
psi_H_D[] <- user()

## Parameters of the ICU_W_R classes
s_ICU_W_R <- user()
gamma_ICU_W_R <- user(0.1)

## Parameters of the ICU_W_D classes
s_ICU_W_D <- user()
gamma_ICU_W_D <- user(0.1)

## Parameters of the I_ICU classes
s_ICU_D <- user()
gamma_ICU_D <- user(0.1)
dim(p_ICU_D_step) <- user()
p_ICU_D_step[] <- user()
psi_ICU_D[] <- user()

## Waning of immunity
waning_rate[] <- user()
dim(waning_rate) <- n_groups

## Parameters of the W_R classes
s_W_R <- user()
gamma_W_R <- user(0.1)

## Parameters of the W_D classes
s_W_D <- user()
gamma_W_D <- user(0.1)

## Parameters of the R_pre classes
gamma_sero_pre_1 <- user(0.1)
gamma_sero_pre_2 <- user(0.1)
gamma_sero_pre[1] <- gamma_sero_pre_1
gamma_sero_pre[2] <- gamma_sero_pre_2
## Governs the mixing - pretty much only makes sense at 0.5
p_sero_pre_1 <- user(0.5)
p_sero_pos[] <- user()

## Parameters of the R_pos classes
s_sero_pos <- user()
gamma_sero_pos <- user(0.1)

## Parameters relating to testing
gamma_U <- user(0.1)
dim(p_star_step) <- user()
p_star_step[] <- user()
psi_star[] <- user()

## Parameters relating to PCR positivity
s_PCR_pre <- user()
gamma_PCR_pre <- user(0.1)
s_PCR_pos <- user()
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
hosp_transmission <- user()
ICU_transmission <- user()
G_D_transmission <- user()
strain_transmission[] <- user()
dim(strain_transmission) <- user() # use length as provided by the user
n_strains <- length(strain_transmission)

## Dimensions of the different "vectors" here vectors stand for
## multi-dimensional arrays

## Vectors handling the S class
dim(S) <- c(n_groups, n_vacc_classes)
dim(new_S) <- c(n_groups, n_vacc_classes)

## Vectors handling the E class
dim(E) <- c(n_groups, n_strains, s_E, n_vacc_classes)
dim(aux_EE) <- c(n_groups, n_strains, s_E, n_vacc_classes)
dim(new_E) <- c(n_groups, n_strains, s_E, n_vacc_classes)
dim(n_EE) <- c(n_groups, n_strains, s_E, n_vacc_classes)

## Vectors handling the I_A class
dim(I_A) <- c(n_groups, n_strains, s_A, n_vacc_classes)
dim(aux_II_A) <- c(n_groups, n_strains, s_A, n_vacc_classes)
dim(new_I_A) <- c(n_groups, n_strains, s_A, n_vacc_classes)
dim(n_II_A) <- c(n_groups, n_strains, s_A, n_vacc_classes)

## Vectors handling the I_C class
dim(I_C) <- c(n_groups, n_strains, s_C, n_vacc_classes)
dim(aux_II_C) <- c(n_groups, n_strains, s_C, n_vacc_classes)
dim(new_I_C) <- c(n_groups, n_strains, s_C, n_vacc_classes)
dim(n_II_C) <- c(n_groups, n_strains, s_C, n_vacc_classes)
dim(prob_H) <- n_groups
dim(psi_H) <- n_groups

## Vectors handling the G_D class
dim(G_D) <- c(n_groups, n_strains, s_G_D, n_vacc_classes)
dim(aux_GG_D) <- c(n_groups, n_strains, s_G_D, n_vacc_classes)
dim(new_G_D) <- c(n_groups, n_strains, s_G_D, n_vacc_classes)
dim(n_GG_D) <- c(n_groups, n_strains, s_G_D, n_vacc_classes)
dim(prob_G_D) <- n_groups
dim(psi_G_D) <- n_groups

## Vectors handling the ICU_pre class
dim(ICU_pre_unconf) <- c(n_groups, n_strains, s_ICU_pre, n_vacc_classes)
dim(aux_ICU_pre_unconf) <- c(n_groups, n_strains, s_ICU_pre, n_vacc_classes)
dim(new_ICU_pre_unconf) <- c(n_groups, n_strains, s_ICU_pre, n_vacc_classes)
dim(n_ICU_pre_unconf) <- c(n_groups, n_strains, s_ICU_pre, n_vacc_classes)
dim(ICU_pre_conf) <- c(n_groups, n_strains, s_ICU_pre, n_vacc_classes)
dim(aux_ICU_pre_conf) <- c(n_groups, n_strains, s_ICU_pre, n_vacc_classes)
dim(new_ICU_pre_conf) <- c(n_groups, n_strains, s_ICU_pre, n_vacc_classes)
dim(n_ICU_pre_conf) <- c(n_groups, n_strains, s_ICU_pre, n_vacc_classes)
dim(n_ICU_pre_unconf_to_conf) <-
  c(n_groups, n_strains, s_ICU_pre, n_vacc_classes)

## Vector handling who progress to ICU
dim(prob_ICU) <- n_groups
dim(psi_ICU) <- n_groups

## Vectors handling the H_R class
dim(H_R_unconf) <- c(n_groups, n_strains, s_H_R, n_vacc_classes)
dim(aux_H_R_unconf) <- c(n_groups, n_strains, s_H_R, n_vacc_classes)
dim(new_H_R_unconf) <- c(n_groups, n_strains, s_H_R, n_vacc_classes)
dim(n_HH_R_unconf) <- c(n_groups, n_strains, s_H_R, n_vacc_classes)
dim(H_R_conf) <- c(n_groups, n_strains, s_H_R, n_vacc_classes)
dim(aux_H_R_conf) <- c(n_groups, n_strains, s_H_R, n_vacc_classes)
dim(new_H_R_conf) <- c(n_groups, n_strains, s_H_R, n_vacc_classes)
dim(n_HH_R_conf) <- c(n_groups, n_strains, s_H_R, n_vacc_classes)
dim(n_H_R_unconf_to_conf) <-
  c(n_groups, n_strains, s_H_R, n_vacc_classes)

## Vectors handling the H_D class
dim(H_D_unconf) <- c(n_groups, n_strains, s_H_D, n_vacc_classes)
dim(aux_H_D_unconf) <- c(n_groups, n_strains, s_H_D, n_vacc_classes)
dim(new_H_D_unconf) <- c(n_groups, n_strains, s_H_D, n_vacc_classes)
dim(n_HH_D_unconf) <- c(n_groups, n_strains, s_H_D, n_vacc_classes)
dim(H_D_conf) <- c(n_groups, n_strains, s_H_D, n_vacc_classes)
dim(aux_H_D_conf) <- c(n_groups, n_strains, s_H_D, n_vacc_classes)
dim(new_H_D_conf) <- c(n_groups, n_strains, s_H_D, n_vacc_classes)
dim(n_HH_D_conf) <- c(n_groups, n_strains, s_H_D, n_vacc_classes)
dim(n_H_D_unconf_to_conf) <-
  c(n_groups, n_strains, s_H_D, n_vacc_classes)

## Vectors handling the ICU_W_R class
dim(ICU_W_R_unconf) <- c(n_groups, n_strains, s_ICU_W_R, n_vacc_classes)
dim(aux_ICU_W_R_unconf) <- c(n_groups, n_strains, s_ICU_W_R, n_vacc_classes)
dim(new_ICU_W_R_unconf) <- c(n_groups, n_strains, s_ICU_W_R, n_vacc_classes)
dim(n_IICU_W_R_unconf) <- c(n_groups, n_strains, s_ICU_W_R, n_vacc_classes)
dim(ICU_W_R_conf) <- c(n_groups, n_strains, s_ICU_W_R, n_vacc_classes)
dim(aux_ICU_W_R_conf) <- c(n_groups, n_strains, s_ICU_W_R, n_vacc_classes)
dim(new_ICU_W_R_conf) <- c(n_groups, n_strains, s_ICU_W_R, n_vacc_classes)
dim(n_IICU_W_R_conf) <- c(n_groups, n_strains, s_ICU_W_R, n_vacc_classes)
dim(n_ICU_W_R_unconf_to_conf) <-
  c(n_groups, n_strains, s_ICU_W_R, n_vacc_classes)

## Vectors handling the ICU_W_D class
dim(ICU_W_D_unconf) <- c(n_groups, n_strains, s_ICU_W_D, n_vacc_classes)
dim(aux_ICU_W_D_unconf) <- c(n_groups, n_strains, s_ICU_W_D, n_vacc_classes)
dim(new_ICU_W_D_unconf) <- c(n_groups, n_strains, s_ICU_W_D, n_vacc_classes)
dim(n_IICU_W_D_unconf) <- c(n_groups, n_strains, s_ICU_W_D, n_vacc_classes)
dim(ICU_W_D_conf) <- c(n_groups, n_strains, s_ICU_W_D, n_vacc_classes)
dim(aux_ICU_W_D_conf) <- c(n_groups, n_strains, s_ICU_W_D, n_vacc_classes)
dim(new_ICU_W_D_conf) <- c(n_groups, n_strains, s_ICU_W_D, n_vacc_classes)
dim(n_IICU_W_D_conf) <- c(n_groups, n_strains, s_ICU_W_D, n_vacc_classes)
dim(n_ICU_W_D_unconf_to_conf) <-
  c(n_groups, n_strains, s_ICU_W_D, n_vacc_classes)

## Vectors handling the ICU_D class
dim(ICU_D_unconf) <- c(n_groups, n_strains, s_ICU_D, n_vacc_classes)
dim(aux_ICU_D_unconf) <- c(n_groups, n_strains, s_ICU_D, n_vacc_classes)
dim(new_ICU_D_unconf) <- c(n_groups, n_strains, s_ICU_D, n_vacc_classes)
dim(n_IICU_D_unconf) <- c(n_groups, n_strains, s_ICU_D, n_vacc_classes)
dim(ICU_D_conf) <- c(n_groups, n_strains, s_ICU_D, n_vacc_classes)
dim(aux_ICU_D_conf) <- c(n_groups, n_strains, s_ICU_D, n_vacc_classes)
dim(new_ICU_D_conf) <- c(n_groups, n_strains, s_ICU_D, n_vacc_classes)
dim(n_IICU_D_conf) <- c(n_groups, n_strains, s_ICU_D, n_vacc_classes)
dim(n_ICU_D_unconf_to_conf) <- c(n_groups, n_strains, s_ICU_D, n_vacc_classes)

## Vectors handling the W_R class
dim(W_R_unconf) <- c(n_groups, n_strains, s_W_R, n_vacc_classes)
dim(aux_W_R_unconf) <-
  c(n_groups, n_strains, s_W_R, n_vacc_classes)
dim(new_W_R_unconf) <-
  c(n_groups, n_strains, s_W_R, n_vacc_classes)
dim(n_WW_R_unconf) <-
  c(n_groups, n_strains, s_W_R, n_vacc_classes)
dim(W_R_conf) <- c(n_groups, n_strains, s_W_R, n_vacc_classes)
dim(aux_W_R_conf) <-
  c(n_groups, n_strains, s_W_R, n_vacc_classes)
dim(new_W_R_conf) <-
  c(n_groups, n_strains, s_W_R, n_vacc_classes)
dim(n_WW_R_conf) <- c(n_groups, n_strains, s_W_R, n_vacc_classes)
dim(n_W_R_unconf_to_conf) <-
  c(n_groups, n_strains, s_W_R, n_vacc_classes)

## Vectors handling the W_D class
dim(W_D_unconf) <- c(n_groups, n_strains, s_W_D, n_vacc_classes)
dim(aux_W_D_unconf) <-
  c(n_groups, n_strains, s_W_D, n_vacc_classes)
dim(new_W_D_unconf) <-
  c(n_groups, n_strains, s_W_D, n_vacc_classes)
dim(n_WW_D_unconf) <-
  c(n_groups, n_strains, s_W_D, n_vacc_classes)
dim(W_D_conf) <-
  c(n_groups, n_strains, s_W_D, n_vacc_classes)
dim(aux_W_D_conf) <-
  c(n_groups, n_strains, s_W_D, n_vacc_classes)
dim(new_W_D_conf) <-
  c(n_groups, n_strains, s_W_D, n_vacc_classes)
dim(n_WW_D_conf) <- c(n_groups, n_strains, s_W_D, n_vacc_classes)
dim(n_W_D_unconf_to_conf) <-
  c(n_groups, n_strains, s_W_D, n_vacc_classes)

## Vectors handling the R_pos class
dim(R) <- c(n_groups, n_strains, n_vacc_classes)
dim(new_R) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the R_pre class and seroconversion
dim(R_pre) <- c(n_groups, n_strains, 2, n_vacc_classes)
dim(new_R_pre) <- c(n_groups, n_strains, 2, n_vacc_classes)
dim(n_R_pre) <- c(n_groups, n_strains, 2, n_vacc_classes)
dim(gamma_sero_pre) <- 2
dim(p_R_pre) <- c(n_groups, n_strains, 2, n_vacc_classes)
dim(p_sero_pos) <- n_groups

## Vectors handling the R_pos class
dim(R_pos) <- c(n_groups, n_strains, s_sero_pos, n_vacc_classes)
dim(n_R_pos) <- c(n_groups, n_strains, s_sero_pos, n_vacc_classes)
dim(new_R_pos) <- c(n_groups, n_strains, s_sero_pos, n_vacc_classes)
dim(n_R_pre_to_R_pos) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the T_sero_neg class
dim(T_sero_neg) <- c(n_groups, n_strains, n_vacc_classes)
dim(new_T_sero_neg) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the D_hosp class
dim(D_hosp) <- n_groups
dim(new_D_hosp) <- n_groups

## Vectors handling the D_comm class
dim(D_comm) <- n_groups
dim(new_D_comm) <- n_groups

## Vectors handling the PCR classes
dim(T_PCR_pre) <- c(n_groups, n_strains, s_PCR_pre, n_vacc_classes)
dim(n_T_PCR_pre) <- c(n_groups, n_strains, s_PCR_pre, n_vacc_classes)
dim(new_T_PCR_pre) <- c(n_groups, n_strains, s_PCR_pre, n_vacc_classes)
dim(T_PCR_pos) <- c(n_groups, n_strains, s_PCR_pos, n_vacc_classes)
dim(n_T_PCR_pos) <- c(n_groups, n_strains, s_PCR_pos, n_vacc_classes)
dim(new_T_PCR_pos) <- c(n_groups, n_strains, s_PCR_pos, n_vacc_classes)
dim(T_PCR_neg) <- c(n_groups, n_strains, n_vacc_classes)
dim(new_T_PCR_neg) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the S->S transitions i.e. moving between vaccination classes
dim(p_S_next_vacc_class) <- c(n_groups, n_vacc_classes)
dim(n_S_next_vacc_class) <- c(n_groups, n_vacc_classes)

dim(p_E_next_vacc_class) <- c(n_groups, n_strains, s_E, n_vacc_classes)
dim(n_E_next_vacc_class) <- c(n_groups, n_strains, s_E, n_vacc_classes)
dim(n_E_progress) <- c(n_groups, n_strains, s_E, n_vacc_classes)
dim(n_EE_next_vacc_class) <- c(n_groups, n_strains, s_E, n_vacc_classes)

dim(p_I_A_next_vacc_class) <-
  c(n_groups, n_strains, s_A, n_vacc_classes)
dim(n_I_A_next_vacc_class) <-
  c(n_groups, n_strains, s_A, n_vacc_classes)
dim(n_I_A_progress) <- c(n_groups, n_strains, s_A, n_vacc_classes)
dim(n_II_A_next_vacc_class) <-
  c(n_groups, n_strains, s_A, n_vacc_classes)

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
dim(n_EI_C) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_EI_A_next_vacc_class) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_EI_C_next_vacc_class) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling I_C to R, G_D transition
dim(n_I_C_to_G_D) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_to_R) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling number of new hospitalisations, ICU admissions and
## recoveries in hospital
dim(n_I_C_to_hosp) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_to_ICU_pre) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_to_ICU_pre_conf) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_hosp_non_ICU) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_to_H_D) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_to_H_D_conf) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_to_H_R) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_to_H_R_conf) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_ICU_pre_unconf_to_ICU_D_unconf) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_ICU_pre_conf_to_ICU_D_conf) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_ICU_pre_unconf_to_ICU_W_R_unconf) <-
  c(n_groups, n_strains, n_vacc_classes)
dim(n_ICU_pre_conf_to_ICU_W_R_conf) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_ICU_pre_unconf_to_ICU_W_D_unconf) <-
  c(n_groups, n_strains, n_vacc_classes)
dim(n_ICU_pre_conf_to_ICU_W_D_conf) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the serology flow
dim(n_com_to_R_pre) <- c(n_groups, n_strains, 2, n_vacc_classes)

## Vectors handling the severity profile
dim(p_C) <- n_groups

## Vectors handling the potential death in hospital (general beds and ICU)
dim(prob_H_D) <- n_groups
dim(psi_H_D) <- n_groups
dim(prob_ICU_D) <- n_groups
dim(psi_ICU_D) <- n_groups
dim(prob_W_D) <- n_groups
dim(psi_W_D) <- n_groups

## Vector handling the probability of being admitted as confirmed
dim(prob_admit_conf) <- n_groups
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
update(N_tot[]) <- sum(S[i, ]) + sum(R[i, , ]) + D_hosp[i] + sum(E[i, , , ]) +
  sum(I_A[i, , , ]) + sum(I_C[i, , , ]) +
  sum(ICU_pre_conf[i, , , ]) + sum(ICU_pre_unconf[i, , , ])  +
  sum(H_R_conf[i, , , ]) + sum(H_R_unconf[i, , , ]) +
  sum(H_D_conf[i, , , ]) + sum(H_D_unconf[i, , , ]) +
  sum(ICU_W_R_conf[i, , , ]) + sum(ICU_W_R_unconf[i, , , ]) +
  sum(ICU_W_D_conf[i, , , ]) + sum(ICU_W_D_unconf[i, , , ]) +
  sum(ICU_D_conf[i, , , ]) + sum(ICU_D_unconf[i, , , ]) +
  sum(W_R_conf[i, , , ]) + sum(W_R_unconf[i, , , ]) +
  sum(W_D_conf[i, , , ]) + sum(W_D_unconf[i, , , ]) +
  sum(G_D[i, , , ]) + D_comm[i]
dim(N_tot) <- n_groups

## Total population calculated with seroconversion flow
initial(N_tot2) <- 0
update(N_tot2) <- sum(S) + sum(R_pre) + sum(R_pos) + sum(T_sero_neg) + sum(E)

## Total population calculated with PCR flow
initial(N_tot3) <- 0
update(N_tot3) <- sum(S) + sum(T_PCR_pre) + sum(T_PCR_pos) + sum(T_PCR_neg)

## Aggregate our reporting statistics by summing across age (simple
## for everything except for seropositivity data, done last)
initial(I_ICU_tot) <- 0
new_I_ICU_tot <- sum(new_ICU_W_R_conf) + sum(new_ICU_W_D_conf) +
  sum(new_ICU_D_conf)
update(I_ICU_tot) <- new_I_ICU_tot

initial(general_tot) <- 0
new_general_tot <- sum(new_ICU_pre_conf) + sum(new_H_R_conf) +
  sum(new_H_D_conf) + sum(new_W_R_conf) +
  sum(new_W_D_conf)
update(general_tot) <- new_general_tot

initial(hosp_tot) <- 0
update(hosp_tot) <- new_I_ICU_tot + new_general_tot

initial(D_hosp_tot) <- 0
new_D_hosp_tot <- sum(new_D_hosp)
update(D_hosp_tot) <- new_D_hosp_tot

initial(D_comm_tot) <- 0
new_D_comm_tot <- sum(new_D_comm)
update(D_comm_tot) <- new_D_comm_tot

initial(D_tot) <- 0
update(D_tot) <- new_D_hosp_tot + new_D_comm_tot

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
update(sero_pos) <- sum(new_R_pos[4:13, , , ])

initial(cum_sympt_cases) <- 0
update(cum_sympt_cases) <- cum_sympt_cases + sum(n_EI_C)

## only over 25s (exclude groups 1 to 5)
initial(cum_sympt_cases_over25) <- 0
update(cum_sympt_cases_over25) <- cum_sympt_cases_over25 +
  sum(n_EI_C[6:n_groups, , ])

## For REACT we exclude the 0-4 (1) and CHR (19) groups
initial(react_pos) <- 0
update(react_pos) <- sum(new_T_PCR_pos[2:18, , , ])

## Vaccination engine

## First, the number of candidates
vaccine_n_candidates[] <- S[i, 1] + sum(E[i, , , 1]) + sum(I_A[i, , , 1]) +
  sum(R[i, , 1])
dim(vaccine_n_candidates) <- n_groups

## The total population reluctant to be vaccinated. Currently modelled
## as a fixed population, rather than as (say) a stratification of
## compartments.
vaccine_population_reluctant[] <- user()
dim(vaccine_population_reluctant) <- n_groups

## We will refuse to vaccine the reluctant population; this is just an
## approximation of that for now.
##
## TODO: this *should* work with
## > max(0, vaccine_n_candidates[i] - vaccine_population_reluctant[i])
## But that is generating invalid code
vaccine_population_possible[] <-
  (if (vaccine_population_reluctant[i] > vaccine_n_candidates[i]) 0
   else vaccine_n_candidates[i] - vaccine_population_reluctant[i])
dim(vaccine_population_possible) <- n_groups

## The number of doses of vaccine available each day:
vaccine_daily_doses <- user(0)

## We'll set this up to treat the first column specially, as that is
## the compartment through which people ge vaccinated, others are
## taken through the vaccination rate
config(include) <- "vaccination.cpp"
vaccine_progression_rate[, 1] <-
  vaccination_schedule(i, vaccine_daily_doses, dt,
                       vaccine_n_candidates, vaccine_population_possible)
vaccine_progression_rate[, 2:n_vacc_classes] <-
  vaccine_progression_rate_base[i, j]

dim(vaccine_progression_rate) <- c(n_groups, n_vacc_classes)
