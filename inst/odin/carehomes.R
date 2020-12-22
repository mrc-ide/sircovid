## E and R stage indexed by i, j, k with
## i for the age group
## j for the progression (not exponential latent and infectious period)
## k for the infectivity group (for I) or vacc. group (for S)

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
  n_S_next_vacc_class[i, j] + n_SE_next_vacc_class[i, j]
dim(cum_n_S_vaccinated) <- c(n_groups, n_vacc_classes)
## vaccinated E
initial(cum_n_E_vaccinated[, ]) <- 0
update(cum_n_E_vaccinated[, ]) <- cum_n_E_vaccinated[i, j] +
  sum(n_E_next_vacc_class[i, , j]) + sum(n_EE_next_vacc_class[i, , j])
dim(cum_n_E_vaccinated) <- c(n_groups, n_vacc_classes)
## vaccinated I_asympt
initial(cum_n_I_asympt_vaccinated[, ]) <- 0
update(cum_n_I_asympt_vaccinated[, ]) <- cum_n_I_asympt_vaccinated[i, j] +
  sum(n_I_asympt_next_vacc_class[i, , j]) +
  sum(n_II_asympt_next_vacc_class[i, , j])
dim(cum_n_I_asympt_vaccinated) <- c(n_groups, n_vacc_classes)
## vaccinated R
initial(cum_n_R_vaccinated[, ]) <- 0
update(cum_n_R_vaccinated[, ]) <- cum_n_R_vaccinated[i, j] +
  n_R_next_vacc_class[i, j] + n_RS_next_vacc_class[i, j]
dim(cum_n_R_vaccinated) <- c(n_groups, n_vacc_classes)

## Total number of vaccinations over S, E, I_asypmt, R for convenience
initial(cum_n_vaccinated[, ]) <- 0
update(cum_n_vaccinated[, ]) <-
  cum_n_S_vaccinated[i, j] +
  cum_n_E_vaccinated[i, j] +
  cum_n_I_asympt_vaccinated[i, j] +
  cum_n_R_vaccinated[i, j]
dim(cum_n_vaccinated) <- c(n_groups, n_vacc_classes)

## Core equations for transitions between compartments:
update(S[, ]) <- new_S[i, j]
update(E[, , ]) <- new_E[i, j, k]
update(I_asympt[, , ]) <- new_I_asympt[i, j, k]
update(I_sympt[, , ]) <- new_I_sympt[i, j, k]
update(I_comm_D[, , ]) <- new_I_comm_D[i, j, k]
update(I_triage_unconf[, , ]) <- new_I_triage_unconf[i, j, k]
update(I_triage_conf[, , ]) <- new_I_triage_conf[i, j, k]
update(I_hosp_R_unconf[, , ]) <- new_I_hosp_R_unconf[i, j, k]
update(I_hosp_R_conf[, , ]) <- new_I_hosp_R_conf[i, j, k]
update(I_hosp_D_unconf[, , ]) <- new_I_hosp_D_unconf[i, j, k]
update(I_hosp_D_conf[, , ]) <- new_I_hosp_D_conf[i, j, k]
update(I_ICU_S_R_unconf[, , ]) <- new_I_ICU_S_R_unconf[i, j, k]
update(I_ICU_S_R_conf[, , ]) <- new_I_ICU_S_R_conf[i, j, k]
update(I_ICU_S_D_unconf[, , ]) <- new_I_ICU_S_D_unconf[i, j, k]
update(I_ICU_S_D_conf[, , ]) <- new_I_ICU_S_D_conf[i, j, k]
update(I_ICU_D_unconf[, , ]) <- new_I_ICU_D_unconf[i, j, k]
update(I_ICU_D_conf[, , ]) <- new_I_ICU_D_conf[i, j, k]
update(R_stepdown_R_unconf[, , ]) <- new_R_stepdown_R_unconf[i, j, k]
update(R_stepdown_R_conf[, , ]) <- new_R_stepdown_R_conf[i, j, k]
update(R_stepdown_D_unconf[, , ]) <- new_R_stepdown_D_unconf[i, j, k]
update(R_stepdown_D_conf[, , ]) <- new_R_stepdown_D_conf[i, j, k]
update(R_pre[, , ]) <- new_R_pre[i, j, k]
update(R_pos[, , ]) <- new_R_pos[i, j, k]
update(R_neg[, ]) <- new_R_neg[i, j]
update(R[, ]) <- new_R[i, j]
update(D_hosp[]) <- new_D_hosp[i]
update(D_comm[]) <- new_D_comm[i]
update(PCR_pre[, , ]) <- new_PCR_pre[i, j, k]
update(PCR_pos[, , ]) <- new_PCR_pos[i, j, k]
update(PCR_neg[, ]) <- new_PCR_neg[i, j]
update(cum_admit_conf) <-
  cum_admit_conf +
  sum(n_sympt_to_hosp_D_conf) +
  sum(n_sympt_to_hosp_R_conf) +
  sum(n_sympt_to_triage_conf)
update(cum_new_conf) <-
  cum_new_conf +
  sum(n_I_hosp_D_unconf_to_conf) +
  sum(n_I_hosp_R_unconf_to_conf) +
  sum(n_I_triage_unconf_to_conf) +
  sum(n_I_ICU_D_unconf_to_conf) +
  sum(n_I_ICU_S_R_unconf_to_conf) +
  sum(n_I_ICU_S_D_unconf_to_conf) +
  sum(n_R_stepdown_R_unconf_to_conf) +
  sum(n_R_stepdown_D_unconf_to_conf)
update(cum_admit_by_age[]) <- cum_admit_by_age[i] + sum(n_sympt_to_hosp[i, ])

## Individual probabilities of transition:

## vaccination
p_S_next_vacc_class[, ] <- 1 - exp(-vaccine_progression_rate[i, j] * dt)
p_E_next_vacc_class[, , ] <- 1 - exp(-vaccine_progression_rate[i, k] * dt)
p_I_asympt_next_vacc_class[, , ] <-
  1 - exp(-vaccine_progression_rate[i, k] * dt)
p_R_next_vacc_class[, ] <-
  1 - exp(-vaccine_progression_rate[i, j] * dt)
## clinical progression
p_SE[, ] <- 1 - exp(-lambda[i] *
                      rel_susceptibility[i, j] * dt) # S to I age/vacc dependent
p_EE <- 1 - exp(-gamma_E * dt) # progression of latent period
p_II_asympt <- 1 - exp(-gamma_asympt * dt) # progression of infectious period
p_II_sympt <- 1 - exp(-gamma_sympt * dt)
p_II_comm_D <- 1 - exp(-gamma_comm_D * dt)
p_II_triage <- 1 - exp(-gamma_triage * dt)
p_II_hosp_R <- 1 - exp(-gamma_hosp_R * dt)
p_II_hosp_D <- 1 - exp(-gamma_hosp_D * dt)
p_II_ICU_S_R <- 1 - exp(-gamma_ICU_S_R * dt)
p_II_ICU_S_D <- 1 - exp(-gamma_ICU_S_D * dt)
p_II_ICU_D <- 1 - exp(-gamma_ICU_D * dt)
p_R_stepdown_R <- 1 - exp(-gamma_stepdown_R * dt)
p_R_stepdown_D <- 1 - exp(-gamma_stepdown_D * dt)
p_R_pre[, , ] <- 1 - exp(-gamma_R_pre[j] * dt)
p_R_pos <- 1 - exp(-gamma_R_pos * dt)
p_test <- 1 - exp(-gamma_test * dt)
p_PCR_pre <- 1 - exp(-gamma_PCR_pre * dt)
p_PCR_pos <- 1 - exp(-gamma_PCR_pos * dt)
p_RS[] <- 1 - exp(-waning_rate[i] * dt) # R to S age dependent

## Work out time-varying probabilities
p_ICU_hosp <- if (step >= length(p_ICU_hosp_step))
  p_ICU_hosp_step[length(p_ICU_hosp_step)] else p_ICU_hosp_step[step + 1]
prob_ICU_hosp[] <- p_ICU_hosp * psi_ICU_hosp[i]

p_hosp_sympt <- if (step >= length(p_hosp_sympt_step))
  p_hosp_sympt_step[length(p_hosp_sympt_step)] else p_hosp_sympt_step[step + 1]
prob_hosp_sympt[] <- p_hosp_sympt * psi_hosp_sympt[i]

p_death_ICU <- if (step >= length(p_death_ICU_step))
  p_death_ICU_step[length(p_death_ICU_step)] else p_death_ICU_step[step + 1]
prob_death_ICU[] <- p_death_ICU * psi_death_ICU[i]

p_death_hosp_D <- if (step >= length(p_death_hosp_D_step))
  p_death_hosp_D_step[length(p_death_hosp_D_step)] else
    p_death_hosp_D_step[step + 1]
prob_death_hosp_D[] <- p_death_hosp_D * psi_death_hosp_D[i]

p_death_stepdown <- if (step >= length(p_death_stepdown_step))
  p_death_stepdown_step[length(p_death_stepdown_step)] else
    p_death_stepdown_step[step + 1]
prob_death_stepdown[] <- p_death_stepdown * psi_death_stepdown[i]

p_death_comm <- if (step >= length(p_death_comm_step))
  p_death_comm_step[length(p_death_comm_step)] else p_death_comm_step[step + 1]
prob_death_comm[] <- p_death_comm * psi_death_comm[i]

p_admit_conf <- if (step >= length(p_admit_conf_step))
  p_admit_conf_step[length(p_admit_conf_step)] else p_admit_conf_step[step + 1]
prob_admit_conf[] <- p_admit_conf * psi_admit_conf[i]

## Draws from binomial distributions for numbers changing between
## compartments:

## modelling infections and vaccine progression, which can happen simultaneously

#### flow out of S ####

## new infections
n_S_progress[, ] <- rbinom(S[i, j], p_SE[i, j])
## of those some can also be vaccinated or progress through vaccination classes
## --> number transitioning from S[j] to E[j+1] (j vaccination class)
n_SE_next_vacc_class[, ] <-
  rbinom(n_S_progress[i, j], p_S_next_vacc_class[i, j])
## resulting transitions from S[j] to E[j]
## (j vaccination class)
n_SE[, 1:n_vacc_classes] <- n_S_progress[i, j] - n_SE_next_vacc_class[i, j]

## vaccine progression
n_S_next_vacc_class[, ] <- rbinom(S[i, j] - n_S_progress[i, j],
                                  p_S_next_vacc_class[i, j])

#### flow out of E ####

n_E_progress[, , ] <- rbinom(E[i, j, k], p_EE)
## of those some can also be vaccinated or progress through vaccination classes
## --> number transitioning from E[j, k] to E[j+1, k+1] (k vaccination class)
n_EE_next_vacc_class[, , ] <-
  rbinom(n_E_progress[i, j, k], p_E_next_vacc_class[i, j, k])
## resulting transitions from E[j, k] to E[j+1, k]
## (k vaccination class)
n_EE[, , ] <- n_E_progress[i, j, k] - n_EE_next_vacc_class[i, j, k]

## vaccine progression
n_E_next_vacc_class[, , ] <- rbinom(E[i, j, k] - n_E_progress[i, j, k],
                                  p_E_next_vacc_class[i, j, k])

#### flow out of I_asympt ####

n_I_asympt_progress[, , ] <- rbinom(I_asympt[i, j, k], p_II_asympt)
## of those some can also be vaccinated or progress through vaccination classes
## --> number transitioning from I_asympt[j, k] to I_asympt[j+1, k+1]
## (k vaccination class)
n_II_asympt_next_vacc_class[, , ] <-
  rbinom(n_I_asympt_progress[i, j, k], p_I_asympt_next_vacc_class[i, j, k])
## resulting transitions from I_asympt[j, k] to I_asympt[j+1, k]
## (k vaccination class)
n_II_asympt[, , ] <- n_I_asympt_progress[i, j, k] -
  n_II_asympt_next_vacc_class[i, j, k]

## vaccine progression
n_I_asympt_next_vacc_class[, , ] <- rbinom(
  I_asympt[i, j, k] - n_I_asympt_progress[i, j, k],
  p_I_asympt_next_vacc_class[i, j, k])


#### flow out of R ####

n_R_progress_tmp[, ] <- rbinom(R[i, j], p_RS[i])
## cap on people who can move out of R based on numbers in R_neg and PCR_neg
n_R_progress_capped[, ] <-
  min(n_R_progress_tmp[i, j], R_neg[i, j], PCR_neg[i, j])
## use cap or not depending on model_pcr_and_serology value
n_R_progress[, ] <- if (model_pcr_and_serology == 1)
  n_R_progress_capped[i, j] else n_R_progress_tmp[i, j]
## of those some can also be vaccinated or progress through vaccination classes
## --> number transitioning from R[j] to S[j+1]
## (j vaccination class)
n_RS_next_vacc_class[, ] <-
  rbinom(n_R_progress[i, j], p_R_next_vacc_class[i, j])
## resulting transitions from R[j] to S[j]
## (j vaccination class)
n_RS[, ] <- n_R_progress[i, j] - n_RS_next_vacc_class[i, j]

## vaccine progression
n_R_next_vacc_class_tmp[, ] <- rbinom(
  R[i, j] - n_R_progress[i, j], p_R_next_vacc_class[i, j])
n_R_next_vacc_class_capped[, ] <- min(n_R_next_vacc_class_tmp[i, j],
  R_neg[i, j] - n_R_progress[i, j], PCR_neg[i, j] - n_R_progress[i, j])
n_R_next_vacc_class[, ] <- if (model_pcr_and_serology == 1)
  n_R_next_vacc_class_capped[i, j] else n_R_next_vacc_class_tmp[i, j]
  
#### other transitions ####

n_II_sympt[, , ] <- rbinom(I_sympt[i, j, k], p_II_sympt)
n_II_comm_D[, , ] <- rbinom(I_comm_D[i, j, k], p_II_comm_D)
n_II_triage_unconf[, , ] <- rbinom(I_triage_unconf[i, j, k], p_II_triage)
n_II_triage_conf[, , ] <- rbinom(I_triage_conf[i, j, k], p_II_triage)
n_II_hosp_R_unconf[, , ] <- rbinom(I_hosp_R_unconf[i, j, k], p_II_hosp_R)
n_II_hosp_R_conf[, , ] <- rbinom(I_hosp_R_conf[i, j, k], p_II_hosp_R)
n_II_hosp_D_unconf[, , ] <- rbinom(I_hosp_D_unconf[i, j, k], p_II_hosp_D)
n_II_hosp_D_conf[, , ] <- rbinom(I_hosp_D_conf[i, j, k], p_II_hosp_D)
n_II_ICU_S_R_unconf[, , ] <- rbinom(I_ICU_S_R_unconf[i, j, k], p_II_ICU_S_R)
n_II_ICU_S_R_conf[, , ] <- rbinom(I_ICU_S_R_conf[i, j, k], p_II_ICU_S_R)
n_II_ICU_S_D_unconf[, , ] <- rbinom(I_ICU_S_D_unconf[i, j, k], p_II_ICU_S_D)
n_II_ICU_S_D_conf[, , ] <- rbinom(I_ICU_S_D_conf[i, j, k], p_II_ICU_S_D)
n_II_ICU_D_unconf[, , ] <- rbinom(I_ICU_D_unconf[i, j, k], p_II_ICU_D)
n_II_ICU_D_conf[, , ] <- rbinom(I_ICU_D_conf[i, j, k], p_II_ICU_D)
n_R_stepdown_R_unconf[, , ] <-
  rbinom(R_stepdown_R_unconf[i, j, k], p_R_stepdown_R)
n_R_stepdown_R_conf[, , ] <- rbinom(R_stepdown_R_conf[i, j, k], p_R_stepdown_R)
n_R_stepdown_D_unconf[, , ] <-
  rbinom(R_stepdown_D_unconf[i, j, k], p_R_stepdown_D)
n_R_stepdown_D_conf[, , ] <- rbinom(R_stepdown_D_conf[i, j, k], p_R_stepdown_D)
n_R_pre[, , ] <- rbinom(R_pre[i, j, k], p_R_pre[i, j, k])
n_R_pos[, , ] <- rbinom(R_pos[i, j, k], p_R_pos)
n_PCR_pre[, , ] <- rbinom(PCR_pre[i, j, k], p_PCR_pre)
n_PCR_pos[, , ] <- rbinom(PCR_pos[i, j, k], p_PCR_pos)

## Cumulative infections, summed over all age groups
initial(cum_infections) <- 0
update(cum_infections) <- cum_infections + sum(n_S_progress)


## Work out the new S (i for age, j for vaccination status)
new_S[, ] <- S[i, j] + n_RS[i, j] - n_S_progress[i, j] -
  n_S_next_vacc_class[i, j]
new_S[, 1] <- new_S[i, 1] + n_S_next_vacc_class[i, n_vacc_classes] +
  n_RS_next_vacc_class[i, n_vacc_classes]
new_S[, 2:n_vacc_classes] <- new_S[i, j] + n_S_next_vacc_class[i, j - 1] +
  n_RS_next_vacc_class[i, j - 1]


## Computes the number of asymptomatic
n_EI_asympt[, ] <- rbinom(n_EE[i, s_E, j],
                          1 - p_sympt[i] * rel_p_sympt[i, j])
n_EI_asympt_next_vacc_class[, ] <-
  rbinom(n_EE_next_vacc_class[i, s_E, j],
         1 - p_sympt[i] * rel_p_sympt[i, j])

## Computes the number of symptomatic cases
n_EI_sympt[, ] <- n_EE[i, s_E, j] - n_EI_asympt[i, j]
n_EI_sympt_next_vacc_class[, ] <- n_EE_next_vacc_class[i, s_E, j] -
  n_EI_asympt_next_vacc_class[i, j]

  
## Work out the S->E and E->E transitions
aux_EE[, 1, ] <- n_SE[i, k]
aux_EE[, 2:s_E, ] <- n_EE[i, j - 1, k]
aux_EE[, , ] <- aux_EE[i, j, k] - n_EE[i, j, k] -
  n_EE_next_vacc_class[i, j, k] -
  n_E_next_vacc_class[i, j, k]
aux_EE[, , 1] <- aux_EE[i, j, 1] +
  n_E_next_vacc_class[i, j, n_vacc_classes]
aux_EE[, , 2:n_vacc_classes] <- aux_EE[i, j, k] +
  n_E_next_vacc_class[i, j, k - 1]
aux_EE[, 1, 1] <- aux_EE[i, 1, 1] +
  n_SE_next_vacc_class[i, n_vacc_classes]
aux_EE[, 1, 2:n_vacc_classes] <- aux_EE[i, j, k] +
  n_SE_next_vacc_class[i, k - 1]
aux_EE[, 2:s_E, 1] <- aux_EE[i, j, k] +
  n_EE_next_vacc_class[i, j - 1, n_vacc_classes]
aux_EE[, 2:s_E, 2:n_vacc_classes] <- aux_EE[i, j, k] +
  n_EE_next_vacc_class[i, j - 1, k - 1]
new_E[, , ] <- E[i, j, k] + aux_EE[i, j, k]

## Work out the I_asympt->I_asympt transitions
aux_II_asympt[, 1, ] <- n_EI_asympt[i, k]
aux_II_asympt[, 2:s_asympt, ] <- n_II_asympt[i, j - 1, k]
aux_II_asympt[, , ] <- aux_II_asympt[i, j, k] - n_II_asympt[i, j, k] -
  n_II_asympt_next_vacc_class[i, j, k] -
  n_I_asympt_next_vacc_class[i, j, k]
aux_II_asympt[, , 1] <- aux_II_asympt[i, j, k]  +
  n_I_asympt_next_vacc_class[i, 1, n_vacc_classes]
aux_II_asympt[, , 2:n_vacc_classes] <- aux_II_asympt[i, j, k] +
  n_I_asympt_next_vacc_class[i, j, k - 1]
aux_II_asympt[, 1, 1] <- aux_II_asympt[i, j, k] +
  n_EI_asympt_next_vacc_class[i, n_vacc_classes]
aux_II_asympt[, 1, 2:n_vacc_classes] <- aux_II_asympt[i, j, k] +
  n_EI_asympt_next_vacc_class[i, k - 1]
aux_II_asympt[, 2:s_asympt, 1] <- aux_II_asympt[i, j, k] +
  n_II_asympt_next_vacc_class[i, j - 1, n_vacc_classes]
aux_II_asympt[, 2:s_asympt, 2:n_vacc_classes] <- aux_II_asympt[i, j, k] +
  n_II_asympt_next_vacc_class[i, j - 1, k - 1]
new_I_asympt[, , ] <- I_asympt[i, j, k] + aux_II_asympt[i, j, k]

## Work out the I_sympt->I_sympt transitions
aux_II_sympt[, 1, 1] <- n_EI_sympt[i, 1] +
  n_EI_sympt_next_vacc_class[i, n_vacc_classes]
aux_II_sympt[, 1, 2:n_vacc_classes] <-
  n_EI_sympt[i, k] + n_EI_sympt_next_vacc_class[i, k - 1]

aux_II_sympt[, 2:s_sympt, ] <- n_II_sympt[i, j - 1, k]
aux_II_sympt[, 1:s_sympt, ] <- aux_II_sympt[i, j, k] - n_II_sympt[i, j, k]
new_I_sympt[, , ] <- I_sympt[i, j, k] + aux_II_sympt[i, j, k]

## Work out the flow from I_sympt -> R, comm_D, hosp
n_sympt_to_R[, ] <- rbinom(n_II_sympt[i, s_sympt, j],
                         1 - prob_hosp_sympt[i] * rel_p_hosp_if_sympt[i, j])
n_sympt_to_comm_D[, ] <-
  rbinom(n_II_sympt[i, s_sympt, j] - n_sympt_to_R[i, j], prob_death_comm[i])
n_sympt_to_hosp[, ] <-
  n_II_sympt[i, s_sympt, j] - n_sympt_to_R[i, j] - n_sympt_to_comm_D[i, j]

## Work out the I_comm_D -> I_comm_D transitions
aux_II_comm_D[, 1, ] <- n_sympt_to_comm_D[i, k]
aux_II_comm_D[, 2:s_comm_D, ] <- n_II_comm_D[i, j - 1, k]
aux_II_comm_D[, 1:s_comm_D, ] <- aux_II_comm_D[i, j, k] - n_II_comm_D[i, j, k]
new_I_comm_D[, , ] <- I_comm_D[i, j, k] + aux_II_comm_D[i, j, k]

## Work out the split in hospitals between hosp_D, hosp_R and triage
n_sympt_to_triage[, ] <- rbinom(n_sympt_to_hosp[i, j], prob_ICU_hosp[i])
n_sympt_to_triage_conf[, ] <- rbinom(n_sympt_to_triage[i, j],
                                   prob_admit_conf[i])
n_hosp_non_ICU[, ] <- n_sympt_to_hosp[i, j] - n_sympt_to_triage[i, j]
n_sympt_to_hosp_D[, ] <- rbinom(n_hosp_non_ICU[i, j], prob_death_hosp_D[i])
n_sympt_to_hosp_D_conf[, ] <- rbinom(n_sympt_to_hosp_D[i, j],
                                     prob_admit_conf[i])
n_sympt_to_hosp_R[, ] <- n_hosp_non_ICU[i, j] - n_sympt_to_hosp_D[i, j]
n_sympt_to_hosp_R_conf[, ] <- rbinom(n_sympt_to_hosp_R[i, j],
                                     prob_admit_conf[i])

## Work out the I_triage -> I_triage transitions
aux_II_triage_unconf[, , ] <- I_triage_unconf[i, j, k]
aux_II_triage_unconf[, 2:s_triage, ] <-
  aux_II_triage_unconf[i, j, k] + n_II_triage_unconf[i, j - 1, k]
aux_II_triage_unconf[, 1:s_triage, ] <-
  aux_II_triage_unconf[i, j, k] - n_II_triage_unconf[i, j, k]
aux_II_triage_conf[, , ] <-
  I_triage_conf[i, j, k]
aux_II_triage_conf[, 2:s_triage, ] <-
  aux_II_triage_conf[i, j, k] + n_II_triage_conf[i, j - 1, k]
aux_II_triage_conf[, 1:s_triage, ] <-
  aux_II_triage_conf[i, j, k] - n_II_triage_conf[i, j, k]
n_I_triage_unconf_to_conf[, , ] <-
  rbinom(aux_II_triage_unconf[i, j, k], p_test)
new_I_triage_unconf[, , ] <-
  aux_II_triage_unconf[i, j, k] - n_I_triage_unconf_to_conf[i, j, k]
new_I_triage_unconf[, 1, ] <-
  new_I_triage_unconf[i, 1, k] + n_sympt_to_triage[i, k] -
  n_sympt_to_triage_conf[i, k]
new_I_triage_conf[, , ] <-
  aux_II_triage_conf[i, j, k] + n_I_triage_unconf_to_conf[i, j, k]
new_I_triage_conf[, 1, ] <-
  new_I_triage_conf[i, 1, k] + n_sympt_to_triage_conf[i, k]

## Work out the I_hosp_R->I_hosp_R transitions
aux_II_hosp_R_unconf[, , ] <- I_hosp_R_unconf[i, j, k]
aux_II_hosp_R_unconf[, 2:s_hosp_R, ] <-
  aux_II_hosp_R_unconf[i, j, k] + n_II_hosp_R_unconf[i, j - 1, k]
aux_II_hosp_R_unconf[, 1:s_hosp_R, ] <-
  aux_II_hosp_R_unconf[i, j, k] - n_II_hosp_R_unconf[i, j, k]
aux_II_hosp_R_conf[, , ] <- I_hosp_R_conf[i, j, k]
aux_II_hosp_R_conf[, 2:s_hosp_R, ] <-
  aux_II_hosp_R_conf[i, j, k] + n_II_hosp_R_conf[i, j - 1, k]
aux_II_hosp_R_conf[, 1:s_hosp_R, ] <-
  aux_II_hosp_R_conf[i, j, k] - n_II_hosp_R_conf[i, j, k]
n_I_hosp_R_unconf_to_conf[, , ] <-
  rbinom(aux_II_hosp_R_unconf[i, j, k], p_test)
new_I_hosp_R_unconf[, , ] <-
  aux_II_hosp_R_unconf[i, j, k] - n_I_hosp_R_unconf_to_conf[i, j, k]
new_I_hosp_R_unconf[, 1, ] <-
  new_I_hosp_R_unconf[i, 1, k] + n_sympt_to_hosp_R[i, k] -
  n_sympt_to_hosp_R_conf[i, k]
new_I_hosp_R_conf[, , ] <-
  aux_II_hosp_R_conf[i, j, k] + n_I_hosp_R_unconf_to_conf[i, j, k]
new_I_hosp_R_conf[, 1, ] <-
  new_I_hosp_R_conf[i, 1, k] + n_sympt_to_hosp_R_conf[i, k]

## Work out the I_hosp_D->I_hosp_D transitions
aux_II_hosp_D_unconf[, , ] <- I_hosp_D_unconf[i, j, k]
aux_II_hosp_D_unconf[, 2:s_hosp_D, ] <-
  aux_II_hosp_D_unconf[i, j, k] + n_II_hosp_D_unconf[i, j - 1, k]
aux_II_hosp_D_unconf[, 1:s_hosp_D, ] <-
  aux_II_hosp_D_unconf[i, j, k] - n_II_hosp_D_unconf[i, j, k]
aux_II_hosp_D_conf[, , ] <- I_hosp_D_conf[i, j, k]
aux_II_hosp_D_conf[, 2:s_hosp_D, ] <-
  aux_II_hosp_D_conf[i, j, k] + n_II_hosp_D_conf[i, j - 1, k]
aux_II_hosp_D_conf[, 1:s_hosp_D, ] <-
  aux_II_hosp_D_conf[i, j, k] - n_II_hosp_D_conf[i, j, k]
n_I_hosp_D_unconf_to_conf[, , ] <-
  rbinom(aux_II_hosp_D_unconf[i, j, k], p_test)
new_I_hosp_D_unconf[, , ] <-
  aux_II_hosp_D_unconf[i, j, k] - n_I_hosp_D_unconf_to_conf[i, j, k]
new_I_hosp_D_unconf[, 1, ] <-
  new_I_hosp_D_unconf[i, 1, k] + n_sympt_to_hosp_D[i, k] -
  n_sympt_to_hosp_D_conf[i, k]
new_I_hosp_D_conf[, , ] <-
  aux_II_hosp_D_conf[i, j, k] + n_I_hosp_D_unconf_to_conf[i, j, k]
new_I_hosp_D_conf[, 1, ] <-
  new_I_hosp_D_conf[i, 1, k] + n_sympt_to_hosp_D_conf[i, k]

## Work out the triage to ICU_D, ICU_S_R and ICU_S_D splits
n_triage_unconf_to_ICU_D_unconf[, ] <-
  rbinom(n_II_triage_unconf[i, s_triage, j], prob_death_ICU[i])
n_triage_conf_to_ICU_D_conf[, ] <-
  rbinom(n_II_triage_conf[i, s_triage, j], prob_death_ICU[i])
n_triage_unconf_to_ICU_S_D_unconf[, ] <-
  rbinom(n_II_triage_unconf[i, s_triage, j] -
           n_triage_unconf_to_ICU_D_unconf[i, j],
         prob_death_stepdown[i])
n_triage_unconf_to_ICU_S_R_unconf[, ] <- n_II_triage_unconf[i, s_triage, j] -
  n_triage_unconf_to_ICU_D_unconf[i, j] -
  n_triage_unconf_to_ICU_S_D_unconf[i, j]
n_triage_conf_to_ICU_S_D_conf[, ] <-
  rbinom(n_II_triage_conf[i, s_triage, j] - n_triage_conf_to_ICU_D_conf[i, j],
         prob_death_stepdown[i])
n_triage_conf_to_ICU_S_R_conf[, ] <- n_II_triage_conf[i, s_triage, j] -
  n_triage_conf_to_ICU_D_conf[i, j] - n_triage_conf_to_ICU_S_D_conf[i, j]


## Work out the I_ICU_S_R->I_ICU_S_R transitions
aux_II_ICU_S_R_unconf[, , ] <- I_ICU_S_R_unconf[i, j, k]
aux_II_ICU_S_R_unconf[, 1, ] <-
  aux_II_ICU_S_R_unconf[i, j, k] + n_triage_unconf_to_ICU_S_R_unconf[i, k]
aux_II_ICU_S_R_unconf[, 2:s_ICU_S_R, ] <-
  aux_II_ICU_S_R_unconf[i, j, k] + n_II_ICU_S_R_unconf[i, j - 1, k]
aux_II_ICU_S_R_unconf[, 1:s_ICU_S_R, ] <-
  aux_II_ICU_S_R_unconf[i, j, k] - n_II_ICU_S_R_unconf[i, j, k]
aux_II_ICU_S_R_conf[, , ] <- I_ICU_S_R_conf[i, j, k]
aux_II_ICU_S_R_conf[, 1, ] <-
  aux_II_ICU_S_R_conf[i, j, k] + n_triage_conf_to_ICU_S_R_conf[i, k]
aux_II_ICU_S_R_conf[, 2:s_ICU_S_R, ] <-
  aux_II_ICU_S_R_conf[i, j, k] + n_II_ICU_S_R_conf[i, j - 1, k]
aux_II_ICU_S_R_conf[, 1:s_ICU_S_R, ] <-
  aux_II_ICU_S_R_conf[i, j, k] - n_II_ICU_S_R_conf[i, j, k]
n_I_ICU_S_R_unconf_to_conf[, , ] <-
  rbinom(aux_II_ICU_S_R_unconf[i, j, k], p_test)
new_I_ICU_S_R_unconf[, , ] <-
  aux_II_ICU_S_R_unconf[i, j, k] - n_I_ICU_S_R_unconf_to_conf[i, j, k]
new_I_ICU_S_R_conf[, , ] <-
  aux_II_ICU_S_R_conf[i, j, k] + n_I_ICU_S_R_unconf_to_conf[i, j, k]

## Work out the I_ICU_S_D->I_ICU_S_D transitions
aux_II_ICU_S_D_unconf[, , ] <- I_ICU_S_D_unconf[i, j, k]
aux_II_ICU_S_D_unconf[, 1, ] <-
  aux_II_ICU_S_D_unconf[i, j, k] + n_triage_unconf_to_ICU_S_D_unconf[i, k]
aux_II_ICU_S_D_unconf[, 2:s_ICU_S_D, ] <-
  aux_II_ICU_S_D_unconf[i, j, k] + n_II_ICU_S_D_unconf[i, j - 1, k]
aux_II_ICU_S_D_unconf[, 1:s_ICU_S_D, ] <-
  aux_II_ICU_S_D_unconf[i, j, k] - n_II_ICU_S_D_unconf[i, j, k]
aux_II_ICU_S_D_conf[, , ] <- I_ICU_S_D_conf[i, j, k]
aux_II_ICU_S_D_conf[, 1, ] <-
  aux_II_ICU_S_D_conf[i, j, k] + n_triage_conf_to_ICU_S_D_conf[i, k]
aux_II_ICU_S_D_conf[, 2:s_ICU_S_D, ] <-
  aux_II_ICU_S_D_conf[i, j, k] + n_II_ICU_S_D_conf[i, j - 1, k]
aux_II_ICU_S_D_conf[, 1:s_ICU_S_D, ] <-
  aux_II_ICU_S_D_conf[i, j, k] - n_II_ICU_S_D_conf[i, j, k]
n_I_ICU_S_D_unconf_to_conf[, , ] <-
  rbinom(aux_II_ICU_S_D_unconf[i, j, k], p_test)
new_I_ICU_S_D_unconf[, , ] <-
  aux_II_ICU_S_D_unconf[i, j, k] - n_I_ICU_S_D_unconf_to_conf[i, j, k]
new_I_ICU_S_D_conf[, , ] <-
  aux_II_ICU_S_D_conf[i, j, k] + n_I_ICU_S_D_unconf_to_conf[i, j, k]

## Work out the I_ICU_D->I_ICU_D transitions
aux_II_ICU_D_unconf[, , ] <- I_ICU_D_unconf[i, j, k]
aux_II_ICU_D_unconf[, 1, ] <-
  aux_II_ICU_D_unconf[i, j, k] + n_triage_unconf_to_ICU_D_unconf[i, k]
aux_II_ICU_D_unconf[, 2:s_ICU_D, ] <-
  aux_II_ICU_D_unconf[i, j, k] + n_II_ICU_D_unconf[i, j - 1, k]
aux_II_ICU_D_unconf[, 1:s_ICU_D, ] <-
  aux_II_ICU_D_unconf[i, j, k] - n_II_ICU_D_unconf[i, j, k]
aux_II_ICU_D_conf[, , ] <- I_ICU_D_conf[i, j, k]
aux_II_ICU_D_conf[, 1, ] <-
  aux_II_ICU_D_conf[i, j, k] + n_triage_conf_to_ICU_D_conf[i, k]
aux_II_ICU_D_conf[, 2:s_ICU_D, ] <-
  aux_II_ICU_D_conf[i, j, k] + n_II_ICU_D_conf[i, j - 1, k]
aux_II_ICU_D_conf[, 1:s_ICU_D, ] <-
  aux_II_ICU_D_conf[i, j, k] - n_II_ICU_D_conf[i, j, k]
n_I_ICU_D_unconf_to_conf[, , ] <-
  rbinom(aux_II_ICU_D_unconf[i, j, k], p_test)
new_I_ICU_D_unconf[, , ] <-
  aux_II_ICU_D_unconf[i, j, k] - n_I_ICU_D_unconf_to_conf[i, j, k]
new_I_ICU_D_conf[, , ] <-
  aux_II_ICU_D_conf[i, j, k] + n_I_ICU_D_unconf_to_conf[i, j, k]

## Work out the R_stepdown_R->R_stepdown_R transitions
aux_R_stepdown_R_unconf[, , ] <- R_stepdown_R_unconf[i, j, k]
aux_R_stepdown_R_unconf[, 1, ] <-
  aux_R_stepdown_R_unconf[i, j, k] + n_II_ICU_S_R_unconf[i, s_ICU_S_R, k]
aux_R_stepdown_R_unconf[, 2:s_stepdown_R, ] <-
  aux_R_stepdown_R_unconf[i, j, k] + n_R_stepdown_R_unconf[i, j - 1, k]
aux_R_stepdown_R_unconf[, 1:s_stepdown_R, ] <-
  aux_R_stepdown_R_unconf[i, j, k] - n_R_stepdown_R_unconf[i, j, k]
aux_R_stepdown_R_conf[, , ] <- R_stepdown_R_conf[i, j, k]
aux_R_stepdown_R_conf[, 1, ] <-
  aux_R_stepdown_R_conf[i, j, k] + n_II_ICU_S_R_conf[i, s_ICU_S_R, k]
aux_R_stepdown_R_conf[, 2:s_stepdown_R, ] <-
  aux_R_stepdown_R_conf[i, j, k] + n_R_stepdown_R_conf[i, j - 1, k]
aux_R_stepdown_R_conf[, 1:s_stepdown_R, ] <-
  aux_R_stepdown_R_conf[i, j, k] - n_R_stepdown_R_conf[i, j, k]
n_R_stepdown_R_unconf_to_conf[, , ] <-
  rbinom(aux_R_stepdown_R_unconf[i, j, k], p_test)
new_R_stepdown_R_unconf[, , ] <-
  aux_R_stepdown_R_unconf[i, j, k] - n_R_stepdown_R_unconf_to_conf[i, j, k]
new_R_stepdown_R_conf[, , ] <-
  aux_R_stepdown_R_conf[i, j, k] + n_R_stepdown_R_unconf_to_conf[i, j, k]

## Work out the R_stepdown_D->R_stepdown_D transitions
aux_R_stepdown_D_unconf[, , ] <- R_stepdown_D_unconf[i, j, k]
aux_R_stepdown_D_unconf[, 1, ] <-
  aux_R_stepdown_D_unconf[i, j, k] + n_II_ICU_S_D_unconf[i, s_ICU_S_D, k]
aux_R_stepdown_D_unconf[, 2:s_stepdown_D, ] <-
  aux_R_stepdown_D_unconf[i, j, k] + n_R_stepdown_D_unconf[i, j - 1, k]
aux_R_stepdown_D_unconf[, 1:s_stepdown_D, ] <-
  aux_R_stepdown_D_unconf[i, j, k] - n_R_stepdown_D_unconf[i, j, k]
aux_R_stepdown_D_conf[, , ] <- R_stepdown_D_conf[i, j, k]
aux_R_stepdown_D_conf[, 1, ] <-
  aux_R_stepdown_D_conf[i, j, k] + n_II_ICU_S_D_conf[i, s_ICU_S_D, k]
aux_R_stepdown_D_conf[, 2:s_stepdown_D, ] <-
  aux_R_stepdown_D_conf[i, j, k] + n_R_stepdown_D_conf[i, j - 1, k]
aux_R_stepdown_D_conf[, 1:s_stepdown_D, ] <-
  aux_R_stepdown_D_conf[i, j, k] - n_R_stepdown_D_conf[i, j, k]
n_R_stepdown_D_unconf_to_conf[, , ] <-
  rbinom(aux_R_stepdown_D_unconf[i, j, k], p_test)
new_R_stepdown_D_unconf[, , ] <-
  aux_R_stepdown_D_unconf[i, j, k] - n_R_stepdown_D_unconf_to_conf[i, j, k]
new_R_stepdown_D_conf[, , ] <-
  aux_R_stepdown_D_conf[i, j, k] + n_R_stepdown_D_unconf_to_conf[i, j, k]

## Work out the number of deaths in hospital
new_D_hosp[] <- D_hosp[i] +
  sum(n_II_hosp_D_unconf[i, s_hosp_D, ]) +
  sum(n_II_hosp_D_conf[i, s_hosp_D, ]) +
  sum(n_II_ICU_D_unconf[i, s_ICU_D, ]) +
  sum(n_II_ICU_D_conf[i, s_ICU_D, ]) +
  sum(n_R_stepdown_D_unconf[i, s_stepdown_D, ]) +
  sum(n_R_stepdown_D_conf[i, s_stepdown_D, ])

## Work out the number of deaths in the community
new_D_comm[] <- D_comm[i] + sum(n_II_comm_D[i, s_comm_D, ])

## Work out the number of people entering the seroconversion flow
n_com_to_R_pre[, 1, 1] <- rbinom(
  n_EE[i, s_E, 1] + n_EE_next_vacc_class[i, s_E, n_vacc_classes], p_R_pre_1)
n_com_to_R_pre[, 1, 2:n_vacc_classes] <- rbinom(
  n_EE[i, s_E, k] + n_EE_next_vacc_class[i, s_E, k - 1], p_R_pre_1)
n_com_to_R_pre[, 2, 1] <- n_EE[i, s_E, 1] +
  n_EE_next_vacc_class[i, s_E, n_vacc_classes] - n_com_to_R_pre[i, 1, 1]
n_com_to_R_pre[, 2, 2:n_vacc_classes] <- n_EE[i, s_E, k] +
  n_EE_next_vacc_class[i, s_E, k - 1] - n_com_to_R_pre[i, 1, k]
new_R_pre[, , ] <- R_pre[i, j, k] + n_com_to_R_pre[i, j, k] - n_R_pre[i, j, k]


## Split the seroconversion flow between people who are going to
## seroconvert and people who are not
n_R_pre_to_R_pos[, ] <- rbinom(sum(n_R_pre[i, , j]), p_seroconversion[i])

new_R_pos[, , ] <- R_pos[i, j, k] - n_R_pos[i, j, k]
new_R_pos[, 1, ] <- new_R_pos[i, 1, k] + n_R_pre_to_R_pos[i, k]
new_R_pos[, 2:s_R_pos, ] <- new_R_pos[i, j, k] + n_R_pos[i, j - 1, k]

new_R_neg[, ] <- R_neg[i, j] + sum(n_R_pre[i, , j]) - n_R_pre_to_R_pos[i, j] +
  n_R_pos[i, s_R_pos, j] - model_pcr_and_serology * n_R_progress[i, j] -
  model_pcr_and_serology * n_R_next_vacc_class[i, j]
new_R_neg[, 1] <- new_R_neg[i, 1] +
  model_pcr_and_serology * n_R_next_vacc_class[i, n_vacc_classes]
new_R_neg[, 2:n_vacc_classes] <- new_R_neg[i, j] +
  model_pcr_and_serology * n_R_next_vacc_class[i, j - 1]


## Work out the total number of recovery
new_R[, ] <- R[i, j] +
  n_II_asympt[i, s_asympt, j] +
  n_sympt_to_R[i, j] +
  n_II_hosp_R_conf[i, s_hosp_R, j] +
  n_II_hosp_R_unconf[i, s_hosp_R, j] +
  n_R_stepdown_R_conf[i, s_stepdown_R, j] +
  n_R_stepdown_R_unconf[i, s_stepdown_R, j] -
  n_R_progress[i, j] -
  n_R_next_vacc_class[i, j]
new_R[, 1] <- new_R[i, 1] +
  n_II_asympt_next_vacc_class[i, s_asympt, n_vacc_classes] +
  n_R_next_vacc_class[i, n_vacc_classes]
new_R[, 2:n_vacc_classes] <- new_R[i, j] +
  n_II_asympt_next_vacc_class[i, s_asympt, j - 1] +
  n_R_next_vacc_class[i, j - 1]



## Work out the PCR positivity
new_PCR_pre[, , ] <- PCR_pre[i, j, k] - n_PCR_pre[i, j, k]
new_PCR_pre[, 1, ] <- new_PCR_pre[i, 1, k] + n_S_progress[i, k]
new_PCR_pre[, 2:s_PCR_pre, ] <- new_PCR_pre[i, j, k] + n_PCR_pre[i, j - 1, k]

new_PCR_pos[, , ] <- PCR_pos[i, j, k] - n_PCR_pos[i, j, k]
new_PCR_pos[, 1, ] <- new_PCR_pos[i, 1, k] + n_PCR_pre[i, s_PCR_pre, k]
new_PCR_pos[, 2:s_PCR_pos, ] <- new_PCR_pos[i, j, k] + n_PCR_pos[i, j - 1, k]

new_PCR_neg[, ] <- PCR_neg[i, j] + n_PCR_pos[i, s_PCR_pos, j] -
  model_pcr_and_serology * n_R_progress[i, j] -
  model_pcr_and_serology * n_R_next_vacc_class[i, j]
new_PCR_neg[, 1] <- new_PCR_neg[i, 1] +
  model_pcr_and_serology * n_R_next_vacc_class[i, n_vacc_classes]
new_PCR_neg[, 2:n_vacc_classes] <- new_PCR_neg[i, j] +
  model_pcr_and_serology * n_R_next_vacc_class[i, j - 1]


## Compute the force of infection
I_with_diff_trans[, ] <-
  (sum(I_asympt[i, , j]) + sum(I_sympt[i, , j]) +
    hosp_transmission * (
      sum(I_triage_unconf[i, , j]) +
      sum(I_triage_conf[i, , j]) +
      sum(I_hosp_R_unconf[i, , j]) +
      sum(I_hosp_R_conf[i, , j]) +
      sum(I_hosp_D_unconf[i, , j]) +
      sum(I_hosp_D_conf[i, , j])) +
    ICU_transmission * (
      sum(I_ICU_S_R_unconf[i, , j]) +
      sum(I_ICU_S_R_conf[i, , j]) +
      sum(I_ICU_S_D_unconf[i, , j]) +
      sum(I_ICU_S_D_conf[i, , j]) +
      sum(I_ICU_D_unconf[i, , j]) +
      sum(I_ICU_D_conf[i, , j])) +
    comm_D_transmission * sum(I_comm_D[i, , j]))


## NOTE: "age groups" 1-17 are age groups, 18 are CHW and 19 CHR. Here we apply
## beta to all contacts *except* within care home contacts
s_ij[, ] <- m[i, j] * sum(I_with_diff_trans[j, ])
s_ij[1:n_age_groups, 1:n_groups] <- beta * s_ij[i, j]
s_ij[(n_age_groups + 1):n_groups, 1:n_age_groups] <- beta * s_ij[i, j]
lambda[] <- sum(s_ij[i, ])

## Initial states are all zerod as we will provide a state vector
## setting S and I based on the seeding model.
initial(S[, ]) <- 0
initial(E[, , ]) <- 0
initial(I_asympt[, , ]) <- 0
initial(I_sympt[, , ]) <- 0
initial(I_comm_D[, , ]) <- 0
initial(I_triage_unconf[, , ]) <- 0
initial(I_triage_conf[, , ]) <- 0
initial(I_hosp_R_unconf[, , ]) <- 0
initial(I_hosp_R_conf[, , ]) <- 0
initial(I_hosp_D_unconf[, , ]) <- 0
initial(I_hosp_D_conf[, , ]) <- 0
initial(I_ICU_S_R_unconf[, , ]) <- 0
initial(I_ICU_S_R_conf[, , ]) <- 0
initial(I_ICU_S_D_unconf[, , ]) <- 0
initial(I_ICU_S_D_conf[, , ]) <- 0
initial(I_ICU_D_unconf[, , ]) <- 0
initial(I_ICU_D_conf[, , ]) <- 0
initial(R_stepdown_R_unconf[, , ]) <- 0
initial(R_stepdown_R_conf[, , ]) <- 0
initial(R_stepdown_D_unconf[, , ]) <- 0
initial(R_stepdown_D_conf[, , ]) <- 0
initial(R_pre[, , ]) <- 0
initial(R_pos[, , ]) <- 0
initial(R_neg[, ]) <- 0
initial(R[, ]) <- 0
initial(D_hosp[]) <- 0
initial(D_comm[]) <- 0
initial(PCR_pre[, , ]) <- 0
initial(PCR_pos[, , ]) <- 0
initial(PCR_neg[, ]) <- 0
initial(cum_admit_conf) <- 0
initial(cum_new_conf) <- 0
initial(cum_admit_by_age[]) <- 0

## User defined parameters - default in parentheses:

## Vaccination parameters
rel_susceptibility[, ] <- user()
dim(rel_susceptibility) <- user() # use length as provided by the user
n_vacc_classes <- dim(rel_susceptibility, 2)
rel_p_sympt[, ] <- user()
dim(rel_p_sympt) <- c(n_groups, n_vacc_classes)
rel_p_hosp_if_sympt[, ] <- user()
dim(rel_p_hosp_if_sympt) <- c(n_groups, n_vacc_classes)

vaccine_progression_rate_base[, ] <- user()
dim(vaccine_progression_rate_base) <- c(n_groups, n_vacc_classes)

## Parameters of the E classes
s_E <- user()
gamma_E <- user(0.1)

## Probability of transitioning from the E to the asymptomatic class,
## the rest go into the symptomatic class
p_sympt[] <- user()

## Parameters of the I_asympt classes
s_asympt <- user()
gamma_asympt <- user(0.1)

## Parameters of the I_sympt classes
s_sympt <- user()
gamma_sympt <- user(0.1)
dim(p_hosp_sympt_step) <- user()
p_hosp_sympt_step[] <- user()
psi_hosp_sympt[] <- user()

## Parameters of the I_comm_D class
s_comm_D <- user()
gamma_comm_D <- user(0.1)
dim(p_death_comm_step) <- user()
p_death_comm_step[] <- user()
psi_death_comm[] <- user()

## Parameters of the I_triage classes
s_triage <- user()
gamma_triage <- user(0.1)

## Proportion of hospital cases progressing to ICU
dim(p_ICU_hosp_step) <- user()
p_ICU_hosp_step[] <- user()
psi_ICU_hosp[] <- user()

## Proportion of stepdown cases dying
dim(p_death_stepdown_step) <- user()
p_death_stepdown_step[] <- user()
psi_death_stepdown[] <- user()

## Parameters of the I_hosp_R classes
s_hosp_R <- user()
gamma_hosp_R <- user(0.1)

## Parameters of the I_hosp_D classes
s_hosp_D <- user()
gamma_hosp_D <- user(0.1)
dim(p_death_hosp_D_step) <- user()
p_death_hosp_D_step[] <- user()
psi_death_hosp_D[] <- user()

## Parameters of the I_ICU_S_R classes
s_ICU_S_R <- user()
gamma_ICU_S_R <- user(0.1)

## Parameters of the I_ICU_S_D classes
s_ICU_S_D <- user()
gamma_ICU_S_D <- user(0.1)

## Parameters of the I_ICU classes
s_ICU_D <- user()
gamma_ICU_D <- user(0.1)
dim(p_death_ICU_step) <- user()
p_death_ICU_step[] <- user()
psi_death_ICU[] <- user()

## Waning of immunity
waning_rate[] <- user()
dim(waning_rate) <- n_groups

## Parameters of the R_stepdown_R classes
s_stepdown_R <- user()
gamma_stepdown_R <- user(0.1)

## Parameters of the R_stepdown_D classes
s_stepdown_D <- user()
gamma_stepdown_D <- user(0.1)

## Parameters of the R_pre classes
gamma_R_pre_1 <- user(0.1)
gamma_R_pre_2 <- user(0.1)
gamma_R_pre[1] <- gamma_R_pre_1
gamma_R_pre[2] <- gamma_R_pre_2
## Governs the mixing - pretty much only makes sense at 0.5
p_R_pre_1 <- user(0.5)
p_seroconversion[] <- user()

# Parameters of the R_pos classes
s_R_pos <- user()
gamma_R_pos <- user(0.1)

## Parameters relating to testing
gamma_test <- user(0.1)
dim(p_admit_conf_step) <- user()
p_admit_conf_step[] <- user()
psi_admit_conf[] <- user()

## Parameters relating to PCR positivity
s_PCR_pre <- user()
gamma_PCR_pre <- user(0.1)
s_PCR_pos <- user()
gamma_PCR_pos <- user(0.1)

## Parameters of the age stratified transmission
beta_step[] <- user()
dim(beta_step) <- user()
## What we really want is min(step + 1, length(beta_step)) but that's not
## supported by odin (it could be made to support this). This code
## does currently create a compiler warning with -Wsign-compare on
## because we have an unsigned/signed integer comparison
beta <- if (step >= length(beta_step))
          beta_step[length(beta_step)] else beta_step[step + 1]

## Useful for debugging
initial(beta_out) <- beta_step[1]
update(beta_out) <- beta

m[, ] <- user()
hosp_transmission <- user()
ICU_transmission <- user()
comm_D_transmission <- user()

## Dimensions of the different "vectors" here vectors stand for
## multi-dimensional arrays

## Vectors handling the S class
dim(S) <- c(n_groups, n_vacc_classes)
dim(new_S) <- c(n_groups, n_vacc_classes)

## Vectors handling the E class
dim(E) <- c(n_groups, s_E, n_vacc_classes)
dim(aux_EE) <- c(n_groups, s_E, n_vacc_classes)
dim(new_E) <- c(n_groups, s_E, n_vacc_classes)
dim(n_EE) <- c(n_groups, s_E, n_vacc_classes)

## Vectors handling the I_asympt class
dim(I_asympt) <- c(n_groups, s_asympt, n_vacc_classes)
dim(aux_II_asympt) <- c(n_groups, s_asympt, n_vacc_classes)
dim(new_I_asympt) <- c(n_groups, s_asympt, n_vacc_classes)
dim(n_II_asympt) <- c(n_groups, s_asympt, n_vacc_classes)

## Vectors handling the I_sympt class
dim(I_sympt) <- c(n_groups, s_sympt, n_vacc_classes)
dim(aux_II_sympt) <- c(n_groups, s_sympt, n_vacc_classes)
dim(new_I_sympt) <- c(n_groups, s_sympt, n_vacc_classes)
dim(n_II_sympt) <- c(n_groups, s_sympt, n_vacc_classes)
dim(prob_hosp_sympt) <- n_groups
dim(psi_hosp_sympt) <- n_groups

## Vectors handling the I_comm_D class
dim(I_comm_D) <- c(n_groups, s_comm_D, n_vacc_classes)
dim(aux_II_comm_D) <- c(n_groups, s_comm_D, n_vacc_classes)
dim(new_I_comm_D) <- c(n_groups, s_comm_D, n_vacc_classes)
dim(n_II_comm_D) <- c(n_groups, s_comm_D, n_vacc_classes)
dim(prob_death_comm) <- n_groups
dim(psi_death_comm) <- n_groups

## Vectors handling the I_triage class
dim(I_triage_unconf) <- c(n_groups, s_triage, n_vacc_classes)
dim(aux_II_triage_unconf) <- c(n_groups, s_triage, n_vacc_classes)
dim(new_I_triage_unconf) <- c(n_groups, s_triage, n_vacc_classes)
dim(n_II_triage_unconf) <- c(n_groups, s_triage, n_vacc_classes)
dim(I_triage_conf) <- c(n_groups, s_triage, n_vacc_classes)
dim(aux_II_triage_conf) <- c(n_groups, s_triage, n_vacc_classes)
dim(new_I_triage_conf) <- c(n_groups, s_triage, n_vacc_classes)
dim(n_II_triage_conf) <- c(n_groups, s_triage, n_vacc_classes)
dim(n_I_triage_unconf_to_conf) <- c(n_groups, s_triage, n_vacc_classes)

## Vector handling who progress to ICU
dim(prob_ICU_hosp) <- n_groups
dim(psi_ICU_hosp) <- n_groups

## Vectors handling the I_hosp_R class
dim(I_hosp_R_unconf) <- c(n_groups, s_hosp_R, n_vacc_classes)
dim(aux_II_hosp_R_unconf) <- c(n_groups, s_hosp_R, n_vacc_classes)
dim(new_I_hosp_R_unconf) <- c(n_groups, s_hosp_R, n_vacc_classes)
dim(n_II_hosp_R_unconf) <- c(n_groups, s_hosp_R, n_vacc_classes)
dim(I_hosp_R_conf) <- c(n_groups, s_hosp_R, n_vacc_classes)
dim(aux_II_hosp_R_conf) <- c(n_groups, s_hosp_R, n_vacc_classes)
dim(new_I_hosp_R_conf) <- c(n_groups, s_hosp_R, n_vacc_classes)
dim(n_II_hosp_R_conf) <- c(n_groups, s_hosp_R, n_vacc_classes)
dim(n_I_hosp_R_unconf_to_conf) <- c(n_groups, s_hosp_R, n_vacc_classes)

## Vectors handling the I_hosp_D class
dim(I_hosp_D_unconf) <- c(n_groups, s_hosp_D, n_vacc_classes)
dim(aux_II_hosp_D_unconf) <- c(n_groups, s_hosp_D, n_vacc_classes)
dim(new_I_hosp_D_unconf) <- c(n_groups, s_hosp_D, n_vacc_classes)
dim(n_II_hosp_D_unconf) <- c(n_groups, s_hosp_D, n_vacc_classes)
dim(I_hosp_D_conf) <- c(n_groups, s_hosp_D, n_vacc_classes)
dim(aux_II_hosp_D_conf) <- c(n_groups, s_hosp_D, n_vacc_classes)
dim(new_I_hosp_D_conf) <- c(n_groups, s_hosp_D, n_vacc_classes)
dim(n_II_hosp_D_conf) <- c(n_groups, s_hosp_D, n_vacc_classes)
dim(n_I_hosp_D_unconf_to_conf) <- c(n_groups, s_hosp_D, n_vacc_classes)

## Vectors handling the I_ICU_S_R class
dim(I_ICU_S_R_unconf) <- c(n_groups, s_ICU_S_R, n_vacc_classes)
dim(aux_II_ICU_S_R_unconf) <- c(n_groups, s_ICU_S_R, n_vacc_classes)
dim(new_I_ICU_S_R_unconf) <- c(n_groups, s_ICU_S_R, n_vacc_classes)
dim(n_II_ICU_S_R_unconf) <- c(n_groups, s_ICU_S_R, n_vacc_classes)
dim(I_ICU_S_R_conf) <- c(n_groups, s_ICU_S_R, n_vacc_classes)
dim(aux_II_ICU_S_R_conf) <- c(n_groups, s_ICU_S_R, n_vacc_classes)
dim(new_I_ICU_S_R_conf) <- c(n_groups, s_ICU_S_R, n_vacc_classes)
dim(n_II_ICU_S_R_conf) <- c(n_groups, s_ICU_S_R, n_vacc_classes)
dim(n_I_ICU_S_R_unconf_to_conf) <- c(n_groups, s_ICU_S_R, n_vacc_classes)

## Vectors handling the I_ICU_S_D class
dim(I_ICU_S_D_unconf) <- c(n_groups, s_ICU_S_D, n_vacc_classes)
dim(aux_II_ICU_S_D_unconf) <- c(n_groups, s_ICU_S_D, n_vacc_classes)
dim(new_I_ICU_S_D_unconf) <- c(n_groups, s_ICU_S_D, n_vacc_classes)
dim(n_II_ICU_S_D_unconf) <- c(n_groups, s_ICU_S_D, n_vacc_classes)
dim(I_ICU_S_D_conf) <- c(n_groups, s_ICU_S_D, n_vacc_classes)
dim(aux_II_ICU_S_D_conf) <- c(n_groups, s_ICU_S_D, n_vacc_classes)
dim(new_I_ICU_S_D_conf) <- c(n_groups, s_ICU_S_D, n_vacc_classes)
dim(n_II_ICU_S_D_conf) <- c(n_groups, s_ICU_S_D, n_vacc_classes)
dim(n_I_ICU_S_D_unconf_to_conf) <- c(n_groups, s_ICU_S_D, n_vacc_classes)

## Vectors handling the I_ICU_D class
dim(I_ICU_D_unconf) <- c(n_groups, s_ICU_D, n_vacc_classes)
dim(aux_II_ICU_D_unconf) <- c(n_groups, s_ICU_D, n_vacc_classes)
dim(new_I_ICU_D_unconf) <- c(n_groups, s_ICU_D, n_vacc_classes)
dim(n_II_ICU_D_unconf) <- c(n_groups, s_ICU_D, n_vacc_classes)
dim(I_ICU_D_conf) <- c(n_groups, s_ICU_D, n_vacc_classes)
dim(aux_II_ICU_D_conf) <- c(n_groups, s_ICU_D, n_vacc_classes)
dim(new_I_ICU_D_conf) <- c(n_groups, s_ICU_D, n_vacc_classes)
dim(n_II_ICU_D_conf) <- c(n_groups, s_ICU_D, n_vacc_classes)
dim(n_I_ICU_D_unconf_to_conf) <- c(n_groups, s_ICU_D, n_vacc_classes)

## Vectors handling the R_stepdown_R class
dim(R_stepdown_R_unconf) <- c(n_groups, s_stepdown_R, n_vacc_classes)
dim(aux_R_stepdown_R_unconf) <- c(n_groups, s_stepdown_R, n_vacc_classes)
dim(new_R_stepdown_R_unconf) <- c(n_groups, s_stepdown_R, n_vacc_classes)
dim(n_R_stepdown_R_unconf) <- c(n_groups, s_stepdown_R, n_vacc_classes)
dim(R_stepdown_R_conf) <- c(n_groups, s_stepdown_R, n_vacc_classes)
dim(aux_R_stepdown_R_conf) <- c(n_groups, s_stepdown_R, n_vacc_classes)
dim(new_R_stepdown_R_conf) <- c(n_groups, s_stepdown_R, n_vacc_classes)
dim(n_R_stepdown_R_conf) <- c(n_groups, s_stepdown_R, n_vacc_classes)
dim(n_R_stepdown_R_unconf_to_conf) <- c(n_groups, s_stepdown_R, n_vacc_classes)

## Vectors handling the R_stepdown_D class
dim(R_stepdown_D_unconf) <- c(n_groups, s_stepdown_D, n_vacc_classes)
dim(aux_R_stepdown_D_unconf) <- c(n_groups, s_stepdown_D, n_vacc_classes)
dim(new_R_stepdown_D_unconf) <- c(n_groups, s_stepdown_D, n_vacc_classes)
dim(n_R_stepdown_D_unconf) <- c(n_groups, s_stepdown_D, n_vacc_classes)
dim(R_stepdown_D_conf) <- c(n_groups, s_stepdown_D, n_vacc_classes)
dim(aux_R_stepdown_D_conf) <- c(n_groups, s_stepdown_D, n_vacc_classes)
dim(new_R_stepdown_D_conf) <- c(n_groups, s_stepdown_D, n_vacc_classes)
dim(n_R_stepdown_D_conf) <- c(n_groups, s_stepdown_D, n_vacc_classes)
dim(n_R_stepdown_D_unconf_to_conf) <- c(n_groups, s_stepdown_D, n_vacc_classes)

## Vectors handling the R_pos class
dim(R) <- c(n_groups, n_vacc_classes)
dim(new_R) <- c(n_groups, n_vacc_classes)

## Vectors handling the R_pre class and seroconversion
dim(R_pre) <- c(n_groups, 2, n_vacc_classes)
dim(new_R_pre) <- c(n_groups, 2, n_vacc_classes)
dim(n_R_pre) <- c(n_groups, 2, n_vacc_classes)
dim(gamma_R_pre) <- 2
dim(p_R_pre) <- c(n_groups, 2, n_vacc_classes)
dim(p_seroconversion) <- n_groups

## Vectors handling the R_pos class
dim(R_pos) <- c(n_groups, s_R_pos, n_vacc_classes)
dim(n_R_pos) <- c(n_groups, s_R_pos, n_vacc_classes)
dim(new_R_pos) <- c(n_groups, s_R_pos, n_vacc_classes)
dim(n_R_pre_to_R_pos) <- c(n_groups, n_vacc_classes)

## Vectors handling the R_neg class
dim(R_neg) <- c(n_groups, n_vacc_classes)
dim(new_R_neg) <- c(n_groups, n_vacc_classes)

## Vectors handling the D_hosp class
dim(D_hosp) <- n_groups
dim(new_D_hosp) <- n_groups

## Vectors handling the D_comm class
dim(D_comm) <- n_groups
dim(new_D_comm) <- n_groups

## Vectors handling the PCR classes
dim(PCR_pre) <- c(n_groups, s_PCR_pre, n_vacc_classes)
dim(n_PCR_pre) <- c(n_groups, s_PCR_pre, n_vacc_classes)
dim(new_PCR_pre) <- c(n_groups, s_PCR_pre, n_vacc_classes)
dim(PCR_pos) <- c(n_groups, s_PCR_pos, n_vacc_classes)
dim(n_PCR_pos) <- c(n_groups, s_PCR_pos, n_vacc_classes)
dim(new_PCR_pos) <- c(n_groups, s_PCR_pos, n_vacc_classes)
dim(PCR_neg) <- c(n_groups, n_vacc_classes)
dim(new_PCR_neg) <- c(n_groups, n_vacc_classes)

## Vectors handling the S->S transitions i.e. moving between vaccination classes
dim(p_S_next_vacc_class) <- c(n_groups, n_vacc_classes)
dim(n_S_next_vacc_class) <- c(n_groups, n_vacc_classes)

dim(p_E_next_vacc_class) <- c(n_groups, s_E, n_vacc_classes)
dim(n_E_next_vacc_class) <- c(n_groups, s_E, n_vacc_classes)
dim(n_E_progress) <- c(n_groups, s_E, n_vacc_classes)
dim(n_EE_next_vacc_class) <- c(n_groups, s_E, n_vacc_classes)

dim(p_I_asympt_next_vacc_class) <- c(n_groups, s_asympt, n_vacc_classes)
dim(n_I_asympt_next_vacc_class) <- c(n_groups, s_asympt, n_vacc_classes)
dim(n_I_asympt_progress) <- c(n_groups, s_asympt, n_vacc_classes)
dim(n_II_asympt_next_vacc_class) <- c(n_groups, s_asympt, n_vacc_classes)

dim(p_R_next_vacc_class) <- c(n_groups, n_vacc_classes)
dim(n_R_next_vacc_class) <- c(n_groups, n_vacc_classes)
dim(n_R_next_vacc_class_capped) <- c(n_groups, n_vacc_classes)
dim(n_R_next_vacc_class_tmp) <- c(n_groups, n_vacc_classes)
dim(n_R_progress) <- c(n_groups, n_vacc_classes)
dim(n_R_progress_capped) <- c(n_groups, n_vacc_classes)
dim(n_R_progress_tmp) <- c(n_groups, n_vacc_classes)
dim(n_RS_next_vacc_class) <- c(n_groups, n_vacc_classes)

## Vectors handling the S->E transition where infected are split
## between level of infectivity
dim(p_SE) <- c(n_groups, n_vacc_classes)
dim(n_SE) <- c(n_groups, n_vacc_classes)
dim(n_S_progress) <- c(n_groups, n_vacc_classes)
dim(n_SE_next_vacc_class) <- c(n_groups, n_vacc_classes)

## Vectors handling the E->I transition where newly infectious cases
## are split between level of severity
dim(n_EI_asympt) <- c(n_groups, n_vacc_classes)
dim(n_EI_sympt) <- c(n_groups, n_vacc_classes)
dim(n_EI_asympt_next_vacc_class) <- c(n_groups, n_vacc_classes)
dim(n_EI_sympt_next_vacc_class) <- c(n_groups, n_vacc_classes)

## Vectors handling I_sympt to R, I_comm_D transition
dim(n_sympt_to_comm_D) <- c(n_groups, n_vacc_classes)
dim(n_sympt_to_R) <- c(n_groups, n_vacc_classes)

## Vectors handling number of new hospitalisations, ICU admissions and
## recoveries in hospital
dim(n_sympt_to_hosp) <- c(n_groups, n_vacc_classes)
dim(n_sympt_to_triage) <- c(n_groups, n_vacc_classes)
dim(n_sympt_to_triage_conf) <- c(n_groups, n_vacc_classes)
dim(n_hosp_non_ICU) <- c(n_groups, n_vacc_classes)
dim(n_sympt_to_hosp_D) <- c(n_groups, n_vacc_classes)
dim(n_sympt_to_hosp_D_conf) <- c(n_groups, n_vacc_classes)
dim(n_sympt_to_hosp_R) <- c(n_groups, n_vacc_classes)
dim(n_sympt_to_hosp_R_conf) <- c(n_groups, n_vacc_classes)
dim(n_triage_unconf_to_ICU_D_unconf) <- c(n_groups, n_vacc_classes)
dim(n_triage_conf_to_ICU_D_conf) <- c(n_groups, n_vacc_classes)
dim(n_triage_unconf_to_ICU_S_R_unconf) <- c(n_groups, n_vacc_classes)
dim(n_triage_conf_to_ICU_S_R_conf) <- c(n_groups, n_vacc_classes)
dim(n_triage_unconf_to_ICU_S_D_unconf) <- c(n_groups, n_vacc_classes)
dim(n_triage_conf_to_ICU_S_D_conf) <- c(n_groups, n_vacc_classes)

## Vectors handling the serology flow
dim(n_com_to_R_pre) <- c(n_groups, 2, n_vacc_classes)

## Vectors handling the severity profile
dim(p_sympt) <- n_groups

## Vectors handling the potential death in hospital (general beds and ICU)
dim(prob_death_hosp_D) <- n_groups
dim(psi_death_hosp_D) <- n_groups
dim(prob_death_ICU) <- n_groups
dim(psi_death_ICU) <- n_groups
dim(prob_death_stepdown) <- n_groups
dim(psi_death_stepdown) <- n_groups

## Vector handling the probability of being admitted as confirmed
dim(prob_admit_conf) <- n_groups
dim(psi_admit_conf) <- n_groups

dim(cum_admit_by_age) <- n_groups

## Vectors handling the age specific heterogeneous transmission process
dim(lambda) <- n_groups
dim(s_ij) <- c(n_groups, n_groups)
dim(m) <- c(n_groups, n_groups)
dim(I_with_diff_trans) <- c(n_groups, n_vacc_classes)

## Vectors handling the loss of immunity
dim(n_RS) <- c(n_groups, n_vacc_classes)
dim(p_RS) <- n_groups

## Total population
initial(N_tot[]) <- 0
update(N_tot[]) <- sum(S[i, ]) + sum(R[i, ]) + D_hosp[i] + sum(E[i, , ]) +
  sum(I_asympt[i, , ]) + sum(I_sympt[i, , ]) +
  sum(I_triage_conf[i, , ]) + sum(I_triage_unconf[i, , ])  +
  sum(I_hosp_R_conf[i, , ]) + sum(I_hosp_R_unconf[i, , ]) +
  sum(I_hosp_D_conf[i, , ]) + sum(I_hosp_D_unconf[i, , ]) +
  sum(I_ICU_S_R_conf[i, , ]) + sum(I_ICU_S_R_unconf[i, , ]) +
  sum(I_ICU_S_D_conf[i, , ]) + sum(I_ICU_S_D_unconf[i, , ]) +
  sum(I_ICU_D_conf[i, , ]) + sum(I_ICU_D_unconf[i, , ]) +
  sum(R_stepdown_R_conf[i, , ]) + sum(R_stepdown_R_unconf[i, , ]) +
  sum(R_stepdown_D_conf[i, , ]) + sum(R_stepdown_D_unconf[i, , ]) +
  sum(I_comm_D[i, , ]) + D_comm[i]
dim(N_tot) <- n_groups

## Total population calculated with seroconversion flow
initial(N_tot2) <- 0
update(N_tot2) <- sum(S) + sum(R_pre) + sum(R_pos) + sum(R_neg) + sum(E)

## Total population calculated with PCR flow
initial(N_tot3) <- 0
update(N_tot3) <- sum(S) + sum(PCR_pre) + sum(PCR_pos) + sum(PCR_neg)

## Aggregate our reporting statistics by summing across age (simple
## for everything except for seropositivity data, done last)
initial(I_ICU_tot) <- 0
new_I_ICU_tot <- sum(new_I_ICU_S_R_conf) + sum(new_I_ICU_S_D_conf) +
  sum(new_I_ICU_D_conf)
update(I_ICU_tot) <- new_I_ICU_tot

initial(general_tot) <- 0
new_general_tot <- sum(new_I_triage_conf) + sum(new_I_hosp_R_conf) +
  sum(new_I_hosp_D_conf) + sum(new_R_stepdown_R_conf) +
  sum(new_R_stepdown_D_conf)
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
update(sero_pos) <- sum(new_R_pos[4:13, , ])

initial(cum_sympt_cases) <- 0
update(cum_sympt_cases) <- cum_sympt_cases + sum(n_EI_sympt)

## only over 25s (exclude groups 1 to 5)
initial(cum_sympt_cases_over25) <- 0
update(cum_sympt_cases_over25) <- cum_sympt_cases_over25 +
  sum(n_EI_sympt[6:n_groups, ])

## For REACT we exclude the 0-4 (1) and CHR (19) groups
initial(react_pos) <- 0
update(react_pos) <- sum(new_PCR_pos[2:18, , ])

## Vaccination engine

## First, the number of candidates
vaccine_n_candidates[] <- S[i, 1] + sum(E[i, , 1]) + sum(I_asympt[i, , 1]) +
  R[i, 1]
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
