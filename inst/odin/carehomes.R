## E and R stage indexed by i, j, k with
## i for the age group
## j for the progression (not exponential latent and infectious period)
## k for the infectivity group (for I) or vacc. group (for S)

## Number of classes (age & vaccination)

## Number of "groups", being the age classes, Carehome workers and
## Carehome residents. This will be 19 in all but experimental uses.
n_age_groups <- user()
n_groups <- user()

## Definition of the time-step and output as "time"
dt <- user()
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
update(S[, 1]) <- S[i, 1] - n_S_next_vacc_class[i, 1] - n_infections[i, 1] +
  n_S_next_vacc_class[i, n_vacc_classes] + n_RS[i, 1] # age, vaccination status
update(S[, 2:n_vacc_classes]) <- S[i, j] - n_S_next_vacc_class[i, j] +
  n_S_next_vacc_class[i, j - 1] - n_infections[i, j] # age, vaccination status
update(E[, , ]) <- new_E[i, j, k]
update(I_asympt[, , ]) <- new_I_asympt[i, j, k]
update(I_mild[, , ]) <- new_I_mild[i, j, k]
update(I_ILI[, , ]) <- new_I_ILI[i, j, k]
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
update(R_neg[, ]) <- new_R_neg[i, j] - n_RS[i, j]
update(R[, ]) <- R[i, j] + delta_R[i, j] - n_RS[i, j]
update(D_hosp[]) <- new_D_hosp[i]
update(D_comm[]) <- new_D_comm[i]
update(PCR_pre[, , ]) <- new_PCR_pre[i, j, k]
update(PCR_pos[, , ]) <- new_PCR_pos[i, j, k]
update(PCR_neg[, ]) <- PCR_neg[i, j] + n_PCR_pos[i, s_PCR_pos, j] - n_RS[i, j]
update(cum_admit_conf) <-
  cum_admit_conf +
  sum(n_ILI_to_hosp_D_conf) +
  sum(n_ILI_to_hosp_R_conf) +
  sum(n_ILI_to_triage_conf)
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
update(cum_admit_by_age[]) <- cum_admit_by_age[i] + sum(n_ILI_to_hosp[i, ])

## Individual probabilities of transition:
p_S_next_vacc_class[, ] <- 1 - exp(-vaccine_progression_rate[i, j] * dt)
p_SE[, ] <- 1 - exp(-lambda[i] *
                      rel_susceptibility[i, j] * dt) # S to I age/vacc dependent
p_EE <- 1 - exp(-gamma_E * dt) # progression of latent period
p_II_asympt <- 1 - exp(-gamma_asympt * dt) # progression of infectious period
p_II_mild <- 1 - exp(-gamma_mild * dt)
p_II_ILI <- 1 - exp(-gamma_ILI * dt)
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

p_hosp_ILI <- if (step >= length(p_hosp_ILI_step))
  p_hosp_ILI_step[length(p_hosp_ILI_step)] else p_hosp_ILI_step[step + 1]
prob_hosp_ILI[] <- p_hosp_ILI * psi_hosp_ILI[i]

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

## new infections
n_infections[, ] <- rbinom(S[i, j], p_SE[i, j])
## of those some can also be vaccinated or progress through vaccination classes
## --> number transitioning from S[j] to E[j+1] (j vaccination class)
n_SE_S_next_vacc_class[, ] <-
  rbinom(n_infections[i, j], p_S_next_vacc_class[i, j])
## resulting transitions from S[j] to E[j]
## (j vaccination class)
n_SE[, 1:n_vacc_classes] <- n_infections[i, j] -
  n_SE_S_next_vacc_class[i, j]

## vaccine progression
n_S_next_vacc_class[, ] <- rbinom(S[i, j] - n_infections[i, j],
                                  p_S_next_vacc_class[i, j])

## other transitions
n_EE[, , ] <- rbinom(E[i, j, k], p_EE)
n_II_asympt[, , ] <- rbinom(I_asympt[i, j, k], p_II_asympt)
n_II_mild[, , ] <- rbinom(I_mild[i, j, k], p_II_mild)
n_II_ILI[, , ] <- rbinom(I_ILI[i, j, k], p_II_ILI)
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
n_RS_tmp[, ] <- rbinom(R[i, j], p_RS[i])
n_RS[, ] <- min(n_RS_tmp[i, j], R_neg[i, j], PCR_neg[i, j])

## Cumulative infections, summed over all age groups
initial(cum_infections) <- 0
update(cum_infections) <- cum_infections + sum(n_infections)

## Computes the number of asymptomatic
n_EI_asympt[, ] <- rbinom(n_EE[i, s_E, j], p_asympt[i])

## Computes the number of mild cases - p_sympt_ILI gives the
## proportion of febrile/ILI cases among the symptomatics
n_EI_mild[, ] <-
  rbinom(n_EE[i, s_E, j] - n_EI_asympt[i, j], 1 - p_sympt_ILI[i])

## Computes the number of ILI cases
n_EI_ILI[, ] <- n_EE[i, s_E, j] - n_EI_asympt[i, j] - n_EI_mild[i, j]

## Work out the S->E and E->E transitions
aux_EE[, 1, 1] <- n_SE[i, 1] + n_SE_S_next_vacc_class[i, n_vacc_classes]
aux_EE[, 1, 2:n_vacc_classes] <-
  n_SE[i, k] + n_SE_S_next_vacc_class[i, k - 1]
aux_EE[, 2:s_E, ] <- n_EE[i, j - 1, k]
aux_EE[, 1:s_E, ] <- aux_EE[i, j, k] - n_EE[i, j, k]
new_E[, , ] <- E[i, j, k] + aux_EE[i, j, k]

## Work out the I_asympt->I_asympt transitions
aux_II_asympt[, 1, ] <- n_EI_asympt[i, k]
aux_II_asympt[, 2:s_asympt, ] <- n_II_asympt[i, j - 1, k]
aux_II_asympt[, 1:s_asympt, ] <- aux_II_asympt[i, j, k] - n_II_asympt[i, j, k]
new_I_asympt[, , ] <- I_asympt[i, j, k] + aux_II_asympt[i, j, k]

## Work out the I_mild->I_mild transitions
aux_II_mild[, 1, ] <- n_EI_mild[i, k]
aux_II_mild[, 2:s_mild, ] <- n_II_mild[i, j - 1, k]
aux_II_mild[, 1:s_mild, ] <- aux_II_mild[i, j, k] - n_II_mild[i, j, k]
new_I_mild[, , ] <- I_mild[i, j, k] + aux_II_mild[i, j, k]

## Work out the I_ILI->I_ILI transitions
aux_II_ILI[, 1, ] <- n_EI_ILI[i, k]
aux_II_ILI[, 2:s_ILI, ] <- n_II_ILI[i, j - 1, k]
aux_II_ILI[, 1:s_ILI, ] <- aux_II_ILI[i, j, k] - n_II_ILI[i, j, k]
new_I_ILI[, , ] <- I_ILI[i, j, k] + aux_II_ILI[i, j, k]

## Work out the flow from I_ILI -> R, comm_D, hosp
n_ILI_to_R[, ] <- rbinom(n_II_ILI[i, s_ILI, j], 1 - prob_hosp_ILI[i])
n_ILI_to_comm_D[, ] <-
  rbinom(n_II_ILI[i, s_ILI, j] - n_ILI_to_R[i, j], prob_death_comm[i])
n_ILI_to_hosp[, ] <-
  n_II_ILI[i, s_ILI, j] - n_ILI_to_R[i, j] - n_ILI_to_comm_D[i, j]

## Work out the I_comm_D -> I_comm_D transitions
aux_II_comm_D[, 1, ] <- n_ILI_to_comm_D[i, k]
aux_II_comm_D[, 2:s_comm_D, ] <- n_II_comm_D[i, j - 1, k]
aux_II_comm_D[, 1:s_comm_D, ] <- aux_II_comm_D[i, j, k] - n_II_comm_D[i, j, k]
new_I_comm_D[, , ] <- I_comm_D[i, j, k] + aux_II_comm_D[i, j, k]

## Work out the split in hospitals between hosp_D, hosp_R and triage
n_ILI_to_triage[, ] <- rbinom(n_ILI_to_hosp[i, j], prob_ICU_hosp[i])
n_ILI_to_triage_conf[, ] <- rbinom(n_ILI_to_triage[i, j],
                                   prob_admit_conf[i])
n_hosp_non_ICU[, ] <- n_ILI_to_hosp[i, j] - n_ILI_to_triage[i, j]
n_ILI_to_hosp_D[, ] <- rbinom(n_hosp_non_ICU[i, j], prob_death_hosp_D[i])
n_ILI_to_hosp_D_conf[, ] <- rbinom(n_ILI_to_hosp_D[i, j], prob_admit_conf[i])
n_ILI_to_hosp_R[, ] <- n_hosp_non_ICU[i, j] - n_ILI_to_hosp_D[i, j]
n_ILI_to_hosp_R_conf[, ] <- rbinom(n_ILI_to_hosp_R[i, j], prob_admit_conf[i])

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
  new_I_triage_unconf[i, 1, k] + n_ILI_to_triage[i, k] -
  n_ILI_to_triage_conf[i, k]
new_I_triage_conf[, , ] <-
  aux_II_triage_conf[i, j, k] + n_I_triage_unconf_to_conf[i, j, k]
new_I_triage_conf[, 1, ] <-
  new_I_triage_conf[i, 1, k] + n_ILI_to_triage_conf[i, k]

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
  new_I_hosp_R_unconf[i, 1, k] + n_ILI_to_hosp_R[i, k] -
  n_ILI_to_hosp_R_conf[i, k]
new_I_hosp_R_conf[, , ] <-
  aux_II_hosp_R_conf[i, j, k] + n_I_hosp_R_unconf_to_conf[i, j, k]
new_I_hosp_R_conf[, 1, ] <-
  new_I_hosp_R_conf[i, 1, k] + n_ILI_to_hosp_R_conf[i, k]

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
  new_I_hosp_D_unconf[i, 1, k] + n_ILI_to_hosp_D[i, k] -
  n_ILI_to_hosp_D_conf[i, k]
new_I_hosp_D_conf[, , ] <-
  aux_II_hosp_D_conf[i, j, k] + n_I_hosp_D_unconf_to_conf[i, j, k]
new_I_hosp_D_conf[, 1, ] <-
  new_I_hosp_D_conf[i, 1, k] + n_ILI_to_hosp_D_conf[i, k]

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
delta_D_hosp[] <-
  sum(n_II_hosp_D_unconf[i, s_hosp_D, ]) +
  sum(n_II_hosp_D_conf[i, s_hosp_D, ]) +
  sum(n_II_ICU_D_unconf[i, s_ICU_D, ]) +
  sum(n_II_ICU_D_conf[i, s_ICU_D, ]) +
  sum(n_R_stepdown_D_unconf[i, s_stepdown_D, ]) +
  sum(n_R_stepdown_D_conf[i, s_stepdown_D, ] )
new_D_hosp[] <- D_hosp[i] + delta_D_hosp[i]

## Work out the number of deaths in the community
delta_D_comm[] <- sum(n_II_comm_D[i, s_comm_D, ])
new_D_comm[] <- D_comm[i] + delta_D_comm[i]

## Work out the number of people entering the seroconversion flow
n_com_to_R_pre[, 1, ] <- rbinom(n_EE[i, s_E, k], p_R_pre_1)
n_com_to_R_pre[, 2, ] <- n_EE[i, s_E, k] - n_com_to_R_pre[i, 1, k]
new_R_pre[, , ] <- R_pre[i, j, k] + n_com_to_R_pre[i, j, k] - n_R_pre[i, j, k]


## Split the seroconversion flow between people who are going to
## seroconvert and people who are not
n_R_pre_to_R_pos[, ] <- rbinom(sum(n_R_pre[i, , j]), p_seroconversion[i])

delta_R_pos[, 1, ] <- n_R_pre_to_R_pos[i, k]
delta_R_pos[, 2:s_R_pos, ] <- n_R_pos[i, j - 1, k]
delta_R_pos[, , ] <- delta_R_pos[i, j, k] - n_R_pos[i, j, k]
new_R_pos[, , ] <- R_pos[i, j, k] + delta_R_pos[i, j, k]

delta_R_neg[, ] <- sum(n_R_pre[i, , j]) - n_R_pre_to_R_pos[i, j] +
  n_R_pos[i, s_R_pos, j]
new_R_neg[, ] <- R_neg[i, j] + delta_R_neg[i, j]

## Work out the total number of recovery
delta_R[, ] <-
  n_II_asympt[i, s_asympt, j] +
  n_II_mild[i, s_mild, j] +
  n_ILI_to_R[i, j] +
  n_II_hosp_R_conf[i, s_hosp_R, j] +
  n_II_hosp_R_unconf[i, s_hosp_R, j] +
  n_R_stepdown_R_conf[i, s_stepdown_R, j] +
  n_R_stepdown_R_unconf[i, s_stepdown_R, j]

## Work out the PCR positivity
delta_PCR_pre[, 1, ] <- n_infections[i, k]
delta_PCR_pre[, 2:s_PCR_pre, ] <- n_PCR_pre[i, j - 1, k]
delta_PCR_pre[, , ] <- delta_PCR_pre[i, j, k] - n_PCR_pre[i, j, k]
new_PCR_pre[, , ] <- PCR_pre[i, j, k] + delta_PCR_pre[i, j, k]

delta_PCR_pos[, 1, ] <- n_PCR_pre[i, s_PCR_pre, k]
delta_PCR_pos[, 2:s_PCR_pos, ] <- n_PCR_pos[i, j - 1, k]
delta_PCR_pos[, , ] <- delta_PCR_pos[i, j, k] - n_PCR_pos[i, j, k]
new_PCR_pos[, , ] <- PCR_pos[i, j, k] + delta_PCR_pos[i, j, k]

## Compute the force of infection
I_with_diff_trans[, ] <-
  (sum(I_asympt[i, , j]) + sum(I_mild[i, , j]) + sum(I_ILI[i, , j]) +
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
initial(I_mild[, , ]) <- 0
initial(I_ILI[, , ]) <- 0
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

## Parameters of the S classes
rel_susceptibility[, ] <- user()
dim(rel_susceptibility) <- user() # use length as provided by the user
n_vacc_classes <- dim(rel_susceptibility, 2)

vaccine_progression_rate[, ] <- user()
dim(vaccine_progression_rate) <- user()

## Parameters of the E classes
s_E <- user()
gamma_E <- user(0.1)

## Probability of transitioning from the E to the asymptomatic,
## febrile classes, the rest goes in mild i.e with mild symptoms
p_asympt[] <- user()
p_sympt_ILI[] <- user()

## Parameters of the I_asympt classes
s_asympt <- user()
gamma_asympt <- user(0.1)

## Parameters of the I_mild classes
s_mild <- user()
gamma_mild <- user(0.1)

## Parameters of the I_ILI classes
s_ILI <- user()
gamma_ILI <- user(0.1)
dim(p_hosp_ILI_step) <- user()
p_hosp_ILI_step[] <- user()
psi_hosp_ILI[] <- user()

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

## Vectors handling the I_mild class
dim(I_mild) <- c(n_groups, s_mild, n_vacc_classes)
dim(aux_II_mild) <- c(n_groups, s_mild, n_vacc_classes)
dim(new_I_mild) <- c(n_groups, s_mild, n_vacc_classes)
dim(n_II_mild) <- c(n_groups, s_mild, n_vacc_classes)

## Vectors handling the I_ILI class
dim(I_ILI) <- c(n_groups, s_ILI, n_vacc_classes)
dim(aux_II_ILI) <- c(n_groups, s_ILI, n_vacc_classes)
dim(new_I_ILI) <- c(n_groups, s_ILI, n_vacc_classes)
dim(n_II_ILI) <- c(n_groups, s_ILI, n_vacc_classes)
dim(prob_hosp_ILI) <- n_groups
dim(psi_hosp_ILI) <- n_groups

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
dim(delta_R) <- c(n_groups, n_vacc_classes)

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
dim(delta_R_pos) <- c(n_groups, s_R_pos, n_vacc_classes)
dim(n_R_pre_to_R_pos) <- c(n_groups, n_vacc_classes)

## Vectors handling the R_neg class
dim(R_neg) <- c(n_groups, n_vacc_classes)
dim(new_R_neg) <- c(n_groups, n_vacc_classes)
dim(delta_R_neg) <- c(n_groups, n_vacc_classes)

## Vectors handling the D_hosp class
dim(D_hosp) <- n_groups
dim(new_D_hosp) <- n_groups
dim(delta_D_hosp) <- n_groups

## Vectors handling the D_comm class
dim(D_comm) <- n_groups
dim(new_D_comm) <- n_groups
dim(delta_D_comm) <- n_groups

## Vectors handling the PCR classes
dim(PCR_pre) <- c(n_groups, s_PCR_pre, n_vacc_classes)
dim(delta_PCR_pre) <- c(n_groups, s_PCR_pre, n_vacc_classes)
dim(n_PCR_pre) <- c(n_groups, s_PCR_pre, n_vacc_classes)
dim(new_PCR_pre) <- c(n_groups, s_PCR_pre, n_vacc_classes)
dim(PCR_pos) <- c(n_groups, s_PCR_pos, n_vacc_classes)
dim(delta_PCR_pos) <- c(n_groups, s_PCR_pos, n_vacc_classes)
dim(n_PCR_pos) <- c(n_groups, s_PCR_pos, n_vacc_classes)
dim(new_PCR_pos) <- c(n_groups, s_PCR_pos, n_vacc_classes)
dim(PCR_neg) <- c(n_groups, n_vacc_classes)

## Vectors handling the S->S transitions i.e. moving between vaccination classes
dim(p_S_next_vacc_class) <- c(n_groups, n_vacc_classes)
dim(n_S_next_vacc_class) <- c(n_groups, n_vacc_classes)

## Vectors handling the S->E transition where infected are split
## between level of infectivity
dim(p_SE) <- c(n_groups, n_vacc_classes)
dim(n_SE) <- c(n_groups, n_vacc_classes)
dim(n_infections) <- c(n_groups, n_vacc_classes)
dim(n_SE_S_next_vacc_class) <- c(n_groups, n_vacc_classes)

## Vectors handling the E->I transition where newly infectious cases
## are split between level of severity
dim(n_EI_asympt) <- c(n_groups, n_vacc_classes)
dim(n_EI_mild) <- c(n_groups, n_vacc_classes)
dim(n_EI_ILI) <- c(n_groups, n_vacc_classes)

## Vectors handling I_ILI to R, I_comm_D transition
dim(n_ILI_to_comm_D) <- c(n_groups, n_vacc_classes)
dim(n_ILI_to_R) <- c(n_groups, n_vacc_classes)

## Vectors handling number of new hospitalisations, ICU admissions and
## recoveries in hospital
dim(n_ILI_to_hosp) <- c(n_groups, n_vacc_classes)
dim(n_ILI_to_triage) <- c(n_groups, n_vacc_classes)
dim(n_ILI_to_triage_conf) <- c(n_groups, n_vacc_classes)
dim(n_hosp_non_ICU) <- c(n_groups, n_vacc_classes)
dim(n_ILI_to_hosp_D) <- c(n_groups, n_vacc_classes)
dim(n_ILI_to_hosp_D_conf) <- c(n_groups, n_vacc_classes)
dim(n_ILI_to_hosp_R) <- c(n_groups, n_vacc_classes)
dim(n_ILI_to_hosp_R_conf) <- c(n_groups, n_vacc_classes)
dim(n_triage_unconf_to_ICU_D_unconf) <- c(n_groups, n_vacc_classes)
dim(n_triage_conf_to_ICU_D_conf) <- c(n_groups, n_vacc_classes)
dim(n_triage_unconf_to_ICU_S_R_unconf) <- c(n_groups, n_vacc_classes)
dim(n_triage_conf_to_ICU_S_R_conf) <- c(n_groups, n_vacc_classes)
dim(n_triage_unconf_to_ICU_S_D_unconf) <- c(n_groups, n_vacc_classes)
dim(n_triage_conf_to_ICU_S_D_conf) <- c(n_groups, n_vacc_classes)

## Vectors handling the serology flow
dim(n_com_to_R_pre) <- c(n_groups, 2, n_vacc_classes)

## Vectors handling the severity profile
dim(p_asympt) <- n_groups
dim(p_sympt_ILI) <- n_groups

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
dim(n_RS_tmp) <- c(n_groups, n_vacc_classes)
dim(n_RS) <- c(n_groups, n_vacc_classes)
dim(p_RS) <- n_groups

## Total population
initial(N_tot[]) <- 0
update(N_tot[]) <- sum(S[i, ]) + sum(R[i, ]) + D_hosp[i] + sum(E[i, , ]) +
  sum(I_asympt[i, , ]) + sum(I_mild[i, , ]) + sum(I_ILI[i, , ]) +
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
update(cum_sympt_cases) <- cum_sympt_cases + sum(n_EI_mild) + sum(n_EI_ILI)

## only over 25s (exclude groups 1 to 5)
initial(cum_sympt_cases_over25) <- 0
update(cum_sympt_cases_over25) <- cum_sympt_cases_over25 +
  sum(n_EI_mild[6:n_groups, ]) + sum(n_EI_ILI[6:n_groups, ])

## For REACT we exclude the 0-4 (1) and CHR (19) groups
initial(react_pos) <- 0
update(react_pos) <- sum(new_PCR_pos[2:18, , ])
