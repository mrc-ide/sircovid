## E and R stage indexed by i, j, k, l with
## i for the age group
## j for the strain
## k for the progression (not exponential latent and infectious period)
## l for the vacc. group

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
## stage without progressing disease stages)
##
## Note that if an individual makes a vaccine skip move, then in the following
## they are counted as making as having moved through each of the intermediate
## classes they skip over, even though they do not spend any time in them.
## Note this should not affect model dynamics, the assumption is purely for
## bookkeeping purposes
n_S_vaccinated[, ] <- n_S_next_vacc_class[i, j] +
  (if (vacc_skipped[j] > 0) n_S_vacc_skip[i, vacc_skipped[j]]
   else 0)
dim(n_S_vaccinated) <- c(n_groups, n_vacc_classes)
n_E_vaccinated[, ] <- sum(n_E_next_vacc_class[i, , , j]) +
  (if (vacc_skipped[j] > 0) sum(n_E_vacc_skip[i, , , vacc_skipped[j]])
   else 0)
dim(n_E_vaccinated) <- c(n_groups, n_vacc_classes)
n_I_A_vaccinated[, ] <- sum(n_I_A_next_vacc_class[i, , , j]) +
  (if (vacc_skipped[j] > 0) sum(n_I_A_vacc_skip[i, , , vacc_skipped[j]])
   else 0)
dim(n_I_A_vaccinated) <- c(n_groups, n_vacc_classes)
n_I_P_vaccinated[, ] <- sum(n_I_P_next_vacc_class[i, , , j]) +
  (if (vacc_skipped[j] > 0) sum(n_I_P_vacc_skip[i, , , vacc_skipped[j]])
   else 0)
dim(n_I_P_vaccinated) <- c(n_groups, n_vacc_classes)
n_R_vaccinated[, ] <- sum(n_R_next_vacc_class[i, , j]) +
  (if (vacc_skipped[j] > 0) sum(n_R_vacc_skip[i, , vacc_skipped[j]])
   else 0)
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


## Output the number of individuals making vaccine skip moves
## Note that these are counted in cum_n_vaccinated etc above, but we output
## the vaccine skip moves separately so we can monitor these relative to non-
## skip moves
initial(cum_n_S_vacc_skip[, ]) <- 0
update(cum_n_S_vacc_skip[, ]) <- cum_n_S_vacc_skip[i, j] + n_S_vacc_skip[i, j]
dim(cum_n_S_vacc_skip) <- c(n_groups, n_vacc_classes)
initial(cum_n_E_vacc_skip[, ]) <- 0
update(cum_n_E_vacc_skip[, ]) <- cum_n_E_vacc_skip[i, j] +
  sum(n_E_vacc_skip[i, , , j])
dim(cum_n_E_vacc_skip) <- c(n_groups, n_vacc_classes)
initial(cum_n_I_A_vacc_skip[, ]) <- 0
update(cum_n_I_A_vacc_skip[, ]) <- cum_n_I_A_vacc_skip[i, j] +
  sum(n_I_A_vacc_skip[i, , , j])
dim(cum_n_I_A_vacc_skip) <- c(n_groups, n_vacc_classes)
initial(cum_n_I_P_vacc_skip[, ]) <- 0
update(cum_n_I_P_vacc_skip[, ]) <- cum_n_I_P_vacc_skip[i, j] +
  sum(n_I_P_vacc_skip[i, , , j])
dim(cum_n_I_P_vacc_skip) <- c(n_groups, n_vacc_classes)
initial(cum_n_R_vacc_skip[, ]) <- 0
update(cum_n_R_vacc_skip[, ]) <- cum_n_R_vacc_skip[i, j] +
  sum(n_R_vacc_skip[i, , j])
dim(cum_n_R_vacc_skip) <- c(n_groups, n_vacc_classes)

n_vacc_skip[, ] <-
  n_S_vacc_skip[i, j] +
  sum(n_E_vacc_skip[i, , , j]) +
  sum(n_I_A_vacc_skip[i, , , j]) +
  sum(n_I_P_vacc_skip[i, , , j]) +
  sum(n_R_vacc_skip[i, , j])
dim(n_vacc_skip) <- c(n_groups, n_vacc_classes)

initial(cum_n_vacc_skip[, ]) <- 0
update(cum_n_vacc_skip[, ]) <- cum_n_vacc_skip[i, j] + n_vacc_skip[i, j]
dim(cum_n_vacc_skip) <- c(n_groups, n_vacc_classes)

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
update(T_sero_pre_1[, , , ]) <- new_T_sero_pre_1[i, j, k, l]
update(T_sero_pos_1[, , , ]) <- new_T_sero_pos_1[i, j, k, l]
update(T_sero_neg_1[, , ]) <- new_T_sero_neg_1[i, j, k]
update(T_sero_pre_2[, , , ]) <- new_T_sero_pre_2[i, j, k, l]
update(T_sero_pos_2[, , , ]) <- new_T_sero_pos_2[i, j, k, l]
update(T_sero_neg_2[, , ]) <- new_T_sero_neg_2[i, j, k]
update(R[, , ]) <- new_R[i, j, k]
update(D_hosp[]) <- D_hosp[i] + delta_D_hosp[i]
update(D_non_hosp[]) <- D_non_hosp[i] + delta_D_non_hosp[i]
update(T_PCR_pre[, , , ]) <- new_T_PCR_pre[i, j, k, l]
update(T_PCR_pos[, , , ]) <- new_T_PCR_pos[i, j, k, l]
update(T_PCR_neg[, , ]) <- new_T_PCR_neg[i, j, k]

# Confirmed admissions by age
delta_admit_conf <-
  sum(n_I_C_2_to_H_D_conf) +
  sum(n_I_C_2_to_H_R_conf) +
  sum(n_I_C_2_to_ICU_pre_conf)
update(cum_admit_conf) <- cum_admit_conf + delta_admit_conf


# New in-hospital diagnoses by age
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


delta_all_admission_0_9_conf <-
  sum(delta_diagnoses_admitted[1:2, ])

delta_all_admission_10_19_conf <-
  sum(delta_diagnoses_admitted[3:4, ])

delta_all_admission_20_29_conf <-
  sum(delta_diagnoses_admitted[5:6, ]) +
  sum(delta_diagnoses_admitted[18, ]) * 1 / 8

delta_all_admission_30_39_conf <-
  sum(delta_diagnoses_admitted[7:8, ]) +
  sum(delta_diagnoses_admitted[18, ]) * 2 / 8

delta_all_admission_40_49_conf <-
  sum(delta_diagnoses_admitted[9:10, ]) +
  sum(delta_diagnoses_admitted[18, ]) * 2 / 8

delta_all_admission_50_59_conf <-
  sum(delta_diagnoses_admitted[11:12, ]) +
  sum(delta_diagnoses_admitted[18, ]) * 2 / 8

delta_all_admission_60_69_conf <-
  sum(delta_diagnoses_admitted[13:14, ]) +
  sum(delta_diagnoses_admitted[18, ]) * 1 / 8 +
  sum(delta_diagnoses_admitted[19, ]) * 0.05

delta_all_admission_70_79_conf <-
  sum(delta_diagnoses_admitted[15:16, ]) +
  sum(delta_diagnoses_admitted[19, ]) * 0.2

delta_all_admission_80_plus_conf <-
  sum(delta_diagnoses_admitted[17, ]) +
  sum(delta_diagnoses_admitted[19, ]) * 0.75

initial(diagnoses_admitted[, ]) <- 0
update(diagnoses_admitted[, ]) <- diagnoses_admitted[i, j] +
  delta_diagnoses_admitted[i, j]
dim(diagnoses_admitted) <- c(n_groups, n_vacc_classes)

delta_diagnoses_admitted[, ] <-
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
dim(delta_diagnoses_admitted) <- c(n_groups, n_vacc_classes)

initial(cum_infections_disag[, ]) <- 0
update(cum_infections_disag[, ]) <- cum_infections_disag[i, j] +
  sum(n_S_progress[i, , j]) +
  sum(n_RE[i, , j])
dim(cum_infections_disag) <- c(n_groups, n_vacc_classes)

initial(admit_conf_inc) <- 0
update(admit_conf_inc) <- if (step %% steps_per_day == 0)
  delta_admit_conf else admit_conf_inc + delta_admit_conf

initial(new_conf_inc) <- 0
update(new_conf_inc) <- if (step %% steps_per_day == 0)
  delta_new_conf else new_conf_inc + delta_new_conf

# Admissions + confirmed by age

initial(all_admission_0_9_conf_inc) <- 0
update(all_admission_0_9_conf_inc) <- if (step %% steps_per_day == 0)
  delta_all_admission_0_9_conf else all_admission_0_9_conf_inc +
  delta_all_admission_0_9_conf

initial(all_admission_10_19_conf_inc) <- 0
update(all_admission_10_19_conf_inc) <- if (step %% steps_per_day == 0)
  delta_all_admission_10_19_conf else all_admission_10_19_conf_inc +
  delta_all_admission_10_19_conf

initial(all_admission_20_29_conf_inc) <- 0
update(all_admission_20_29_conf_inc) <- if (step %% steps_per_day == 0)
  delta_all_admission_20_29_conf else all_admission_20_29_conf_inc +
  delta_all_admission_20_29_conf

initial(all_admission_30_39_conf_inc) <- 0
update(all_admission_30_39_conf_inc) <- if (step %% steps_per_day == 0)
  delta_all_admission_30_39_conf else all_admission_30_39_conf_inc +
  delta_all_admission_30_39_conf

initial(all_admission_40_49_conf_inc) <- 0
update(all_admission_40_49_conf_inc) <- if (step %% steps_per_day == 0)
  delta_all_admission_40_49_conf else all_admission_40_49_conf_inc +
  delta_all_admission_40_49_conf

initial(all_admission_50_59_conf_inc) <- 0
update(all_admission_50_59_conf_inc) <- if (step %% steps_per_day == 0)
  delta_all_admission_50_59_conf else all_admission_50_59_conf_inc +
  delta_all_admission_50_59_conf

initial(all_admission_60_69_conf_inc) <- 0
update(all_admission_60_69_conf_inc) <- if (step %% steps_per_day == 0)
  delta_all_admission_60_69_conf else all_admission_60_69_conf_inc +
  delta_all_admission_60_69_conf

initial(all_admission_70_79_conf_inc) <- 0
update(all_admission_70_79_conf_inc) <- if (step %% steps_per_day == 0)
  delta_all_admission_70_79_conf else all_admission_70_79_conf_inc +
  delta_all_admission_70_79_conf

initial(all_admission_80_plus_conf_inc) <- 0
update(all_admission_80_plus_conf_inc) <- if (step %% steps_per_day == 0)
  delta_all_admission_80_plus_conf else all_admission_80_plus_conf_inc +
  delta_all_admission_80_plus_conf

update(cum_admit_by_age[]) <- cum_admit_by_age[i] + sum(n_I_C_2_to_hosp[i, , ])

## Individual probabilities of transition:

## vaccination
p_S_next_vacc_class[, ] <- vaccine_probability[i, j]
p_E_next_vacc_class[, , , ] <- vaccine_probability[i, l]
p_I_A_next_vacc_class[, , , ] <- vaccine_probability[i, l]
p_I_P_next_vacc_class[, , , ] <- vaccine_probability[i, l]
p_R_next_vacc_class[, , ] <- vaccine_probability[i, k]

p_S_vacc_skip[, ] <- vacc_skip_probability[i, j]
p_E_vacc_skip[, , , ] <- vacc_skip_probability[i, l]
p_I_A_vacc_skip[, , , ] <- vacc_skip_probability[i, l]
p_I_P_vacc_skip[, , , ] <- vacc_skip_probability[i, l]
p_R_vacc_skip[, , ] <- vacc_skip_probability[i, k]

## clinical progression
p_SE[, ] <- 1 - exp(- sum(lambda_susc[i, , j]) * dt) # S to I age/vacc dependent
p_E_progress[] <- 1 - exp(-gamma_E[i] * dt) # progression of latent period
p_I_A_progress[] <- 1 - exp(-gamma_A[i] * dt) # progression of infectious period
p_I_P_progress[] <- 1 - exp(-gamma_P[i] * dt)
p_I_C_1_progress[] <- 1 - exp(-gamma_C_1[i] * dt)
p_I_C_2_progress[] <- 1 - exp(-gamma_C_2[i] * dt)
p_G_D_progress[] <- 1 - exp(-gamma_G_D[i] * dt)
p_ICU_pre_progress[] <- 1 - exp(-gamma_ICU_pre[i] * dt)
p_H_R_progress[] <- 1 - exp(-gamma_H_R[i] * dt)
p_H_D_progress[] <- 1 - exp(-gamma_H_D[i] * dt)
p_ICU_W_R_progress[] <- 1 - exp(-gamma_ICU_W_R[i] * dt)
p_ICU_W_D_progress[] <- 1 - exp(-gamma_ICU_W_D[i] * dt)
p_ICU_D_progress[] <- 1 - exp(-gamma_ICU_D[i] * dt)
p_W_R_progress[] <- 1 - exp(-gamma_W_R[i] * dt)
p_W_D_progress[] <- 1 - exp(-gamma_W_D[i] * dt)
p_T_sero_pre_1_progress <- 1 - exp(-gamma_sero_pre_1 * dt)
p_T_sero_pos_1_progress <- 1 - exp(-gamma_sero_pos_1 * dt)
p_T_sero_pre_2_progress <- 1 - exp(-gamma_sero_pre_2 * dt)
p_T_sero_pos_2_progress <- 1 - exp(-gamma_sero_pos_2 * dt)
p_test <- 1 - exp(-gamma_U * dt)
p_T_PCR_pre_progress <- 1 - exp(-gamma_PCR_pre * dt)
p_T_PCR_pos_progress <- 1 - exp(-gamma_PCR_pos * dt)

dim(p_E_progress) <- n_strains
dim(p_I_A_progress) <- n_strains
dim(p_I_P_progress) <- n_strains
dim(p_I_C_1_progress) <- n_strains
dim(p_I_C_2_progress) <- n_strains
dim(p_G_D_progress) <- n_strains
dim(p_ICU_pre_progress) <- n_strains
dim(p_H_R_progress) <- n_strains
dim(p_H_D_progress) <- n_strains
dim(p_ICU_W_R_progress) <- n_strains
dim(p_ICU_W_D_progress) <- n_strains
dim(p_ICU_D_progress) <- n_strains
dim(p_W_R_progress) <- n_strains
dim(p_W_D_progress) <- n_strains


## Work out time-varying probabilities
p_C[, , ] <- if (as.integer(step) >= n_p_C_steps)
  min(p_C_step[n_p_C_steps, i] * rel_p_sympt[i, j, k] *
        strain_rel_p_sympt[j], as.numeric(1)) else
    min(p_C_step[step + 1, i] * rel_p_sympt[i, j, k] *
          strain_rel_p_sympt[j], as.numeric(1))

p_H[, , ] <- if (as.integer(step) >= n_p_H_steps)
  min(p_H_step[n_p_H_steps, i] * rel_p_hosp_if_sympt[i, j, k] *
        strain_rel_p_hosp_if_sympt[j],
      as.numeric(1)) else
        min(p_H_step[step + 1, i] * rel_p_hosp_if_sympt[i, j, k] *
              strain_rel_p_hosp_if_sympt[j],
            as.numeric(1))

p_ICU[, , ] <- if (as.integer(step) >= n_p_ICU_steps)
  min(p_ICU_step[n_p_ICU_steps, i] * rel_p_ICU[i, j, k] *
        strain_rel_p_icu[j], as.numeric(1)) else
    min(p_ICU_step[step + 1, i] * rel_p_ICU[i, j, k] *
          strain_rel_p_icu[j], as.numeric(1))

p_ICU_D[, , ] <- if (as.integer(step) >= n_p_ICU_D_steps)
  min(p_ICU_D_step[n_p_ICU_D_steps, i] * rel_p_ICU_D[i, j, k] *
        strain_rel_p_ICU_D[j],
      as.numeric(1)) else
        min(p_ICU_D_step[step + 1, i] * rel_p_ICU_D[i, j, k] *
              strain_rel_p_ICU_D[j], as.numeric(1))

p_H_D[, , ] <- if (as.integer(step) >= n_p_H_D_steps)
  min(p_H_D_step[n_p_H_D_steps, i] * rel_p_H_D[i, j, k] *
        strain_rel_p_H_D[j], as.numeric(1)) else
          min(p_H_D_step[step + 1, i] * rel_p_H_D[i, j, k] *
                strain_rel_p_H_D[j], as.numeric(1))

p_W_D[, , ] <- if (as.integer(step) >= n_p_W_D_steps)
  min(p_W_D_step[n_p_W_D_steps, i] * rel_p_W_D[i, j, k] *
        strain_rel_p_W_D[j], as.numeric(1)) else
          min(p_W_D_step[step + 1, i] * rel_p_W_D[i, j, k] *
                strain_rel_p_W_D[j], as.numeric(1))

p_G_D[, , ] <- if (as.integer(step) >= n_p_G_D_steps)
  min(p_G_D_step[n_p_G_D_steps, i] * rel_p_G_D[i, j, k] *
        strain_rel_p_G_D[j], as.numeric(1)) else
          min(p_G_D_step[step + 1, i] * rel_p_G_D[i, j, k] *
                strain_rel_p_G_D[j], as.numeric(1))

p_R[, , ] <- if (as.integer(step) >= n_p_R_steps)
  min(p_R_step[n_p_R_steps, i] * rel_p_R[i, j, k], as.numeric(1)) else
    min(p_R_step[step + 1, i] * rel_p_R[i, j, k], as.numeric(1))

p_star[] <- if (as.integer(step) >= n_p_star_steps)
  p_star_step[n_p_star_steps, i] else p_star_step[step + 1, i]

## Work out time-varying gammas
gamma_E[] <- if (as.integer(step) >= n_gamma_E_steps)
  gamma_E_step[n_gamma_E_steps] * rel_gamma_E[i] else
    gamma_E_step[step + 1] * rel_gamma_E[i]

gamma_A[] <- if (as.integer(step) >= n_gamma_A_steps)
  gamma_A_step[n_gamma_A_steps] * rel_gamma_A[i] else
    gamma_A_step[step + 1] * rel_gamma_A[i]

gamma_P[] <- if (as.integer(step) >= n_gamma_P_steps)
  gamma_P_step[n_gamma_P_steps] * rel_gamma_P[i] else
    gamma_P_step[step + 1] * rel_gamma_P[i]

gamma_C_1[] <- if (as.integer(step) >= n_gamma_C_1_steps)
  gamma_C_1_step[n_gamma_C_1_steps] * rel_gamma_C_1[i] else
    gamma_C_1_step[step + 1] * rel_gamma_C_1[i]

gamma_C_2[] <- if (as.integer(step) >= n_gamma_C_2_steps)
  gamma_C_2_step[n_gamma_C_2_steps] * rel_gamma_C_2[i] else
    gamma_C_2_step[step + 1] * rel_gamma_C_2[i]

gamma_G_D[] <- if (as.integer(step) >= n_gamma_G_D_steps)
  gamma_G_D_step[n_gamma_G_D_steps] * rel_gamma_G_D[i] else
    gamma_G_D_step[step + 1] * rel_gamma_G_D[i]

gamma_ICU_pre[] <- if (as.integer(step) >= n_gamma_ICU_pre_steps)
  gamma_ICU_pre_step[n_gamma_ICU_pre_steps] * rel_gamma_ICU_pre[i] else
    gamma_ICU_pre_step[step + 1] * rel_gamma_ICU_pre[i]

gamma_H_R[] <- if (as.integer(step) >= n_gamma_H_R_steps)
  gamma_H_R_step[n_gamma_H_R_steps] * rel_gamma_H_R[i] else
    gamma_H_R_step[step + 1] * rel_gamma_H_R[i]

gamma_H_D[] <- if (as.integer(step) >= n_gamma_H_D_steps)
  gamma_H_D_step[n_gamma_H_D_steps] * rel_gamma_H_D[i] else
    gamma_H_D_step[step + 1] * rel_gamma_H_D[i]

gamma_ICU_W_R[] <- if (as.integer(step) >= n_gamma_ICU_W_R_steps)
  gamma_ICU_W_R_step[n_gamma_ICU_W_R_steps] * rel_gamma_ICU_W_R[i] else
    gamma_ICU_W_R_step[step + 1] * rel_gamma_ICU_W_R[i]

gamma_ICU_W_D[] <- if (as.integer(step) >= n_gamma_ICU_W_D_steps)
  gamma_ICU_W_D_step[n_gamma_ICU_W_D_steps] * rel_gamma_ICU_W_D[i] else
    gamma_ICU_W_D_step[step + 1] * rel_gamma_ICU_W_D[i]

gamma_ICU_D[] <- if (as.integer(step) >= n_gamma_ICU_D_steps)
  gamma_ICU_D_step[n_gamma_ICU_D_steps] * rel_gamma_ICU_D[i] else
    gamma_ICU_D_step[step + 1] * rel_gamma_ICU_D[i]

gamma_W_R[] <- if (as.integer(step) >= n_gamma_W_R_steps)
  gamma_W_R_step[n_gamma_W_R_steps] * rel_gamma_W_R[i] else
    gamma_W_R_step[step + 1] * rel_gamma_W_R[i]

gamma_W_D[] <- if (as.integer(step) >= n_gamma_W_D_steps)
  gamma_W_D_step[n_gamma_W_D_steps] * rel_gamma_W_D[i] else
    gamma_W_D_step[step + 1] * rel_gamma_W_D[i]

gamma_U <- if (as.integer(step) >= n_gamma_U_steps)
  gamma_U_step[n_gamma_U_steps] else
    gamma_U_step[step + 1]

## Draws from binomial distributions for numbers changing between
## compartments:

## modelling infections and vaccine progression, which can happen simultaneously

#### flow out of S ####

## new infections

## Compute the new infections with multiple strains using nested binomials
## No one can move from S to E3 or E4
n_S_progress_tot[, ] <- rbinom(S[i, j], p_SE[i, j])
n_S_progress[, , ] <- if (j == 1 || n_real_strains == 1)
  rbinom(n_S_progress_tot[i, k], rel_foi_strain[i, j, k]) else
    (if (j == 2) n_S_progress_tot[i, k] - n_S_progress[i, 1, k] else 0)

## Seeding of first wave, we seed in group 4, strain 1, vaccine stratum 1
##
seed_step_end <- seed_step_start + length(seed_value)
seed_rate <- if (step >= seed_step_start && step < seed_step_end)
  seed_value[as.integer(step - seed_step_start + 1)] else 0
seed <- rpois(seed_rate)
seed_age_band <- as.integer(4) # 15-19y band

seed_step_start <- user()
seed_value[] <- user()
dim(seed_value) <- user()

n_S_progress[seed_age_band, 1, 1] <-
  n_S_progress[seed_age_band, 1, 1] +
  min(S[seed_age_band, 1] - n_S_progress_tot[seed_age_band, 1], seed)

## Introduction of new strains. n_S_progress is arranged as:
##
## [age, strain infected with, vaccine stage]
##
## As in the model initialisation we will use the teenager category,
## and only infect *unvaccinated* people. For now we will model only
## movement into the second compartment as that represents our "new"
## strain.
strain_seed_step_end <- strain_seed_step_start + length(strain_seed_value)
strain_seed_rate <-
  if (step >= strain_seed_step_start && step < strain_seed_step_end)
    strain_seed_value[as.integer(step - strain_seed_step_start + 1)] else 0
strain_seed <- rpois(strain_seed_rate)

strain_seed_step_start <- user()
strain_seed_value[] <- user()
dim(strain_seed_value) <- user()

## We must never try to move more individuals from this S category
## than are available, so need to do this with a min()
##
## NOTE: We *must* use the range 2:n_strains here even though only one
## strain variant is allowed exist, otherwise the generated code leads
## us to write out-of-bounds when running with a single strain.
##
## After setting up progress remove all transitions from S to
## strain 3 (1.2) and 4 (2.1) (this is a safety it's handled above).
n_S_progress[4, 2:n_strains, 1] <-
  if (j < 3) min(n_S_progress[i, j, k] + strain_seed,
                 n_S_progress[i, j, k] + S[i, k] -
                   sum(n_S_progress[i, , k])) else 0

## vaccine progression
n_S_next_vacc_class[, ] <- rbinom(S[i, j] - sum(n_S_progress[i, , j]),
                                  p_S_next_vacc_class[i, j])

n_S_vacc_skip[, ] <-
  rbinom(S[i, j] - sum(n_S_progress[i, , j]) -
           n_S_next_vacc_class[i, j],
         p_S_vacc_skip[i, j])

#### flow out of E ####

n_E_progress[, , , ] <- rbinom(E[i, j, k, l], p_E_progress[j])

## vaccine progression
n_E_next_vacc_class[, , , ] <- rbinom(E[i, j, k, l] - n_E_progress[i, j, k, l],
                                      p_E_next_vacc_class[i, j, k, l])

n_E_vacc_skip[, , , ] <-
  rbinom(E[i, j, k, l] - n_E_progress[i, j, k, l] -
           n_E_next_vacc_class[i, j, k, l],
         p_E_vacc_skip[i, j, k, l])

#### flow out of I_A ####

n_I_A_progress[, , , ] <- rbinom(I_A[i, j, k, l], p_I_A_progress[j])

## vaccine progression
n_I_A_next_vacc_class[, , , ] <- rbinom(
  I_A[i, j, k, l] - n_I_A_progress[i, j, k, l],
  p_I_A_next_vacc_class[i, j, k, l])


n_I_A_vacc_skip[, , , ] <-
  rbinom(I_A[i, j, k, l] -
           n_I_A_progress[i, j, k, l] -
           n_I_A_next_vacc_class[i, j, k, l],
         p_I_A_vacc_skip[i, j, k, l])

#### flow out of I_P ####

n_I_P_progress[, , , ] <- rbinom(I_P[i, j, k, l], p_I_P_progress[j])

## vaccine progression
n_I_P_next_vacc_class[, , , ] <- rbinom(
  I_P[i, j, k, l] - n_I_P_progress[i, j, k, l],
  p_I_P_next_vacc_class[i, j, k, l]
)

n_I_P_vacc_skip[, , , ] <-
  rbinom(I_P[i, j, k, l] -
           n_I_P_progress[i, j, k, l] -
           n_I_P_next_vacc_class[i, j, k, l],
         p_I_P_vacc_skip[i, j, k, l])

#### flow out of R ####

## rate of progressing from R w/o superinfection is just waning_rate
## if n_strains is 1 then can only go to S w.r. waning_rate
## if in R3 or R4 then can only go to S w.r. waning_rate
## otherwise:
##  R1 and R2 can progress to S (w.r. waning_rate[i]) or
##  R1 can progress to E3 (w.r. strain 2 (3 - 1))
##  R2 can progress tp E4 (w.r. strain 1 (3 - 2))
##
## Note that (if n_real_strains == 2)
## cross_immunity[1] is the cross immunity of strain 1 against strain 2
## cross_immunity[2] is the cross immunity of strain 2 against strain 1
rate_R_progress[, , ] <- waning_rate[i] +
  if (n_strains == 1 || j > 2) 0 else
    lambda_susc[i, 3 - j, k] * (1 - cross_immunity[j])

p_R_progress[, , ] <- 1 - exp(-rate_R_progress[i, j, k] * dt)

## n_R_progress is total number who either:
##  - leave R for S or E and stay in same vacc class
##  - leave R for S or E and change same vacc class
n_R_progress[, , ] <- rbinom(R[i, j, k], p_R_progress[i, j, k])

## Number going from R to S
## In one-strain model, all progress to S only
## In multi-strain model, R3 and R4 progress to S only
## If not modelling super infection then all progress to S only
##
## In multi-strain model, number of R1 and R2 to S is binomial w.p. waning over
##  waning plus prob strain
## TODO (RS): waning_rate should eventually be variant varying
p_RS[, , ] <- if (n_strains == 1 || j > 2) 1 else
  (if (waning_rate[i] == 0) 0 else
    (waning_rate[i] /
       (waning_rate[i] + lambda_susc[i, 3 - j, k] * (1 - cross_immunity[j]))))
n_RS[, , ] <- rbinom(n_R_progress[i, j, k], p_RS[i, j, k])


## n_RE[i, j, k] is the number going from
## R[age i, strain j, vacc class k] to
## E[age i, strain j + 2, vacc class k]
## Movement possible to E3 and E4 are from R1 and R2 respectively
n_RE[, , ] <- n_R_progress[i, j, k] - n_RS[i, j, k]

## R -> R vaccine progression
n_R_tmp[, , ] <- R[i, j, k] - n_R_progress[i, j, k]
n_R_next_vacc_class[, , ] <- rbinom(n_R_tmp[i, j, k],
                                    p_R_next_vacc_class[i, j, k])

n_R_vacc_skip[, , ] <-
  rbinom(n_R_tmp[i, j, k] -
           n_R_next_vacc_class[i, j, k],
         p_R_vacc_skip[i, j, k])

#### other transitions ####

n_I_C_1_progress[, , , ] <- rbinom(I_C_1[i, j, k, l], p_I_C_1_progress[j])
n_I_C_2_progress[, , , ] <- rbinom(I_C_2[i, j, k, l], p_I_C_2_progress[j])
n_G_D_progress[, , , ] <- rbinom(G_D[i, j, k, l], p_G_D_progress[j])
n_ICU_pre_unconf_progress[, , , ] <-
  rbinom(ICU_pre_unconf[i, j, k, l], p_ICU_pre_progress[j])
n_ICU_pre_conf_progress[, , , ] <-
  rbinom(ICU_pre_conf[i, j, k, l], p_ICU_pre_progress[j])
n_H_R_unconf_progress[, , , ] <-
  rbinom(H_R_unconf[i, j, k, l], p_H_R_progress[j])
n_H_R_conf_progress[, , , ] <- rbinom(H_R_conf[i, j, k, l], p_H_R_progress[j])
n_H_D_unconf_progress[, , , ] <-
  rbinom(H_D_unconf[i, j, k, l], p_H_D_progress[j])
n_H_D_conf_progress[, , , ] <- rbinom(H_D_conf[i, j, k, l], p_H_D_progress[j])
n_ICU_W_R_unconf_progress[, , , ] <-
  rbinom(ICU_W_R_unconf[i, j, k, l], p_ICU_W_R_progress[j])
n_ICU_W_R_conf_progress[, , , ] <-
  rbinom(ICU_W_R_conf[i, j, k, l], p_ICU_W_R_progress[j])
n_ICU_W_D_unconf_progress[, , , ] <-
  rbinom(ICU_W_D_unconf[i, j, k, l], p_ICU_W_D_progress[j])
n_ICU_W_D_conf_progress[, , , ] <-
  rbinom(ICU_W_D_conf[i, j, k, l], p_ICU_W_D_progress[j])
n_ICU_D_unconf_progress[, , , ] <-
  rbinom(ICU_D_unconf[i, j, k, l], p_ICU_D_progress[j])
n_ICU_D_conf_progress[, , , ] <-
  rbinom(ICU_D_conf[i, j, k, l], p_ICU_D_progress[j])
n_W_R_unconf_progress[, , , ] <-
  rbinom(W_R_unconf[i, j, k, l], p_W_R_progress[j])
n_W_R_conf_progress[, , , ] <- rbinom(W_R_conf[i, j, k, l], p_W_R_progress[j])
n_W_D_unconf_progress[, , , ] <-
  rbinom(W_D_unconf[i, j, k, l], p_W_D_progress[j])
n_W_D_conf_progress[, , , ] <- rbinom(W_D_conf[i, j, k, l], p_W_D_progress[j])
n_T_sero_pre_1_progress[, , , ] <-
  rbinom(T_sero_pre_1[i, j, k, l], p_T_sero_pre_1_progress)
n_T_sero_pos_1_progress[, , , ] <-
  rbinom(T_sero_pos_1[i, j, k, l], p_T_sero_pos_1_progress)
n_T_sero_pre_2_progress[, , , ] <-
  rbinom(T_sero_pre_2[i, j, k, l], p_T_sero_pre_2_progress)
n_T_sero_pos_2_progress[, , , ] <-
  rbinom(T_sero_pos_2[i, j, k, l], p_T_sero_pos_2_progress)
n_T_PCR_pre_progress[, , , ] <-
  rbinom(T_PCR_pre[i, j, k, l], p_T_PCR_pre_progress)
n_T_PCR_pos_progress[, , , ] <-
  rbinom(T_PCR_pos[i, j, k, l], p_T_PCR_pos_progress)

## Cumulative infections, summed over all age groups
initial(cum_infections) <- 0
delta_infections <- sum(n_S_progress) + sum(n_RE)
update(cum_infections) <- cum_infections + delta_infections

initial(infections_inc) <- 0
update(infections_inc) <- if (step %% steps_per_day == 0)
  delta_infections else infections_inc + delta_infections

initial(cum_infections_per_strain[]) <- 0
delta_infections_per_strain[] <-
  sum(n_S_progress[, i, ]) +
  (if (i > 2)
    (sum(n_RE[, i - 2, ]))
   else
     0)
update(cum_infections_per_strain[]) <-
  cum_infections_per_strain[i] + delta_infections_per_strain[i]
dim(delta_infections_per_strain) <- n_strains
dim(cum_infections_per_strain) <- n_strains

initial(infections_inc_per_strain[]) <- 0
update(infections_inc_per_strain[]) <-
  if (step %% steps_per_day == 0)
    delta_infections_per_strain[i] else
      infections_inc_per_strain[i] + delta_infections_per_strain[i]
dim(infections_inc_per_strain) <- n_strains

## Work out the new S (i for age, j for vaccination status)
new_S[, ] <- S[i, j] + sum(n_RS[i, , j]) + sum(n_infected_to_S[i, , j]) -
  sum(n_S_progress[i, , j]) - n_S_next_vacc_class[i, j]
new_S[, ] <- new_S[i, j] +
  (if (j == 1) n_S_next_vacc_class[i, n_vacc_classes] else
    n_S_next_vacc_class[i, j - 1]) -
  (if (vacc_skip_to[j] > 0) n_S_vacc_skip[i, j] else 0) +
  (if (vacc_skip_from[j] > 0) n_S_vacc_skip[i, vacc_skip_from[j]] else 0)

## Computes the number of asymptomatic
n_EI_A[, , ] <- rbinom(n_E_progress[i, j, k_E, k], 1 - p_C[i, j, k])

## Computes the number of symptomatic cases
n_EI_P[, , ] <- n_E_progress[i, j, k_E, k] - n_EI_A[i, j, k]

## Work out the S->E and E->E transitions
aux_E[, , , ] <- (if (k == 1) n_S_progress[i, j, l] +
                    (if (j > 2) n_RE[i, j - 2, l] else 0)
                  else n_E_progress[i, j, k - 1, l]) -
  n_E_progress[i, j, k, l] -
  n_E_next_vacc_class[i, j, k, l] +
  (if (l == 1)
    n_E_next_vacc_class[i, j, k, n_vacc_classes]
   else
     n_E_next_vacc_class[i, j, k, l - 1]) -
  (if (vacc_skip_to[l] > 0) n_E_vacc_skip[i, j, k, l] else 0) +
  (if (vacc_skip_from[l] > 0) n_E_vacc_skip[i, j, k, vacc_skip_from[l]] else 0)

new_E[, , , ] <- E[i, j, k, l] + aux_E[i, j, k, l]

## Work out the I_A->I_A transitions
aux_I_A[, , , ] <- (if (k == 1) n_EI_A[i, j, l] else
  n_I_A_progress[i, j, k - 1, l]) -
  n_I_A_progress[i, j, k, l] -
  n_I_A_next_vacc_class[i, j, k, l] +
  (if (l == 1) n_I_A_next_vacc_class[i, j, k, n_vacc_classes] else
    n_I_A_next_vacc_class[i, j, k, l - 1]) -
  (if (vacc_skip_to[l] > 0) n_I_A_vacc_skip[i, j, k, l] else 0) +
  (if (vacc_skip_from[l] > 0)
    n_I_A_vacc_skip[i, j, k, vacc_skip_from[l]] else 0)

new_I_A[, , , ] <- I_A[i, j, k, l] + aux_I_A[i, j, k, l]

## Work out the I_P->I_P transitions
aux_I_P[, , , ] <- (if (k == 1) n_EI_P[i, j, l] else
  n_I_P_progress[i, j, k - 1, l]) -
  n_I_P_progress[i, j, k, l] -
  n_I_P_next_vacc_class[i, j, k, l] +
  (if (l == 1) n_I_P_next_vacc_class[i, j, k, n_vacc_classes] else
    n_I_P_next_vacc_class[i, j, k, l - 1]) -
  (if (vacc_skip_to[l] > 0) n_I_P_vacc_skip[i, j, k, l] else 0) +
  (if (vacc_skip_from[l] > 0)
    n_I_P_vacc_skip[i, j, k, vacc_skip_from[l]] else 0)

new_I_P[, , , ] <- I_P[i, j, k, l] + aux_I_P[i, j, k, l]

## Work out the I_C_1->I_C_1 transitions
aux_I_C_1[, , , ] <- (if (k == 1) n_I_P_progress[i, j, k_P, l] else
  n_I_C_1_progress[i, j, k - 1, l]) - n_I_C_1_progress[i, j, k, l]

new_I_C_1[, , , ] <- I_C_1[i, j, k, l] + aux_I_C_1[i, j, k, l]

## Work out the I_C_2->I_C_2 transitions
aux_I_C_2[, , , ] <- (if (k == 1) n_I_C_1_progress[i, j, k_C_1, l] else
  n_I_C_2_progress[i, j, k - 1, l]) - n_I_C_2_progress[i, j, k, l]

new_I_C_2[, , , ] <- I_C_2[i, j, k, l] + aux_I_C_2[i, j, k, l]

## Work out the flow from I_C_2 -> R, G_D, hosp
n_I_C_2_to_RS[, , ] <-
  rbinom(n_I_C_2_progress[i, j, k_C_2, k], 1 - p_H[i, j, k])
n_I_C_2_to_G_D[, , ] <-
  rbinom(n_I_C_2_progress[i, j, k_C_2, k] - n_I_C_2_to_RS[i, j, k],
         p_G_D[i, j, k])
n_I_C_2_to_hosp[, , ] <- n_I_C_2_progress[i, j, k_C_2, k] -
  n_I_C_2_to_RS[i, j, k] - n_I_C_2_to_G_D[i, j, k]

## Work out the G_D -> G_D transitions
aux_G_D[, , , ] <- (if (k == 1) n_I_C_2_to_G_D[i, j, l] else
  n_G_D_progress[i, j, k - 1, l]) - n_G_D_progress[i, j, k, l]

new_G_D[, , , ] <- G_D[i, j, k, l] + aux_G_D[i, j, k, l]

## Work out the split in hospitals between H_D, H_R and ICU_pre
n_I_C_2_to_ICU_pre[, , ] <- rbinom(n_I_C_2_to_hosp[i, j, k], p_ICU[i, j, k])
n_I_C_2_to_ICU_pre_conf[, , ] <- rbinom(n_I_C_2_to_ICU_pre[i, j, k],
                                        p_star[i])
n_hosp_non_ICU[, , ] <- n_I_C_2_to_hosp[i, j, k] - n_I_C_2_to_ICU_pre[i, j, k]
n_I_C_2_to_H_D[, , ] <- rbinom(n_hosp_non_ICU[i, j, k], p_H_D[i, j, k])
n_I_C_2_to_H_D_conf[, , ] <- rbinom(n_I_C_2_to_H_D[i, j, k],
                                    p_star[i])
n_I_C_2_to_H_R[, , ] <- n_hosp_non_ICU[i, j, k] - n_I_C_2_to_H_D[i, j, k]
n_I_C_2_to_H_R_conf[, , ] <- rbinom(n_I_C_2_to_H_R[i, j, k],
                                    p_star[i])

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
  rbinom(n_ICU_pre_unconf_progress[i, j, k_ICU_pre, k], p_ICU_D[i, j, k])
n_ICU_pre_conf_to_ICU_D_conf[, , ] <-
  rbinom(n_ICU_pre_conf_progress[i, j, k_ICU_pre, k], p_ICU_D[i, j, k])
n_ICU_pre_unconf_to_ICU_W_D_unconf[, , ] <-
  rbinom(n_ICU_pre_unconf_progress[i, j, k_ICU_pre, k] -
           n_ICU_pre_unconf_to_ICU_D_unconf[i, j, k],
         p_W_D[i, j, k])
n_ICU_pre_unconf_to_ICU_W_R_unconf[, , ] <-
  n_ICU_pre_unconf_progress[i, j, k_ICU_pre, k] -
  n_ICU_pre_unconf_to_ICU_D_unconf[i, j, k] -
  n_ICU_pre_unconf_to_ICU_W_D_unconf[i, j, k]
n_ICU_pre_conf_to_ICU_W_D_conf[, , ] <-
  rbinom(n_ICU_pre_conf_progress[i, j, k_ICU_pre, k] -
           n_ICU_pre_conf_to_ICU_D_conf[i, j, k], p_W_D[i, j, k])
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

delta_D_hosp[] <- sum(delta_D_hosp_disag[i, ])

## Work out the number of deaths in the community
delta_D_non_hosp[] <- sum(delta_D_non_hosp_disag[i, ])

## Work out the number of people entering the seroconversion flow
n_com_to_T_sero_pre[, , ] <- n_E_progress[i, j, k_E, k]

## Calculate for sero flow 1
new_T_sero_pre_1[, , , ] <- T_sero_pre_1[i, j, k, l] -
  n_T_sero_pre_1_progress[i, j, k, l] +
  (if (k == 1) n_com_to_T_sero_pre[i, j, l] else
    n_T_sero_pre_1_progress[i, j, k - 1, l])


## Split the seroconversion flow between people who are going to
## seroconvert and people who are not
n_T_sero_pre_1_to_T_sero_pos_1[, , ] <-
  rbinom(n_T_sero_pre_1_progress[i, j, k_sero_pre_1, k], p_sero_pos_1[i])

new_T_sero_pos_1[, , , ] <- T_sero_pos_1[i, j, k, l] -
  n_T_sero_pos_1_progress[i, j, k, l] +
  (if (k == 1)  n_T_sero_pre_1_to_T_sero_pos_1[i, j, l] else
    n_T_sero_pos_1_progress[i, j, k - 1, l])

new_T_sero_neg_1[, , ] <- T_sero_neg_1[i, j, k] +
  n_T_sero_pre_1_progress[i, j, k_sero_pre_1, k] -
  n_T_sero_pre_1_to_T_sero_pos_1[i, j, k] +
  n_T_sero_pos_1_progress[i, j, k_sero_pos_1, k]


## Calculate for sero flow 2
new_T_sero_pre_2[, , , ] <- T_sero_pre_2[i, j, k, l] -
  n_T_sero_pre_2_progress[i, j, k, l] +
  (if (k == 1) n_com_to_T_sero_pre[i, j, l] else
    n_T_sero_pre_2_progress[i, j, k - 1, l])


## Split the seroconversion flow between people who are going to
## seroconvert and people who are not
n_T_sero_pre_2_to_T_sero_pos_2[, , ] <-
  rbinom(n_T_sero_pre_2_progress[i, j, k_sero_pre_2, k], p_sero_pos_2[i])

new_T_sero_pos_2[, , , ] <- T_sero_pos_2[i, j, k, l] -
  n_T_sero_pos_2_progress[i, j, k, l] +
  (if (k == 1)  n_T_sero_pre_2_to_T_sero_pos_2[i, j, l] else
    n_T_sero_pos_2_progress[i, j, k - 1, l])

new_T_sero_neg_2[, , ] <- T_sero_neg_2[i, j, k] +
  n_T_sero_pre_2_progress[i, j, k_sero_pre_2, k] -
  n_T_sero_pre_2_to_T_sero_pos_2[i, j, k] +
  n_T_sero_pos_2_progress[i, j, k_sero_pos_2, k]

n_infection_end[, , ] <- n_I_A_progress[i, j, k_A, k] +
  n_I_C_2_to_RS[i, j, k] +
  n_H_R_conf_progress[i, j, k_H_R, k] +
  n_H_R_unconf_progress[i, j, k_H_R, k] +
  n_W_R_conf_progress[i, j, k_W_R, k] +
  n_W_R_unconf_progress[i, j, k_W_R, k]

n_infected_to_R[, , ] <- rbinom(n_infection_end[i, j, k], p_R[i, j, k])

n_infected_to_S[, , ] <- n_infection_end[i, j, k] - n_infected_to_R[i, j, k]

## Work out the total number of recovery
new_R[, , ] <- R[i, j, k] -
  n_R_progress[i, j, k] -
  n_R_next_vacc_class[i, j, k] +
  n_infected_to_R[i, j, k] +
  (if (k == 1) n_R_next_vacc_class[i, j, n_vacc_classes] else
    n_R_next_vacc_class[i, j, k - 1])  -
  (if (vacc_skip_to[k] > 0) n_R_vacc_skip[i, j, k] else 0) +
  (if (vacc_skip_from[k] > 0) n_R_vacc_skip[i, j, vacc_skip_from[k]] else 0)

## Work out the PCR positivity
new_T_PCR_pre[, , , ] <- T_PCR_pre[i, j, k, l] -
  n_T_PCR_pre_progress[i, j, k, l] +
  (if (k == 1) n_S_progress[i, j, l] else
    n_T_PCR_pre_progress[i, j, k - 1, l])

new_T_PCR_pos[, , , ] <- T_PCR_pos[i, j, k, l] -
  n_T_PCR_pos_progress[i, j, k, l] +
  (if (k == 1) n_T_PCR_pre_progress[i, j, k_PCR_pre, l] +
     n_RE[i, j, l] else
       n_T_PCR_pos_progress[i, j, k - 1, l])

new_T_PCR_neg[, , ] <- T_PCR_neg[i, j, k] +
  n_T_PCR_pos_progress[i, j, k_PCR_pos, k]

## Compute the force of infection

I_with_diff_trans[, , ] <-
  rel_infectivity[i, j, k] * strain_transmission[j] * (
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
## P(Strain = 1) := P(Strain = Only 1) + P(Strain = 2->1), same for Strain = 2
lambda[, ] <- if (n_real_strains == 1) sum(s_ij[i, , 1]) else
  (if (j == 1) sum(s_ij[i, , 1]) + sum(s_ij[i, , 4]) else
    sum(s_ij[i, , 2]) + sum(s_ij[i, , 3]))
lambda_susc[, , ] <- lambda[i, j] * rel_susceptibility[i, j, k]

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
initial(T_sero_pre_1[, , , ]) <- 0
initial(T_sero_pos_1[, , , ]) <- 0
initial(T_sero_neg_1[, , ]) <- 0
initial(T_sero_pre_2[, , , ]) <- 0
initial(T_sero_pos_2[, , , ]) <- 0
initial(T_sero_neg_2[, , ]) <- 0
initial(R[, , ]) <- 0
initial(D_hosp[]) <- 0
initial(D_non_hosp[]) <- 0
initial(T_PCR_pre[, , , ]) <- 0
initial(T_PCR_pos[, , , ]) <- 0
initial(T_PCR_neg[, , ]) <- 0
initial(cum_admit_conf) <- 0
initial(cum_new_conf) <- 0
initial(cum_admit_by_age[]) <- 0

## User defined parameters - default in parentheses:

## Vaccination/strain effect parameters
n_vacc_classes <- user()
rel_susceptibility[, , ] <- user()
dim(rel_susceptibility) <- c(n_groups, n_strains, n_vacc_classes)
rel_p_sympt[, , ] <- user()
dim(rel_p_sympt) <- c(n_groups, n_strains, n_vacc_classes)
strain_rel_p_sympt[] <- user()
dim(strain_rel_p_sympt) <- n_strains
rel_p_hosp_if_sympt[, , ] <- user()
dim(rel_p_hosp_if_sympt) <- c(n_groups, n_strains, n_vacc_classes)
strain_rel_p_hosp_if_sympt[] <- user()
dim(strain_rel_p_hosp_if_sympt) <- n_strains
rel_p_ICU[, , ] <- user()
dim(rel_p_ICU) <- c(n_groups, n_strains, n_vacc_classes)
strain_rel_p_icu[] <- user()
dim(strain_rel_p_icu) <- n_strains
rel_p_ICU_D[, , ] <- user()
dim(rel_p_ICU_D) <- c(n_groups, n_strains, n_vacc_classes)
rel_p_H_D[, , ] <- user()
dim(rel_p_H_D) <- c(n_groups, n_strains, n_vacc_classes)
rel_p_W_D[, , ] <- user()
dim(rel_p_W_D) <- c(n_groups, n_strains, n_vacc_classes)
rel_p_G_D[, , ] <- user()
dim(rel_p_G_D) <- c(n_groups, n_strains, n_vacc_classes)
strain_rel_p_ICU_D[] <- user()
dim(strain_rel_p_ICU_D) <- n_strains
strain_rel_p_H_D[] <- user()
dim(strain_rel_p_H_D) <- n_strains
strain_rel_p_W_D[] <- user()
dim(strain_rel_p_W_D) <- n_strains
strain_rel_p_G_D[] <- user()
dim(strain_rel_p_G_D) <- n_strains
rel_p_R[, , ] <- user()
dim(rel_p_R) <- c(n_groups, n_strains, n_vacc_classes)
rel_infectivity[, , ] <- user()
dim(rel_infectivity) <- c(n_groups, n_strains, n_vacc_classes)

vaccine_progression_rate_base[, ] <- user()
dim(vaccine_progression_rate_base) <- c(n_groups, n_vacc_classes)

## Parameters of the E classes
k_E <- user()
dim(gamma_E) <- n_strains
gamma_E_step[] <- user()
n_gamma_E_steps <- user()
dim(gamma_E_step) <- n_gamma_E_steps
rel_gamma_E[] <- user()
dim(rel_gamma_E) <- n_strains

## Probability of transitioning from the E to the symptomatic class,
## the rest go into the asymptomatic class
p_C_step[, ] <- user()
n_p_C_steps <- user()
dim(p_C) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_C_step) <- c(n_p_C_steps, n_groups)

## Parameters of the I_A classes
k_A <- user()
dim(gamma_A) <- n_strains
gamma_A_step[] <- user()
n_gamma_A_steps <- user()
dim(gamma_A_step) <- n_gamma_A_steps
rel_gamma_A[] <- user()
dim(rel_gamma_A) <- n_strains

## Parameters of the I_P classes
k_P <- user()
dim(gamma_P) <- n_strains
gamma_P_step[] <- user()
n_gamma_P_steps <- user()
dim(gamma_P_step) <- n_gamma_P_steps
rel_gamma_P[] <- user()
dim(rel_gamma_P) <- n_strains

## Parameters of the I_C_1 classes
k_C_1 <- user()
dim(gamma_C_1) <- n_strains
gamma_C_1_step[] <- user()
n_gamma_C_1_steps <- user()
dim(gamma_C_1_step) <- n_gamma_C_1_steps
rel_gamma_C_1[] <- user()
dim(rel_gamma_C_1) <- n_strains

## Parameters of the I_C_2 classes
k_C_2 <- user()
dim(gamma_C_2) <- n_strains
gamma_C_2_step[] <- user()
n_gamma_C_2_steps <- user()
dim(gamma_C_2_step) <- n_gamma_C_2_steps
rel_gamma_C_2[] <- user()
dim(rel_gamma_C_2) <- n_strains

## Proportion of cases requiring hospitalisation
p_H_step[, ] <- user()
n_p_H_steps <- user()
dim(p_H) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_H_step) <- c(n_p_H_steps, n_groups)

## Parameters of the G_D class
k_G_D <- user()
dim(gamma_G_D) <- n_strains
gamma_G_D_step[] <- user()
n_gamma_G_D_steps <- user()
dim(gamma_G_D_step) <- n_gamma_G_D_steps
rel_gamma_G_D[] <- user()
dim(rel_gamma_G_D) <- n_strains
p_G_D_step[, ] <- user()
n_p_G_D_steps <- user()
dim(p_G_D) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_G_D_step) <- c(n_p_G_D_steps, n_groups)

## Parameters of the ICU_pre classes
k_ICU_pre <- user()
dim(gamma_ICU_pre) <- n_strains
gamma_ICU_pre_step[] <- user()
n_gamma_ICU_pre_steps <- user()
dim(gamma_ICU_pre_step) <- n_gamma_ICU_pre_steps
rel_gamma_ICU_pre[] <- user()
dim(rel_gamma_ICU_pre) <- n_strains

## Proportion of hospital cases progressing to ICU
p_ICU_step[, ] <- user()
n_p_ICU_steps <- user()
dim(p_ICU) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_ICU_step) <- c(n_p_ICU_steps, n_groups)

## Proportion of stepdown cases dying
p_W_D_step[, ] <- user()
n_p_W_D_steps <- user()
dim(p_W_D) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_W_D_step) <- c(n_p_W_D_steps, n_groups)

## Parameters of the H_R classes
k_H_R <- user()
dim(gamma_H_R) <- n_strains
gamma_H_R_step[] <- user()
n_gamma_H_R_steps <- user()
dim(gamma_H_R_step) <- n_gamma_H_R_steps
rel_gamma_H_R[] <- user()
dim(rel_gamma_H_R) <- n_strains

## Parameters of the H_D classes
k_H_D <- user()
dim(gamma_H_D) <- n_strains
gamma_H_D_step[] <- user()
n_gamma_H_D_steps <- user()
dim(gamma_H_D_step) <- n_gamma_H_D_steps
rel_gamma_H_D[] <- user()
dim(rel_gamma_H_D) <- n_strains
p_H_D_step[, ] <- user()
n_p_H_D_steps <- user()
dim(p_H_D) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_H_D_step) <- c(n_p_H_D_steps, n_groups)

## Parameters of the ICU_W_R classes
k_ICU_W_R <- user()
dim(gamma_ICU_W_R) <- n_strains
gamma_ICU_W_R_step[] <- user()
n_gamma_ICU_W_R_steps <- user()
dim(gamma_ICU_W_R_step) <- n_gamma_ICU_W_R_steps
rel_gamma_ICU_W_R[] <- user()
dim(rel_gamma_ICU_W_R) <- n_strains

## Parameters of the ICU_W_D classes
k_ICU_W_D <- user()
dim(gamma_ICU_W_D) <- n_strains
gamma_ICU_W_D_step[] <- user()
n_gamma_ICU_W_D_steps <- user()
dim(gamma_ICU_W_D_step) <- n_gamma_ICU_W_D_steps
rel_gamma_ICU_W_D[] <- user()
dim(rel_gamma_ICU_W_D) <- n_strains

## Parameters of the ICU_D classes
k_ICU_D <- user()
dim(gamma_ICU_D) <- n_strains
gamma_ICU_D_step[] <- user()
n_gamma_ICU_D_steps <- user()
dim(gamma_ICU_D_step) <- n_gamma_ICU_D_steps
rel_gamma_ICU_D[] <- user()
dim(rel_gamma_ICU_D) <- n_strains
p_ICU_D_step[, ] <- user()
n_p_ICU_D_steps <- user()
dim(p_ICU_D) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_ICU_D_step) <- c(n_p_ICU_D_steps, n_groups)

## Parameters of the W_R classes
k_W_R <- user()
dim(gamma_W_R) <- n_strains
gamma_W_R_step[] <- user()
n_gamma_W_R_steps <- user()
dim(gamma_W_R_step) <- n_gamma_W_R_steps
rel_gamma_W_R[] <- user()
dim(rel_gamma_W_R) <- n_strains

## Parameters of the W_D classes
k_W_D <- user()
dim(gamma_W_D) <- n_strains
gamma_W_D_step[] <- user()
n_gamma_W_D_steps <- user()
dim(gamma_W_D_step) <- n_gamma_W_D_steps
rel_gamma_W_D[] <- user()
dim(rel_gamma_W_D) <- n_strains

## Parameters of the T_sero_pre_1 classes
k_sero_pre_1 <- user()
gamma_sero_pre_1 <- user(0.1)
p_sero_pos_1[] <- user()

## Parameters of the T_sero_pos_1 classes
k_sero_pos_1 <- user()
gamma_sero_pos_1 <- user(0.1)

## Parameters of the T_sero_pre_2 classes
k_sero_pre_2 <- user()
gamma_sero_pre_2 <- user(0.1)
p_sero_pos_2[] <- user()

## Parameters of the T_sero_pos_2 classes
k_sero_pos_2 <- user()
gamma_sero_pos_2 <- user(0.1)

## Parameters relating to testing
gamma_U_step[] <- user()
n_gamma_U_steps <- user()
dim(gamma_U_step) <- n_gamma_U_steps
p_star_step[, ] <- user()
n_p_star_steps <- user()
dim(p_star) <- n_groups
dim(p_star_step) <- c(n_p_star_steps, n_groups)

## Parameters relating to PCR positivity
k_PCR_pre <- user()
gamma_PCR_pre <- user(0.1)
k_PCR_pos <- user()
gamma_PCR_pos <- user(0.1)

## Waning of immunity
waning_rate[] <- user()
dim(waning_rate) <- n_groups

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
n_real_strains <- if (n_strains == 4) 2 else 1

## Dimensions of the different "vectors" here vectors stand for
## multi-dimensional arrays

## Vectors handling the S class
dim(S) <- c(n_groups, n_vacc_classes)
dim(new_S) <- c(n_groups, n_vacc_classes)

## Vectors handling the E class
dim(E) <- c(n_groups, n_strains, k_E, n_vacc_classes)
dim(aux_E) <- c(n_groups, n_strains, k_E, n_vacc_classes)
dim(new_E) <- c(n_groups, n_strains, k_E, n_vacc_classes)

## Vectors handling the I_A class
dim(I_A) <- c(n_groups, n_strains, k_A, n_vacc_classes)
dim(aux_I_A) <- c(n_groups, n_strains, k_A, n_vacc_classes)
dim(new_I_A) <- c(n_groups, n_strains, k_A, n_vacc_classes)

## Vectors handling the I_P class
dim(I_P) <- c(n_groups, n_strains, k_P, n_vacc_classes)
dim(aux_I_P) <- c(n_groups, n_strains, k_P, n_vacc_classes)
dim(new_I_P) <- c(n_groups, n_strains, k_P, n_vacc_classes)

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


## Vectors handling the G_D class
dim(G_D) <- c(n_groups, n_strains, k_G_D, n_vacc_classes)
dim(aux_G_D) <- c(n_groups, n_strains, k_G_D, n_vacc_classes)
dim(new_G_D) <- c(n_groups, n_strains, k_G_D, n_vacc_classes)
dim(n_G_D_progress) <- c(n_groups, n_strains, k_G_D, n_vacc_classes)

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

## Vectors handling the R class
dim(R) <- c(n_groups, n_strains, n_vacc_classes)
dim(new_R) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the T_sero_pre_1 class and seroconversion
dim(T_sero_pre_1) <- c(n_groups, n_strains, k_sero_pre_1, n_vacc_classes)
dim(new_T_sero_pre_1) <- c(n_groups, n_strains, k_sero_pre_1, n_vacc_classes)
dim(n_T_sero_pre_1_progress) <-
  c(n_groups, n_strains, k_sero_pre_1, n_vacc_classes)
dim(p_sero_pos_1) <- n_groups

## Vectors handling the T_sero_pos_1 class
dim(T_sero_pos_1) <- c(n_groups, n_strains, k_sero_pos_1, n_vacc_classes)
dim(n_T_sero_pos_1_progress) <-
  c(n_groups, n_strains, k_sero_pos_1, n_vacc_classes)
dim(new_T_sero_pos_1) <- c(n_groups, n_strains, k_sero_pos_1, n_vacc_classes)
dim(n_T_sero_pre_1_to_T_sero_pos_1) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the T_sero_neg_1 class
dim(T_sero_neg_1) <- c(n_groups, n_strains, n_vacc_classes)
dim(new_T_sero_neg_1) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the T_sero_pre_2 class and seroconversion
dim(T_sero_pre_2) <- c(n_groups, n_strains, k_sero_pre_2, n_vacc_classes)
dim(new_T_sero_pre_2) <- c(n_groups, n_strains, k_sero_pre_2, n_vacc_classes)
dim(n_T_sero_pre_2_progress) <-
  c(n_groups, n_strains, k_sero_pre_2, n_vacc_classes)
dim(p_sero_pos_2) <- n_groups

## Vectors handling the T_sero_pos_2 class
dim(T_sero_pos_2) <- c(n_groups, n_strains, k_sero_pos_2, n_vacc_classes)
dim(n_T_sero_pos_2_progress) <-
  c(n_groups, n_strains, k_sero_pos_2, n_vacc_classes)
dim(new_T_sero_pos_2) <- c(n_groups, n_strains, k_sero_pos_2, n_vacc_classes)
dim(n_T_sero_pre_2_to_T_sero_pos_2) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the T_sero_neg_2 class
dim(T_sero_neg_2) <- c(n_groups, n_strains, n_vacc_classes)
dim(new_T_sero_neg_2) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling the D_hosp class
dim(D_hosp) <- n_groups
dim(delta_D_hosp) <- n_groups

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
dim(p_S_vacc_skip) <- c(n_groups, n_vacc_classes)
dim(n_S_vacc_skip) <- c(n_groups, n_vacc_classes)

dim(p_E_next_vacc_class) <- c(n_groups, n_strains, k_E, n_vacc_classes)
dim(n_E_next_vacc_class) <- c(n_groups, n_strains, k_E, n_vacc_classes)
dim(p_E_vacc_skip) <- c(n_groups, n_strains, k_E, n_vacc_classes)
dim(n_E_vacc_skip) <- c(n_groups, n_strains, k_E, n_vacc_classes)
dim(n_E_progress) <- c(n_groups, n_strains, k_E, n_vacc_classes)

dim(p_I_A_next_vacc_class) <-
  c(n_groups, n_strains, k_A, n_vacc_classes)
dim(n_I_A_next_vacc_class) <-
  c(n_groups, n_strains, k_A, n_vacc_classes)
dim(p_I_A_vacc_skip) <- c(n_groups, n_strains, k_A, n_vacc_classes)
dim(n_I_A_vacc_skip) <- c(n_groups, n_strains, k_A, n_vacc_classes)
dim(n_I_A_progress) <- c(n_groups, n_strains, k_A, n_vacc_classes)

dim(p_I_P_next_vacc_class) <-
  c(n_groups, n_strains, k_P, n_vacc_classes)
dim(n_I_P_next_vacc_class) <-
  c(n_groups, n_strains, k_P, n_vacc_classes)
dim(p_I_P_vacc_skip) <- c(n_groups, n_strains, k_P, n_vacc_classes)
dim(n_I_P_vacc_skip) <- c(n_groups, n_strains, k_P, n_vacc_classes)
dim(n_I_P_progress) <- c(n_groups, n_strains, k_P, n_vacc_classes)

## Vectors handling the S->E transition where infected are split
## between level of infectivity
dim(p_SE) <- c(n_groups, n_vacc_classes)
dim(n_S_progress) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_S_progress_tot) <- c(n_groups, n_vacc_classes)

## Vectors handling the E->I transition where newly infectious cases
## are split between level of severity
dim(n_EI_A) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_EI_P) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling I_C_2 to R, G_D transition
dim(n_I_C_2_to_G_D) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_I_C_2_to_RS) <- c(n_groups, n_strains, n_vacc_classes)

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

## Numbers transitioning from infected compartments to R or S
dim(n_infection_end) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_infected_to_R) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_infected_to_S) <- c(n_groups, n_strains, n_vacc_classes)
p_R_step[, ] <- user()
n_p_R_steps <- user()
dim(p_R) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_R_step) <- c(n_p_R_steps, n_groups)

## Vectors handling the serology flow
dim(n_com_to_T_sero_pre) <- c(n_groups, n_strains, n_vacc_classes)

dim(cum_admit_by_age) <- n_groups

## Vectors handling the age specific heterogeneous transmission process
dim(lambda) <- c(n_groups, n_real_strains)
dim(lambda_susc) <- c(n_groups, n_real_strains, n_vacc_classes)
dim(s_ij) <- c(n_groups, n_groups, n_strains)
dim(m) <- c(n_groups, n_groups)
dim(I_with_diff_trans) <- c(n_groups, n_strains, n_vacc_classes)

## Vectors handling progress from R
dim(p_R_next_vacc_class) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_R_next_vacc_class) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_R_vacc_skip) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_R_vacc_skip) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_R_progress) <- c(n_groups, n_strains, n_vacc_classes)

dim(n_R_tmp) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_RS) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_RS) <- c(n_groups, n_strains, n_vacc_classes)
dim(n_RE) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_R_progress) <- c(n_groups, n_strains, n_vacc_classes)
dim(rate_R_progress) <- c(n_groups, n_strains, n_vacc_classes)

dim(cross_immunity) <- n_real_strains
cross_immunity[] <- user()

## Total population
initial(N_tot[]) <- 0
update(N_tot[]) <- sum(S[i, ]) + sum(R[i, , ]) + D_hosp[i] + sum(E[i, , , ]) +
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
initial(N_tot_sero_1) <- 0
update(N_tot_sero_1) <- sum(S) + sum(T_sero_pre_1) +
  sum(T_sero_pos_1) + sum(T_sero_neg_1) + sum(E)

initial(N_tot_sero_2) <- 0
update(N_tot_sero_2) <- sum(S) + sum(T_sero_pre_2) +
  sum(T_sero_pos_2) + sum(T_sero_neg_2) + sum(E)

## Total population calculated with PCR flow
initial(N_tot_PCR) <- 0
update(N_tot_PCR) <- sum(S) + sum(T_PCR_pre) + sum(T_PCR_pos) + sum(T_PCR_neg)

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

## cumulative deaths in hospital by age
initial(D_hosp_tot) <- 0
delta_D_hosp_tot <- sum(delta_D_hosp)
update(D_hosp_tot) <- D_hosp_tot + delta_D_hosp_tot

initial(D_hosp_0_49_tot) <- 0
delta_D_hosp_0_49_tot <- sum(delta_D_hosp[1:10]) + delta_D_hosp[18] * 3 / 8
update(D_hosp_0_49_tot) <- D_hosp_0_49_tot + delta_D_hosp_0_49_tot

initial(D_hosp_50_54_tot) <- 0
delta_D_hosp_50_54_tot <- delta_D_hosp[11] + delta_D_hosp[18] * 1 / 8
update(D_hosp_50_54_tot) <- D_hosp_50_54_tot + delta_D_hosp_50_54_tot

initial(D_hosp_55_59_tot) <- 0
delta_D_hosp_55_59_tot <- delta_D_hosp[12] + delta_D_hosp[18] * 2 / 8
update(D_hosp_55_59_tot) <- D_hosp_55_59_tot + delta_D_hosp_55_59_tot

initial(D_hosp_60_64_tot) <- 0
delta_D_hosp_60_64_tot <- delta_D_hosp[13] + delta_D_hosp[18] * 2 / 8
update(D_hosp_60_64_tot) <- D_hosp_60_64_tot + delta_D_hosp_60_64_tot

initial(D_hosp_65_69_tot) <- 0
delta_D_hosp_65_69_tot <- delta_D_hosp[14] + delta_D_hosp[19] * 0.05
update(D_hosp_65_69_tot) <- D_hosp_65_69_tot + delta_D_hosp_65_69_tot

initial(D_hosp_70_74_tot) <- 0
delta_D_hosp_70_74_tot <- delta_D_hosp[15] + delta_D_hosp[19] * 0.05
update(D_hosp_70_74_tot) <- D_hosp_70_74_tot + delta_D_hosp_70_74_tot

initial(D_hosp_75_79_tot) <- 0
delta_D_hosp_75_79_tot <- delta_D_hosp[16] + delta_D_hosp[19] * 0.15
update(D_hosp_75_79_tot) <- D_hosp_75_79_tot + delta_D_hosp_75_79_tot

initial(D_hosp_80_plus_tot) <- 0
delta_D_hosp_80_plus_tot <- delta_D_hosp[17] + delta_D_hosp[19] * 0.75
update(D_hosp_80_plus_tot) <- D_hosp_80_plus_tot + delta_D_hosp_80_plus_tot

## community deaths are non-hospital deaths in groups 1 to 18
initial(D_comm_tot) <- 0
delta_D_comm_tot <- sum(delta_D_non_hosp[1:18])
update(D_comm_tot) <- D_comm_tot + delta_D_comm_tot

initial(D_comm_inc) <- 0
update(D_comm_inc) <- if (step %% steps_per_day == 0)
  delta_D_comm_tot else D_comm_inc + delta_D_comm_tot

initial(D_comm_0_49_inc) <- 0
delta_D_comm_0_49 <- sum(delta_D_non_hosp[1:10]) + delta_D_non_hosp[18] * 3 / 8
update(D_comm_0_49_inc) <- if (step %% steps_per_day == 0)
  delta_D_comm_0_49 else D_comm_0_49_inc + delta_D_comm_0_49

initial(D_comm_50_54_inc) <- 0
delta_D_comm_50_54 <- delta_D_non_hosp[11] + delta_D_non_hosp[18] * 1 / 8
update(D_comm_50_54_inc) <- if (step %% steps_per_day == 0)
  delta_D_comm_50_54 else D_comm_50_54_inc + delta_D_comm_50_54

initial(D_comm_55_59_inc) <- 0
delta_D_comm_55_59 <- delta_D_non_hosp[12] + delta_D_non_hosp[18] * 2 / 8
update(D_comm_55_59_inc) <- if (step %% steps_per_day == 0)
  delta_D_comm_55_59 else D_comm_55_59_inc + delta_D_comm_55_59

initial(D_comm_60_64_inc) <- 0
delta_D_comm_60_64 <- delta_D_non_hosp[13] + delta_D_non_hosp[18] * 2 / 8
update(D_comm_60_64_inc) <- if (step %% steps_per_day == 0)
  delta_D_comm_60_64 else D_comm_60_64_inc + delta_D_comm_60_64

initial(D_comm_65_69_inc) <- 0
delta_D_comm_65_69 <- delta_D_non_hosp[14]
update(D_comm_65_69_inc) <- if (step %% steps_per_day == 0)
  delta_D_comm_65_69 else D_comm_65_69_inc + delta_D_comm_65_69

initial(D_comm_70_74_inc) <- 0
delta_D_comm_70_74 <- delta_D_non_hosp[15]
update(D_comm_70_74_inc) <- if (step %% steps_per_day == 0)
  delta_D_comm_70_74 else D_comm_70_74_inc + delta_D_comm_70_74

initial(D_comm_75_79_inc) <- 0
delta_D_comm_75_79 <- delta_D_non_hosp[16]
update(D_comm_75_79_inc) <- if (step %% steps_per_day == 0)
  delta_D_comm_75_79 else D_comm_75_79_inc + delta_D_comm_75_79

initial(D_comm_80_plus_inc) <- 0
delta_D_comm_80_plus <- delta_D_non_hosp[17]
update(D_comm_80_plus_inc) <- if (step %% steps_per_day == 0)
  delta_D_comm_80_plus else D_comm_80_plus_inc + delta_D_comm_80_plus


## carehome deaths are non-hospital deaths in group 19
initial(D_carehomes_tot) <- 0
delta_D_carehomes_tot <- delta_D_non_hosp[19]
update(D_carehomes_tot) <- D_carehomes_tot + delta_D_carehomes_tot

initial(D_carehomes_inc) <- 0
update(D_carehomes_inc) <- if (step %% steps_per_day == 0)
  delta_D_carehomes_tot else D_carehomes_inc + delta_D_carehomes_tot

initial(D_tot) <- 0
delta_D_tot <- delta_D_hosp_tot + delta_D_comm_tot + delta_D_carehomes_tot
update(D_tot) <- D_tot + delta_D_tot

initial(D_inc) <- 0
update(D_inc) <- if (step %% steps_per_day == 0)
  delta_D_tot else D_inc + delta_D_tot

## Incident deaths in hospital overall and then by age
initial(D_hosp_inc) <- 0
update(D_hosp_inc) <- if (step %% steps_per_day == 0)
  delta_D_hosp_tot else D_hosp_inc + delta_D_hosp_tot

initial(D_hosp_0_49_inc) <- 0
update(D_hosp_0_49_inc) <- if (step %% steps_per_day == 0)
  delta_D_hosp_0_49_tot else D_hosp_0_49_inc + delta_D_hosp_0_49_tot

initial(D_hosp_50_54_inc) <- 0
update(D_hosp_50_54_inc) <- if (step %% steps_per_day == 0)
  delta_D_hosp_50_54_tot else D_hosp_50_54_inc + delta_D_hosp_50_54_tot

initial(D_hosp_55_59_inc) <- 0
update(D_hosp_55_59_inc) <- if (step %% steps_per_day == 0)
  delta_D_hosp_55_59_tot else D_hosp_55_59_inc + delta_D_hosp_55_59_tot

initial(D_hosp_60_64_inc) <- 0
update(D_hosp_60_64_inc) <- if (step %% steps_per_day == 0)
  delta_D_hosp_60_64_tot else D_hosp_60_64_inc + delta_D_hosp_60_64_tot

initial(D_hosp_65_69_inc) <- 0
update(D_hosp_65_69_inc) <- if (step %% steps_per_day == 0)
  delta_D_hosp_65_69_tot else D_hosp_65_69_inc + delta_D_hosp_65_69_tot

initial(D_hosp_70_74_inc) <- 0
update(D_hosp_70_74_inc) <- if (step %% steps_per_day == 0)
  delta_D_hosp_70_74_tot else D_hosp_70_74_inc + delta_D_hosp_70_74_tot

initial(D_hosp_75_79_inc) <- 0
update(D_hosp_75_79_inc) <- if (step %% steps_per_day == 0)
  delta_D_hosp_75_79_tot else D_hosp_75_79_inc + delta_D_hosp_75_79_tot

initial(D_hosp_80_plus_inc) <- 0
update(D_hosp_80_plus_inc) <- if (step %% steps_per_day == 0)
  delta_D_hosp_80_plus_tot else D_hosp_80_plus_inc + delta_D_hosp_80_plus_tot


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
initial(sero_pos_1) <- 0
update(sero_pos_1) <- sum(new_T_sero_pos_1[4:13, , , ])

initial(sero_pos_2) <- 0
update(sero_pos_2) <- sum(new_T_sero_pos_2[4:13, , , ])

initial(cum_sympt_cases) <- 0
new_sympt_cases <- sum(n_EI_P)
update(cum_sympt_cases) <- cum_sympt_cases + new_sympt_cases

initial(cum_sympt_cases_non_variant) <- 0
new_sympt_cases_non_variant <-
  sum(n_EI_P[, 1, ]) +
  (if (n_real_strains == 2) sum(n_EI_P[, 4, ]) else 0)
update(cum_sympt_cases_non_variant) <-
  cum_sympt_cases_non_variant + new_sympt_cases_non_variant

## only over 25s (exclude groups 1 to 5)
initial(cum_sympt_cases_over25) <- 0
new_sympt_cases_over25 <- sum(n_EI_P[6:n_groups, , ])
update(cum_sympt_cases_over25) <- cum_sympt_cases_over25 +
  new_sympt_cases_over25

initial(cum_sympt_cases_non_variant_over25) <- 0
new_sympt_cases_non_variant_over25 <-
  sum(n_EI_P[6:n_groups, 1, ]) +
  (if (n_real_strains == 2) sum(n_EI_P[6:n_groups, 4, ]) else 0)
update(cum_sympt_cases_non_variant_over25) <-
  cum_sympt_cases_non_variant_over25 + new_sympt_cases_non_variant_over25

## new pillar 2 by age
initial(cum_sympt_cases_under15) <- 0
new_sympt_cases_under15 <- sum(n_EI_P[1:3, , ])
update(cum_sympt_cases_under15) <- cum_sympt_cases_under15 +
  new_sympt_cases_under15

initial(cum_sympt_cases_15_24) <- 0
new_sympt_cases_15_24 <- sum(n_EI_P[4:5, , ])
update(cum_sympt_cases_15_24) <- cum_sympt_cases_15_24 +
  new_sympt_cases_15_24

## assume CHW [18] are equally distributed amongst 25-64 age bands
initial(cum_sympt_cases_25_49) <- 0
new_sympt_cases_25_49 <- sum(n_EI_P[6:10, , ]) + (sum(n_EI_P[18, , ]) / 8) * 5
update(cum_sympt_cases_25_49) <- cum_sympt_cases_25_49 +
  new_sympt_cases_25_49

initial(cum_sympt_cases_50_64) <- 0
new_sympt_cases_50_64 <- sum(n_EI_P[11:13, , ]) + (sum(n_EI_P[18, , ]) / 8) * 3
update(cum_sympt_cases_50_64) <- cum_sympt_cases_50_64 +
  new_sympt_cases_50_64

## assume CHR [19] are 1/4 aged 65-79 and 3/4 80 plus
initial(cum_sympt_cases_65_79) <- 0
new_sympt_cases_65_79 <- sum(n_EI_P[14:16, , ]) + (sum(n_EI_P[19, , ]) * 0.25)
update(cum_sympt_cases_65_79) <- cum_sympt_cases_65_79 +
  new_sympt_cases_65_79

initial(cum_sympt_cases_80_plus) <- 0
new_sympt_cases_80_plus <- sum(n_EI_P[17, , ]) + (sum(n_EI_P[19, , ]) * 0.75)
update(cum_sympt_cases_80_plus) <- cum_sympt_cases_80_plus +
  new_sympt_cases_80_plus

## And incidence:
initial(sympt_cases_inc) <- 0
update(sympt_cases_inc) <- (
  if (step %% steps_per_day == 0) new_sympt_cases
  else sympt_cases_inc + new_sympt_cases)

initial(sympt_cases_non_variant_inc) <- 0
update(sympt_cases_non_variant_inc) <- (
  if (step %% steps_per_day == 0) new_sympt_cases_non_variant
  else sympt_cases_non_variant_inc + new_sympt_cases_non_variant)

initial(sympt_cases_over25_inc) <- 0
update(sympt_cases_over25_inc) <- (
  if (step %% steps_per_day == 0) new_sympt_cases_over25
  else sympt_cases_over25_inc + new_sympt_cases_over25)

initial(sympt_cases_non_variant_over25_inc) <- 0
update(sympt_cases_non_variant_over25_inc) <- (
  if (step %% steps_per_day == 0) new_sympt_cases_non_variant_over25
  else sympt_cases_non_variant_over25_inc + new_sympt_cases_non_variant_over25)

initial(sympt_cases_under15_inc) <- 0
update(sympt_cases_under15_inc) <- (
  if (step %% steps_per_day == 0) new_sympt_cases_under15
  else sympt_cases_under15_inc + new_sympt_cases_under15)

initial(sympt_cases_15_24_inc) <- 0
update(sympt_cases_15_24_inc) <- (
  if (step %% steps_per_day == 0) new_sympt_cases_15_24
  else sympt_cases_15_24_inc + new_sympt_cases_15_24)

initial(sympt_cases_25_49_inc) <- 0
update(sympt_cases_25_49_inc) <- (
  if (step %% steps_per_day == 0) new_sympt_cases_25_49
  else sympt_cases_25_49_inc + new_sympt_cases_25_49)

initial(sympt_cases_50_64_inc) <- 0
update(sympt_cases_50_64_inc) <- (
  if (step %% steps_per_day == 0) new_sympt_cases_50_64
  else sympt_cases_50_64_inc + new_sympt_cases_50_64)

initial(sympt_cases_65_79_inc) <- 0
update(sympt_cases_65_79_inc) <- (
  if (step %% steps_per_day == 0) new_sympt_cases_65_79
  else sympt_cases_65_79_inc + new_sympt_cases_65_79)

initial(sympt_cases_80_plus_inc) <- 0
update(sympt_cases_80_plus_inc) <- (
  if (step %% steps_per_day == 0) new_sympt_cases_80_plus
  else sympt_cases_80_plus_inc + new_sympt_cases_80_plus)

## For REACT we exclude the 0-4 (1) and CHR (19) groups
initial(react_pos) <- 0
update(react_pos) <- sum(new_T_PCR_pos[2:18, , , ])

initial(react_5_24_pos) <- 0
update(react_5_24_pos) <- sum(new_T_PCR_pos[2:5, , , ])

initial(react_25_34_pos) <- 0
update(react_25_34_pos) <- sum(new_T_PCR_pos[6:7, , , ]) +
  sum(new_T_PCR_pos[18, , , ]) * 2 / 8

initial(react_35_44_pos) <- 0
update(react_35_44_pos) <- sum(new_T_PCR_pos[8:9, , , ]) +
  sum(new_T_PCR_pos[18, , , ]) * 2 / 8

initial(react_45_54_pos) <- 0
update(react_45_54_pos) <- sum(new_T_PCR_pos[10:11, , , ]) +
  sum(new_T_PCR_pos[18, , , ]) * 2 / 8

initial(react_55_64_pos) <- 0
update(react_55_64_pos) <- sum(new_T_PCR_pos[12:13, , , ]) +
  sum(new_T_PCR_pos[18, , , ]) * 2 / 8

initial(react_65_plus_pos) <- 0
update(react_65_plus_pos) <- sum(new_T_PCR_pos[14:17, , , ])


## rel_foi_strain is probability of an infection in group i, vaccination class k
## being of strain j
##
## NOTE: the min(x / sum(x), 1) is required here because with floats,
## on a GPU, and with fast math, the sum can include sufficient
## rounding error that x / sum(x) can be > 1 by a very small amount;
## this keeps us bounded correctly.
rel_foi_strain[, , ] <-
  (if (sum(lambda_susc[i, , k]) == 0)
    (if (j == 1) 1 else 0) else
      min(lambda_susc[i, j, k] / sum(lambda_susc[i, , k]),
          as.numeric(1)))
dim(rel_foi_strain) <- c(n_groups, n_real_strains, n_vacc_classes)

## I_weighted used in IFR calculation
dim(new_I_weighted) <- c(n_groups, n_strains, n_vacc_classes)
new_I_weighted[, , ] <-
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
  G_D_transmission * sum(new_G_D[i, j, , k])
sum_new_I_weighted <- sum(new_I_weighted)
initial(I_weighted[, , ]) <- 0
dim(I_weighted) <- c(n_groups, n_strains, n_vacc_classes)
## If there are zero infectives we will just default to putting weight
## in group 4/strain 1/vaccine stratum 1. This will avoid NAs in IFR
update(I_weighted[, , ]) <-
  (if (sum_new_I_weighted == 0)
    (if (i == seed_age_band && j == 1 && k == 1) 1 else 0)
   else new_I_weighted[i, j, k])

## prob_strain is proportion of total I_weighted in each strain
## If there are zero infectives, we default to full weight on strain 1
## to avoid NAs in Rt
prob_strain_1 <- if (n_real_strains == 1 || sum_new_I_weighted == 0) 1 else
  (sum(new_I_weighted[, 1, ]) + sum(new_I_weighted[, 4, ])) /
  sum_new_I_weighted
initial(prob_strain[1:n_real_strains]) <- 0
initial(prob_strain[1]) <- 1
update(prob_strain[]) <- if (i == 1) prob_strain_1 else 1 - prob_strain_1
dim(prob_strain) <- n_real_strains

## Calculate effective susceptibles to each strain
## Weight each person in S/R by their relative susceptibility
## Note that for those in R we further account for cross immunity
## to strains. Those recovered from strain 3 - j will be (partially)
## susceptible to strain j
dim(eff_S) <- c(n_groups, n_real_strains, n_vacc_classes)
eff_S[, , ] <- if (n_real_strains == 1)
  S[i, k] * rel_susceptibility[i, j, k] else
    (S[i, k] + (1 - cross_immunity[3 - j]) * R[i, 3 - j, k]) *
  rel_susceptibility[i, j, k]

initial(effective_susceptible[]) <- 0
update(effective_susceptible[]) <- sum(eff_S[, i, ])
dim(effective_susceptible) <- n_real_strains

## Vaccination engine
n_doses <- user()
index_dose[] <- user(integer = TRUE)
dim(index_dose) <- n_doses

index_dose_inverse[] <- user(integer = TRUE)
dim(index_dose_inverse) <- n_vacc_classes

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

vacc_skip_n_candidates[, ] <-
  (if (vacc_skip_dose[j] > 0)
    S[i, vacc_skip_dose[j]] +
    sum(E[i, , , vacc_skip_dose[j]]) +
    sum(I_A[i, , , vacc_skip_dose[j]]) +
    sum(I_P[i, , , vacc_skip_dose[j]]) +
    sum(R[i, , vacc_skip_dose[j]])
   else 0)
dim(vacc_skip_n_candidates) <- c(n_groups, n_doses)

## Work out the vaccination probability via doses, driven by the
## schedule
vaccine_probability_doses[, ] <- min(
  if (vaccine_n_candidates[i, j] > 0)
    vaccine_attempted_doses[i, j] / vaccine_n_candidates[i, j] else 0,
  as.numeric(1))
dim(vaccine_probability_doses) <- c(n_groups, n_doses)

## Work out the total attempted doses
total_attempted_doses[, ] <- vaccine_missed_doses[i, j] + (
  if (as.integer(step) >= dim(vaccine_dose_step, 3)) 0
  else vaccine_dose_step[i, j, step + 1])
dim(total_attempted_doses) <- c(n_groups, n_doses)

## Now we work out the split of the total attempted doses, firstly for the
## next vaccine class moves, then for the vaccine skip moves.
##
## Note attempted doses for next vaccine class moves competing with a vaccine
## skip are weighted here, with remaining doses made available to the vaccine
## skip move below.
##
## Note that since we require vacc_skip_dose_weight <= 1, it is not possible for
## there to be an excess of doses for vaccine skip moves while having not
## enough doses for the next vaccine class candidates.
vaccine_attempted_doses[, ] <-
  (if (vaccine_n_candidates[i, j] == 0) 0
     else
       min(vaccine_n_candidates[i, j] /
             (vaccine_n_candidates[i, j] +
                vacc_skip_dose_weight[j] * vacc_skip_n_candidates[i, j])
           * total_attempted_doses[i, j],
           vaccine_n_candidates[i, j]))
dim(vaccine_attempted_doses) <- c(n_groups, n_doses)

vacc_skip_attempted_doses[, ] <-
  (if (vacc_skip_dose_weight[j] > 0)
    (if (vacc_skip_dose[j] > 0)
      total_attempted_doses[i, j] -
       vaccine_attempted_doses[i, j]
     else 0)
   else 0)
dim(vacc_skip_attempted_doses) <- c(n_groups, n_doses)

initial(vaccine_missed_doses[, ]) <- 0
update(vaccine_missed_doses[, ]) <-
  vaccine_catchup_fraction *
  max(total_attempted_doses[i, j] - n_vaccinated[i, index_dose[j]],
      as.numeric(0))
dim(vaccine_missed_doses) <- c(n_groups, n_doses)

vaccine_catchup_fraction <- user(0)


## Then either fix everything based on progression at a constant rate,
## or take from the supplied time-varying probabilities.
vaccine_probability[, ] <- (
  if (index_dose_inverse[j] > 0)
    vaccine_probability_doses[i, index_dose_inverse[j]]
  else
    1 - exp(-vaccine_progression_rate_base[i, j] * dt))
dim(vaccine_probability) <- c(n_groups, n_vacc_classes)

initial(tmp_vaccine_n_candidates[, ]) <- 0
update(tmp_vaccine_n_candidates[, ]) <- vaccine_n_candidates[i, j]
dim(tmp_vaccine_n_candidates) <- c(n_groups, n_doses)

initial(tmp_vaccine_probability[, ]) <- 0
update(tmp_vaccine_probability[, ]) <- vaccine_probability[i, j]
dim(tmp_vaccine_probability) <- c(n_groups, n_vacc_classes)

vacc_skip_probability[, ] <- (
  if (vacc_skip_dose_inverse[j] > 0) (
    if (vacc_skip_n_candidates[i, vacc_skip_dose_inverse[j]] > 0)
      min(vacc_skip_attempted_doses[i, vacc_skip_dose_inverse[j]] /
            vacc_skip_n_candidates[i, vacc_skip_dose_inverse[j]],
          as.numeric(1))
    else 0)
  else
    1 - exp(-vacc_skip_progression_rate_base[j] * dt))
dim(vacc_skip_probability) <- c(n_groups, n_vacc_classes)

## Vaccine skip inputs
## 1. vacc_skip_to[j] is the vaccine stratum that the vaccine skip move
##    from stratum j goes to (0 represents no vaccine skip move from j)
## 2. vacc_skip_from[j] is the vaccine stratum that the vaccine skip move
##    to stratum j comes from (0 represents no vaccine skip move to j)
## 3. vacc_skip_progression_rate_base[j] is the progression rate used for
##    the vaccine skip from stratum j (unless the vaccine skip move is
##    controlled by doses)
## 4. vacc_skip_dose[j] is the vaccine stratum that the vaccine skip move comes
##    moves from that is controlled by dose j (0 represents no vaccine skip
##    move controlled by dose j)
## 5. vacc_skip_dose_inverse[j] is the dose that the vaccine skip move from
##    stratum j is controlled by (0 represents that no dose controls the vaccine
##    skip move)
## 6. vacc_skip_dose_weight[j] represents how much vaccine skip candidates are
##    weighted for distribution of dose j relative to standard vaccine
##    candidates for that dose
## 7. vacc_skipped is used for bookkeeping - vacc_skipped[j] is the vaccine
##    stratum a vaccine skip move goes from that either starts at j or skips
##    over j (if there is no such move then the vacc_skipped[j] is 0)

vacc_skip_to[] <- user(integer = TRUE)
dim(vacc_skip_to) <- n_vacc_classes
vacc_skip_from[] <- user(integer = TRUE)
dim(vacc_skip_from) <- n_vacc_classes
vacc_skip_progression_rate_base[] <- user()
dim(vacc_skip_progression_rate_base) <- n_vacc_classes
vacc_skip_dose[] <- user(integer = TRUE)
dim(vacc_skip_dose) <- n_doses
vacc_skip_dose_inverse[] <- user(integer = TRUE)
dim(vacc_skip_dose_inverse) <- n_vacc_classes
vacc_skip_dose_weight[] <- user()
dim(vacc_skip_dose_weight) <- n_doses
vacc_skipped[] <- user(integer = TRUE)
dim(vacc_skipped) <- n_vacc_classes

## Severity outputs by age - vacc class - infection class
dim(new_inf) <- c(n_groups, n_strains, n_vacc_classes)
dim(IHR_disag) <- c(n_groups, n_strains, n_vacc_classes)
dim(IHR_disag_weighted) <- c(n_groups, n_strains, n_vacc_classes)
dim(HFR_disag) <- c(n_groups, n_strains, n_vacc_classes)
dim(HFR_disag_weighted) <- c(n_groups, n_strains, n_vacc_classes)
dim(IFR_disag) <- c(n_groups, n_strains, n_vacc_classes)
dim(IFR_disag_weighted) <- c(n_groups, n_strains, n_vacc_classes)

new_inf[, , ] <- n_S_progress[i, j, k] +
  (if (j > 2) n_RE[i, j - 2, k] else 0)

IHR_disag[, , ] <- p_C[i, j, k] * p_H[i, j, k] * (1 - p_G_D[i, j, k])
IHR_disag_weighted[, , ] <- IHR_disag[i, j, k] * new_inf[i, j, k]
initial(ihr) <- 0
update(ihr) <- sum(IHR_disag_weighted) / sum(new_inf)

HFR_disag[, , ] <- (1 - p_ICU[i, j, k]) * p_H_D[i, j, k] +
  p_ICU[i, j, k] * p_ICU_D[i, j, k] +
  p_ICU[i, j, k] * (1 - p_ICU_D[i, j, k]) * p_W_D[i, j, k]
HFR_disag_weighted[, , ] <- HFR_disag[i, j, k] * n_I_C_2_to_hosp[i, j, k]
initial(hfr) <- 0
update(hfr) <- sum(HFR_disag_weighted) / sum(n_I_C_2_to_hosp)

IFR_disag[, , ] <- IHR_disag[i, j, k] * HFR_disag[i, j, k] +
  p_C[i, j, k] * p_H[i, j, k] * p_G_D[i, j, k]
IFR_disag_weighted[, , ] <- IFR_disag[i, j, k] * new_inf[i, j, k]
initial(ifr) <- 0
update(ifr) <- sum(IFR_disag_weighted) / sum(new_inf)

## By strain
initial(ifr_strain[]) <- 0
update(ifr_strain[]) <- if (n_real_strains == 1)
  sum(IFR_disag_weighted[, 1, ]) / sum(new_inf[, 1, ]) else
    (sum(IFR_disag_weighted[, i, ]) + sum(IFR_disag_weighted[, 5 - i, ])) /
  (sum(new_inf[, i, ]) + sum(new_inf[, 5 - i, ]))
dim(ifr_strain) <- n_real_strains

initial(ihr_strain[]) <- 0
update(ihr_strain[]) <- if (n_real_strains == 1)
  sum(IHR_disag_weighted[, 1, ]) / sum(new_inf[, 1, ]) else
    (sum(IHR_disag_weighted[, i, ]) + sum(IHR_disag_weighted[, 5 - i, ])) /
  (sum(new_inf[, i, ]) + sum(new_inf[, 5 - i, ]))
dim(ihr_strain) <- n_real_strains

initial(hfr_strain[]) <- 0
update(hfr_strain[]) <- if (n_real_strains == 1)
  sum(HFR_disag_weighted[, 1, ]) / sum(n_I_C_2_to_hosp[, 1, ]) else
    (sum(HFR_disag_weighted[, i, ]) + sum(HFR_disag_weighted[, 5 - i, ])) /
  (sum(n_I_C_2_to_hosp[, i, ]) + sum(n_I_C_2_to_hosp[, 5 - i, ]))
dim(hfr_strain) <- n_real_strains


## By age and vaccination class
initial(ifr_age_vacc_class[, ]) <- 0
update(ifr_age_vacc_class[, ]) <- sum(IFR_disag_weighted[i, , j]) /
  sum(new_inf[i, , j])
dim(ifr_age_vacc_class) <- c(n_groups, n_vacc_classes)

initial(ihr_age_vacc_class[, ]) <- 0
update(ihr_age_vacc_class[, ]) <- sum(IHR_disag_weighted[i, , j]) /
  sum(new_inf[i, , j])
dim(ihr_age_vacc_class) <- c(n_groups, n_vacc_classes)

initial(hfr_age_vacc_class[, ]) <- 0
update(hfr_age_vacc_class[, ]) <- sum(HFR_disag_weighted[i, , j]) /
  sum(n_I_C_2_to_hosp[i, , j])
dim(hfr_age_vacc_class) <- c(n_groups, n_vacc_classes)


config(compare) <- "compare_lancelot.cpp"
## Parameters and code to support the compare function. Because these
## do not appear in any odin equation we mark them as "ignore.unused"
## so that odin doesn't complain that they appear redundant. This
## means that unwanted things can accumulate here so keep on top of
## it.
N_tot_all <- user() # ignore.unused
N_tot_over25 <- user() # ignore.unused
N_tot_under15 <- user() # ignore.unused
N_tot_15_24 <- user() # ignore.unused
N_tot_25_49 <- user() # ignore.unused
N_tot_50_64 <- user() # ignore.unused
N_tot_65_79 <- user() # ignore.unused
N_tot_80_plus <- user() # ignore.unused
N_tot_react <- user() # ignore.unused
N_5_24_react <- user() # ignore.unused
N_25_34_react <- user() # ignore.unused
N_35_44_react <- user() # ignore.unused
N_45_54_react <- user() # ignore.unused
N_55_64_react <- user() # ignore.unused
N_65_plus_react <- user() # ignore.unused
N_tot_15_64 <- user() # ignore.unused

p_NC_under15 <- user() # ignore.unused
p_NC_weekend_under15 <- user() # ignore.unused
p_NC_15_24 <- user() # ignore.unused
p_NC_weekend_15_24 <- user() # ignore.unused
p_NC_25_49 <- user() # ignore.unused
p_NC_weekend_25_49 <- user() # ignore.unused
p_NC_50_64 <- user() # ignore.unused
p_NC_weekend_50_64 <- user() # ignore.unused
p_NC_65_79 <- user() # ignore.unused
p_NC_weekend_65_79 <- user() # ignore.unused
p_NC_80_plus <- user() # ignore.unused
p_NC_weekend_80_plus <- user() # ignore.unused
pillar2_sensitivity <- user() # ignore.unused
pillar2_specificity <- user() # ignore.unused
react_sensitivity <- user() # ignore.unused
react_specificity <- user() # ignore.unused
sero_sensitivity_1 <- user() # ignore.unused
sero_specificity_1 <- user() # ignore.unused
sero_sensitivity_2 <- user() # ignore.unused
sero_specificity_2 <- user() # ignore.unused
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
kappa_pillar2_cases <- user() # ignore.unused
phi_pillar2_cases_under15 <- user() # ignore.unused
phi_pillar2_cases_weekend_under15 <- user() # ignore.unused
phi_pillar2_cases_15_24 <- user() # ignore.unused
phi_pillar2_cases_weekend_15_24 <- user() # ignore.unused
phi_pillar2_cases_25_49 <- user() # ignore.unused
phi_pillar2_cases_weekend_25_49 <- user() # ignore.unused
phi_pillar2_cases_50_64 <- user() # ignore.unused
phi_pillar2_cases_weekend_50_64 <- user() # ignore.unused
phi_pillar2_cases_65_79 <- user() # ignore.unused
phi_pillar2_cases_weekend_65_79 <- user() # ignore.unused
phi_pillar2_cases_80_plus <- user() # ignore.unused
phi_pillar2_cases_weekend_80_plus <- user() # ignore.unused
