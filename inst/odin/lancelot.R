## E and R stage indexed by i, j, k, l with
## i for the age group
## j for the strain
## k for the progression (not exponential latent and infectious period)
## l for the vacc. group

## Number of classes (age & vaccination)

## Number of "groups", being the age classes, Carehome workers and
## Carehome residents. This will be 19 in all but experimental uses.
n_age_groups <- parameter()
n_groups <- parameter()
has_carehomes <- parameter(type = "integer")

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
update(D_hosp[, ]) <- D_hosp[i, j] + delta_D_hosp_disag[i, j]
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
  (if (has_carehomes == 1) sum(delta_diagnoses_admitted[18, ]) * 1 / 8
   else 0)

delta_all_admission_30_39_conf <-
  sum(delta_diagnoses_admitted[7:8, ]) +
  (if (has_carehomes == 1) sum(delta_diagnoses_admitted[18, ]) * 2 / 8
   else 0)

delta_all_admission_40_49_conf <-
  sum(delta_diagnoses_admitted[9:10, ]) +
  (if (has_carehomes == 1) sum(delta_diagnoses_admitted[18, ]) * 2 / 8
   else 0)

delta_all_admission_50_59_conf <-
  sum(delta_diagnoses_admitted[11:12, ]) +
  (if (has_carehomes == 1) sum(delta_diagnoses_admitted[18, ]) * 2 / 8
   else 0)

delta_all_admission_60_69_conf <-
  sum(delta_diagnoses_admitted[13:14, ]) +
  (if (has_carehomes == 1) sum(delta_diagnoses_admitted[18, ]) * 1 / 8 +
  sum(delta_diagnoses_admitted[19, ]) * 0.05
  else 0)

delta_all_admission_70_79_conf <-
  sum(delta_diagnoses_admitted[15:16, ]) +
  (if (has_carehomes == 1) sum(delta_diagnoses_admitted[19, ]) * 0.2
   else 0)

delta_all_admission_80_plus_conf <-
  sum(delta_diagnoses_admitted[17, ]) +
  (if (has_carehomes == 1) sum(delta_diagnoses_admitted[19, ]) * 0.75
   else 0)

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

delta_infections[, , ] <-
  n_S_progress[i, j, k] +
  (if (j > 2)
    n_RE[i, j, k]
   else
     0)
dim(delta_infections) <- c(n_groups, n_strains, n_vacc_classes)

initial(cum_infections_disag[, ]) <- 0
update(cum_infections_disag[, ]) <- cum_infections_disag[i, j] +
  sum(delta_infections[i, , j])
dim(cum_infections_disag) <- c(n_groups, n_vacc_classes)

initial(admit_conf_inc, zero_every = 1) <- 0
update(admit_conf_inc) <- admit_conf_inc + delta_admit_conf

initial(new_conf_inc, zero_every = 1) <- 0
update(new_conf_inc) <- new_conf_inc + delta_new_conf

# Admissions + confirmed by age

initial(all_admission_0_9_conf_inc, zero_every = 1) <- 0
update(all_admission_0_9_conf_inc) <- all_admission_0_9_conf_inc +
  delta_all_admission_0_9_conf

initial(all_admission_10_19_conf_inc, zero_every = 1) <- 0
update(all_admission_10_19_conf_inc) <- all_admission_10_19_conf_inc +
  delta_all_admission_10_19_conf

initial(all_admission_20_29_conf_inc, zero_every = 1) <- 0
update(all_admission_20_29_conf_inc) <- all_admission_20_29_conf_inc +
  delta_all_admission_20_29_conf

initial(all_admission_30_39_conf_inc, zero_every = 1) <- 0
update(all_admission_30_39_conf_inc) <- all_admission_30_39_conf_inc +
  delta_all_admission_30_39_conf

initial(all_admission_40_49_conf_inc, zero_every = 1) <- 0
update(all_admission_40_49_conf_inc) <- all_admission_40_49_conf_inc +
  delta_all_admission_40_49_conf

initial(all_admission_50_59_conf_inc, zero_every = 1) <- 0
update(all_admission_50_59_conf_inc) <- all_admission_50_59_conf_inc +
  delta_all_admission_50_59_conf

initial(all_admission_60_69_conf_inc, zero_every = 1) <- 0
update(all_admission_60_69_conf_inc) <- all_admission_60_69_conf_inc +
  delta_all_admission_60_69_conf

initial(all_admission_70_79_conf_inc, zero_every = 1) <- 0
update(all_admission_70_79_conf_inc) <- all_admission_70_79_conf_inc +
  delta_all_admission_70_79_conf

initial(all_admission_80_plus_conf_inc, zero_every = 1) <- 0
update(all_admission_80_plus_conf_inc) <- all_admission_80_plus_conf_inc +
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
p_T_PCR_pre_progress[] <- 1 - exp(-gamma_PCR_pre[i] * dt)
p_T_PCR_pos_progress[] <- 1 - exp(-gamma_PCR_pos[i] * dt)

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
dim(p_T_PCR_pre_progress) <- n_strains
dim(p_T_PCR_pos_progress) <- n_strains


## Work out time-varying probabilities
p_C_t <- interpolate(p_C_time, p_C_value, "linear")
p_C[, , ] <- min(p_C_t[i] * rel_p_sympt[i, j, k] * strain_rel_p_sympt[j],
                 as.numeric(1))

p_H_t <- interpolate(p_H_time, p_H_value, "linear")
p_H[, , ] <- min(p_H_t[i] * rel_p_hosp_if_sympt[i, j, k] *
                   strain_rel_p_hosp_if_sympt[j], as.numeric(1))

p_ICU_t <- interpolate(p_ICU_time, p_ICU_value, "linear")
p_ICU[, , ] <- min(p_ICU_t[i] * rel_p_ICU[i, j, k] * strain_rel_p_icu[j],
                   as.numeric(1))

p_ICU_D_t <- interpolate(p_ICU_D_time, p_ICU_D_value, "linear")
p_ICU_D[, , ] <- min(p_ICU_D_t[i] * rel_p_ICU_D[i, j, k] *
                       strain_rel_p_ICU_D[j], as.numeric(1))

p_H_D_t <- interpolate(p_H_D_time, p_H_D_value, "linear")
p_H_D[, , ] <- min(p_H_D_t[i] * rel_p_H_D[i, j, k] * strain_rel_p_H_D[j],
                   as.numeric(1))

p_W_D_t <- interpolate(p_W_D_time, p_W_D_value, "linear")
p_W_D[, , ] <- min(p_W_D_t[i] * rel_p_W_D[i, j, k] * strain_rel_p_W_D[j],
                   as.numeric(1))

p_G_D_t <- interpolate(p_G_D_time, p_G_D_value, "linear")
p_G_D[, , ] <- min(p_G_D_t[i] * rel_p_G_D[i, j, k] * strain_rel_p_G_D[j],
                   as.numeric(1))

p_R_t <- interpolate(p_R_time, p_R_value, "linear")
p_R[, , ] <- min(p_R_t[i] * rel_p_R[i, j, k], as.numeric(1))

p_star <- interpolate(p_star_time, p_star_value, "linear")

## Work out time-varying gammas
gamma_E_t <- interpolate(gamma_E_time, gamma_E_value, "linear")
gamma_E[] <- gamma_E_t * rel_gamma_E[i]

gamma_A_t <- interpolate(gamma_A_time, gamma_A_value, "linear")
gamma_A[] <- gamma_A_t * rel_gamma_A[i]

gamma_P_t <- interpolate(gamma_P_time, gamma_P_value, "linear")
gamma_P[] <- gamma_P_t * rel_gamma_P[i]

gamma_C_1_t <- interpolate(gamma_C_1_time, gamma_C_1_value, "linear")
gamma_C_1[] <- gamma_C_1_t * rel_gamma_C_1[i]

gamma_C_2_t <- interpolate(gamma_C_2_time, gamma_C_2_value, "linear")
gamma_C_2[] <- gamma_C_2_t * rel_gamma_C_2[i]

gamma_G_D_t <- interpolate(gamma_G_D_time, gamma_G_D_value, "linear")
gamma_G_D[] <- gamma_G_D_t * rel_gamma_G_D[i]

gamma_ICU_pre_t <-
  interpolate(gamma_ICU_pre_time, gamma_ICU_pre_value, "linear")
gamma_ICU_pre[] <- gamma_ICU_pre_t * rel_gamma_ICU_pre[i]

gamma_H_R_t <- interpolate(gamma_H_R_time, gamma_H_R_value, "linear")
gamma_H_R[] <- gamma_H_R_t * rel_gamma_H_R[i]

gamma_H_D_t <- interpolate(gamma_H_D_time, gamma_H_D_value, "linear")
gamma_H_D[] <- gamma_H_D_t * rel_gamma_H_D[i]

gamma_ICU_W_R_t <-
  interpolate(gamma_ICU_W_R_time, gamma_ICU_W_R_value, "linear")
gamma_ICU_W_R[] <- gamma_ICU_W_R_t * rel_gamma_ICU_W_R[i]

gamma_ICU_W_D_t <-
  interpolate(gamma_ICU_W_D_time, gamma_ICU_W_D_value, "linear")
gamma_ICU_W_D[] <- gamma_ICU_W_D_t * rel_gamma_ICU_W_D[i]

gamma_ICU_D_t <- interpolate(gamma_ICU_D_time, gamma_ICU_D_value, "linear")
gamma_ICU_D[] <- gamma_ICU_D_t * rel_gamma_ICU_D[i]

gamma_W_R_t <- interpolate(gamma_W_R_time, gamma_W_R_value, "linear")
gamma_W_R[] <- gamma_W_R_t * rel_gamma_W_R[i]

gamma_W_D_t <- interpolate(gamma_W_D_time, gamma_W_D_value, "linear")
gamma_W_D[] <- gamma_W_D_t * rel_gamma_W_D[i]

gamma_PCR_pre_t <-
  interpolate(gamma_PCR_pre_time, gamma_PCR_pre_value, "linear")
gamma_PCR_pre[] <- gamma_PCR_pre_t * rel_gamma_PCR_pre[i]

gamma_PCR_pos_t <-
  interpolate(gamma_PCR_pos_time, gamma_PCR_pos_value, "linear")
gamma_PCR_pos[] <- gamma_PCR_pos_t * rel_gamma_PCR_pos[i]

gamma_U <- interpolate(gamma_U_time, gamma_U_value, "linear")

## Draws from binomial distributions for numbers changing between
## compartments:

## modelling infections and vaccine progression, which can happen simultaneously

#### flow out of S ####

## new infections

## Seeding of first wave, we seed in group 4, strain 1, vaccine stratum 1
##
seed_rate <- interpolate(seed_time, seed_value, "constant")
seed <- Poisson(seed_rate)
seed_age_band <- as.integer(4) # 15-19y band

seed_time <- parameter()
seed_value <- parameter()
dim(seed_time) <- parameter(rank = 1)
dim(seed_value) <- parameter(rank = 1)

## Introduction of new strains. n_S_progress is arranged as:
##
## [age, strain infected with, vaccine stage]
##
## As in the model initialisation we will use the teenager category,
## and only infect *unvaccinated* people. For now we will model only
## movement into the second compartment as that represents our "new"
## strain.
strain_seed_rate <- interpolate(strain_seed_time, strain_seed_value, "constant")
strain_seed <- Poisson(strain_seed_rate)

strain_seed_time <- parameter()
strain_seed_value <- parameter()
dim(strain_seed_time) <- parameter(rank = 1)
dim(strain_seed_value) <- parameter(rank = 1)

## Compute the new infections with multiple strains using nested binomials
## No one can move from S to E3 or E4
n_S_progress_tot[, ] <- Binomial(S[i, j], p_SE[i, j])
n_S_progress[, , ] <-
  if (j == 1 || n_real_strains == 1)
    Binomial(n_S_progress_tot[i, k], rel_foi_strain[i, j, k]) else
      (if (j == 2) n_S_progress_tot[i, k] - n_S_progress[i, 1, k] else 0)

n_S_progress[seed_age_band, 1, 1] <-
  n_S_progress[seed_age_band, 1, 1] +
  min(S[seed_age_band, 1] - n_S_progress_tot[seed_age_band, 1], seed)

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
n_S_next_vacc_class[, ] <-
  Binomial(S[i, j] - sum(n_S_progress[i, , j]), p_S_next_vacc_class[i, j])

n_S_vacc_skip[, ] <-
  Binomial(S[i, j] - sum(n_S_progress[i, , j]) -
             n_S_next_vacc_class[i, j], p_S_vacc_skip[i, j])

#### flow out of E ####

n_E_progress[, , , ] <- Binomial(E[i, j, k, l], p_E_progress[j])

## vaccine progression
n_E_next_vacc_class[, , , ] <-
  Binomial(E[i, j, k, l] - n_E_progress[i, j, k, l],
           p_E_next_vacc_class[i, j, k, l])

n_E_vacc_skip[, , , ] <-
  Binomial(E[i, j, k, l] - n_E_progress[i, j, k, l] -
             n_E_next_vacc_class[i, j, k, l], p_E_vacc_skip[i, j, k, l])

#### flow out of I_A ####

n_I_A_progress[, , , ] <- Binomial(I_A[i, j, k, l], p_I_A_progress[j])

## vaccine progression
n_I_A_next_vacc_class[, , , ] <-
  Binomial(I_A[i, j, k, l] - n_I_A_progress[i, j, k, l],
           p_I_A_next_vacc_class[i, j, k, l])


n_I_A_vacc_skip[, , , ] <-
  Binomial(I_A[i, j, k, l] - n_I_A_progress[i, j, k, l] -
             n_I_A_next_vacc_class[i, j, k, l], p_I_A_vacc_skip[i, j, k, l])

#### flow out of I_P ####

n_I_P_progress[, , , ] <- Binomial(I_P[i, j, k, l], p_I_P_progress[j])

## vaccine progression
n_I_P_next_vacc_class[, , , ] <-
  Binomial(I_P[i, j, k, l] - n_I_P_progress[i, j, k, l],
           p_I_P_next_vacc_class[i, j, k, l])

n_I_P_vacc_skip[, , , ] <-
  Binomial(I_P[i, j, k, l] - n_I_P_progress[i, j, k, l] -
             n_I_P_next_vacc_class[i, j, k, l], p_I_P_vacc_skip[i, j, k, l])

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
rate_RE_progress[, , ] <- if (n_strains == 1) 0 else
  lambda_susc[i, 3 - j, k] * (1 - cross_immunity[j])
dim(rate_RE_progress) <- c(n_groups, n_real_strains, n_vacc_classes)

rate_R_progress[, , ] <- waning_rate[i] +
  if (n_strains == 1) 0 else
    if (j == 1 || j == 4)
    rate_RE_progress[i, 1, k] else
      if (j == 5) sum(rate_RE_progress[i, , k]) else 0

p_R_progress[, , ] <- 1 - exp(-rate_R_progress[i, j, k] * dt)

## n_R_progress is total number who either:
##  - leave R for S or E and stay in same vacc class
##  - leave R for S or E and change same vacc class
n_R_progress[, , ] <- Binomial(R[i, j, k], p_R_progress[i, j, k])

## Number going from R to S
## In one-strain model, all progress to S only
## In multi-strain model, R3 and R4 progress to S only
## If not modelling super infection then all progress to S only
##
## In multi-strain model, number of R1 and R2 to S is binomial w.p. waning over
##  waning plus prob strain
## TODO (RS): waning_rate should eventually be variant varying
p_RS[, , ] <- if (n_strains == 1) 1 else
  (if (waning_rate[i] == 0) 0 else
    waning_rate[i] / rate_R_progress[i, j, k])
n_RS[, , ] <- Binomial(n_R_progress[i, j, k], p_RS[i, j, k])

p_R5_to_E3[, ] <- if (n_strains == 1) 1 else
  if (rate_RE_progress[i, 1, j] == 0) 0 else
    (rate_RE_progress[i, 1, j] /
       sum(rate_RE_progress[i, , j]))

n_R5_to_E3[, ] <- if (n_strains == 1) 0 else
  Binomial(n_R_progress[i, 5, j] - n_RS[i, 5, j], p_R5_to_E3[i, j])

dim(p_R5_to_E3) <- c(n_groups, n_vacc_classes)
dim(n_R5_to_E3) <- c(n_groups, n_vacc_classes)

## n_RE[i, j, k] is the number going from
## R[age i, , vacc class k] to
## E[age i, strain j, vacc class k]
## R1 and R4 can go to E3, R5 can go to E3 or E4
n_RE[, , ] <- if (n_strains == 1 || j < 3) 0 else
  if (j == 3) (n_R_progress[i, 1, k] - n_RS[i, 1, k] +
                 n_R_progress[i, 4, k] - n_RS[i, 4, k] + n_R5_to_E3[i, k]) else
                   (n_R_progress[i, 5, k] - n_RS[i, 5, k] - n_R5_to_E3[i, k])

## R -> R vaccine progression
n_R_tmp[, , ] <- R[i, j, k] - n_R_progress[i, j, k]
n_R_next_vacc_class[, , ] <-
  Binomial(n_R_tmp[i, j, k], p_R_next_vacc_class[i, j, k])

n_R_vacc_skip[, , ] <- Binomial(n_R_tmp[i, j, k] - n_R_next_vacc_class[i, j, k],
                                p_R_vacc_skip[i, j, k])

#### other transitions ####

n_I_C_1_progress[, , , ] <- Binomial(I_C_1[i, j, k, l], p_I_C_1_progress[j])
n_I_C_2_progress[, , , ] <- Binomial(I_C_2[i, j, k, l], p_I_C_2_progress[j])
n_G_D_progress[, , , ] <- Binomial(G_D[i, j, k, l], p_G_D_progress[j])
n_ICU_pre_unconf_progress[, , , ] <-
  Binomial(ICU_pre_unconf[i, j, k, l], p_ICU_pre_progress[j])
n_ICU_pre_conf_progress[, , , ] <-
  Binomial(ICU_pre_conf[i, j, k, l], p_ICU_pre_progress[j])
n_H_R_unconf_progress[, , , ] <-
  Binomial(H_R_unconf[i, j, k, l], p_H_R_progress[j])
n_H_R_conf_progress[, , , ] <-
  Binomial(H_R_conf[i, j, k, l], p_H_R_progress[j])
n_H_D_unconf_progress[, , , ] <-
  Binomial(H_D_unconf[i, j, k, l], p_H_D_progress[j])
n_H_D_conf_progress[, , , ] <-
  Binomial(H_D_conf[i, j, k, l], p_H_D_progress[j])
n_ICU_W_R_unconf_progress[, , , ] <-
  Binomial(ICU_W_R_unconf[i, j, k, l], p_ICU_W_R_progress[j])
n_ICU_W_R_conf_progress[, , , ] <-
  Binomial(ICU_W_R_conf[i, j, k, l], p_ICU_W_R_progress[j])
n_ICU_W_D_unconf_progress[, , , ] <-
  Binomial(ICU_W_D_unconf[i, j, k, l], p_ICU_W_D_progress[j])
n_ICU_W_D_conf_progress[, , , ] <-
  Binomial(ICU_W_D_conf[i, j, k, l], p_ICU_W_D_progress[j])
n_ICU_D_unconf_progress[, , , ] <-
  Binomial(ICU_D_unconf[i, j, k, l], p_ICU_D_progress[j])
n_ICU_D_conf_progress[, , , ] <-
  Binomial(ICU_D_conf[i, j, k, l], p_ICU_D_progress[j])
n_W_R_unconf_progress[, , , ] <-
  Binomial(W_R_unconf[i, j, k, l], p_W_R_progress[j])
n_W_R_conf_progress[, , , ] <-
  Binomial(W_R_conf[i, j, k, l], p_W_R_progress[j])
n_W_D_unconf_progress[, , , ] <-
  Binomial(W_D_unconf[i, j, k, l], p_W_D_progress[j])
n_W_D_conf_progress[, , , ] <-
  Binomial(W_D_conf[i, j, k, l], p_W_D_progress[j])
n_T_sero_pre_1_progress[, , , ] <-
  Binomial(T_sero_pre_1[i, j, k, l], p_T_sero_pre_1_progress)
n_T_sero_pos_1_progress[, , , ] <-
  Binomial(T_sero_pos_1[i, j, k, l], p_T_sero_pos_1_progress)
n_T_sero_pre_2_progress[, , , ] <-
  Binomial(T_sero_pre_2[i, j, k, l], p_T_sero_pre_2_progress)
n_T_sero_pos_2_progress[, , , ] <-
  Binomial(T_sero_pos_2[i, j, k, l], p_T_sero_pos_2_progress)
n_T_PCR_pre_progress[, , , ] <-
  Binomial(T_PCR_pre[i, j, k, l], p_T_PCR_pre_progress[j])
n_T_PCR_pos_progress[, , , ] <-
  Binomial(T_PCR_pos[i, j, k, l], p_T_PCR_pos_progress[j])

## Cumulative infections, summed over all age groups
initial(cum_infections) <- 0
delta_infections_total <- sum(delta_infections)
update(cum_infections) <- cum_infections + delta_infections_total

initial(infections_inc, zero_every = 1) <- 0
new_infections_inc <- infections_inc + delta_infections_total
update(infections_inc) <- new_infections_inc

initial(cum_infections_strain[]) <- 0
delta_infections_strain[] <- sum(delta_infections[, i, ])
update(cum_infections_strain[]) <-
  cum_infections_strain[i] + delta_infections_strain[i]
dim(delta_infections_strain) <- n_strains
dim(cum_infections_strain) <- n_strains

initial(infections_inc_strain[], zero_every = 1) <- 0
new_infections_inc_strain[] <-
  infections_inc_strain[i] + delta_infections_strain[i]
update(infections_inc_strain[]) <- new_infections_inc_strain[i]
dim(new_infections_inc_strain) <- n_strains
dim(infections_inc_strain) <- n_strains

initial(infections_inc_age[], zero_every = 1) <- 0
delta_infections_age[] <- sum(delta_infections[i, , ])
new_infections_inc_age[] <- infections_inc_age[i] + delta_infections_age[i]
update(infections_inc_age[]) <- new_infections_inc_age[i]
dim(delta_infections_age) <- n_groups
dim(new_infections_inc_age) <- n_groups
dim(infections_inc_age) <- n_groups

## Hospitalisations
initial(hospitalisations_inc, zero_every = 1) <- 0
delta_hospitalisations_total <- sum(n_I_C_2_to_hosp)
new_hospitalisations_inc <- hospitalisations_inc + delta_hospitalisations_total
update(hospitalisations_inc) <- new_hospitalisations_inc

initial(hospitalisations_inc_strain[], zero_every = 1) <- 0
delta_hospitalisations_strain[] <- sum(n_I_C_2_to_hosp[, i, ])
new_hospitalisations_inc_strain[] <-
  hospitalisations_inc_strain[i] + delta_hospitalisations_strain[i]
update(hospitalisations_inc_strain[]) <-
  new_hospitalisations_inc_strain[i]
dim(delta_hospitalisations_strain) <- n_strains
dim(new_hospitalisations_inc_strain) <- n_strains
dim(hospitalisations_inc_strain) <- n_strains

initial(hospitalisations_inc_age[], zero_every = 1) <- 0
delta_hospitalisations_age[] <- sum(n_I_C_2_to_hosp[i, , ])
new_hospitalisations_inc_age[] <-
  hospitalisations_inc_age[i] + delta_hospitalisations_age[i]
update(hospitalisations_inc_age[]) <- new_hospitalisations_inc_age[i]
dim(delta_hospitalisations_age) <- n_groups
dim(new_hospitalisations_inc_age) <- n_groups
dim(hospitalisations_inc_age) <- n_groups

## Work out the new S (i for age, j for vaccination status)
new_S[, ] <- S[i, j] + sum(n_RS[i, , j]) + sum(n_infected_to_S[i, , j]) -
  sum(n_S_progress[i, , j]) - n_S_next_vacc_class[i, j]
new_S[, ] <- new_S[i, j] +
  (if (j == 1) n_S_next_vacc_class[i, n_vacc_classes] else
    n_S_next_vacc_class[i, j - 1]) -
  (if (vacc_skip_to[j] > 0) n_S_vacc_skip[i, j] else 0) +
  (if (vacc_skip_from[j] > 0) n_S_vacc_skip[i, vacc_skip_from[j]] else 0)

## Computes the number of asymptomatic
n_EI_A[, , ] <- Binomial(n_E_progress[i, j, k_E, k], 1 - p_C[i, j, k])

## Computes the number of symptomatic cases
n_EI_P[, , ] <- n_E_progress[i, j, k_E, k] - n_EI_A[i, j, k]

## Work out the S->E and E->E transitions
new_E[, , , ] <- E[i, j, k, l] +
  (if (k == 1) n_S_progress[i, j, l] +
     (if (j > 2) n_RE[i, j, l] else 0)
   else n_E_progress[i, j, k - 1, l]) -
  n_E_progress[i, j, k, l] -
  n_E_next_vacc_class[i, j, k, l] +
  (if (l == 1)
    n_E_next_vacc_class[i, j, k, n_vacc_classes]
   else
     n_E_next_vacc_class[i, j, k, l - 1]) -
  (if (vacc_skip_to[l] > 0) n_E_vacc_skip[i, j, k, l] else 0) +
  (if (vacc_skip_from[l] > 0) n_E_vacc_skip[i, j, k, vacc_skip_from[l]] else 0)

## Work out the I_A->I_A transitions
new_I_A[, , , ] <- I_A[i, j, k, l] +
  (if (k == 1) n_EI_A[i, j, l] else
    n_I_A_progress[i, j, k - 1, l]) -
    n_I_A_progress[i, j, k, l] -
    n_I_A_next_vacc_class[i, j, k, l] +
  (if (l == 1) n_I_A_next_vacc_class[i, j, k, n_vacc_classes] else
    n_I_A_next_vacc_class[i, j, k, l - 1]) -
  (if (vacc_skip_to[l] > 0) n_I_A_vacc_skip[i, j, k, l] else 0) +
  (if (vacc_skip_from[l] > 0)
    n_I_A_vacc_skip[i, j, k, vacc_skip_from[l]] else 0)

## Work out the I_P->I_P transitions
new_I_P[, , , ] <- I_P[i, j, k, l] +
  (if (k == 1) n_EI_P[i, j, l] else
    n_I_P_progress[i, j, k - 1, l]) -
  n_I_P_progress[i, j, k, l] -
  n_I_P_next_vacc_class[i, j, k, l] +
  (if (l == 1) n_I_P_next_vacc_class[i, j, k, n_vacc_classes] else
    n_I_P_next_vacc_class[i, j, k, l - 1]) -
  (if (vacc_skip_to[l] > 0) n_I_P_vacc_skip[i, j, k, l] else 0) +
  (if (vacc_skip_from[l] > 0)
    n_I_P_vacc_skip[i, j, k, vacc_skip_from[l]] else 0)

## Work out the I_C_1->I_C_1 transitions
new_I_C_1[, , , ] <- I_C_1[i, j, k, l] +
  (if (k == 1) n_I_P_progress[i, j, k_P, l] else
    n_I_C_1_progress[i, j, k - 1, l]) - n_I_C_1_progress[i, j, k, l]

## Work out the I_C_2->I_C_2 transitions
new_I_C_2[, , , ] <- I_C_2[i, j, k, l] +
  (if (k == 1) n_I_C_1_progress[i, j, k_C_1, l] else
    n_I_C_2_progress[i, j, k - 1, l]) - n_I_C_2_progress[i, j, k, l]

## Work out the flow from I_C_2 -> R, G_D, hosp
n_I_C_2_to_RS[, , ] <-
  Binomial(n_I_C_2_progress[i, j, k_C_2, k], 1 - p_H[i, j, k])
n_I_C_2_to_G_D[, , ] <- Binomial(n_I_C_2_progress[i, j, k_C_2, k] -
                                   n_I_C_2_to_RS[i, j, k], p_G_D[i, j, k])
n_I_C_2_to_hosp[, , ] <- n_I_C_2_progress[i, j, k_C_2, k] -
  n_I_C_2_to_RS[i, j, k] - n_I_C_2_to_G_D[i, j, k]

## Work out the G_D -> G_D transitions
new_G_D[, , , ] <- G_D[i, j, k, l] +
  (if (k == 1) n_I_C_2_to_G_D[i, j, l] else
    n_G_D_progress[i, j, k - 1, l]) - n_G_D_progress[i, j, k, l]

## Work out the split in hospitals between H_D, H_R and ICU_pre
n_I_C_2_to_ICU_pre[, , ] <- Binomial(n_I_C_2_to_hosp[i, j, k], p_ICU[i, j, k])
n_I_C_2_to_ICU_pre_conf[, , ] <-
  Binomial(n_I_C_2_to_ICU_pre[i, j, k], p_star[i])
n_hosp_non_ICU[, , ] <- n_I_C_2_to_hosp[i, j, k] - n_I_C_2_to_ICU_pre[i, j, k]
n_I_C_2_to_H_D[, , ] <- Binomial(n_hosp_non_ICU[i, j, k], p_H_D[i, j, k])
n_I_C_2_to_H_D_conf[, , ] <- Binomial(n_I_C_2_to_H_D[i, j, k], p_star[i])
n_I_C_2_to_H_R[, , ] <- n_hosp_non_ICU[i, j, k] - n_I_C_2_to_H_D[i, j, k]
n_I_C_2_to_H_R_conf[, , ] <- Binomial(n_I_C_2_to_H_R[i, j, k], p_star[i])

## Work out the ICU_pre -> ICU_pre transitions
aux_ICU_pre_unconf[, , , ] <- ICU_pre_unconf[i, j, k, l] +
  (if (k > 1) n_ICU_pre_unconf_progress[i, j, k - 1, l] else 0) -
  n_ICU_pre_unconf_progress[i, j, k, l]
aux_ICU_pre_conf[, , , ] <- ICU_pre_conf[i, j, k, l] +
  (if (k > 1) n_ICU_pre_conf_progress[i, j, k - 1, l] else 0) -
  n_ICU_pre_conf_progress[i, j, k, l]

n_ICU_pre_unconf_to_conf[, , , ] <-
  Binomial(aux_ICU_pre_unconf[i, j, k, l], p_test)

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

n_H_R_unconf_to_conf[, , , ] <- Binomial(aux_H_R_unconf[i, j, k, l], p_test)

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

n_H_D_unconf_to_conf[, , , ] <- Binomial(aux_H_D_unconf[i, j, k, l], p_test)

new_H_D_unconf[, , , ] <-
  aux_H_D_unconf[i, j, k, l] - n_H_D_unconf_to_conf[i, j, k, l] +
  (if (k == 1) n_I_C_2_to_H_D[i, j, l] - n_I_C_2_to_H_D_conf[i, j, l] else 0)
new_H_D_conf[, , , ] <-
  aux_H_D_conf[i, j, k, l] + n_H_D_unconf_to_conf[i, j, k, l] +
  (if (k == 1) n_I_C_2_to_H_D_conf[i, j, l] else 0)

## Work out the ICU_pre to ICU_D, ICU_W_R and ICU_W_D splits
n_ICU_pre_unconf_to_ICU_D_unconf[, , ] <-
  Binomial(n_ICU_pre_unconf_progress[i, j, k_ICU_pre, k], p_ICU_D[i, j, k])
n_ICU_pre_conf_to_ICU_D_conf[, , ] <-
  Binomial(n_ICU_pre_conf_progress[i, j, k_ICU_pre, k], p_ICU_D[i, j, k])
n_ICU_pre_unconf_to_ICU_W_D_unconf[, , ] <-
  Binomial(n_ICU_pre_unconf_progress[i, j, k_ICU_pre, k] -
             n_ICU_pre_unconf_to_ICU_D_unconf[i, j, k], p_W_D[i, j, k])
n_ICU_pre_unconf_to_ICU_W_R_unconf[, , ] <-
  n_ICU_pre_unconf_progress[i, j, k_ICU_pre, k] -
  n_ICU_pre_unconf_to_ICU_D_unconf[i, j, k] -
  n_ICU_pre_unconf_to_ICU_W_D_unconf[i, j, k]
n_ICU_pre_conf_to_ICU_W_D_conf[, , ] <-
  Binomial(n_ICU_pre_conf_progress[i, j, k_ICU_pre, k] -
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
  Binomial(aux_ICU_W_R_unconf[i, j, k, l], p_test)
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
  Binomial(aux_ICU_W_D_unconf[i, j, k, l], p_test)
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

n_ICU_D_unconf_to_conf[, , , ] <- Binomial(aux_ICU_D_unconf[i, j, k, l], p_test)
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

n_W_R_unconf_to_conf[, , , ] <- Binomial(aux_W_R_unconf[i, j, k, l], p_test)
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

n_W_D_unconf_to_conf[, , , ] <- Binomial(aux_W_D_unconf[i, j, k, l], p_test)
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
  Binomial(n_T_sero_pre_1_progress[i, j, k_sero_pre_1, k], p_sero_pos_1[i])

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
  Binomial(n_T_sero_pre_2_progress[i, j, k_sero_pre_2, k], p_sero_pos_2[i])

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

n_infected_to_R[, , ] <- Binomial(n_infection_end[i, j, k], p_R[i, j, k])

n_infected_to_S[, , ] <- n_infection_end[i, j, k] - n_infected_to_R[i, j, k]

## Work out the total number of recovery
new_R[, , ] <- R[i, j, k] -
  n_R_progress[i, j, k] -
  n_R_next_vacc_class[i, j, k] +
  (if (n_strains == 1 || j < 5) n_infected_to_R[i, j, k] else 0) +
  (if (k == 1) n_R_next_vacc_class[i, j, n_vacc_classes] else
    n_R_next_vacc_class[i, j, k - 1])  -
  (if (vacc_skip_to[k] > 0) n_R_vacc_skip[i, j, k] else 0) +
  (if (vacc_skip_from[k] > 0) n_R_vacc_skip[i, j, vacc_skip_from[k]] else 0)

## Work out the PCR positivity
new_T_PCR_pre[, , , ] <- T_PCR_pre[i, j, k, l] -
  n_T_PCR_pre_progress[i, j, k, l] +
  (if (k == 1) n_S_progress[i, j, l] + n_RE[i, j, l] else
    n_T_PCR_pre_progress[i, j, k - 1, l])

new_T_PCR_pos[, , , ] <- T_PCR_pos[i, j, k, l] -
  n_T_PCR_pos_progress[i, j, k, l] +
  (if (k == 1) n_T_PCR_pre_progress[i, j, k_PCR_pre, l] else
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
initial(D_hosp[, ]) <- 0
initial(D_non_hosp[]) <- 0
initial(T_PCR_pre[, , , ]) <- 0
initial(T_PCR_pos[, , , ]) <- 0
initial(T_PCR_neg[, , ]) <- 0
initial(cum_admit_conf) <- 0
initial(cum_new_conf) <- 0
initial(cum_admit_by_age[]) <- 0

## User defined parameters - default in parentheses:

## Vaccination/strain effect parameters
n_vacc_classes <- parameter()
rel_susceptibility <- parameter()
dim(rel_susceptibility) <- c(n_groups, n_strains, n_vacc_classes)
rel_p_sympt <- parameter()
dim(rel_p_sympt) <- c(n_groups, n_strains, n_vacc_classes)
strain_rel_p_sympt <- parameter()
dim(strain_rel_p_sympt) <- n_strains
rel_p_hosp_if_sympt <- parameter()
dim(rel_p_hosp_if_sympt) <- c(n_groups, n_strains, n_vacc_classes)
strain_rel_p_hosp_if_sympt <- parameter()
dim(strain_rel_p_hosp_if_sympt) <- n_strains
rel_p_ICU <- parameter()
dim(rel_p_ICU) <- c(n_groups, n_strains, n_vacc_classes)
strain_rel_p_icu <- parameter()
dim(strain_rel_p_icu) <- n_strains
rel_p_ICU_D <- parameter()
dim(rel_p_ICU_D) <- c(n_groups, n_strains, n_vacc_classes)
rel_p_H_D <- parameter()
dim(rel_p_H_D) <- c(n_groups, n_strains, n_vacc_classes)
rel_p_W_D <- parameter()
dim(rel_p_W_D) <- c(n_groups, n_strains, n_vacc_classes)
rel_p_G_D <- parameter()
dim(rel_p_G_D) <- c(n_groups, n_strains, n_vacc_classes)
strain_rel_p_ICU_D <- parameter()
dim(strain_rel_p_ICU_D) <- n_strains
strain_rel_p_H_D <- parameter()
dim(strain_rel_p_H_D) <- n_strains
strain_rel_p_W_D <- parameter()
dim(strain_rel_p_W_D) <- n_strains
strain_rel_p_G_D <- parameter()
dim(strain_rel_p_G_D) <- n_strains
rel_p_R <- parameter()
dim(rel_p_R) <- c(n_groups, n_strains, n_vacc_classes)
rel_infectivity <- parameter()
dim(rel_infectivity) <- c(n_groups, n_strains, n_vacc_classes)

vaccine_progression_rate_base <- parameter()
dim(vaccine_progression_rate_base) <- c(n_groups, n_vacc_classes)

## Parameters of the E classes
k_E <- parameter()
dim(gamma_E) <- n_strains
gamma_E_value <- parameter()
gamma_E_time <- parameter()
n_gamma_E_time <- parameter()
dim(gamma_E_value) <- n_gamma_E_time
dim(gamma_E_time) <- n_gamma_E_time
rel_gamma_E <- parameter()
dim(rel_gamma_E) <- n_strains

## Probability of transitioning from the E to the symptomatic class,
## the rest go into the asymptomatic class
p_C_value <- parameter()
n_p_C_time <- parameter()
p_C_time <- parameter()
dim(p_C) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_C_t) <- n_groups
dim(p_C_value) <- c(n_groups, n_p_C_time)
dim(p_C_time) <- n_p_C_time

## Parameters of the I_A classes
k_A <- parameter()
dim(gamma_A) <- n_strains
gamma_A_value <- parameter()
gamma_A_time <- parameter()
n_gamma_A_time <- parameter()
dim(gamma_A_value) <- n_gamma_A_time
dim(gamma_A_time) <- n_gamma_A_time
rel_gamma_A <- parameter()
dim(rel_gamma_A) <- n_strains

## Parameters of the I_P classes
k_P <- parameter()
dim(gamma_P) <- n_strains
gamma_P_value <- parameter()
gamma_P_time <- parameter()
n_gamma_P_time <- parameter()
dim(gamma_P_value) <- n_gamma_P_time
dim(gamma_P_time) <- n_gamma_P_time
rel_gamma_P <- parameter()
dim(rel_gamma_P) <- n_strains

## Parameters of the I_C_1 classes
k_C_1 <- parameter()
dim(gamma_C_1) <- n_strains
gamma_C_1_value <- parameter()
gamma_C_1_time <- parameter()
n_gamma_C_1_time <- parameter()
dim(gamma_C_1_value) <- n_gamma_C_1_time
dim(gamma_C_1_time) <- n_gamma_C_1_time
rel_gamma_C_1 <- parameter()
dim(rel_gamma_C_1) <- n_strains

## Parameters of the I_C_2 classes
k_C_2 <- parameter()
dim(gamma_C_2) <- n_strains
gamma_C_2_value <- parameter()
gamma_C_2_time <- parameter()
n_gamma_C_2_time <- parameter()
dim(gamma_C_2_value) <- n_gamma_C_2_time
dim(gamma_C_2_time) <- n_gamma_C_2_time
rel_gamma_C_2 <- parameter()
dim(rel_gamma_C_2) <- n_strains

## Proportion of cases requiring hospitalisation
p_H_value <- parameter()
n_p_H_time <- parameter()
p_H_time <- parameter()
dim(p_H) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_H_t) <- n_groups
dim(p_H_value) <- c(n_groups, n_p_H_time)
dim(p_H_time) <- n_p_H_time

## Parameters of the G_D class
k_G_D <- parameter()
dim(gamma_G_D) <- n_strains
gamma_G_D_value <- parameter()
gamma_G_D_time <- parameter()
n_gamma_G_D_time <- parameter()
dim(gamma_G_D_value) <- n_gamma_G_D_time
dim(gamma_G_D_time) <- n_gamma_G_D_time
rel_gamma_G_D <- parameter()
dim(rel_gamma_G_D) <- n_strains
p_G_D_value <- parameter()
n_p_G_D_time <- parameter()
p_G_D_time <- parameter()
dim(p_G_D) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_G_D_t) <- n_groups
dim(p_G_D_value) <- c(n_groups, n_p_G_D_time)
dim(p_G_D_time) <- n_p_G_D_time

## Parameters of the ICU_pre classes
k_ICU_pre <- parameter()
dim(gamma_ICU_pre) <- n_strains
gamma_ICU_pre_value <- parameter()
gamma_ICU_pre_time <- parameter()
n_gamma_ICU_pre_time <- parameter()
dim(gamma_ICU_pre_value) <- n_gamma_ICU_pre_time
dim(gamma_ICU_pre_time) <- n_gamma_ICU_pre_time
rel_gamma_ICU_pre <- parameter()
dim(rel_gamma_ICU_pre) <- n_strains

## Proportion of hospital cases progressing to ICU
p_ICU_value <- parameter()
n_p_ICU_time <- parameter()
p_ICU_time <- parameter()
dim(p_ICU) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_ICU_t) <- n_groups
dim(p_ICU_value) <- c(n_groups, n_p_ICU_time)
dim(p_ICU_time) <- n_p_ICU_time

## Proportion of stepdown cases dying
p_W_D_value <- parameter()
n_p_W_D_time <- parameter()
p_W_D_time <- parameter()
dim(p_W_D) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_W_D_t) <- n_groups
dim(p_W_D_value) <- c(n_groups, n_p_W_D_time)
dim(p_W_D_time) <- n_p_W_D_time

## Parameters of the H_R classes
k_H_R <- parameter()
dim(gamma_H_R) <- n_strains
gamma_H_R_value <- parameter()
gamma_H_R_time <- parameter()
n_gamma_H_R_time <- parameter()
dim(gamma_H_R_value) <- n_gamma_H_R_time
dim(gamma_H_R_time) <- n_gamma_H_R_time
rel_gamma_H_R <- parameter()
dim(rel_gamma_H_R) <- n_strains

## Parameters of the H_D classes
k_H_D <- parameter()
dim(gamma_H_D) <- n_strains
gamma_H_D_value <- parameter()
gamma_H_D_time <- parameter()
n_gamma_H_D_time <- parameter()
dim(gamma_H_D_value) <- n_gamma_H_D_time
dim(gamma_H_D_time) <- n_gamma_H_D_time
rel_gamma_H_D <- parameter()
dim(rel_gamma_H_D) <- n_strains
p_H_D_value <- parameter()
n_p_H_D_time <- parameter()
p_H_D_time <- parameter()
dim(p_H_D) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_H_D_t) <- n_groups
dim(p_H_D_value) <- c(n_groups, n_p_H_D_time)
dim(p_H_D_time) <- n_p_H_D_time

## Parameters of the ICU_W_R classes
k_ICU_W_R <- parameter()
dim(gamma_ICU_W_R) <- n_strains
gamma_ICU_W_R_value <- parameter()
gamma_ICU_W_R_time <- parameter()
n_gamma_ICU_W_R_time <- parameter()
dim(gamma_ICU_W_R_value) <- n_gamma_ICU_W_R_time
dim(gamma_ICU_W_R_time) <- n_gamma_ICU_W_R_time
rel_gamma_ICU_W_R <- parameter()
dim(rel_gamma_ICU_W_R) <- n_strains

## Parameters of the ICU_W_D classes
k_ICU_W_D <- parameter()
dim(gamma_ICU_W_D) <- n_strains
gamma_ICU_W_D_value <- parameter()
gamma_ICU_W_D_time <- parameter()
n_gamma_ICU_W_D_time <- parameter()
dim(gamma_ICU_W_D_value) <- n_gamma_ICU_W_D_time
dim(gamma_ICU_W_D_time) <- n_gamma_ICU_W_D_time
rel_gamma_ICU_W_D <- parameter()
dim(rel_gamma_ICU_W_D) <- n_strains

## Parameters of the ICU_D classes
k_ICU_D <- parameter()
dim(gamma_ICU_D) <- n_strains
gamma_ICU_D_value <- parameter()
gamma_ICU_D_time <- parameter()
n_gamma_ICU_D_time <- parameter()
dim(gamma_ICU_D_value) <- n_gamma_ICU_D_time
dim(gamma_ICU_D_time) <- n_gamma_ICU_D_time
rel_gamma_ICU_D <- parameter()
dim(rel_gamma_ICU_D) <- n_strains
p_ICU_D_value <- parameter()
n_p_ICU_D_time <- parameter()
p_ICU_D_time <- parameter()
dim(p_ICU_D) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_ICU_D_t) <- n_groups
dim(p_ICU_D_value) <- c(n_groups, n_p_ICU_D_time)
dim(p_ICU_D_time) <- n_p_ICU_D_time

## Parameters of the W_R classes
k_W_R <- parameter()
dim(gamma_W_R) <- n_strains
gamma_W_R_value <- parameter()
gamma_W_R_time <- parameter()
n_gamma_W_R_time <- parameter()
dim(gamma_W_R_value) <- n_gamma_W_R_time
dim(gamma_W_R_time) <- n_gamma_W_R_time
rel_gamma_W_R <- parameter()
dim(rel_gamma_W_R) <- n_strains

## Parameters of the W_D classes
k_W_D <- parameter()
dim(gamma_W_D) <- n_strains
gamma_W_D_value <- parameter()
gamma_W_D_time <- parameter()
n_gamma_W_D_time <- parameter()
dim(gamma_W_D_value) <- n_gamma_W_D_time
dim(gamma_W_D_time) <- n_gamma_W_D_time
rel_gamma_W_D <- parameter()
dim(rel_gamma_W_D) <- n_strains

## Parameters of the T_sero_pre_1 classes
k_sero_pre_1 <- parameter()
gamma_sero_pre_1 <- parameter(0.1)
p_sero_pos_1 <- parameter()

## Parameters of the T_sero_pos_1 classes
k_sero_pos_1 <- parameter()
gamma_sero_pos_1 <- parameter(0.1)

## Parameters of the T_sero_pre_2 classes
k_sero_pre_2 <- parameter()
gamma_sero_pre_2 <- parameter(0.1)
p_sero_pos_2 <- parameter()

## Parameters of the T_sero_pos_2 classes
k_sero_pos_2 <- parameter()
gamma_sero_pos_2 <- parameter(0.1)

## Parameters relating to testing
gamma_U_value <- parameter()
gamma_U_time <- parameter()
n_gamma_U_time <- parameter()
dim(gamma_U_value) <- n_gamma_U_time
dim(gamma_U_time) <- n_gamma_U_time
p_star_value <- parameter()
p_star_time <- parameter()
n_p_star_time <- parameter()
dim(p_star_value) <- c(n_groups, n_p_star_time)
dim(p_star_time) <- n_p_star_time
dim(p_star) <- n_groups

## Parameters relating to PCR positivity
k_PCR_pre <- parameter()
dim(gamma_PCR_pre) <- n_strains
gamma_PCR_pre_value <- parameter()
gamma_PCR_pre_time <- parameter()
n_gamma_PCR_pre_time <- parameter()
dim(gamma_PCR_pre_value) <- n_gamma_PCR_pre_time
dim(gamma_PCR_pre_time) <- n_gamma_PCR_pre_time
rel_gamma_PCR_pre <- parameter()
dim(rel_gamma_PCR_pre) <- n_strains

k_PCR_pos <- parameter()
dim(gamma_PCR_pos) <- n_strains
gamma_PCR_pos_value <- parameter()
gamma_PCR_pos_time <- parameter()
n_gamma_PCR_pos_time <- parameter()
dim(gamma_PCR_pos_value) <- n_gamma_PCR_pos_time
dim(gamma_PCR_pos_time) <- n_gamma_PCR_pos_time
rel_gamma_PCR_pos <- parameter()
dim(rel_gamma_PCR_pos) <- n_strains

## Waning of immunity
waning_rate <- parameter()
dim(waning_rate) <- n_groups

## Parameters of the age stratified transmission
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
I_A_transmission <- parameter()
I_P_transmission <- parameter()
I_C_1_transmission <- parameter()
I_C_2_transmission <- parameter()
hosp_transmission <- parameter()
ICU_transmission <- parameter()
G_D_transmission <- parameter()
strain_transmission <- parameter()
dim(strain_transmission) <- n_strains
n_strains <- parameter()
n_strains_R <- parameter()
n_real_strains <- if (n_strains == 4) 2 else 1

## Dimensions of the different "vectors" here vectors stand for
## multi-dimensional arrays

## Vectors handling the S class
dim(S) <- c(n_groups, n_vacc_classes)
dim(new_S) <- c(n_groups, n_vacc_classes)

## Vectors handling the E class
dim(E) <- c(n_groups, n_strains, k_E, n_vacc_classes)
dim(new_E) <- c(n_groups, n_strains, k_E, n_vacc_classes)

## Vectors handling the I_A class
dim(I_A) <- c(n_groups, n_strains, k_A, n_vacc_classes)
dim(new_I_A) <- c(n_groups, n_strains, k_A, n_vacc_classes)

## Vectors handling the I_P class
dim(I_P) <- c(n_groups, n_strains, k_P, n_vacc_classes)
dim(new_I_P) <- c(n_groups, n_strains, k_P, n_vacc_classes)

## Vectors handling the I_C_2 class
dim(I_C_1) <- c(n_groups, n_strains, k_C_1, n_vacc_classes)
dim(new_I_C_1) <- c(n_groups, n_strains, k_C_1, n_vacc_classes)
dim(n_I_C_1_progress) <- c(n_groups, n_strains, k_C_1, n_vacc_classes)

## Vectors handling the I_C_2 class
dim(I_C_2) <- c(n_groups, n_strains, k_C_2, n_vacc_classes)
dim(new_I_C_2) <- c(n_groups, n_strains, k_C_2, n_vacc_classes)
dim(n_I_C_2_progress) <- c(n_groups, n_strains, k_C_2, n_vacc_classes)


## Vectors handling the G_D class
dim(G_D) <- c(n_groups, n_strains, k_G_D, n_vacc_classes)
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
dim(R) <- c(n_groups, n_strains_R, n_vacc_classes)
dim(new_R) <- c(n_groups, n_strains_R, n_vacc_classes)

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
dim(D_hosp) <- c(n_groups, n_vacc_classes)
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
p_R_value <- parameter()
n_p_R_time <- parameter()
p_R_time <- parameter()
dim(p_R) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_R_t) <- n_groups
dim(p_R_value) <- c(n_groups, n_p_R_time)
dim(p_R_time) <- n_p_R_time

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
dim(p_R_next_vacc_class) <- c(n_groups, n_strains_R, n_vacc_classes)
dim(n_R_next_vacc_class) <- c(n_groups, n_strains_R, n_vacc_classes)
dim(p_R_vacc_skip) <- c(n_groups, n_strains_R, n_vacc_classes)
dim(n_R_vacc_skip) <- c(n_groups, n_strains_R, n_vacc_classes)
dim(n_R_progress) <- c(n_groups, n_strains_R, n_vacc_classes)

dim(n_R_tmp) <- c(n_groups, n_strains_R, n_vacc_classes)
dim(n_RS) <- c(n_groups, n_strains_R, n_vacc_classes)
dim(p_RS) <- c(n_groups, n_strains_R, n_vacc_classes)
dim(n_RE) <- c(n_groups, n_strains, n_vacc_classes)
dim(p_R_progress) <- c(n_groups, n_strains_R, n_vacc_classes)
dim(rate_R_progress) <- c(n_groups, n_strains_R, n_vacc_classes)

dim(cross_immunity) <- n_real_strains
cross_immunity <- parameter()

## Total population
initial(N_tot[]) <- 0
update(N_tot[]) <- sum(S[i, ]) + sum(R[i, , ]) + sum(D_hosp[i, ]) +
  sum(E[i, , , ]) + sum(I_A[i, , , ]) + sum(I_P[i, , , ]) +
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
delta_D_hosp_0_49_tot <- sum(delta_D_hosp[1:10]) +
  (if (has_carehomes == 1) delta_D_hosp[18] * 3 / 8 else 0)
update(D_hosp_0_49_tot) <- D_hosp_0_49_tot + delta_D_hosp_0_49_tot

initial(D_hosp_50_54_tot) <- 0
delta_D_hosp_50_54_tot <- delta_D_hosp[11] +
  (if (has_carehomes == 1) delta_D_hosp[18] * 1 / 8 else 0)
update(D_hosp_50_54_tot) <- D_hosp_50_54_tot + delta_D_hosp_50_54_tot

initial(D_hosp_55_59_tot) <- 0
delta_D_hosp_55_59_tot <- delta_D_hosp[12] +
  (if (has_carehomes == 1) delta_D_hosp[18] * 2 / 8 else 0)
update(D_hosp_55_59_tot) <- D_hosp_55_59_tot + delta_D_hosp_55_59_tot

initial(D_hosp_60_64_tot) <- 0
delta_D_hosp_60_64_tot <- delta_D_hosp[13] +
  (if (has_carehomes == 1) delta_D_hosp[18] * 2 / 8 else 0)
update(D_hosp_60_64_tot) <- D_hosp_60_64_tot + delta_D_hosp_60_64_tot

initial(D_hosp_65_69_tot) <- 0
delta_D_hosp_65_69_tot <- delta_D_hosp[14] +
  (if (has_carehomes == 1) delta_D_hosp[19] * 0.05 else 0)
update(D_hosp_65_69_tot) <- D_hosp_65_69_tot + delta_D_hosp_65_69_tot

initial(D_hosp_70_74_tot) <- 0
delta_D_hosp_70_74_tot <- delta_D_hosp[15] +
  (if (has_carehomes == 1) delta_D_hosp[19] * 0.05 else 0)
update(D_hosp_70_74_tot) <- D_hosp_70_74_tot + delta_D_hosp_70_74_tot

initial(D_hosp_75_79_tot) <- 0
delta_D_hosp_75_79_tot <- delta_D_hosp[16] +
  (if (has_carehomes == 1) delta_D_hosp[19] * 0.15 else 0)
update(D_hosp_75_79_tot) <- D_hosp_75_79_tot + delta_D_hosp_75_79_tot

initial(D_hosp_80_plus_tot) <- 0
delta_D_hosp_80_plus_tot <- delta_D_hosp[17] +
  (if (has_carehomes == 1) delta_D_hosp[19] * 0.75 else 0)
update(D_hosp_80_plus_tot) <- D_hosp_80_plus_tot + delta_D_hosp_80_plus_tot

## community deaths are non-hospital deaths in groups 1 to 18
initial(D_comm_tot) <- 0
delta_D_comm_tot <- sum(delta_D_non_hosp[1:18])
update(D_comm_tot) <- D_comm_tot + delta_D_comm_tot

initial(D_comm_inc, zero_every = 1) <- 0
update(D_comm_inc) <- D_comm_inc + delta_D_comm_tot

initial(D_comm_0_49_inc, zero_every = 1) <- 0
delta_D_comm_0_49 <- sum(delta_D_non_hosp[1:10]) +
  (if (has_carehomes == 1) delta_D_non_hosp[18] * 3 / 8 else 0)
update(D_comm_0_49_inc) <- D_comm_0_49_inc + delta_D_comm_0_49

initial(D_comm_50_54_inc, zero_every = 1) <- 0
delta_D_comm_50_54 <- delta_D_non_hosp[11] +
  (if (has_carehomes == 1) delta_D_non_hosp[18] * 1 / 8 else 0)
update(D_comm_50_54_inc) <- D_comm_50_54_inc + delta_D_comm_50_54

initial(D_comm_55_59_inc, zero_every = 1) <- 0
delta_D_comm_55_59 <- delta_D_non_hosp[12] +
  (if (has_carehomes == 1) delta_D_non_hosp[18] * 2 / 8 else 0)
update(D_comm_55_59_inc) <- D_comm_55_59_inc + delta_D_comm_55_59

initial(D_comm_60_64_inc, zero_every = 1) <- 0
delta_D_comm_60_64 <- delta_D_non_hosp[13] +
  (if (has_carehomes == 1) delta_D_non_hosp[18] * 2 / 8 else 0)
update(D_comm_60_64_inc) <- D_comm_60_64_inc + delta_D_comm_60_64

initial(D_comm_65_69_inc, zero_every = 1) <- 0
delta_D_comm_65_69 <- delta_D_non_hosp[14]
update(D_comm_65_69_inc) <- D_comm_65_69_inc + delta_D_comm_65_69

initial(D_comm_70_74_inc, zero_every = 1) <- 0
delta_D_comm_70_74 <- delta_D_non_hosp[15]
update(D_comm_70_74_inc) <- D_comm_70_74_inc + delta_D_comm_70_74

initial(D_comm_75_79_inc, zero_every = 1) <- 0
delta_D_comm_75_79 <- delta_D_non_hosp[16]
update(D_comm_75_79_inc) <- D_comm_75_79_inc + delta_D_comm_75_79

initial(D_comm_80_plus_inc, zero_every = 1) <- 0
delta_D_comm_80_plus <- delta_D_non_hosp[17]
update(D_comm_80_plus_inc) <- D_comm_80_plus_inc + delta_D_comm_80_plus


## carehome deaths are non-hospital deaths in group 19
initial(D_carehomes_tot) <- 0
delta_D_carehomes_tot <- if (has_carehomes == 1) delta_D_non_hosp[19] else 0
update(D_carehomes_tot) <- D_carehomes_tot + delta_D_carehomes_tot

initial(D_carehomes_inc, zero_every = 1) <- 0
update(D_carehomes_inc) <- D_carehomes_inc + delta_D_carehomes_tot

initial(D_tot) <- 0
delta_D_tot <- delta_D_hosp_tot + delta_D_comm_tot + delta_D_carehomes_tot
update(D_tot) <- D_tot + delta_D_tot

initial(D_inc, zero_every = 1) <- 0
update(D_inc) <- D_inc + delta_D_tot

## Incident deaths in hospital overall and then by age
initial(D_hosp_inc, zero_every = 1) <- 0
update(D_hosp_inc) <- D_hosp_inc + delta_D_hosp_tot

initial(D_hosp_0_49_inc, zero_every = 1) <- 0
update(D_hosp_0_49_inc) <- D_hosp_0_49_inc + delta_D_hosp_0_49_tot

initial(D_hosp_50_54_inc, zero_every = 1) <- 0
update(D_hosp_50_54_inc) <- D_hosp_50_54_inc + delta_D_hosp_50_54_tot

initial(D_hosp_55_59_inc, zero_every = 1) <- 0
update(D_hosp_55_59_inc) <- D_hosp_55_59_inc + delta_D_hosp_55_59_tot

initial(D_hosp_60_64_inc, zero_every = 1) <- 0
update(D_hosp_60_64_inc) <- D_hosp_60_64_inc + delta_D_hosp_60_64_tot

initial(D_hosp_65_69_inc, zero_every = 1) <- 0
update(D_hosp_65_69_inc) <- D_hosp_65_69_inc + delta_D_hosp_65_69_tot

initial(D_hosp_70_74_inc, zero_every = 1) <- 0
update(D_hosp_70_74_inc) <- D_hosp_70_74_inc + delta_D_hosp_70_74_tot

initial(D_hosp_75_79_inc, zero_every = 1) <- 0
update(D_hosp_75_79_inc) <- D_hosp_75_79_inc + delta_D_hosp_75_79_tot

initial(D_hosp_80_plus_inc, zero_every = 1) <- 0
update(D_hosp_80_plus_inc) <- D_hosp_80_plus_inc + delta_D_hosp_80_plus_tot


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
new_sympt_cases_25_49 <- sum(n_EI_P[6:10, , ]) +
  (if (has_carehomes == 1) (sum(n_EI_P[18, , ]) / 8) * 5 else 0)
update(cum_sympt_cases_25_49) <- cum_sympt_cases_25_49 +
  new_sympt_cases_25_49

initial(cum_sympt_cases_50_64) <- 0
new_sympt_cases_50_64 <- sum(n_EI_P[11:13, , ]) +
  (if (has_carehomes == 1) (sum(n_EI_P[18, , ]) / 8) * 3 else 0)
update(cum_sympt_cases_50_64) <- cum_sympt_cases_50_64 +
  new_sympt_cases_50_64

## assume CHR [19] are 1/4 aged 65-79 and 3/4 80 plus
initial(cum_sympt_cases_65_79) <- 0
new_sympt_cases_65_79 <- sum(n_EI_P[14:16, , ]) +
  (if (has_carehomes == 1) (sum(n_EI_P[19, , ]) * 0.25) else 0)
update(cum_sympt_cases_65_79) <- cum_sympt_cases_65_79 +
  new_sympt_cases_65_79

initial(cum_sympt_cases_80_plus) <- 0
new_sympt_cases_80_plus <- sum(n_EI_P[17, , ]) +
  (if (has_carehomes == 1) (sum(n_EI_P[19, , ]) * 0.75) else 0)
update(cum_sympt_cases_80_plus) <- cum_sympt_cases_80_plus +
  new_sympt_cases_80_plus

## And incidence:
initial(sympt_cases_inc, zero_every = 1) <- 0
update(sympt_cases_inc) <- sympt_cases_inc + new_sympt_cases

initial(sympt_cases_non_variant_inc, zero_every = 1) <- 0
update(sympt_cases_non_variant_inc) <-
  sympt_cases_non_variant_inc + new_sympt_cases_non_variant

initial(sympt_cases_over25_inc, zero_every = 1) <- 0
update(sympt_cases_over25_inc) <-
  sympt_cases_over25_inc + new_sympt_cases_over25

initial(sympt_cases_non_variant_over25_inc, zero_every = 1) <- 0
update(sympt_cases_non_variant_over25_inc) <-
  sympt_cases_non_variant_over25_inc + new_sympt_cases_non_variant_over25

initial(sympt_cases_under15_inc, zero_every = 1) <- 0
update(sympt_cases_under15_inc) <-
  sympt_cases_under15_inc + new_sympt_cases_under15

initial(sympt_cases_15_24_inc, zero_every = 1) <- 0
update(sympt_cases_15_24_inc) <- sympt_cases_15_24_inc + new_sympt_cases_15_24

initial(sympt_cases_25_49_inc, zero_every = 1) <- 0
update(sympt_cases_25_49_inc) <- sympt_cases_25_49_inc + new_sympt_cases_25_49

initial(sympt_cases_50_64_inc, zero_every = 1) <- 0
update(sympt_cases_50_64_inc) <- sympt_cases_50_64_inc + new_sympt_cases_50_64

initial(sympt_cases_65_79_inc, zero_every = 1) <- 0
update(sympt_cases_65_79_inc) <- sympt_cases_65_79_inc + new_sympt_cases_65_79

initial(sympt_cases_80_plus_inc, zero_every = 1) <- 0
update(sympt_cases_80_plus_inc) <-
  sympt_cases_80_plus_inc + new_sympt_cases_80_plus

## For ONS we exclude the 0-1 (40% of 1) and CHR (19) groups
initial(ons_positives) <- 0
update(ons_positives) <- sum(new_T_PCR_pos[1, , , ]) * 3 / 5 +
  (if (has_carehomes == 1) sum(new_T_PCR_pos[2:18, , , ]) else
    sum(new_T_PCR_pos[2:17, , , ]))

## For REACT we exclude the 0-4 (1) and CHR (19) groups
initial(react_positives) <- 0
update(react_positives) <-
  (if (has_carehomes == 1) sum(new_T_PCR_pos[2:18, , , ]) else
    sum(new_T_PCR_pos[2:17, , , ]))

initial(react_5_24_positives) <- 0
update(react_5_24_positives) <- sum(new_T_PCR_pos[2:5, , , ])


initial(react_25_34_positives) <- 0
update(react_25_34_positives) <- sum(new_T_PCR_pos[6:7, , , ]) +
  (if (has_carehomes == 1) sum(new_T_PCR_pos[18, , , ]) * 2 / 8 else 0)

initial(react_35_44_positives) <- 0
update(react_35_44_positives) <- sum(new_T_PCR_pos[8:9, , , ]) +
  (if (has_carehomes == 1) sum(new_T_PCR_pos[18, , , ]) * 2 / 8 else 0)

initial(react_45_54_positives) <- 0
update(react_45_54_positives) <- sum(new_T_PCR_pos[10:11, , , ]) +
  (if (has_carehomes == 1) sum(new_T_PCR_pos[18, , , ]) * 2 / 8 else 0)

initial(react_55_64_positives) <- 0
update(react_55_64_positives) <- sum(new_T_PCR_pos[12:13, , , ]) +
  (if (has_carehomes == 1) sum(new_T_PCR_pos[18, , , ]) * 2 / 8 else 0)

initial(react_65_plus_positives) <- 0
update(react_65_plus_positives) <- sum(new_T_PCR_pos[14:17, , , ])


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
## to strains. Those in R1, R4 and R5 will be (partially)
## susceptible to strain 2, those in R5 will also be (partially)
## susceptible to strain 1)
dim(eff_sus_S) <- c(n_groups, n_real_strains, n_vacc_classes)
dim(eff_sus_R) <- c(n_groups, n_real_strains, n_vacc_classes)
eff_sus_S[, , ] <- new_S[i, k] * rel_susceptibility[i, j, k]
eff_sus_R[, , ] <- if (n_real_strains == 1) 0 else
  (new_R[i, 5, k] + if (j == 2) new_R[i, 1, k] + new_R[i, 4, k] else 0) *
  (1 - cross_immunity[3 - j]) * rel_susceptibility[i, j, k]

initial(effective_susceptible[]) <- 0
update(effective_susceptible[]) <- sum(eff_sus_S[, i, ]) + sum(eff_sus_R[, i, ])
dim(effective_susceptible) <- n_real_strains

## Calculate the total number of susceptibles, and recovered by strain
initial(susceptible) <- 0
update(susceptible) <- sum(new_S)
initial(recovered[]) <- 0
update(recovered[]) <- sum(new_R[, i, ])
dim(recovered) <- n_strains_R


## Calculate the (weighted) number of individuals protected against infection
## to each strain in S and R
initial(protected_S_vaccinated[]) <- 0
initial(protected_R_unvaccinated[]) <- 0
initial(protected_R_vaccinated[]) <- 0
update(protected_S_vaccinated[]) <- sum(new_S) - sum(eff_sus_S[, i, ])
update(protected_R_unvaccinated[]) <- sum(new_R[, , 1]) - sum(eff_sus_R[, i, 1])
update(protected_R_vaccinated[]) <- sum(new_R) - sum(new_R[, , 1]) -
  (sum(eff_sus_R[, i, ]) - sum(eff_sus_R[, i, 1]))
dim(protected_S_vaccinated) <- n_real_strains
dim(protected_R_unvaccinated) <- n_real_strains
dim(protected_R_vaccinated) <- n_real_strains

## Vaccination engine
n_doses <- parameter()
index_dose <- parameter(type = "integer")
dim(index_dose) <- n_doses

index_dose_inverse <- parameter(type = "integer")
dim(index_dose_inverse) <- n_vacc_classes

vaccine_dose_value <- parameter()
dim(vaccine_dose_value) <- parameter(rank = 3)
vaccine_dose_time <- parameter()
dim(vaccine_dose_time) <- parameter(rank = 1)

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
vaccine_dose_t <- interpolate(vaccine_dose_time, vaccine_dose_value, "constant")
dim(vaccine_dose_t) <- c(n_groups, n_doses)
total_attempted_doses[, ] <- vaccine_missed_doses[i, j] + vaccine_dose_t[i, j]
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

vaccine_catchup_fraction <- parameter(0)


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

vacc_skip_to <- parameter(type = "integer")
dim(vacc_skip_to) <- n_vacc_classes
vacc_skip_from <- parameter(type = "integer")
dim(vacc_skip_from) <- n_vacc_classes
vacc_skip_progression_rate_base <- parameter()
dim(vacc_skip_progression_rate_base) <- n_vacc_classes
vacc_skip_dose <- parameter(type = "integer")
dim(vacc_skip_dose) <- n_doses
vacc_skip_dose_inverse <- parameter(type = "integer")
dim(vacc_skip_dose_inverse) <- n_vacc_classes
vacc_skip_dose_weight <- parameter()
dim(vacc_skip_dose_weight) <- n_doses
vacc_skipped <- parameter(type = "integer")
dim(vacc_skipped) <- n_vacc_classes

## Severity outputs by age - vacc class - infection class
dim(IHR_disag) <- c(n_groups, n_strains, n_vacc_classes)
dim(IHR_disag_weighted_inc) <- c(n_groups, n_strains, n_vacc_classes)
dim(new_IHR_disag_weighted_inc) <- c(n_groups, n_strains, n_vacc_classes)
dim(HFR_disag) <- c(n_groups, n_strains, n_vacc_classes)
dim(HFR_disag_weighted_inc) <- c(n_groups, n_strains, n_vacc_classes)
dim(new_HFR_disag_weighted_inc) <- c(n_groups, n_strains, n_vacc_classes)
dim(IFR_disag) <- c(n_groups, n_strains, n_vacc_classes)
dim(IFR_disag_weighted_inc) <- c(n_groups, n_strains, n_vacc_classes)
dim(new_IFR_disag_weighted_inc) <- c(n_groups, n_strains, n_vacc_classes)

IHR_disag[, , ] <- p_C[i, j, k] * p_H[i, j, k] * (1 - p_G_D[i, j, k])
new_IHR_disag_weighted_inc[, , ] <- IHR_disag_weighted_inc[i, j, k] +
  IHR_disag[i, j, k] * delta_infections[i, j, k]
initial(IHR_disag_weighted_inc[, , ], zero_every = 1) <- 0
update(IHR_disag_weighted_inc[, , ]) <- new_IHR_disag_weighted_inc[i, j, k]
initial(ihr) <- 0
update(ihr) <- sum(new_IHR_disag_weighted_inc) / new_infections_inc

HFR_disag[, , ] <- (1 - p_ICU[i, j, k]) * p_H_D[i, j, k] +
  p_ICU[i, j, k] * p_ICU_D[i, j, k] +
  p_ICU[i, j, k] * (1 - p_ICU_D[i, j, k]) * p_W_D[i, j, k]
new_HFR_disag_weighted_inc[, , ] <- HFR_disag_weighted_inc[i, j, k] +
  HFR_disag[i, j, k] * n_I_C_2_to_hosp[i, j, k]
initial(HFR_disag_weighted_inc[, , ], zero_every = 1) <- 0
update(HFR_disag_weighted_inc[, , ]) <- new_HFR_disag_weighted_inc[i, j, k]
initial(hfr) <- 0
update(hfr) <- sum(new_HFR_disag_weighted_inc) / new_hospitalisations_inc


IFR_disag[, , ] <- IHR_disag[i, j, k] * HFR_disag[i, j, k] +
  p_C[i, j, k] * p_H[i, j, k] * p_G_D[i, j, k]
new_IFR_disag_weighted_inc[, , ] <- IFR_disag_weighted_inc[i, j, k] +
  IFR_disag[i, j, k] * delta_infections[i, j, k]
initial(IFR_disag_weighted_inc[, , ], zero_every = 1) <- 0
update(IFR_disag_weighted_inc[, , ]) <- new_IFR_disag_weighted_inc[i, j, k]
initial(ifr) <- 0
update(ifr) <- sum(new_IFR_disag_weighted_inc) / new_infections_inc

## By strain
initial(ifr_strain[]) <- 0
update(ifr_strain[]) <- if (n_real_strains == 1)
  sum(new_IFR_disag_weighted_inc[, 1, ]) / new_infections_inc_strain[1] else
    (sum(new_IFR_disag_weighted_inc[, i, ]) +
       sum(new_IFR_disag_weighted_inc[, 5 - i, ])) /
  (new_infections_inc_strain[i] + new_infections_inc_strain[5 - i])
dim(ifr_strain) <- n_real_strains

initial(ihr_strain[]) <- 0
update(ihr_strain[]) <- if (n_real_strains == 1)
  sum(new_IHR_disag_weighted_inc[, 1, ]) / new_infections_inc_strain[1] else
    (sum(new_IHR_disag_weighted_inc[, i, ]) +
       sum(new_IHR_disag_weighted_inc[, 5 - i, ])) /
  (new_infections_inc_strain[i] + new_infections_inc_strain[5 - i])
dim(ihr_strain) <- n_real_strains

initial(hfr_strain[]) <- 0
update(hfr_strain[]) <- if (n_real_strains == 1)
  sum(new_HFR_disag_weighted_inc[, 1, ]) /
  new_hospitalisations_inc_strain[1] else
    (sum(new_HFR_disag_weighted_inc[, i, ]) +
       sum(new_HFR_disag_weighted_inc[, 5 - i, ])) /
  (new_hospitalisations_inc_strain[i] + new_hospitalisations_inc_strain[5 - i])
dim(hfr_strain) <- n_real_strains

## By age
dim(ifr_age) <- n_groups
initial(ifr_age[]) <- 0
update(ifr_age[]) <- sum(new_IFR_disag_weighted_inc[i, , ]) /
  new_infections_inc_age[i]

dim(ihr_age) <- n_groups
initial(ihr_age[]) <- 0
update(ihr_age[]) <- sum(new_IHR_disag_weighted_inc[i, , ]) /
  new_infections_inc_age[i]

dim(hfr_age) <- n_groups
initial(hfr_age[]) <- 0
update(hfr_age[]) <- sum(new_HFR_disag_weighted_inc[i, , ]) /
  new_hospitalisations_inc_age[i]



## COMPARE
exp_noise <- parameter()

## Hospital bed prevalences
icu <- data()
phi_ICU <- parameter()
kappa_ICU <- parameter()
icu_with_noise <- phi_ICU * ICU_tot + Exponential(exp_noise)
icu ~ NegativeBinomial(kappa_ICU, mu = icu_with_noise)

general <- data()
phi_general <- parameter()
kappa_general <- parameter()
general_with_noise <- phi_general * general_tot + Exponential(exp_noise)
general ~ NegativeBinomial(kappa_general, mu = general_with_noise)

hosp <- data()
phi_hosp <- parameter()
kappa_hosp <- parameter()
hosp_with_noise <- phi_hosp * hosp_tot + Exponential(exp_noise)
hosp ~ NegativeBinomial(kappa_hosp, mu = hosp_with_noise)

## Hospital deaths
phi_death_hosp <- parameter()
kappa_death_hosp <- parameter()

deaths_hosp <- data()
D_hosp_inc_with_noise <- phi_death_hosp * D_hosp_inc + Exponential(exp_noise)
deaths_hosp ~ NegativeBinomial(kappa_death_hosp, mu = D_hosp_inc_with_noise)

deaths_hosp_0_49 <- data()
D_hosp_0_49_inc_with_noise <-
  phi_death_hosp * D_hosp_0_49_inc + Exponential(exp_noise)
deaths_hosp_0_49 ~
  NegativeBinomial(kappa_death_hosp, mu = D_hosp_0_49_inc_with_noise)

deaths_hosp_50_54 <- data()
D_hosp_50_54_inc_with_noise <-
  phi_death_hosp * D_hosp_50_54_inc + Exponential(exp_noise)
deaths_hosp_50_54 ~
  NegativeBinomial(kappa_death_hosp, mu = D_hosp_50_54_inc_with_noise)

deaths_hosp_55_59 <- data()
D_hosp_55_59_inc_with_noise <-
  phi_death_hosp * D_hosp_55_59_inc + Exponential(exp_noise)
deaths_hosp_55_59 ~
  NegativeBinomial(kappa_death_hosp, mu = D_hosp_55_59_inc_with_noise)

deaths_hosp_60_64 <- data()
D_hosp_60_64_inc_with_noise <-
  phi_death_hosp * D_hosp_60_64_inc + Exponential(exp_noise)
deaths_hosp_60_64 ~
  NegativeBinomial(kappa_death_hosp, mu = D_hosp_60_64_inc_with_noise)

deaths_hosp_65_69 <- data()
D_hosp_65_69_inc_with_noise <-
  phi_death_hosp * D_hosp_65_69_inc + Exponential(exp_noise)
deaths_hosp_65_69 ~
  NegativeBinomial(kappa_death_hosp, mu = D_hosp_65_69_inc_with_noise)

deaths_hosp_70_74 <- data()
D_hosp_70_74_inc_with_noise <-
  phi_death_hosp * D_hosp_70_74_inc + Exponential(exp_noise)
deaths_hosp_70_74 ~
  NegativeBinomial(kappa_death_hosp, mu = D_hosp_70_74_inc_with_noise)

deaths_hosp_75_79 <- data()
D_hosp_75_79_inc_with_noise <-
  phi_death_hosp * D_hosp_75_79_inc + Exponential(exp_noise)
deaths_hosp_75_79 ~
  NegativeBinomial(kappa_death_hosp, mu = D_hosp_75_79_inc_with_noise)

deaths_hosp_80_plus <- data()
D_hosp_80_plus_inc_with_noise <-
  phi_death_hosp * D_hosp_80_plus_inc + Exponential(exp_noise)
deaths_hosp_80_plus ~
  NegativeBinomial(kappa_death_hosp, mu = D_hosp_80_plus_inc_with_noise)

## Community deaths
phi_death_comm <- parameter()
kappa_death_comm <- parameter()

deaths_comm <- data()
D_comm_inc_with_noise <- phi_death_comm * D_comm_inc + Exponential(exp_noise)
deaths_comm ~ NegativeBinomial(kappa_death_comm, mu = D_comm_inc_with_noise)

deaths_comm_0_49 <- data()
D_comm_0_49_inc_with_noise <-
  phi_death_comm * D_comm_0_49_inc + Exponential(exp_noise)
deaths_comm_0_49 ~
  NegativeBinomial(kappa_death_comm, mu = D_comm_0_49_inc_with_noise)

deaths_comm_50_54 <- data()
D_comm_50_54_inc_with_noise <-
  phi_death_comm * D_comm_50_54_inc + Exponential(exp_noise)
deaths_comm_50_54 ~
  NegativeBinomial(kappa_death_comm, mu = D_comm_50_54_inc_with_noise)

deaths_comm_55_59 <- data()
D_comm_55_59_inc_with_noise <-
  phi_death_comm * D_comm_55_59_inc + Exponential(exp_noise)
deaths_comm_55_59 ~
  NegativeBinomial(kappa_death_comm, mu = D_comm_55_59_inc_with_noise)

deaths_comm_60_64 <- data()
D_comm_60_64_inc_with_noise <-
  phi_death_comm * D_comm_60_64_inc + Exponential(exp_noise)
deaths_comm_60_64 ~
  NegativeBinomial(kappa_death_comm, mu = D_comm_60_64_inc_with_noise)

deaths_comm_65_69 <- data()
D_comm_65_69_inc_with_noise <-
  phi_death_comm * D_comm_65_69_inc + Exponential(exp_noise)
deaths_comm_65_69 ~
  NegativeBinomial(kappa_death_comm, mu = D_comm_65_69_inc_with_noise)

deaths_comm_70_74 <- data()
D_comm_70_74_inc_with_noise <-
  phi_death_comm * D_comm_70_74_inc + Exponential(exp_noise)
deaths_comm_70_74 ~
  NegativeBinomial(kappa_death_comm, mu = D_comm_70_74_inc_with_noise)

deaths_comm_75_79 <- data()
D_comm_75_79_inc_with_noise <-
  phi_death_comm * D_comm_75_79_inc + Exponential(exp_noise)
deaths_comm_75_79 ~
  NegativeBinomial(kappa_death_comm, mu = D_comm_75_79_inc_with_noise)

deaths_comm_80_plus <- data()
D_comm_80_plus_inc_with_noise <-
  phi_death_comm * D_comm_80_plus_inc + Exponential(exp_noise)
deaths_comm_80_plus ~
  NegativeBinomial(kappa_death_comm, mu = D_comm_80_plus_inc_with_noise)

## Other death datastreams
deaths_carehomes <- data()
phi_death_carehomes <- parameter()
kappa_death_carehomes <- parameter()
D_carehomes_inc_with_noise <- phi_death_carehomes * D_carehomes_inc +
  Exponential(exp_noise)
deaths_carehomes ~
  NegativeBinomial(kappa_death_carehomes, mu = D_carehomes_inc_with_noise)

deaths_non_hosp <- data()
kappa_death_non_hosp <- parameter()
D_non_hosp_inc_with_noise <- phi_death_carehomes * D_carehomes_inc +
  phi_death_comm * D_comm_inc + Exponential(exp_noise)
deaths_non_hosp ~
  NegativeBinomial(kappa_death_non_hosp, mu = D_non_hosp_inc_with_noise)

deaths <- data()
kappa_death <- parameter()
D_inc_with_noise <- phi_death_carehomes * D_carehomes_inc +
  phi_death_comm * D_comm_inc + phi_death_hosp * D_hosp_inc +
  Exponential(exp_noise)
deaths ~ NegativeBinomial(kappa_death, mu = D_inc_with_noise)

## Hospital admissions
admitted <- data()
phi_admitted <- parameter()
kappa_admitted <- parameter()
admitted_with_noise <- phi_admitted * admit_conf_inc + Exponential(exp_noise)
admitted ~ NegativeBinomial(kappa_admitted, mu = admitted_with_noise)

diagnoses <- data()
phi_diagnoses <- parameter()
kappa_diagnoses <- parameter()
diagnoses_with_noise <- phi_diagnoses * new_conf_inc + Exponential(exp_noise)
diagnoses ~ NegativeBinomial(kappa_diagnoses, mu = diagnoses_with_noise)

all_admission <- data()
phi_all_admission <- parameter()
kappa_all_admission <- parameter()
all_admission_with_noise <-
  phi_all_admission * (admit_conf_inc + new_conf_inc) + Exponential(exp_noise)
all_admission ~
  NegativeBinomial(kappa_all_admission, mu = all_admission_with_noise)

all_admission_0_9 <- data()
all_admission_0_9_with_noise <-
  phi_all_admission * all_admission_0_9_conf_inc + Exponential(exp_noise)
all_admission_0_9 ~
  NegativeBinomial(kappa_all_admission, mu = all_admission_0_9_with_noise)

all_admission_10_19 <- data()
all_admission_10_19_with_noise <-
  phi_all_admission * all_admission_10_19_conf_inc + Exponential(exp_noise)
all_admission_10_19 ~
  NegativeBinomial(kappa_all_admission, mu = all_admission_10_19_with_noise)

all_admission_20_29 <- data()
all_admission_20_29_with_noise <-
  phi_all_admission * all_admission_20_29_conf_inc + Exponential(exp_noise)
all_admission_20_29 ~
  NegativeBinomial(kappa_all_admission, mu = all_admission_20_29_with_noise)

all_admission_30_39 <- data()
all_admission_30_39_with_noise <-
  phi_all_admission * all_admission_30_39_conf_inc + Exponential(exp_noise)
all_admission_30_39 ~
  NegativeBinomial(kappa_all_admission, mu = all_admission_30_39_with_noise)

all_admission_40_49 <- data()
all_admission_40_49_with_noise <-
  phi_all_admission * all_admission_40_49_conf_inc + Exponential(exp_noise)
all_admission_40_49 ~
  NegativeBinomial(kappa_all_admission, mu = all_admission_40_49_with_noise)

all_admission_50_59 <- data()
all_admission_50_59_with_noise <-
  phi_all_admission * all_admission_50_59_conf_inc + Exponential(exp_noise)
all_admission_50_59 ~
  NegativeBinomial(kappa_all_admission, mu = all_admission_50_59_with_noise)

all_admission_60_69 <- data()
all_admission_60_69_with_noise <-
  phi_all_admission * all_admission_60_69_conf_inc + Exponential(exp_noise)
all_admission_60_69 ~
  NegativeBinomial(kappa_all_admission, mu = all_admission_60_69_with_noise)

all_admission_70_79 <- data()
all_admission_70_79_with_noise <-
  phi_all_admission * all_admission_70_79_conf_inc + Exponential(exp_noise)
all_admission_70_79 ~
  NegativeBinomial(kappa_all_admission, mu = all_admission_70_79_with_noise)

all_admission_80_plus <- data()
all_admission_80_plus_with_noise <-
  phi_all_admission * all_admission_80_plus_conf_inc + Exponential(exp_noise)
all_admission_80_plus ~
  NegativeBinomial(kappa_all_admission, mu = all_admission_80_plus_with_noise)

## Pillar 2 positivity
pillar2_sensitivity <- parameter()
pillar2_specificity <- parameter()
rho_pillar2_tests <- parameter()

is_weekend <- (time + 3) %% 7 < 2

pillar2_under15_pos <- data()
pillar2_under15_tot <- data()
N_tot_under15 <- parameter()
p_NC_under15 <- parameter()
p_NC_weekend_under15 <- parameter()
p_NC_today_under15 <- if (is_weekend) p_NC_weekend_under15 else p_NC_under15
mod_pillar2_under15_pos <- sympt_cases_under15_inc + Exponential(exp_noise)
mod_pillar2_under15_neg <-
  p_NC_today_under15 * (N_tot_under15 - sympt_cases_under15_inc) +
  Exponential(exp_noise)
mod_pillar2_under15_prob_pos <-
  (pillar2_sensitivity * mod_pillar2_under15_pos +
     (1 - pillar2_specificity) * mod_pillar2_under15_neg) /
  (mod_pillar2_under15_pos + mod_pillar2_under15_neg)
pillar2_under15_pos ~ BetaBinomial(pillar2_under15_tot,
                                   prob = mod_pillar2_under15_prob_pos,
                                   rho = rho_pillar2_tests)

pillar2_15_24_pos <- data()
pillar2_15_24_tot <- data()
N_tot_15_24 <- parameter()
p_NC_15_24 <- parameter()
p_NC_weekend_15_24 <- parameter()
p_NC_today_15_24 <- if (is_weekend) p_NC_weekend_15_24 else p_NC_15_24
mod_pillar2_15_24_pos <- sympt_cases_15_24_inc + Exponential(exp_noise)
mod_pillar2_15_24_neg <-
  p_NC_today_15_24 * (N_tot_15_24 - sympt_cases_15_24_inc) +
  Exponential(exp_noise)
mod_pillar2_15_24_prob_pos <-
  (pillar2_sensitivity * mod_pillar2_15_24_pos +
     (1 - pillar2_specificity) * mod_pillar2_15_24_neg) /
  (mod_pillar2_15_24_pos + mod_pillar2_15_24_neg)
pillar2_15_24_pos ~ BetaBinomial(pillar2_15_24_tot,
                                 prob = mod_pillar2_15_24_prob_pos,
                                 rho = rho_pillar2_tests)

pillar2_25_49_pos <- data()
pillar2_25_49_tot <- data()
N_tot_25_49 <- parameter()
p_NC_25_49 <- parameter()
p_NC_weekend_25_49 <- parameter()
p_NC_today_25_49 <- if (is_weekend) p_NC_weekend_25_49 else p_NC_25_49
mod_pillar2_25_49_pos <- sympt_cases_25_49_inc + Exponential(exp_noise)
mod_pillar2_25_49_neg <-
  p_NC_today_25_49 * (N_tot_25_49 - sympt_cases_25_49_inc) +
  Exponential(exp_noise)
mod_pillar2_25_49_prob_pos <-
  (pillar2_sensitivity * mod_pillar2_25_49_pos +
     (1 - pillar2_specificity) * mod_pillar2_25_49_neg) /
  (mod_pillar2_25_49_pos + mod_pillar2_25_49_neg)
pillar2_25_49_pos ~ BetaBinomial(pillar2_25_49_tot,
                                 prob = mod_pillar2_25_49_prob_pos,
                                 rho = rho_pillar2_tests)

pillar2_50_64_pos <- data()
pillar2_50_64_tot <- data()
N_tot_50_64 <- parameter()
p_NC_50_64 <- parameter()
p_NC_weekend_50_64 <- parameter()
p_NC_today_50_64 <- if (is_weekend) p_NC_weekend_50_64 else p_NC_50_64
mod_pillar2_50_64_pos <- sympt_cases_50_64_inc + Exponential(exp_noise)
mod_pillar2_50_64_neg <-
  p_NC_today_50_64 * (N_tot_50_64 - sympt_cases_50_64_inc) +
  Exponential(exp_noise)
mod_pillar2_50_64_prob_pos <-
  (pillar2_sensitivity * mod_pillar2_50_64_pos +
     (1 - pillar2_specificity) * mod_pillar2_50_64_neg) /
  (mod_pillar2_50_64_pos + mod_pillar2_50_64_neg)
pillar2_50_64_pos ~ BetaBinomial(pillar2_50_64_tot,
                                 prob = mod_pillar2_50_64_prob_pos,
                                 rho = rho_pillar2_tests)

pillar2_65_79_pos <- data()
pillar2_65_79_tot <- data()
N_tot_65_79 <- parameter()
p_NC_65_79 <- parameter()
p_NC_weekend_65_79 <- parameter()
p_NC_today_65_79 <- if (is_weekend) p_NC_weekend_65_79 else p_NC_65_79
mod_pillar2_65_79_pos <- sympt_cases_65_79_inc + Exponential(exp_noise)
mod_pillar2_65_79_neg <-
  p_NC_today_65_79 * (N_tot_65_79 - sympt_cases_65_79_inc) +
  Exponential(exp_noise)
mod_pillar2_65_79_prob_pos <-
  (pillar2_sensitivity * mod_pillar2_65_79_pos +
     (1 - pillar2_specificity) * mod_pillar2_65_79_neg) /
  (mod_pillar2_65_79_pos + mod_pillar2_65_79_neg)
pillar2_65_79_pos ~ BetaBinomial(pillar2_65_79_tot,
                                 prob = mod_pillar2_65_79_prob_pos,
                                 rho = rho_pillar2_tests)

pillar2_80_plus_pos <- data()
pillar2_80_plus_tot <- data()
N_tot_80_plus <- parameter()
p_NC_80_plus <- parameter()
p_NC_weekend_80_plus <- parameter()
p_NC_today_80_plus <- if (is_weekend) p_NC_weekend_80_plus else p_NC_80_plus
mod_pillar2_80_plus_pos <- sympt_cases_80_plus_inc + Exponential(exp_noise)
mod_pillar2_80_plus_neg <-
  p_NC_today_80_plus * (N_tot_80_plus - sympt_cases_80_plus_inc) +
  Exponential(exp_noise)
mod_pillar2_80_plus_prob_pos <-
  (pillar2_sensitivity * mod_pillar2_80_plus_pos +
     (1 - pillar2_specificity) * mod_pillar2_80_plus_neg) /
  (mod_pillar2_80_plus_pos + mod_pillar2_80_plus_neg)
pillar2_80_plus_pos ~ BetaBinomial(pillar2_80_plus_tot,
                                   prob = mod_pillar2_80_plus_prob_pos,
                                   rho = rho_pillar2_tests)

pillar2_over25_pos <- data()
pillar2_over25_tot <- data()
mod_pillar2_over25_pos <- sympt_cases_over25_inc + Exponential(exp_noise)
mod_pillar2_over25_neg <-
  p_NC_today_25_49 * (N_tot_25_49 - sympt_cases_25_49_inc) +
  p_NC_today_50_64 * (N_tot_50_64 - sympt_cases_50_64_inc) +
  p_NC_today_65_79 * (N_tot_65_79 - sympt_cases_65_79_inc) +
  p_NC_today_80_plus * (N_tot_80_plus - sympt_cases_80_plus_inc) +
  Exponential(exp_noise)
mod_pillar2_over25_prob_pos <-
  (pillar2_sensitivity * mod_pillar2_over25_pos +
     (1 - pillar2_specificity) * mod_pillar2_over25_neg) /
  (mod_pillar2_over25_pos + mod_pillar2_over25_neg)
pillar2_over25_pos ~ BetaBinomial(pillar2_over25_tot,
                                  prob = mod_pillar2_over25_prob_pos,
                                  rho = rho_pillar2_tests)

pillar2_pos <- data()
pillar2_tot <- data()
mod_pillar2_pos <- sympt_cases_inc + Exponential(exp_noise)
mod_pillar2_neg <-
  p_NC_today_under15 * (N_tot_under15 - sympt_cases_under15_inc) +
  p_NC_today_15_24 * (N_tot_15_24 - sympt_cases_15_24_inc) +
  p_NC_today_25_49 * (N_tot_25_49 - sympt_cases_25_49_inc) +
  p_NC_today_50_64 * (N_tot_50_64 - sympt_cases_50_64_inc) +
  p_NC_today_65_79 * (N_tot_65_79 - sympt_cases_65_79_inc) +
  p_NC_today_80_plus * (N_tot_80_plus - sympt_cases_80_plus_inc) +
  Exponential(exp_noise)
mod_pillar2_prob_pos <-
  (pillar2_sensitivity * mod_pillar2_pos +
     (1 - pillar2_specificity) * mod_pillar2_neg) /
  (mod_pillar2_pos + mod_pillar2_neg)
pillar2_pos ~ BetaBinomial(pillar2_tot, prob = mod_pillar2_prob_pos,
                           rho = rho_pillar2_tests)

## Pillar 2 cases
kappa_pillar2_cases <- parameter()

pillar2_under15_cases <- data()
phi_pillar2_cases_under15 <- parameter()
phi_pillar2_cases_weekend_under15 <- parameter()
phi_pillar2_cases_today_under15 <-
  if (is_weekend) phi_pillar2_cases_weekend_under15 else
    phi_pillar2_cases_under15
mod_pillar2_cases_under15 <-
  phi_pillar2_cases_today_under15 * sympt_cases_under15_inc +
  Exponential(exp_noise)
pillar2_under15_cases ~
  NegativeBinomial(kappa_pillar2_cases, mu = mod_pillar2_cases_under15)

pillar2_15_24_cases <- data()
phi_pillar2_cases_15_24 <- parameter()
phi_pillar2_cases_weekend_15_24 <- parameter()
phi_pillar2_cases_today_15_24 <-
  if (is_weekend) phi_pillar2_cases_weekend_15_24 else phi_pillar2_cases_15_24
mod_pillar2_cases_15_24 <-
  phi_pillar2_cases_today_15_24 * sympt_cases_15_24_inc + Exponential(exp_noise)
pillar2_15_24_cases ~
  NegativeBinomial(kappa_pillar2_cases, mu = mod_pillar2_cases_15_24)

pillar2_25_49_cases <- data()
phi_pillar2_cases_25_49 <- parameter()
phi_pillar2_cases_weekend_25_49 <- parameter()
phi_pillar2_cases_today_25_49 <-
  if (is_weekend) phi_pillar2_cases_weekend_25_49 else phi_pillar2_cases_25_49
mod_pillar2_cases_25_49 <-
  phi_pillar2_cases_today_25_49 * sympt_cases_25_49_inc + Exponential(exp_noise)
pillar2_25_49_cases ~
  NegativeBinomial(kappa_pillar2_cases, mu = mod_pillar2_cases_25_49)

pillar2_50_64_cases <- data()
phi_pillar2_cases_50_64 <- parameter()
phi_pillar2_cases_weekend_50_64 <- parameter()
phi_pillar2_cases_today_50_64 <-
  if (is_weekend) phi_pillar2_cases_weekend_50_64 else phi_pillar2_cases_50_64
mod_pillar2_cases_50_64 <-
  phi_pillar2_cases_today_50_64 * sympt_cases_50_64_inc + Exponential(exp_noise)
pillar2_50_64_cases ~
  NegativeBinomial(kappa_pillar2_cases, mu = mod_pillar2_cases_50_64)

pillar2_65_79_cases <- data()
phi_pillar2_cases_65_79 <- parameter()
phi_pillar2_cases_weekend_65_79 <- parameter()
phi_pillar2_cases_today_65_79 <-
  if (is_weekend) phi_pillar2_cases_weekend_65_79 else phi_pillar2_cases_65_79
mod_pillar2_cases_65_79 <-
  phi_pillar2_cases_today_65_79 * sympt_cases_65_79_inc + Exponential(exp_noise)
pillar2_65_79_cases ~
  NegativeBinomial(kappa_pillar2_cases, mu = mod_pillar2_cases_65_79)

pillar2_80_plus_cases <- data()
phi_pillar2_cases_80_plus <- parameter()
phi_pillar2_cases_weekend_80_plus <- parameter()
phi_pillar2_cases_today_80_plus <-
  if (is_weekend) phi_pillar2_cases_weekend_80_plus else
    phi_pillar2_cases_80_plus
mod_pillar2_cases_80_plus <-
  phi_pillar2_cases_today_80_plus * sympt_cases_80_plus_inc +
  Exponential(exp_noise)
pillar2_80_plus_cases ~
  NegativeBinomial(kappa_pillar2_cases, mu = mod_pillar2_cases_80_plus)

pillar2_over25_cases <- data()
mod_pillar2_cases_over25 <-
  phi_pillar2_cases_today_25_49 * sympt_cases_25_49_inc +
  phi_pillar2_cases_today_50_64 * sympt_cases_50_64_inc +
  phi_pillar2_cases_today_65_79 * sympt_cases_65_79_inc +
  phi_pillar2_cases_today_80_plus * sympt_cases_80_plus_inc +
  Exponential(exp_noise)
pillar2_over25_cases ~
  NegativeBinomial(kappa_pillar2_cases, mu = mod_pillar2_cases_over25)

pillar2_cases <- data()
mod_pillar2_cases <- phi_pillar2_cases_today_under15 * sympt_cases_under15_inc +
  phi_pillar2_cases_today_15_24 * sympt_cases_15_24_inc +
  phi_pillar2_cases_today_25_49 * sympt_cases_25_49_inc +
  phi_pillar2_cases_today_50_64 * sympt_cases_50_64_inc +
  phi_pillar2_cases_today_65_79 * sympt_cases_65_79_inc +
  phi_pillar2_cases_today_80_plus * sympt_cases_80_plus_inc +
  Exponential(exp_noise)
pillar2_cases ~ NegativeBinomial(kappa_pillar2_cases, mu = mod_pillar2_cases)

## Seropositivity
N_tot_15_64 <- parameter()

sero_pos_15_64_1 <- data()
sero_tot_15_64_1 <- data()
sero_sensitivity_1 <- parameter()
sero_specificity_1 <- parameter()
sero_pos_1_capped <- min(sero_pos_1, N_tot_15_64)
mod_sero_pos_1 <- sero_pos_1_capped + Exponential(exp_noise)
mod_sero_neg_1 <- N_tot_15_64 - sero_pos_1_capped + Exponential(exp_noise)
mod_sero_prob_pos_1 <-
  (sero_sensitivity_1 * mod_sero_pos_1 +
     (1 - sero_specificity_1) * mod_sero_neg_1) /
  (mod_sero_pos_1 + mod_sero_neg_1)
sero_pos_15_64_1 ~ Binomial(sero_tot_15_64_1, mod_sero_prob_pos_1)

sero_pos_15_64_2 <- data()
sero_tot_15_64_2 <- data()
sero_sensitivity_2 <- parameter()
sero_specificity_2 <- parameter()
sero_pos_2_capped <- min(sero_pos_2, N_tot_15_64)
mod_sero_pos_2 <- sero_pos_2_capped + Exponential(exp_noise)
mod_sero_neg_2 <- N_tot_15_64 - sero_pos_2_capped + Exponential(exp_noise)
mod_sero_prob_pos_2 <-
  (sero_sensitivity_2 * mod_sero_pos_2 +
     (1 - sero_specificity_2) * mod_sero_neg_2) /
  (mod_sero_pos_2 + mod_sero_neg_2)
sero_pos_15_64_2 ~ Binomial(sero_tot_15_64_2, mod_sero_prob_pos_2)

## ONS positivity
ons_pos <- data()
ons_tot <- data()
N_tot_ons <- parameter()
ons_sensitivity <- parameter()
ons_specificity <- parameter()
ons_positives_capped <- min(ons_positives, N_tot_ons)
mod_ons_pos <- ons_positives_capped + Exponential(exp_noise)
mod_ons_neg <- N_tot_ons - ons_positives_capped + Exponential(exp_noise)
mod_ons_prob_pos <-
  (ons_sensitivity * mod_ons_pos + (1 - ons_specificity) * mod_ons_neg) /
  (mod_ons_pos + mod_ons_neg)
ons_pos ~ Binomial(ons_tot, mod_ons_prob_pos)

## REACT positivity
react_sensitivity <- parameter()
react_specificity <- parameter()

react_pos <- data()
react_tot <- data()
N_tot_react <- parameter()
react_positives_capped <- min(react_positives, N_tot_react)
mod_react_pos <- react_positives_capped + Exponential(exp_noise)
mod_react_neg <- N_tot_react - react_positives_capped +
  Exponential(exp_noise)
mod_react_prob_pos <-
  (react_sensitivity * mod_react_pos +
     (1 - react_specificity) * mod_react_neg) /
  (mod_react_pos + mod_react_neg)
react_pos ~ Binomial(react_tot, mod_react_prob_pos)

react_5_24_pos <- data()
react_5_24_tot <- data()
N_5_24_react <- parameter()
react_5_24_positives_capped <- min(react_5_24_positives, N_5_24_react)
mod_react_5_24_pos <- react_5_24_positives_capped + Exponential(exp_noise)
mod_react_5_24_neg <- N_5_24_react - react_5_24_positives_capped +
  Exponential(exp_noise)
mod_react_5_24_prob_pos <-
  (react_sensitivity * mod_react_5_24_pos +
     (1 - react_specificity) * mod_react_5_24_neg) /
  (mod_react_5_24_pos + mod_react_5_24_neg)
react_5_24_pos ~ Binomial(react_5_24_tot, mod_react_5_24_prob_pos)

react_25_34_pos <- data()
react_25_34_tot <- data()
N_25_34_react <- parameter()
react_25_34_positives_capped <- min(react_25_34_positives, N_25_34_react)
mod_react_25_34_pos <- react_25_34_positives_capped + Exponential(exp_noise)
mod_react_25_34_neg <- N_25_34_react - react_25_34_positives_capped +
  Exponential(exp_noise)
mod_react_25_34_prob_pos <-
  (react_sensitivity * mod_react_25_34_pos +
     (1 - react_specificity) * mod_react_25_34_neg) /
  (mod_react_25_34_pos + mod_react_25_34_neg)
react_25_34_pos ~ Binomial(react_25_34_tot, mod_react_25_34_prob_pos)

react_35_44_pos <- data()
react_35_44_tot <- data()
N_35_44_react <- parameter()
react_35_44_positives_capped <- min(react_35_44_positives, N_35_44_react)
mod_react_35_44_pos <- react_35_44_positives_capped + Exponential(exp_noise)
mod_react_35_44_neg <- N_35_44_react - react_35_44_positives_capped +
  Exponential(exp_noise)
mod_react_35_44_prob_pos <-
  (react_sensitivity * mod_react_35_44_pos +
     (1 - react_specificity) * mod_react_35_44_neg) /
  (mod_react_35_44_pos + mod_react_35_44_neg)
react_35_44_pos ~ Binomial(react_35_44_tot, mod_react_35_44_prob_pos)

react_45_54_pos <- data()
react_45_54_tot <- data()
N_45_54_react <- parameter()
react_45_54_positives_capped <- min(react_45_54_positives, N_45_54_react)
mod_react_45_54_pos <- react_45_54_positives_capped + Exponential(exp_noise)
mod_react_45_54_neg <- N_45_54_react - react_45_54_positives_capped +
  Exponential(exp_noise)
mod_react_45_54_prob_pos <-
  (react_sensitivity * mod_react_45_54_pos +
     (1 - react_specificity) * mod_react_45_54_neg) /
  (mod_react_45_54_pos + mod_react_45_54_neg)
react_45_54_pos ~ Binomial(react_45_54_tot, mod_react_45_54_prob_pos)

react_55_64_pos <- data()
react_55_64_tot <- data()
N_55_64_react <- parameter()
react_55_64_positives_capped <- min(react_55_64_positives, N_55_64_react)
mod_react_55_64_pos <- react_55_64_positives_capped + Exponential(exp_noise)
mod_react_55_64_neg <- N_55_64_react - react_55_64_positives_capped +
  Exponential(exp_noise)
mod_react_55_64_prob_pos <-
  (react_sensitivity * mod_react_55_64_pos +
     (1 - react_specificity) * mod_react_55_64_neg) /
  (mod_react_55_64_pos + mod_react_55_64_neg)
react_55_64_pos ~ Binomial(react_55_64_tot, mod_react_55_64_prob_pos)

react_65_plus_pos <- data()
react_65_plus_tot <- data()
N_65_plus_react <- parameter()
react_65_plus_positives_capped <- min(react_65_plus_positives, N_65_plus_react)
mod_react_65_plus_pos <- react_65_plus_positives_capped + Exponential(exp_noise)
mod_react_65_plus_neg <- N_65_plus_react - react_65_plus_positives_capped +
  Exponential(exp_noise)
mod_react_65_plus_prob_pos <-
  (react_sensitivity * mod_react_65_plus_pos +
     (1 - react_specificity) * mod_react_65_plus_neg) /
  (mod_react_65_plus_pos + mod_react_65_plus_neg)
react_65_plus_pos ~ Binomial(react_65_plus_tot, mod_react_65_plus_prob_pos)

## Strains
strain_non_variant <- data()
strain_tot <- data()
mod_strain_non_variant <- sympt_cases_non_variant_inc + Exponential(exp_noise)
mod_strain_variant <- sympt_cases_inc - sympt_cases_non_variant_inc +
  Exponential(exp_noise)
mod_strain_prob_non_variant <- mod_strain_non_variant /
  (mod_strain_non_variant + mod_strain_variant)
strain_non_variant ~ Binomial(strain_tot, mod_strain_prob_non_variant)

strain_over25_non_variant <- data()
strain_over25_tot <- data()
mod_strain_over25_non_variant <-
  sympt_cases_non_variant_over25_inc + Exponential(exp_noise)
mod_strain_over25_variant <- sympt_cases_over25_inc -
  sympt_cases_non_variant_over25_inc + Exponential(exp_noise)
mod_strain_over25_prob_non_variant <- mod_strain_over25_non_variant /
  (mod_strain_over25_non_variant + mod_strain_over25_variant)
strain_over25_non_variant ~
  Binomial(strain_over25_tot, mod_strain_over25_prob_non_variant)
