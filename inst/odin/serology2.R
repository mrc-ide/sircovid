#E and R stage indexed by i,j,k with
#i for the age group
#j for the progression (not exponential latent and infectious period)
#k for the infectivity group

#Number of age classes & number of transmissibility classes
N_age <- user()
trans_classes <- user()

#definition of the time-step and output as "time"
dt <- user()
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
update(S[]) <- S[i] - n_SE[i]
update(E[,,]) <- E[i,j,k] + delta_E[i,j,k]
update(I_asympt[,,]) <- I_asympt[i,j,k] + delta_I_asympt[i,j,k]
update(I_mild[,,]) <- I_mild[i,j,k] + delta_I_mild[i,j,k]
update(I_ILI[,,]) <- I_ILI[i,j,k] + delta_I_ILI[i,j,k]
update(I_comm_D[,,]) <- I_comm_D[i,j,k] + delta_I_comm_D[i,j,k]
update(I_triage_R_unconf[,,]) <- new_I_triage_R_unconf[i,j,k]
update(I_triage_R_conf[,,]) <- new_I_triage_R_conf[i,j,k]
update(I_triage_D_unconf[,,]) <- new_I_triage_D_unconf[i,j,k]
update(I_triage_D_conf[,,]) <- new_I_triage_D_conf[i,j,k]
update(I_hosp_R_unconf[,,]) <- new_I_hosp_R_unconf[i,j,k]
update(I_hosp_R_conf[,,]) <- new_I_hosp_R_conf[i,j,k]
update(I_hosp_D_unconf[,,]) <- new_I_hosp_D_unconf[i,j,k]
update(I_hosp_D_conf[,,]) <- new_I_hosp_D_conf[i,j,k]
update(I_ICU_R_unconf[,,]) <- new_I_ICU_R_unconf[i,j,k]
update(I_ICU_R_conf[,,]) <- new_I_ICU_R_conf[i,j,k]
update(I_ICU_D_unconf[,,]) <- new_I_ICU_D_unconf[i,j,k]
update(I_ICU_D_conf[,,]) <- new_I_ICU_D_conf[i,j,k]
update(R_stepdown_unconf[,]) <- new_R_stepdown_unconf[i,j]
update(R_stepdown_conf[,]) <- new_R_stepdown_conf[i,j]
update(R_pre[,]) <- R_pre[i,j] + n_com_to_R_pre[i,j] - n_R_pre[i,j]
update(R_pos[]) <- R_pos[i] + delta_R_pos[i]
update(R_neg[]) <- R_neg[i] + delta_R_neg[i]
update(R[]) <- R[i] + delta_R[i]
update(D_hosp[]) <- D_hosp[i] + delta_D_hosp[i]
update(D_comm[]) <- D_comm[i] + delta_D_comm[i]
update(PCR_pos[,]) <- PCR_pos[i,j] + delta_PCR_pos[i,j]
update(cum_admit_conf) <- cum_admit_conf + sum(n_ILI_to_hosp_D_conf) + sum(n_ILI_to_hosp_R_conf) + sum(n_ILI_to_triage_D_conf) + sum(n_ILI_to_triage_R_conf)
update(cum_new_conf) <- cum_new_conf + sum(n_I_hosp_D_unconf_to_conf) + sum(n_I_hosp_R_unconf_to_conf) + sum(n_I_triage_D_unconf_to_conf) + sum(n_I_triage_R_unconf_to_conf) + sum(n_I_ICU_D_unconf_to_conf) + sum(n_I_ICU_R_unconf_to_conf) + sum(n_R_stepdown_unconf_to_conf)
update(cum_admit_by_age[]) <- cum_admit_by_age[i] + sum(n_ILI_to_hosp[i, ])

## Individual probabilities of transition:
p_SE[] <- 1 - exp(-lambda[i]*dt) # S to I - age dependent
p_EE <- 1 - exp(-gamma_E*dt) # progression of latent period
p_II_asympt <- 1 - exp(-gamma_asympt*dt) # progression of infectious period
p_II_mild <- 1 - exp(-gamma_mild*dt)
p_II_ILI <- 1 - exp(-gamma_ILI*dt)
p_II_comm_D <- 1 - exp(-gamma_comm_D*dt)
p_II_triage <- 1 - exp(-gamma_triage*dt)
p_II_hosp_R <- 1 - exp(-gamma_hosp_R*dt)
p_II_hosp_D <- 1 - exp(-gamma_hosp_D*dt)
p_II_ICU_R <- 1 - exp(-gamma_ICU_R*dt)
p_II_ICU_D <- 1 - exp(-gamma_ICU_D*dt)
p_R_stepdown <- 1 - exp(-gamma_stepdown*dt)
p_R_pre[,] <- 1 - exp(-gamma_R_pre[j]*dt)
p_test <- 1- exp(-gamma_test*dt)
p_PCR_pos <- 1 - exp(-gamma_PCR_pos*dt)

## Draws from binomial distributions for numbers changing between
## compartments:
n_SE[] <- rbinom(S[i], p_SE[i])
n_EE[,,] <- rbinom(E[i,j,k], p_EE)
n_II_asympt[,,] <- rbinom(I_asympt[i,j,k], p_II_asympt)
n_II_mild[,,] <- rbinom(I_mild[i,j,k], p_II_mild)
n_II_ILI[,,] <- rbinom(I_ILI[i,j,k], p_II_ILI)
n_II_comm_D[,,] <- rbinom(I_comm_D[i,j,k], p_II_comm_D)
n_II_triage_R_unconf[,,] <- rbinom(I_triage_R_unconf[i,j,k], p_II_triage)
n_II_triage_R_conf[,,] <- rbinom(I_triage_R_conf[i,j,k], p_II_triage)
n_II_triage_D_unconf[,,] <- rbinom(I_triage_D_unconf[i,j,k], p_II_triage)
n_II_triage_D_conf[,,] <- rbinom(I_triage_D_conf[i,j,k], p_II_triage)
n_II_hosp_R_unconf[,,] <- rbinom(I_hosp_R_unconf[i,j,k], p_II_hosp_R)
n_II_hosp_R_conf[,,] <- rbinom(I_hosp_R_conf[i,j,k], p_II_hosp_R)
n_II_hosp_D_unconf[,,] <- rbinom(I_hosp_D_unconf[i,j,k], p_II_hosp_D)
n_II_hosp_D_conf[,,] <- rbinom(I_hosp_D_conf[i,j,k], p_II_hosp_D)
n_II_ICU_R_unconf[,,] <- rbinom(I_ICU_R_unconf[i,j,k], p_II_ICU_R)
n_II_ICU_R_conf[,,] <- rbinom(I_ICU_R_conf[i,j,k], p_II_ICU_R)
n_II_ICU_D_unconf[,,] <- rbinom(I_ICU_D_unconf[i,j,k], p_II_ICU_D)
n_II_ICU_D_conf[,,] <- rbinom(I_ICU_D_conf[i,j,k], p_II_ICU_D)
n_R_stepdown_unconf[,] <- rbinom(R_stepdown_unconf[i,j], p_R_stepdown)
n_R_stepdown_conf[,] <- rbinom(R_stepdown_conf[i,j], p_R_stepdown)
n_R_pre[,] <- rbinom(R_pre[i,j], p_R_pre[i,j])
n_PCR_pos[,] <- rbinom(PCR_pos[i,j], p_PCR_pos)


#Computes the number of asymptomatic
n_EI_asympt[,] <- rbinom(n_EE[i,s_E,j], p_asympt[i])

#Computes the number of mild cases - p_sympt_ILI gives the proportion of febrile/ILI cases among the symptomatics
n_EI_mild[,] <- rbinom(n_EE[i,s_E,j]-n_EI_asympt[i,j], 1-p_sympt_ILI[i])

#Computes the number of ILI cases
n_EI_ILI[,] <- n_EE[i,s_E,j]-n_EI_asympt[i,j]-n_EI_mild[i,j]

#Compute the aux_p_bin matrix of binom nested coeff
aux_p_bin[,1] <- trans_profile[i,1]
aux_p_bin[,2:(trans_classes-1)] <- trans_profile[i,j]/sum(trans_profile[i,j:trans_classes])

#Implementation of multinom via nested binomial
aux_EE[,1,1] <- rbinom(n_SE[i],aux_p_bin[i,1])
aux_EE[,1,2:(trans_classes - 1)] <- rbinom(n_SE[i] - sum(aux_EE[i,1,1:(k-1)]), aux_p_bin[i,k])
aux_EE[,1,trans_classes] <- n_SE[i] - sum(aux_EE[i,1,1:(trans_classes - 1)])

#Work out the E->E transitions
aux_EE[,2:s_E,] <- n_EE[i,j-1,k]
aux_EE[,1:s_E,] <- aux_EE[i,j,k] - n_EE[i,j,k]
delta_E[,,] <- aux_EE[i,j,k]

#Work out the I_asympt->I_asympt transitions
aux_II_asympt[,1,] <- n_EI_asympt[i,k]
aux_II_asympt[,2:s_asympt,] <- n_II_asympt[i,j-1,k]
aux_II_asympt[,1:s_asympt,] <- aux_II_asympt[i,j,k] - n_II_asympt[i,j,k]
delta_I_asympt[,,] <- aux_II_asympt[i,j,k]

#Work out the I_mild->I_mild transitions
aux_II_mild[,1,] <- n_EI_mild[i,k]
aux_II_mild[,2:s_mild,] <- n_II_mild[i,j-1,k]
aux_II_mild[,1:s_mild,] <- aux_II_mild[i,j,k] - n_II_mild[i,j,k]
delta_I_mild[,,] <- aux_II_mild[i,j,k]

#Work out the I_ILI->I_ILI transitions
aux_II_ILI[,1,] <- n_EI_ILI[i,k]
aux_II_ILI[,2:s_ILI,] <- n_II_ILI[i,j-1,k]
aux_II_ILI[,1:s_ILI,] <- aux_II_ILI[i,j,k] - n_II_ILI[i,j,k]
delta_I_ILI[,,] <- aux_II_ILI[i,j,k]

#Work out the flow from I_ILI -> R, comm_D, hosp
n_ILI_to_R[,] <- rbinom(n_II_ILI[i,s_ILI,j],1-p_hosp_ILI[i])
n_ILI_to_comm_D[,] <- rbinom(n_II_ILI[i,s_ILI,j] - n_ILI_to_R[i,j],p_death_comm[i])
n_ILI_to_hosp[,] <- n_II_ILI[i,s_ILI,j] - n_ILI_to_R[i,j] - n_ILI_to_comm_D[i,j]

#Work out the I_comm_D -> I_comm_D transitions
aux_II_comm_D[,1,] <- n_ILI_to_comm_D[i,k]
aux_II_comm_D[,2:s_comm_D,] <- n_II_comm_D[i,j-1,k]
aux_II_comm_D[,1:s_comm_D,] <- aux_II_comm_D[i,j,k] - n_II_comm_D[i,j,k]
delta_I_comm_D[,,] <- aux_II_comm_D[i,j,k]

#Work out the split in hospitals between hosp_D, hosp_R, triage_R and triage_D
n_ILI_to_triage[,] <- rbinom(n_ILI_to_hosp[i,j],p_ICU_hosp[i])
n_hosp_non_ICU[,] <- n_ILI_to_hosp[i,j] - n_ILI_to_triage[i,j]
n_ILI_to_hosp_D[,] <- rbinom(n_hosp_non_ICU[i,j], p_death_hosp_D[i])
n_ILI_to_hosp_D_conf[,] <- rbinom(n_ILI_to_hosp_D[i,j],p_admit_conf[i])
n_ILI_to_hosp_R[,] <- n_hosp_non_ICU[i,j] - n_ILI_to_hosp_D[i,j]
n_ILI_to_hosp_R_conf[,] <- rbinom(n_ILI_to_hosp_R[i,j],p_admit_conf[i])
n_ILI_to_triage_D[,] <- rbinom(n_ILI_to_triage[i,j], p_death_ICU[i])
n_ILI_to_triage_D_conf[,] <- rbinom(n_ILI_to_triage_D[i,j],p_admit_conf[i])
n_ILI_to_triage_R[,] <- n_ILI_to_triage[i,j] - n_ILI_to_triage_D[i,j]
n_ILI_to_triage_R_conf[,] <- rbinom(n_ILI_to_triage_R[i,j],p_admit_conf[i])

#Work out the I_triage_R -> I_triage_R transitions
aux_II_triage_R_unconf[,,] <- I_triage_R_unconf[i,j,k]
aux_II_triage_R_unconf[,2:s_triage,] <- aux_II_triage_R_unconf[i,j,k] + n_II_triage_R_unconf[i,j-1,k]
aux_II_triage_R_unconf[,1:s_triage,] <- aux_II_triage_R_unconf[i,j,k] - n_II_triage_R_unconf[i,j,k]
aux_II_triage_R_conf[,,] <- I_triage_R_conf[i,j,k]
aux_II_triage_R_conf[,2:s_triage,] <- aux_II_triage_R_conf[i,j,k] + n_II_triage_R_conf[i,j-1,k]
aux_II_triage_R_conf[,1:s_triage,] <- aux_II_triage_R_conf[i,j,k] - n_II_triage_R_conf[i,j,k]
n_I_triage_R_unconf_to_conf[,,] <- rbinom(aux_II_triage_R_unconf[i,j,k],p_test)
new_I_triage_R_unconf[,,] <- aux_II_triage_R_unconf[i,j,k] - n_I_triage_R_unconf_to_conf[i,j,k]
new_I_triage_R_unconf[,1,] <- new_I_triage_R_unconf[i,1,k] + n_ILI_to_triage_R[i,k] - n_ILI_to_triage_R_conf[i,k]
new_I_triage_R_conf[,,] <- aux_II_triage_R_conf[i,j,k] + n_I_triage_R_unconf_to_conf[i,j,k]
new_I_triage_R_conf[,1,] <- new_I_triage_R_conf[i,1,k] + n_ILI_to_triage_R_conf[i,k]

#Work out the I_triage_D -> I_triage_D transitions
aux_II_triage_D_unconf[,,] <- I_triage_D_unconf[i,j,k]
aux_II_triage_D_unconf[,2:s_triage,] <- aux_II_triage_D_unconf[i,j,k] + n_II_triage_D_unconf[i,j-1,k]
aux_II_triage_D_unconf[,1:s_triage,] <- aux_II_triage_D_unconf[i,j,k] - n_II_triage_D_unconf[i,j,k]
aux_II_triage_D_conf[,,] <- I_triage_D_conf[i,j,k]
aux_II_triage_D_conf[,2:s_triage,] <- aux_II_triage_D_conf[i,j,k] + n_II_triage_D_conf[i,j-1,k]
aux_II_triage_D_conf[,1:s_triage,] <- aux_II_triage_D_conf[i,j,k] - n_II_triage_D_conf[i,j,k]
n_I_triage_D_unconf_to_conf[,,] <- rbinom(aux_II_triage_D_unconf[i,j,k],p_test)
new_I_triage_D_unconf[,,] <- aux_II_triage_D_unconf[i,j,k] - n_I_triage_D_unconf_to_conf[i,j,k]
new_I_triage_D_unconf[,1,] <- new_I_triage_D_unconf[i,1,k] + n_ILI_to_triage_D[i,k] - n_ILI_to_triage_D_conf[i,k]
new_I_triage_D_conf[,,] <- aux_II_triage_D_conf[i,j,k] + n_I_triage_D_unconf_to_conf[i,j,k]
new_I_triage_D_conf[,1,] <- new_I_triage_D_conf[i,1,k] + n_ILI_to_triage_D_conf[i,k]

#Work out the I_hosp_R->I_hosp_R transitions
aux_II_hosp_R_unconf[,,] <- I_hosp_R_unconf[i,j,k]
aux_II_hosp_R_unconf[,2:s_hosp_R,] <- aux_II_hosp_R_unconf[i,j,k] + n_II_hosp_R_unconf[i,j-1,k]
aux_II_hosp_R_unconf[,1:s_hosp_R,] <- aux_II_hosp_R_unconf[i,j,k] - n_II_hosp_R_unconf[i,j,k]
aux_II_hosp_R_conf[,,] <- I_hosp_R_conf[i,j,k]
aux_II_hosp_R_conf[,2:s_hosp_R,] <- aux_II_hosp_R_conf[i,j,k] + n_II_hosp_R_conf[i,j-1,k]
aux_II_hosp_R_conf[,1:s_hosp_R,] <- aux_II_hosp_R_conf[i,j,k] - n_II_hosp_R_conf[i,j,k]
n_I_hosp_R_unconf_to_conf[,,] <- rbinom(aux_II_hosp_R_unconf[i,j,k],p_test)
new_I_hosp_R_unconf[,,] <- aux_II_hosp_R_unconf[i,j,k] - n_I_hosp_R_unconf_to_conf[i,j,k]
new_I_hosp_R_unconf[,1,] <- new_I_hosp_R_unconf[i,1,k] + n_ILI_to_hosp_R[i,k] - n_ILI_to_hosp_R_conf[i,k]
new_I_hosp_R_conf[,,] <- aux_II_hosp_R_conf[i,j,k] + n_I_hosp_R_unconf_to_conf[i,j,k]
new_I_hosp_R_conf[,1,] <- new_I_hosp_R_conf[i,1,k] + n_ILI_to_hosp_R_conf[i,k]

#Work out the I_hosp_D->I_hosp_D transitions
aux_II_hosp_D_unconf[,,] <- I_hosp_D_unconf[i,j,k]
aux_II_hosp_D_unconf[,2:s_hosp_D,] <- aux_II_hosp_D_unconf[i,j,k] + n_II_hosp_D_unconf[i,j-1,k]
aux_II_hosp_D_unconf[,1:s_hosp_D,] <- aux_II_hosp_D_unconf[i,j,k] - n_II_hosp_D_unconf[i,j,k]
aux_II_hosp_D_conf[,,] <- I_hosp_D_conf[i,j,k]
aux_II_hosp_D_conf[,2:s_hosp_D,] <- aux_II_hosp_D_conf[i,j,k] + n_II_hosp_D_conf[i,j-1,k]
aux_II_hosp_D_conf[,1:s_hosp_D,] <- aux_II_hosp_D_conf[i,j,k] - n_II_hosp_D_conf[i,j,k]
n_I_hosp_D_unconf_to_conf[,,] <- rbinom(aux_II_hosp_D_unconf[i,j,k],p_test)
new_I_hosp_D_unconf[,,] <- aux_II_hosp_D_unconf[i,j,k] - n_I_hosp_D_unconf_to_conf[i,j,k]
new_I_hosp_D_unconf[,1,] <- new_I_hosp_D_unconf[i,1,k] + n_ILI_to_hosp_D[i,k] - n_ILI_to_hosp_D_conf[i,k]
new_I_hosp_D_conf[,,] <- aux_II_hosp_D_conf[i,j,k] + n_I_hosp_D_unconf_to_conf[i,j,k]
new_I_hosp_D_conf[,1,] <- new_I_hosp_D_conf[i,1,k] + n_ILI_to_hosp_D_conf[i,k]

#Work out the I_ICU_R->I_ICU_R transitions
aux_II_ICU_R_unconf[,,] <- I_ICU_R_unconf[i,j,k]
aux_II_ICU_R_unconf[,1,] <- aux_II_ICU_R_unconf[i,j,k] + n_II_triage_R_unconf[i,s_triage,k]
aux_II_ICU_R_unconf[,2:s_ICU_R,] <- aux_II_ICU_R_unconf[i,j,k] + n_II_ICU_R_unconf[i,j-1,k]
aux_II_ICU_R_unconf[,1:s_ICU_R,] <- aux_II_ICU_R_unconf[i,j,k] - n_II_ICU_R_unconf[i,j,k]
aux_II_ICU_R_conf[,,] <- I_ICU_R_conf[i,j,k]
aux_II_ICU_R_conf[,1,] <- aux_II_ICU_R_conf[i,j,k] + n_II_triage_R_conf[i,s_triage,k]
aux_II_ICU_R_conf[,2:s_ICU_R,] <- aux_II_ICU_R_conf[i,j,k] + n_II_ICU_R_conf[i,j-1,k]
aux_II_ICU_R_conf[,1:s_ICU_R,] <- aux_II_ICU_R_conf[i,j,k] - n_II_ICU_R_conf[i,j,k]
n_I_ICU_R_unconf_to_conf[,,] <- rbinom(aux_II_ICU_R_unconf[i,j,k],p_test)
new_I_ICU_R_unconf[,,] <- aux_II_ICU_R_unconf[i,j,k] - n_I_ICU_R_unconf_to_conf[i,j,k]
new_I_ICU_R_conf[,,] <- aux_II_ICU_R_conf[i,j,k] + n_I_ICU_R_unconf_to_conf[i,j,k]

#Work out the I_ICU_D->I_ICU_D transitions
aux_II_ICU_D_unconf[,,] <- I_ICU_D_unconf[i,j,k]
aux_II_ICU_D_unconf[,1,] <- aux_II_ICU_D_unconf[i,j,k] + n_II_triage_D_unconf[i,s_triage,k]
aux_II_ICU_D_unconf[,2:s_ICU_D,] <- aux_II_ICU_D_unconf[i,j,k] + n_II_ICU_D_unconf[i,j-1,k]
aux_II_ICU_D_unconf[,1:s_ICU_D,] <- aux_II_ICU_D_unconf[i,j,k] - n_II_ICU_D_unconf[i,j,k]
aux_II_ICU_D_conf[,,] <- I_ICU_D_conf[i,j,k]
aux_II_ICU_D_conf[,1,] <- aux_II_ICU_D_conf[i,j,k] + n_II_triage_D_conf[i,s_triage,k]
aux_II_ICU_D_conf[,2:s_ICU_D,] <- aux_II_ICU_D_conf[i,j,k] + n_II_ICU_D_conf[i,j-1,k]
aux_II_ICU_D_conf[,1:s_ICU_D,] <- aux_II_ICU_D_conf[i,j,k] - n_II_ICU_D_conf[i,j,k]
n_I_ICU_D_unconf_to_conf[,,] <- rbinom(aux_II_ICU_D_unconf[i,j,k],p_test)
new_I_ICU_D_unconf[,,] <- aux_II_ICU_D_unconf[i,j,k] - n_I_ICU_D_unconf_to_conf[i,j,k]
new_I_ICU_D_conf[,,] <- aux_II_ICU_D_conf[i,j,k] + n_I_ICU_D_unconf_to_conf[i,j,k]

#Work out the R_stepdown->R_stepdown transitions
aux_R_stepdown_unconf[,] <- R_stepdown_unconf[i,j]
aux_R_stepdown_unconf[,1] <- aux_R_stepdown_unconf[i,j] + sum(n_II_ICU_R_unconf[i,s_ICU_R,])
aux_R_stepdown_unconf[,2:s_stepdown] <- aux_R_stepdown_unconf[i,j] + n_R_stepdown_unconf[i,j-1]
aux_R_stepdown_unconf[,1:s_stepdown] <- aux_R_stepdown_unconf[i,j] - n_R_stepdown_unconf[i,j]
aux_R_stepdown_conf[,] <- R_stepdown_conf[i,j]
aux_R_stepdown_conf[,1] <- aux_R_stepdown_conf[i,j] + sum(n_II_ICU_R_conf[i,s_ICU_R,])
aux_R_stepdown_conf[,2:s_stepdown] <- aux_R_stepdown_conf[i,j] + n_R_stepdown_conf[i,j-1]
aux_R_stepdown_conf[,1:s_stepdown] <- aux_R_stepdown_conf[i,j] - n_R_stepdown_conf[i,j]
n_R_stepdown_unconf_to_conf[,] <- rbinom(aux_R_stepdown_unconf[i,j],p_test)
new_R_stepdown_unconf[,] <- aux_R_stepdown_unconf[i,j] - n_R_stepdown_unconf_to_conf[i,j]
new_R_stepdown_conf[,] <- aux_R_stepdown_conf[i,j] + n_R_stepdown_unconf_to_conf[i,j]

#Work out the number of deaths in hospital
delta_D_hosp[] <- sum(n_II_hosp_D_unconf[i,s_hosp_D,]) + sum(n_II_hosp_D_conf[i,s_hosp_D,]) + sum(n_II_ICU_D_unconf[i,s_ICU_D,]) + sum(n_II_ICU_D_conf[i,s_ICU_D,])

#Work out the number of deaths in the community
delta_D_comm[] <- sum(n_II_comm_D[i,s_comm_D,])

#Work out the number of people entering the seroconversion flow
n_com_to_R_pre[,1] <- rbinom(sum(n_EE[i,s_E,]),p_R_pre_1)
n_com_to_R_pre[,2] <- sum(n_EE[i,s_E,]) - n_com_to_R_pre[i,1]

#Split the seroconversion flow between people who are going to seroconvert and people who are not
delta_R_pos[] <- rbinom(sum(n_R_pre[i,]),p_seroconversion[i])
delta_R_neg[] <- sum(n_R_pre[i,]) - delta_R_pos[i]

#Work out the total number of recovery
delta_R[] <- sum(n_II_asympt[i,s_asympt,]) + sum(n_II_mild[i,s_mild,]) + sum(n_ILI_to_R[i,]) +
  sum(n_II_hosp_R_conf[i,s_hosp_R,]) + sum(n_II_hosp_R_unconf[i,s_hosp_R,]) +
  sum(n_R_stepdown_conf[i,s_stepdown]) + sum(n_R_stepdown_unconf[i,s_stepdown])

#Work out the PCR positivity
delta_PCR_pos[,1] <- n_SE[i]
delta_PCR_pos[,2:s_PCR_pos] <- n_PCR_pos[i,j-1]
delta_PCR_pos[,] <- delta_PCR_pos[i,j] - n_PCR_pos[i,j]

#Compute the force of infection
I_with_diff_trans[,] <- trans_increase[i,j]*(sum(I_asympt[i,,j])+
                                               sum(I_mild[i,,j])+sum(I_ILI[i,,j])+
                                               hosp_transmission*
                                               (sum(I_triage_R_unconf[i,,j])+sum(I_triage_R_conf[i,,j])+
                                                  sum(I_triage_D_unconf[i,,j])+sum(I_triage_D_conf[i,,j])+
                                                  sum(I_hosp_R_unconf[i,,j])+sum(I_hosp_R_conf[i,,j])+
                                                  sum(I_hosp_D_unconf[i,,j])+sum(I_hosp_D_conf[i,,j]))+
                                               ICU_transmission*
                                                 (sum(I_ICU_R_unconf[i,,j])+sum(I_ICU_R_conf[i,,j])+
                                                    sum(I_ICU_D_unconf[i,,j])+sum(I_ICU_D_conf[i,,j]))+
                                               comm_D_transmission*
                                                 sum(I_comm_D[i,,j]))
s_ij[,] <- m[i,j] * sum(I_with_diff_trans[j,])
lambda[] <- beta*sum(s_ij[i,]) + psi*importation/N_tot2

## Initial states:
initial(S[]) <- S0[i] # will be user-defined
initial(E[,,]) <- E0[i,j,k] # will be user-defined
initial(I_asympt[,,]) <- I0_asympt[i,j,k] # will be user-defined
initial(I_mild[,,]) <- I0_mild[i,j,k] # will be user-defined
initial(I_ILI[,,]) <- I0_ILI[i,j,k] # will be user-defined
initial(I_comm_D[,,]) <- I0_comm_D[i,j,k]
initial(I_triage_R_unconf[,,]) <- I0_triage_R_unconf[i,j,k]
initial(I_triage_R_conf[,,]) <- I0_triage_R_conf[i,j,k]
initial(I_triage_D_unconf[,,]) <- I0_triage_D_unconf[i,j,k]
initial(I_triage_D_conf[,,]) <- I0_triage_D_conf[i,j,k]
initial(I_hosp_R_unconf[,,]) <- I0_hosp_R_unconf[i,j,k]
initial(I_hosp_R_conf[,,]) <- I0_hosp_R_conf[i,j,k]
initial(I_hosp_D_unconf[,,]) <- I0_hosp_D_unconf[i,j,k]
initial(I_hosp_D_conf[,,]) <- I0_hosp_D_conf[i,j,k]
initial(I_ICU_R_unconf[,,]) <- I0_ICU_R_unconf[i,j,k]
initial(I_ICU_R_conf[,,]) <- I0_ICU_R_conf[i,j,k]
initial(I_ICU_D_unconf[,,]) <- I0_ICU_D_unconf[i,j,k]
initial(I_ICU_D_conf[,,]) <- I0_ICU_D_conf[i,j,k]
initial(R_stepdown_unconf[,]) <- R0_stepdown_unconf[i,j]
initial(R_stepdown_conf[,]) <- R0_stepdown_conf[i,j]
initial(R_pre[,]) <- R0_pre[i,j]
initial(R_pos[]) <- R0_pos[i]
initial(R_neg[]) <- R0_neg[i]
initial(R[]) <- R0[i]
initial(D_hosp[]) <- D0_hosp[i]
initial(D_comm[]) <- D0_comm[i]
initial(PCR_pos[,]) <- PCR0_pos[i,j]
initial(cum_admit_conf) <- cum0_admit_conf
initial(cum_new_conf) <- cum0_new_conf
initial(cum_admit_by_age[]) <- cum0_admit_by_age[i]

## User defined parameters - default in parentheses:

#Initial vectors
S0[] <- user()
E0[,,] <- user()
I0_asympt[,,] <- user()
I0_mild[,,] <- user()
I0_ILI[,,] <- user()
I0_comm_D[,,] <- user()
I0_triage_R_unconf[,,] <- user()
I0_triage_R_conf[,,] <- user()
I0_triage_D_unconf[,,] <- user()
I0_triage_D_conf[,,] <- user()
I0_hosp_R_unconf[,,] <- user()
I0_hosp_R_conf[,,] <- user()
I0_hosp_D_unconf[,,] <- user()
I0_hosp_D_conf[,,] <- user()
I0_ICU_R_unconf[,,] <- user()
I0_ICU_R_conf[,,] <- user()
I0_ICU_D_unconf[,,] <- user()
I0_ICU_D_conf[,,] <- user()
R0_stepdown_unconf[,] <- user()
R0_stepdown_conf[,] <- user()
R0[] <- user()
R0_pre[,] <- user()
R0_neg[] <- user()
R0_pos[] <- user()
D0_hosp[] <- user()
D0_comm[] <- user()
PCR0_pos[,] <- user()
cum0_admit_conf <- user()
cum0_new_conf <- user()
cum0_admit_by_age[] <- user()

#Parameters of the E classes
s_E <- user()
gamma_E <- user(0.1)

#Probability of transitioning from the E to the asymptomatic, febrile classes, the rest goes in mild i.e with mild symptoms
p_asympt[] <- user()
p_sympt_ILI[] <- user()

#Parameters of the I_asympt classes
s_asympt <- user()
gamma_asympt <- user(0.1)

#Parameters of the I_mild classes
s_mild <- user()
gamma_mild <- user(0.1)

#Parameters of the I_ILI classes
s_ILI <- user()
gamma_ILI <- user(0.1)
p_hosp_ILI[] <- user()

#Parameters of the I_comm_D class
s_comm_D <- user()
gamma_comm_D <- user(0.1)
p_death_comm[] <- user()

#Parameters of the I_triage classes
s_triage <- user()
gamma_triage <- user(0.1)

#Proportion of hospital cases progressing to ICU
p_ICU_hosp[] <- user()

#Parameters of the I_hosp_R classes
s_hosp_R <- user()
gamma_hosp_R <- user(0.1)

#Parameters of the I_hosp_D classes
s_hosp_D <- user()
gamma_hosp_D <- user(0.1)
p_death_hosp_D[] <- user()

#Parameters of the I_ICU_R classes
s_ICU_R <- user()
gamma_ICU_R <- user(0.1)

#Parameters of the I_ICU classes
s_ICU_D <- user()
gamma_ICU_D <- user(0.1)
p_death_ICU[] <- user()

#Parameters of the R_stepdown classes
s_stepdown <- user()
gamma_stepdown <- user(0.1)

#Parameters of the R_pre classes
gamma_R_pre_1 <- user(0.1)
gamma_R_pre_2 <- user(0.1)
gamma_R_pre[1] <- gamma_R_pre_1
gamma_R_pre[2] <- gamma_R_pre_2
p_R_pre_1 <- user()
p_seroconversion[] <- user()

#Parameters relating to testing
gamma_test <- user(0.1)
p_admit_conf[] <- user()

#Parameters relating to PCR positivity
s_PCR_pos <- user()
gamma_PCR_pos <- user(0.1)

#Parameters of the age stratified transmission
beta_step[] <- user()
dim(beta_step) <- user()
## What we really want is min(step + 1, length(beta_step)) but that's not
## supported by odin (it could be made to support this). This code
## does currently create a compiler warning with -Wsign-compare on
## because we have an unsigned/signed integer comparison
beta <- if (step >= length(beta_step)) beta_step[length(beta_step)] else beta_step[step + 1]

# Useful for debugging
initial(beta_out) <- beta_step[1]
update(beta_out) <- beta

#Parameters of the external infection rate
importation_step[] <- user()
dim(importation_step) <- user()
importation <- if (step >= length(importation_step)) importation_step[length(importation_step)] else importation_step[step + 1]

psi <- user(0.1)

m[,] <- user()
trans_profile[,] <- user()
trans_increase[,] <- user()
hosp_transmission <- user()
ICU_transmission <- user()
comm_D_transmission <- user()

##Dimensions of the different "vectors" here vectors stand for multi-dimensional arrays

#Vectors handling the S class
dim(S) <- N_age
dim(S0) <- N_age

#Vectors handling the E class
dim(E) <- c(N_age,s_E,trans_classes)
dim(E0) <- c(N_age,s_E,trans_classes)
dim(aux_EE) <- c(N_age,s_E,trans_classes)
dim(delta_E) <- c(N_age,s_E,trans_classes)
dim(n_EE) <- c(N_age,s_E,trans_classes)

#Vectors handling the I_asympt class
dim(I_asympt) <- c(N_age,s_asympt,trans_classes)
dim(I0_asympt) <- c(N_age,s_asympt,trans_classes)
dim(aux_II_asympt) <- c(N_age,s_asympt,trans_classes)
dim(delta_I_asympt) <- c(N_age,s_asympt,trans_classes)
dim(n_II_asympt) <- c(N_age,s_asympt,trans_classes)

#Vectors handling the I_mild class
dim(I_mild) <- c(N_age,s_mild,trans_classes)
dim(I0_mild) <- c(N_age,s_mild,trans_classes)
dim(aux_II_mild) <- c(N_age,s_mild,trans_classes)
dim(delta_I_mild) <- c(N_age,s_mild,trans_classes)
dim(n_II_mild) <- c(N_age,s_mild,trans_classes)

#Vectors handling the I_ILI class
dim(I_ILI) <- c(N_age,s_ILI,trans_classes)
dim(I0_ILI) <- c(N_age,s_ILI,trans_classes)
dim(aux_II_ILI) <- c(N_age,s_ILI,trans_classes)
dim(delta_I_ILI) <- c(N_age,s_ILI,trans_classes)
dim(n_II_ILI) <- c(N_age,s_ILI,trans_classes)
dim(p_hosp_ILI) <- c(N_age)

#Vectors handling the I_comm_D class
dim(I_comm_D) <- c(N_age,s_comm_D,trans_classes)
dim(I0_comm_D) <- c(N_age,s_comm_D,trans_classes)
dim(aux_II_comm_D) <- c(N_age,s_comm_D,trans_classes)
dim(delta_I_comm_D) <- c(N_age,s_comm_D,trans_classes)
dim(n_II_comm_D) <- c(N_age,s_comm_D,trans_classes)
dim(p_death_comm) <- c(N_age)

#Vectors handling the I_triage_R class
dim(I_triage_R_unconf) <- c(N_age,s_triage,trans_classes)
dim(I0_triage_R_unconf) <- c(N_age,s_triage,trans_classes)
dim(aux_II_triage_R_unconf) <- c(N_age,s_triage,trans_classes)
dim(new_I_triage_R_unconf) <- c(N_age,s_triage,trans_classes)
dim(n_II_triage_R_unconf) <- c(N_age,s_triage,trans_classes)
dim(I_triage_R_conf) <- c(N_age,s_triage,trans_classes)
dim(I0_triage_R_conf) <- c(N_age,s_triage,trans_classes)
dim(aux_II_triage_R_conf) <- c(N_age,s_triage,trans_classes)
dim(new_I_triage_R_conf) <- c(N_age,s_triage,trans_classes)
dim(n_II_triage_R_conf) <- c(N_age,s_triage,trans_classes)
dim(n_I_triage_R_unconf_to_conf) <- c(N_age,s_triage,trans_classes)

#Vectors handling the I_triage_D class
dim(I_triage_D_unconf) <- c(N_age,s_triage,trans_classes)
dim(I0_triage_D_unconf) <- c(N_age,s_triage,trans_classes)
dim(aux_II_triage_D_unconf) <- c(N_age,s_triage,trans_classes)
dim(new_I_triage_D_unconf) <- c(N_age,s_triage,trans_classes)
dim(n_II_triage_D_unconf) <- c(N_age,s_triage,trans_classes)
dim(I_triage_D_conf) <- c(N_age,s_triage,trans_classes)
dim(I0_triage_D_conf) <- c(N_age,s_triage,trans_classes)
dim(aux_II_triage_D_conf) <- c(N_age,s_triage,trans_classes)
dim(new_I_triage_D_conf) <- c(N_age,s_triage,trans_classes)
dim(n_II_triage_D_conf) <- c(N_age,s_triage,trans_classes)
dim(n_I_triage_D_unconf_to_conf) <- c(N_age,s_triage,trans_classes)

#Vector handling who progress to ICU
dim(p_ICU_hosp) <- c(N_age)

#Vectors handling the I_hosp_R class
dim(I_hosp_R_unconf) <- c(N_age,s_hosp_R,trans_classes)
dim(I0_hosp_R_unconf) <- c(N_age,s_hosp_R,trans_classes)
dim(aux_II_hosp_R_unconf) <- c(N_age,s_hosp_R,trans_classes)
dim(new_I_hosp_R_unconf) <- c(N_age,s_hosp_R,trans_classes)
dim(n_II_hosp_R_unconf) <- c(N_age,s_hosp_R,trans_classes)
dim(I_hosp_R_conf) <- c(N_age,s_hosp_R,trans_classes)
dim(I0_hosp_R_conf) <- c(N_age,s_hosp_R,trans_classes)
dim(aux_II_hosp_R_conf) <- c(N_age,s_hosp_R,trans_classes)
dim(new_I_hosp_R_conf) <- c(N_age,s_hosp_R,trans_classes)
dim(n_II_hosp_R_conf) <- c(N_age,s_hosp_R,trans_classes)
dim(n_I_hosp_R_unconf_to_conf) <- c(N_age,s_hosp_R,trans_classes)

#Vectors handling the I_hosp_D class
dim(I_hosp_D_unconf) <- c(N_age,s_hosp_D,trans_classes)
dim(I0_hosp_D_unconf) <- c(N_age,s_hosp_D,trans_classes)
dim(aux_II_hosp_D_unconf) <- c(N_age,s_hosp_D,trans_classes)
dim(new_I_hosp_D_unconf) <- c(N_age,s_hosp_D,trans_classes)
dim(n_II_hosp_D_unconf) <- c(N_age,s_hosp_D,trans_classes)
dim(I_hosp_D_conf) <- c(N_age,s_hosp_D,trans_classes)
dim(I0_hosp_D_conf) <- c(N_age,s_hosp_D,trans_classes)
dim(aux_II_hosp_D_conf) <- c(N_age,s_hosp_D,trans_classes)
dim(new_I_hosp_D_conf) <- c(N_age,s_hosp_D,trans_classes)
dim(n_II_hosp_D_conf) <- c(N_age,s_hosp_D,trans_classes)
dim(n_I_hosp_D_unconf_to_conf) <- c(N_age,s_hosp_D,trans_classes)

#Vectors handling the I_ICU_R class
dim(I_ICU_R_unconf) <- c(N_age,s_ICU_R,trans_classes)
dim(I0_ICU_R_unconf) <- c(N_age,s_ICU_R,trans_classes)
dim(aux_II_ICU_R_unconf) <- c(N_age,s_ICU_R,trans_classes)
dim(new_I_ICU_R_unconf) <- c(N_age,s_ICU_R,trans_classes)
dim(n_II_ICU_R_unconf) <- c(N_age,s_ICU_R,trans_classes)
dim(I_ICU_R_conf) <- c(N_age,s_ICU_R,trans_classes)
dim(I0_ICU_R_conf) <- c(N_age,s_ICU_R,trans_classes)
dim(aux_II_ICU_R_conf) <- c(N_age,s_ICU_R,trans_classes)
dim(new_I_ICU_R_conf) <- c(N_age,s_ICU_R,trans_classes)
dim(n_II_ICU_R_conf) <- c(N_age,s_ICU_R,trans_classes)
dim(n_I_ICU_R_unconf_to_conf) <- c(N_age,s_ICU_R,trans_classes)

#Vectors handling the I_ICU_D class
dim(I_ICU_D_unconf) <- c(N_age,s_ICU_D,trans_classes)
dim(I0_ICU_D_unconf) <- c(N_age,s_ICU_D,trans_classes)
dim(aux_II_ICU_D_unconf) <- c(N_age,s_ICU_D,trans_classes)
dim(new_I_ICU_D_unconf) <- c(N_age,s_ICU_D,trans_classes)
dim(n_II_ICU_D_unconf) <- c(N_age,s_ICU_D,trans_classes)
dim(I_ICU_D_conf) <- c(N_age,s_ICU_D,trans_classes)
dim(I0_ICU_D_conf) <- c(N_age,s_ICU_D,trans_classes)
dim(aux_II_ICU_D_conf) <- c(N_age,s_ICU_D,trans_classes)
dim(new_I_ICU_D_conf) <- c(N_age,s_ICU_D,trans_classes)
dim(n_II_ICU_D_conf) <- c(N_age,s_ICU_D,trans_classes)
dim(n_I_ICU_D_unconf_to_conf) <- c(N_age,s_ICU_D,trans_classes)

#Vectors handling the R_stepdown class
dim(R_stepdown_unconf) <- c(N_age,s_stepdown)
dim(R0_stepdown_unconf) <- c(N_age,s_stepdown)
dim(aux_R_stepdown_unconf) <- c(N_age,s_stepdown)
dim(new_R_stepdown_unconf) <- c(N_age,s_stepdown)
dim(n_R_stepdown_unconf) <- c(N_age,s_stepdown)
dim(R_stepdown_conf) <- c(N_age,s_stepdown)
dim(R0_stepdown_conf) <- c(N_age,s_stepdown)
dim(aux_R_stepdown_conf) <- c(N_age,s_stepdown)
dim(new_R_stepdown_conf) <- c(N_age,s_stepdown)
dim(n_R_stepdown_conf) <- c(N_age,s_stepdown)
dim(n_R_stepdown_unconf_to_conf) <- c(N_age,s_stepdown)

#Vectors handling the R_pos class
dim(R) <- c(N_age)
dim(R0) <- c(N_age)
dim(delta_R) <- c(N_age)

#Vectors handling the R_pre class and seroconversion
dim(R_pre) <- c(N_age, 2)
dim(R0_pre) <- c(N_age, 2)
dim(n_R_pre) <- c(N_age, 2)
dim(gamma_R_pre) <- c(2)
dim(p_R_pre) <- c(N_age,2)
dim(p_seroconversion) <- c(N_age)

#Vectors handling the R_pos class
dim(R_pos) <- c(N_age)
dim(R0_pos) <- c(N_age)
dim(delta_R_pos) <- c(N_age)

#Vectors handling the R_neg class
dim(R_neg) <- c(N_age)
dim(R0_neg) <- c(N_age)
dim(delta_R_neg) <- c(N_age)

#Vectors handling the D_hosp class
dim(D_hosp) <- c(N_age)
dim(D0_hosp) <- c(N_age)
dim(delta_D_hosp) <- c(N_age)

#Vectors handling the D_comm class
dim(D_comm) <- c(N_age)
dim(D0_comm) <- c(N_age)
dim(delta_D_comm) <- c(N_age)

#Vectors handling the R_pos class
dim(PCR_pos) <- c(N_age,s_PCR_pos)
dim(PCR0_pos) <- c(N_age,s_PCR_pos)
dim(delta_PCR_pos) <- c(N_age,s_PCR_pos)
dim(n_PCR_pos) <- c(N_age,s_PCR_pos)

#Vectors handling the S->E transition where infected are split between level of infectivity
dim(p_SE) <- N_age
dim(n_SE) <- N_age
dim(aux_p_bin) <- c(N_age,trans_classes)

#Vectors handling the E->I transition where newly infectious cases are split between level of severity
dim(n_EI_asympt) <- c(N_age,trans_classes)
dim(n_EI_mild) <- c(N_age,trans_classes)
dim(n_EI_ILI) <- c(N_age,trans_classes)

#vectors handling I_ILI to R, I_comm_D transition
dim(n_ILI_to_comm_D) <- c(N_age,trans_classes)
dim(n_ILI_to_R) <- c(N_age,trans_classes)

#Vectors handling number of new hospitalisations, ICU admissions and recoveries in hospital
dim(n_ILI_to_hosp) <- c(N_age,trans_classes)
dim(n_ILI_to_triage) <- c(N_age,trans_classes)
dim(n_hosp_non_ICU) <- c(N_age,trans_classes)
dim(n_ILI_to_hosp_D) <- c(N_age,trans_classes)
dim(n_ILI_to_hosp_D_conf) <- c(N_age,trans_classes)
dim(n_ILI_to_hosp_R) <- c(N_age,trans_classes)
dim(n_ILI_to_hosp_R_conf) <- c(N_age,trans_classes)
dim(n_ILI_to_triage_R) <- c(N_age,trans_classes)
dim(n_ILI_to_triage_R_conf) <- c(N_age,trans_classes)
dim(n_ILI_to_triage_D) <- c(N_age,trans_classes)
dim(n_ILI_to_triage_D_conf) <- c(N_age,trans_classes)

#Vectors handling the serology flow
dim(n_com_to_R_pre) <- c(N_age,2)

#Vectors handling the severity profile
dim(p_asympt) <- c(N_age)
dim(p_sympt_ILI) <- c(N_age)

#Vectors handling the potential death in hospital (general beds and ICU)
dim(p_death_hosp_D) <- c(N_age)
dim(p_death_ICU) <- c(N_age)

#Vector handling the probability of being admitted as confirmed
dim(p_admit_conf) <- c(N_age)

dim(cum_admit_by_age) <- c(N_age)
dim(cum0_admit_by_age) <- c(N_age)

#Vectors handling the age specific heterogeneous transmission process
dim(lambda) <- N_age
dim(s_ij) <- c(N_age,N_age)
dim(m) <- c(N_age,N_age)
dim(trans_profile) <- c(N_age,trans_classes)
dim(trans_increase) <- c(N_age,trans_classes)
dim(I_with_diff_trans) <- c(N_age,trans_classes)


#Total population
initial(N_tot[]) <- N0_tot[i]
N0_tot[] <- user()
update(N_tot[]) <- S[i] + R[i] + D_hosp[i] + sum(E[i,,]) + sum(I_asympt[i,,]) + sum(I_mild[i,,]) + sum(I_ILI[i,,]) +
  sum(I_triage_D_conf[i,,]) + sum(I_triage_D_unconf[i,,]) + sum(I_triage_R_conf[i,,]) + sum(I_triage_R_unconf[i,,])  + sum(I_hosp_R_conf[i,,]) + sum(I_hosp_R_unconf[i,,]) +
  sum(I_hosp_D_conf[i,,]) + sum(I_hosp_D_unconf[i,,]) + sum(I_ICU_R_conf[i,,]) + sum(I_ICU_R_unconf[i,,]) + sum(I_ICU_D_conf[i,,]) + sum(I_ICU_D_unconf[i,,]) +
  sum(R_stepdown_conf[i,]) + sum(R_stepdown_unconf[i,]) + sum(I_comm_D[i,,]) + D_comm[i]
dim(N_tot) <- c(N_age)
dim(N0_tot) <- c(N_age)

#Total population calculated with seroconversion flow, exclude triage_R, ICU_R, hosp_R and stepdown, to avoid double counting with R's
initial(N_tot2) <- N0_tot2
N0_tot2 <- user()
update(N_tot2) <- sum(S) + sum(R_pre) + sum(R_pos) + sum(R_neg) + sum(E)

#Tracker of population size
#dim(N) <- N_age
