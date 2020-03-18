## f <- function(i) {
library(sircovid)

survey_pop <- readRDS(file="~/Documents/survey_pop.Rds")

pars_model <- parameters(
  beta = 0.042,
  survey_pop = survey_pop)

if(1==0)
  pars_model$p_asympt <- rep(1,9)  #all asymptomatic

if(1==0){
  pars_model$p_asympt <- rep(0,9)  #all symptomatic
  pars_model$p_sympt_ILI <- rep(0,9)  #all no fever
}

if(1==0){  
  pars_model$p_asympt <- rep(0,9)  #all symptomatic
  pars_model$p_sympt_ILI <- rep(1,9)  #all with fever
}

if(1==0){  
  pars_model$p_asympt <- rep(0,9)  #all symptomatic
  pars_model$p_sympt_ILI <- rep(1,9)  #all with fever
}

if(1==0){  
  pars_model$p_asympt <- rep(0,9)  #all symptomatic
  pars_model$p_sympt_ILI <- rep(1,9)  #all with fever
  pars_model$p_recov_ILI <- rep(1,9) #all recover from ILI, no hospitalisation
}

if(1==1){  
  pars_model$p_asympt <- rep(0,9)  #all symptomatic
  pars_model$p_sympt_ILI <- rep(1,9)  #all with fever
  pars_model$p_recov_ILI <- rep(0,9) #none recover from ILI, all cases hospitalised
  pars_model$p_recov_hosp <- rep(1,9)
}

set.seed(1)
mod <- basic(user = pars_model)
t <- seq(from=1, to=2, by = 1)
set.seed(1)
tmp <- mod$run(t)
results <- mod$transform_variables(tmp)
tt <- mod$contents()
temp_tt <- tt$delta_E + tt$delta_I_asympt + tt$delta_I_hosp + tt$delta_I_ICU + tt$delta_I_mild +
  tt$delta_I_ILI + tt$delta_R_hosp
(delta_N <- - tt$n_SE + tt$delta_D + tt$delta_R + apply(X=temp_tt, MARGIN=1, FUN=sum))

v <- c("S", "R", "D", "E", "I_asympt", "I_mild", "I_ILI", "I_hosp", "I_ICU", "R_hosp")
res <- lapply(v, function(x) apply(results[[x]], 1, sum))
names(res) <- v
res <- data.frame(res)

print(delta_N)
print(res)
