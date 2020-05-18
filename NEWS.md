# sircovid 0.2.33

- 3rd model added to package - serology_model

# sircovid 0.2.33

- Fix bug in miscalculation of R compartment in new_hospital_model

# sircovid 0.2.32

- Fix bug in update_beta()

# sircovid 0.2.31

- Add beta post_lockdown

# sircovid 0.2.30

- Allow beta_end > beta_start or beta_reduction > 1 in generate_beta

# sircovid 0.2.29

- generalised pMCMC to sample gammas
- rename cov_mat -> proposal_kernel
- Working sample_pmcmc (bug fix)

# sircovid 0.2.28

- remove pars_obs and model_params as default arguments
- fix axis naming bug in plot.pmcmc_list

# sircovid 0.2.27

- Sampler/forecasts for pMCMC results

# sircovid 0.2.26

- plot method for parallel pMCMC

# sircovid 0.2.25

- Three parameter pMCMC

# sircovid 0.2.24

- Working pMCMC (bug fix)

# sircovid 0.2.23

- Add function for variable seeding of the epidemic

# sircovid 0.2.22

- Check grid search boundary is set correctly

# sircovid 0.2.21

- Add date stamp to plots
- Add prior to beta grid search
- Add summary function to sample grid results

# sircovid 0.2.20

- Multivariate normal proposals for pMCMC

# sircovid 0.2.19

- Fix name error with p_death_hosp/_D in new odin model

# sircovid 0.2.18

- Fix dates used in pMCMC

# sircovid 0.2.17

- Options for using different models reimplemented
- Changing beta generalised

# sircovid 0.2.16

- Generate parameter function for new_hospital_model using fitted parameters from hospital model

# sircovid 0.2.15

- Generate parameter function for new_hospital_model

# sircovid 0.2.14

- 2nd model added to package - new_hospital_model

# sircovid 0.2.13

- Adds a save_end_states option to particle_filter

# sircovid 0.2.12

- When doing grid sample and projections plot deaths counts rather than cumulative

# sircovid 0.2.11

- Adds pMCMC for beta and start date inference

# sircovid 0.2.10

- Use a time-varying beta in the grid search
- Fix an issue with large memory usage in the grid search and sampling
- Make the run_particle_filter() function more portable

# sircovid 0.2.9

- Fixes a bug with time scaling of beta macthing coding in odin

# sircovid 0.2.8

- sample_grid_scan can now accept forecast_days

# sircovid 0.2.7

- Adds sample_grid_scan for sampling from a grid to produce trajectories
- Plotting function for sample_grid_search output

# sircovid 0.2.6

- Adds a new main, exported, sircovid function, which fits to data using a time-varying beta

# sircovid 0.2.5

- Adds the 'last' option to save_particles 

# sircovid 0.2.4

- Fixes a small bug with non-integer population size

# sircovid 0.2.3

- Add data, model_params and pars_obs to grid search return

# sircovid 0.2.2

- Fixing issues in generate_parameters and changed default parameters

# sircovid 0.2.1

- Adds time varying beta
- Adds parameter grid search function

# sircovid 0.2.0

* Refactor the particle filter and parameter handling

# sircovid 0.1.3

* Facilitate passing the default population structure into the model

# sircovid 0.1.2

* allows for chosing a small time-step in the model, argument passed using "dt" in days

# sircovid 0.1.1

* corrects a bug in the number of people leaving the ICU compartment

# sircovid 0.1.0

* initial version of the model

