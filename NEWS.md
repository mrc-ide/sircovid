# sircovid 0.10.53

* carehomes_Rt now only returns NA when prob_strain is NA and weight_Rt is TRUE

# sircovid 0.10.52

* Can now input `population` and `carehome_beds` to `carehomes_parameters` to allow use for regions where these are not included in the package

# sircovid 0.10.51

* rel vacc params can now be less than one in unvaccinated classes for second infections (strains 3 and 4)

# sircovid 0.10.50

* Fix Rt calculations for multistrain models with cross-immunity

# sircovid 0.10.49

* Add `lag_groups` and `lag_days` parameters to `vaccine_schedule_future` to manually add lags for given age groups to the vaccine schedule

# sircovid 0.10.47

* New function `sircovid::inflate_state_strains` for adding empty strain compartments to a model state that was run without multiple strains

# sircovid 0.10.46

* Added `carehomes_check_severity` to check all severity probabilities lie in [0, 1] after transformation and all relative parameters are mirrored as expected


# sircovid 0.10.45

* Deaths and admissions by age and vaccination class compartments exported as
predicted samples

# sircovid 0.10.44

* Improve efficiency of Rt calculation for single strains when supplying multistrain parameters

# sircovid 0.10.43

* Fix multistrain Rt calculation for single particle case

# sircovid 0.10.41

* Fix IFR_t calculation for multi-strain models

# sircovid 0.10.40

* New function `sircovid::inflate_state_strains` for adding empty strain compartments to a model state that was run without multiple strains

# sircovid 0.10.39

* R compartment now exported as part of predicted samples

# sircovid 0.10.39

* Added ability to fit 2nd serology assay with different serology parameters

# sircovid 0.10.37

* Added ability to directly set time-varying severity and progression parameters in `carehomes_parameters_progression` and `carehomes_parameters_severity`.
* `sircovid_parameters_beta` renamed to `sircovid_parameters_piecewise_linear` and `sircovid_parameters_beta_expand` renamed to `sircovid_parameters_expand_step`

# sircovid 0.10.36

* Fixed Rt calculation for multistrain so that Rt is NA only in the steps where prob_strain is NA (previously all Rt set to NA if any prob_strain was NA)

# sircovid 0.10.34

* Fixed `carehomes_rt` function to have flexibility for multiple dimensions, given strains, age categories and vaccinations classes.
* test-support.R updated accordingly.

# sircovid 0.10.33

* Can now input a list of observation parameters into the `carehomes_parameters` function instead of taking the default observation parameters.


# sircovid 0.10.31

* Add `cross_immunity` to control the amount of immunity conferred by one strain that decreases the probability of immediately being exposed to the other strain (ignored in single strain model).
* Removed `model_super_infection` as `model_super_infection = 0` is equivalent to `cross_immunity = 1` (default)

# sircovid 0.10.30

* Parameters describing vaccine efficacy (`rel_susceptibility`, `rel_p_sympt`, `rel_p_hosp_if_sympt` and `rel_infectivity`) can now take array values for
varying across age groups, vaccination classes, and pathogen strains.

# sircovid 0.10.29

* Multi-strain model now allows for 'super-infections', after recovering from strain 1 or 2, one can now be immediately infected with the 'other' strain. This can be turned on with `model_super_infection`.

# sircovid 0.10.28

* Relaxed restrictions on waning immunity - the number of recovered individuals losing immunity is no longer capped by the number of individuals in the `T_sero_neg` or `T_PCR_neg` compartments.

# sircovid 0.10.26

* Add parameter `strain_rel_severity` for varying the probabilities of `p_G_D`, `p_H_D`, `p_W_D` and `p_ICU_D` by strain.

# sircovid 0.10.25

* Add additional disaggregated outputs for deaths and diagnosed admissions

# sircovid 0.10.24

* New set of parameter `strain_rel_gamma_` for varying `gamma_A`, `gamma_P`, `gamma_C_1`, and `gamma_C_2` by strain.

# sircovid 0.10.23

* Support to allow catch up of vaccine doses that could not be delivered
(e.g. when people were infected at the time they should have been vaccinated)

# sircovid 0.10.22

* Added support for time-varying hospital durations

# sircovid 0.10.20

* Output cumulative numbers vaccinated

# sircovid 0.10.19

* Changed fixed `strain_seed_value` to stochastic `strain_seed_rate` vector
* Changed `strain_seed_date` from vector of start and end seeding dates to a vector of dates corresponding to rates set by `strain_seed_rate`

# sircovid 0.10.18

* EpiEstim now allows age-varying p_C

# sircovid 0.10.17

* New function `sircovid::compile_gpu` which compiles sircovid to run on a GPU (#237)

# sircovid 0.10.16

* New functions to support reordering samples in mcstate_pmcmc objects


# sircovid 0.10.12

* Fix corner case in vaccination exposed by recent dust updates

# sircovid 0.10.8

* New functions to estimate effective Rt using EpiEstim

# sircovid 0.10.7

* Added utils for nested pmcmc

# sircovid 0.10.5

* Can control the number of initial infected individuals used in seeding

# sircovid 0.10.3

* Predefined vaccination schedule implemented

# sircovid 0.10.1

* Rt corrected for restructured model

# sircovid 0.10.0

* Carehomes model restructured
* Symptomatic pathway now includes presymptomatic compartment and two symptomatic compartments
* Cases for pillar 2 testing arise in the movement from presymptomatic to the first symptomatic compartment
* NHS regions populations updated to ONS mid-2019 estimates

# sircovid 0.9.28

* Output deaths by age

# sircovid 0.9.27

* Can now fit non-hospital deaths as two separate streams: community deaths and care home deaths

# sircovid 0.9.23

* Allow vector of absolute values in `future_Rt()`

# sircovid 0.9.22

* Faster Rt calculation (#209)

# sircovid 0.9.17

* IFR calculation added for carehomes model

# sircovid 0.9.14

* Allow a subset of Rt types to be calculated, perhaps in parallel

# sircovid 0.9.12

* Enables daily vaccine doses to be time-varying

# sircovid 0.9.11

* Fixed Rt calculation for multistrain so that Rt is NA when prob_strain is NA

# sircovid 0.9.7

* Allow specification of Rt calculation in `add_future_betas()`

# sircovid 0.9.4

* Add over 25 strain data stream to likelihood

# sircovid 0.9.3

* Corrected calculation of NGM for multistrain model

# sircovid 0.9.2

* Added different infectivity levels depending on vaccination stage
* Corrected calculation of NGM with vaccination (but not for multistrain model)

# sircovid 0.9.1

* Rename compartments and parameters to match paper notation

# sircovid 0.9.0

* Multistrain model: ability to model 2 strains (and flexibility for more in the future) with strain-specific transmissibility. 2nd strain can be seeded over a number of pre-specified days at a pre-specified seeding level.

# sircovid 0.8.2

* Added functions `add_future_betas` and `future_Rt` for constructing future transmission scenarios
* Added new helper functions `add_trajectory_incidence`, `drop_trajectory_incidence` and `drop_trajectory_predicted`

# sircovid 0.8.0

* Move vaccination process inside odin code

# sircovid 0.7.4

* Now drawing equal numbers of care home workers across worker age groups.

# sircovid 0.7.3

* Correct Rt calculation accounting for vaccination

# sircovid 0.7.2

* I_mild compartment removed from both models and I_ILI renamed to I_sympt

# sircovid 0.7.0

* Extended model structure for vaccination modelling

# sircovid 0.6.10

* Deaths in stepdown now incorporated in carehomes model
* Triage compartments merged in carehomes model

# sircovid 0.6.9

* Implement trajectory and Rt aggregation (e.g., across regions)

# sircovid 0.6.8

* Support for vaccination modelling, with vaccination affecting only levels of susceptibility (age-specific)

# sircovid 0.6.5

* Parameter `N_age` was split into `n_age_groups` (both models) and `n_groups` (carehomes)
* Implement preliminary support for waning immunity

# sircovid 0.6.4

* Implement preliminary support for seroreversion

# sircovid 0.6.3

* Implement preliminary support for vaccination modelling

# sircovid 0.6.2

* Support for probability of positive calculation directly in the package
* Pillar 2 data can now be fitted to over 25s only, as well as all ages
* Sensitivity is explicitly considered in the sero probability calculation

# sircovid 0.6.1

* Add REACT to log-likelihood
