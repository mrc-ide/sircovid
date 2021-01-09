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
