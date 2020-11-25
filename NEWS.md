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
