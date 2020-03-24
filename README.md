## sircovid

* Make changes to the model in `inst/odin`
* Run `odin::odin_package(".")` from the root directory, which will generate updated files `R/odin.R` and `src/odin.c` (along with `inst/odin/<modelname>.json`, which you can ignore)

### versions

#### 0.1.2

allows for chosing a small time-step in the model, argument passed using "dt" in days

#### 0.1.1

corrects a bug in the number of people leaving the ICU compartment

#### 0.1.0

initial version of the model

