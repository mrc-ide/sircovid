# sircovid <img src='man/figures/logo.png' align="right" height="138.5" />

<!-- badges: start -->
[![Build Status](https://travis-ci.com/mrc-ide/sircovid.svg?branch=master)](https://travis-ci.com/mrc-ide/sircovid)
[![Codecov test coverage](https://codecov.io/gh/mrc-ide/sircovid/branch/master/graph/badge.svg)](https://codecov.io/gh/mrc-ide/sircovid?branch=master)
<!-- badges: end -->

* Make changes to the model in `inst/odin`
* Run `odin::odin_package(".")` from the root directory, which will generate updated files `R/odin.R` and `src/odin.c` (along with `inst/odin/<modelname>.json`, which you can ignore)

sircovid implements a series of models to help modelling the transmission of the SARS-Cov-2 virus using a stochastic compartmental model. Sircovid also provides some tools to perfom evidence synthesis from several surveillance data streams in order infer transmission parameters.

<img src="man/figures/sircovid_diagram.png" align="center" style = "border: none; float: center;" width = "800px">

## Installation

Install from the ncov drat:

```
drat:::add("ncov-ic")
install.packages("sircovid")
```

or install directly from GitHub with:

```r
remotes::install_github("mrc-ide/sircovid", build_vignettes = TRUE)
```

You will need the most recent version of `dde`, not yet on CRAN, to use the package.  That will be installed automatically if you install `sircovid` from drat, or manually with:

```r
drat:::add("ncov-ic")
install.packages("dde")
```

## Documentation
Run:
```
browseVignettes("sircovid")
```

Current vignettes:
- Adding a new model