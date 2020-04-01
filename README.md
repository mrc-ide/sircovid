## sircovid

<!-- badges: start -->
[![Build Status](https://travis-ci.com/mrc-ide/sircovid.svg?branch=master)](https://travis-ci.com/mrc-ide/sircovid)
[![Codecov test coverage](https://codecov.io/gh/mrc-ide/sircovid/branch/master/graph/badge.svg)](https://codecov.io/gh/mrc-ide/sircovid?branch=master)
<!-- badges: end -->

* Make changes to the model in `inst/odin`
* Run `odin::odin_package(".")` from the root directory, which will generate updated files `R/odin.R` and `src/odin.c` (along with `inst/odin/<modelname>.json`, which you can ignore)

## Installation

Install from the ncov drat:

```
drat:::add("ncov-ic")
install.packages("sircovid")
```

or install directly from GitHub with:

```r
remotes::install_github("mrc-ide/sircovid")
```

You will need the most recent version of `dde`, not yet on CRAN, to use the package.  That will be installed automatically if you install `sircovid` from drat, or manually with:

```r
drat:::add("ncov-ic")
install.packages("dde")
```
