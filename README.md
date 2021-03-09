# sircovid <img src='man/figures/logo.png' align="right" height="138.5" />

<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R build status](https://github.com/mrc-ide/sircovid/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/sircovid/actions)
[![CodeFactor](https://www.codefactor.io/repository/github/mrc-ide/sircovid/badge)](https://www.codefactor.io/repository/github/mrc-ide/sircovid)
[![Codecov test coverage](https://codecov.io/gh/mrc-ide/sircovid/branch/master/graph/badge.svg)](https://codecov.io/gh/mrc-ide/sircovid?branch=master)
[![Docker build status](https://badge.buildkite.com/8b8c5742874fc1dc137e5c085f107a1e4346e9cdf65d72934b.svg?branch=master)](https://buildkite.com/mrc-ide/sircovid)
<!-- badges: end -->

`sircovid` implements a series of mechanistic models to help modelling the transmission of the SARS-Cov-2 virus using stochastic compartmental models. `sircovid` also provides some tools to perfom Bayesian evidence synthesis from several surveillance data streams through the estimation of transmission parameters.

<img src="man/figures/sircovid_diagram.png" align="center" style = "border: none; float: center;" width = "800px">

## Important information for users

Please note that whilst this code is free to use and adapt, Imperial College London does not endorse the outputs, results, or conclusions drawn from the implementation of this model to other settings. While we encourage the use and modification of our model for research and scientific purposes please do not refer to such results as the "Imperial model" or similar unless referring to specific use in publications by Imperial College researchers.

## Installation

Install from the ncov drat:

```r
drat:::add("ncov-ic")
install.packages("sircovid")
```

or install directly from GitHub with:

```r
remotes::install_github("mrc-ide/sircovid", upgrade = FALSE)
```

You will need the most recent version of [`dust`](https://mrc-ide.github.io/dust) and [`mcstate`](https://mrc-ide.github.io/mcstate) to use the package.  These will be installed automatically if you install `sircovid` from drat, or manually with:

```r
drat:::add("ncov-ic")
install.packages(c("dust", "mcstate"))
```


We use OpenMP for parallelism, and this may not be available on your system. If not then compilation will fail with an error like:

```r
clang: error: unsupported option '-fopenmp'
make[1]: *** [basic.o] Error 1
ERROR: compilation failed for package ‘sircovid’
```

You can either install OpenMP support, or edit your personal `Makevars` file to tell R that you do not have it. To do this, you can run

```r
usethis::edit_r_makevars("user")
```

and add the lines

```r
SHLIB_OPENMP_CFLAGS=
SHLIB_OPENMP_CXXFLAGS=
```

after which compilation will succeed, but the model will only run on one core.

## Development

In addition to the above you need to install [`odin`](https://mrc-ide.github.io/odin) and [`odin.dust`](https://mrc-ide.github.io/odin.dust/)

```r
drat:::add("ncov-ic")
install.packages(c("odin", "odin.dust"))
```

* Make changes to the models in `inst/odin`
* Run `odin.dust::odin_dust_package(here::here())` from the root directory, which will generate updated files `R/dust.R` and `src/basic.cpp`, `src/carehomes.cpp`, along with `R/cpp11.R` and `src/cpp11.cpp`

Alternatively, run `./scripts/generate_odin`

## License

MIT © Imperial College of Science, Technology and Medicine
