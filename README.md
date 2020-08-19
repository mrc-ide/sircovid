# sircovid2 <img src='man/figures/logo.png' align="right" height="138.5" />

<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.com/mrc-ide/sircovid2.svg?branch=master)](https://travis-ci.com/mrc-ide/sircovid2)
[![CodeFactor](https://www.codefactor.io/repository/github/mrc-ide/sircovid2/badge)](https://www.codefactor.io/repository/github/mrc-ide/sircovid2)
[![codecov.io](https://codecov.io/github/mrc-ide/sircovid2/coverage.svg?branch=master)](https://codecov.io/github/mrc-ide/sircovid2?branch=master)
<!-- badges: end -->

`sircovid2` implements a series of models to help modelling the transmission of the SARS-Cov-2 virus using a stochastic compartmental model. `sircovid2` also provides some tools to perfom evidence synthesis from several surveillance data streams in order infer transmission parameters.

<img src="man/figures/sircovid_diagram.png" align="center" style = "border: none; float: center;" width = "800px">

## Installation

Install from the ncov drat:

```r
drat:::add("ncov-ic")
install.packages("sircovid2")
```

or install directly from GitHub with:

```r
remotes::install_github("mrc-ide/sircovid2", upgrade = FALSE)
```

You will need the most recent version of [`dust`](https://mrc-ide.github.io/dust) and [`mcstate`](https://mrc-ide.github.io/mcstate) to use the package.  These will be installed automatically if you install `sircovid2` from drat, or manually with:

```r
drat:::add("ncov-ic")
install.packages(c("dust", "mcstate"))
```

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
