# sircovid2

<!-- badges: start -->
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![Build Status](https://travis-ci.com/mrc-ide/sircovid2.svg?branch=master)](https://travis-ci.com/mrc-ide/sircovid2)
[![CodeFactor](https://www.codefactor.io/repository/github/mrc-ide/sircovid2/badge)](https://www.codefactor.io/repository/github/mrc-ide/sircovid2)
[![codecov.io](https://codecov.io/github/mrc-ide/sircovid2/coverage.svg?branch=master)](https://codecov.io/github/mrc-ide/sircovid2?branch=master)
<!-- badges: end -->

## Development

* Make changes to the models in `inst/odin`
* Run `odin.dust::odin_dust_package(here::here())` from the root directory, which will generate updated files `R/dust.R` and `src/basic.cpp`, `src/carehomes.cpp`, along with `R/cpp11.R` and `src/cpp11.cpp`

## License

MIT © Imperial College of Science, Technology and Medicine
