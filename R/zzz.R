##' @useDynLib sircovid, .registration = TRUE
NULL

cache <- new.env(parent = emptyenv())

clear_cache <- function() {
  rm(list = ls(all.names = TRUE, envir = cache), envir = cache)
}
