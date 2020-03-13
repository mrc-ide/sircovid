`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}


sircovid_file <- function(path) {
  system.file(path, package = "sircovid", mustWork = TRUE)
}
