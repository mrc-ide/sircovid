sircovid_model <- function(x) {
  if (inherits(x, "sircovid_model")) {
    return(x)
  }
  if (is.null(x) || x == "basic") {
    model <- basic
    compare <- compare_icu
  } else {
    path <- system.file("odin", package = "sircovid", mustWork = TRUE)
    possible <- sub("\\.json$", "", dir(path, pattern = "\\.json$"))
    if (x %in% possible) {
      env <- asNamespace("sircovid")
      model <- get(x, envir = env, mode = "function", inherits = FALSE)
      compare <- get(paste0("compare_", x), envir = env,
                     mode = "function", inherits = FALSE)
    } else {
      stop("Unknown model: ", model)
    }
  }
  ret <- list(model = model, compare = compare)
  class(ret) <- "sircovid_model"
  ret
}
