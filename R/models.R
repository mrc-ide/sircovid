sircovid_model <- function(name) {
  if (is.null(model) || model == "basic") {
    model <- basic
    compare <- compare_icu
  } else {
    path <- system.file("odin", package = "sircovid", mustWork = TRUE)
    possible <- sub("\\.json$", "", dir(path, pattern = "\\.json$"))
    if (name %in% possible) {
      env <- asNamespace("sircovid")
      model <- get(name, envir = env, mode = "function", inherits = FALSE)
      compare <- get(paste0("compare_", name), envir = env,
                     mode = "function", inherits = FALSE)
    } else {
      stop("Unknown model: ", model)
    }
  }
  list(model = model, compare = compare)
}
