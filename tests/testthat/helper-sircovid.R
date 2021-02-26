dput_str <- function(x, name, prefix, width) {
  txt <- trimws(deparse(x, width.cutoff = width))
  if (last(txt) == ")") {
    txt[length(txt) - 1] <- paste0(txt[length(txt) - 1], ")")
    txt <- txt[-length(txt)]
  }
  pad1 <- strrep(" ", prefix - nchar(name))
  txt[[1]] <- sprintf("%s%s = %s", name, pad1, txt[[1]])
  if (length(txt) > 1) {
    txt[-1] <- paste0(strrep(" ", prefix + 5), txt[-1])
  }
  paste0(strrep(" ", 10), txt, collapse = "\n")
}


dput_named_matrix <- function(m, width = 55) {
  prefix1 <- "   rbind("
  prefix2 <- "         "
  nms <- rownames(m)
  width <- width - max(nchar(nms))
  contents <- paste(vcapply(seq_len(nrow(m)), function(i)
    dput_str(m[i, ], nms[i], max(nchar(nms)), width)),
    collapse = ",\n")
  res <- sprintf("    rbind(%s)", trimws(contents, "left"))
  writeLines(res)
  invisible(res)
}


on_ci <- function() {
  isTRUE(as.logical(Sys.getenv("CI")))
}


on_windows <- function() {
  tolower(Sys.info()[["sysname"]]) == "windows"
}


on_mac <- function() {
  tolower(Sys.info()[["sysname"]]) == "darwin"
}


skip_on_windows_gha <- function() {
  if (on_ci() && on_windows()) {
    testthat::skip("On Windows Github Actions")
  }
}


skip_on_mac_gha <- function() {
  if (on_ci() && on_mac()) {
    testthat::skip("On macOS Github Actions")
  }
}
