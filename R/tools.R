##' Upgrade model state
##' @title Upgrade model state after change
##'
##' @param state_orig Old model info
##'
##' @param info_orig Old model info
##'
##' @param info_new New model info
##'
##' @param allowed Character vector of allowed states to add; new
##'   values will be initialised as zeros.
##'
##' @return A new copy of the state
##' @export
upgrade_state <- function(state_orig, info_orig, info_new, allowed = NULL) {
  if (!is.matrix(state_orig)) {
    stop("Expected a matrix for 'state_orig'")
  }

  extra <- setdiff(names(info_orig$index), names(info_new$index))
  if (length(extra) > 0) {
    stop(sprintf("Can't downgrade state (previously had variables %s)",
                 paste(squote(extra), collapse = ", ")))
  }

  msg <- setdiff(names(info_new$index), names(info_orig$index))
  ## We can tolerate any vaccine-related variable that can start
  ## zero'd. This will be the case through brief windows of upgrading
  ## sircovid only (e.g., between sircovid 0.7.2 and 0.8.0)
  allowed <- allowed %||% c("D", "diagnoses_admitted")
  err <- setdiff(msg, allowed)
  if (length(err) > 0) {
    stop(sprintf("Can't remap state (can't add variables %s)",
                 paste(squote(err), collapse = ", ")))
  }

  state_new <- matrix(0.0, info_new$len, ncol(state_orig))
  for (nm in names(info_orig$index)) {
    i_orig <- info_orig$index[[nm]]
    i_new <- info_new$index[[nm]]
    if (length(i_orig) != length(i_new)) {
      stop(sprintf("States are incompatible lengths for '%s'", nm))
    }
    state_new[i_new, ] <- state_orig[i_orig, ]
  }

  state_new
}
