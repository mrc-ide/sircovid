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
  if (nrow(state_orig) != info_orig$len) {
    stop(sprintf("Expected a matrix with %d rows for 'state_orig'",
                 info_orig$len))
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
  allowed <- allowed %||% c("cum_n_vaccinated", "D",
                            "cum_infections_disag", "diagnoses_admitted")
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


##' Inflate model state, so that additional strain compartments are
##' created but zerod. Use this to run the model for a while, then set
##' it up to work with new strains.
##'
##' @title Inflate model state
##' @param state1 The state matrix with one strain
##'
##' @param info1 The model info with one strain
##'
##' @param info2 The model info with two strains
##'
##' @return An expanded model state with two strains
##'
##' @author Richard Fitzjohn
inflate_state_strains <- function(state1, info1, info2) {
  if (!is.matrix(state1)) {
    stop("Expected a matrix for 'state1'")
  }
  if (nrow(state1) != info1$len) {
    stop(sprintf("Expected a matrix with %d rows for 'state1'",
                 info1$len))
  }

  if (!setequal(names(info1$index), names(info2$index))) {
    stop("Can't inflate state (try upgrading first)")
  }

  ny <- ncol(state1)
  state2 <- matrix(0.0, info2$len, ny)
  for (nm in names(info1$index)) {
    d1 <- info1$dim[[nm]]
    d2 <- info2$dim[[nm]]
    i1 <- info1$index[[nm]]
    i2 <- info2$index[[nm]]

    if (!identical(d1, d2)) {
      x1 <- state1[i1, ]
      dim(x1) <- c(d1, ny)
      x2 <- array(0, c(d2, ny))

      if (length(d1) == 1) {
        x2[1, ] <- x1
      } else if (length(d2) == 3) {
        x2[, 1, , ] <- x1
      } else if (length(d2) == 4) {
        x2[, 1, , , ] <- x1
      } else {
        ## This will trigger if someone adds a strain-including
        ## compartment with rank 2 or 5+
        stop("Unexpected dimension of output") # nocov
      }
      state2[i2, ] <- x2
    } else {
      state2[i2, ] <- state1[i1, ]
    }
  }

  state2
}
