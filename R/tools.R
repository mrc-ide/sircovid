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
##' @title Inflate strains in model state
##' @param state1 The state matrix with one strain
##'
##' @param info1 The model info with one strain
##'
##' @param info2 The model info with two strains
##'
##' @return An expanded model state with two strains
##'
##' @author Richard Fitzjohn
##' @export
inflate_state_strains <- function(state1, info1, info2) {

  fn <- function(state1, info1, info2, d1, d2, i1, i2, ny, x1, x2, nm) {
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

    x2
  }

  inflate_state(state1, info1, info2, fn)
}


##' Inflate model state, so that additional vacc class compartments are
##' created but zerod. Use this to run the model for a while, then set
##' it up to work with new vacc classes (e.g. booster doses).
##'
##' @title Inflate vacc classes in model state
##'
##' @param state1 The state matrix with X vacc classes
##' @param info1 The model info with X vacc classes
##' @param info2 The model info with Y > X vacc classes
##'
##' @return An expanded model state with Y > X vacc classes
##'
##' @export
inflate_state_vacc_classes <- function(state1, info1, info2) {

  fn <- function(state1, info1, info2, d1, d2, i1, i2, ny, x1, x2, nm) {
    if (length(d2) == 2) {
      x2[, seq_len(d1[2L]), ] <- x1
    } else if (length(d2) == 3) {
      x2[, , seq(d1[3L]), ] <- x1
    } else if (length(d2) == 4) {
      x2[, , , seq(d1[4L]), ] <- x1
    } else {
      ## This will trigger if someone adds a vacc-including
      ## compartment with rank 1 or 4+
      stop("Unexpected dimension of output") # nocov
    }

    x2
  }

  inflate_state(state1, info1, info2, fn)
}


inflate_state <- function(state1, info1, info2, fn) {
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
      state2[i2, ] <- fn(state1, info1, info2, d1, d2, i1, i2, ny, x1, x2, nm)
    } else {
      state2[i2, ] <- state1[i1, ]
    }
  }

  state2
}

##' List the models available in this package
##'
##' @title Available sircovid Models
##' @export
sircovid_models <- function() {
  c("basic", "lancelot")
}


##' Check that a model type is valid (i.e. in [sircovid_models])
##'
##' @title Check a model type
##' @param sircovid_model String
##' @export
check_sircovid_model <- function(sircovid_model) {
  assert_scalar(sircovid_model)
  if (!(sircovid_model %in% sircovid_models())) {
    stop(sprintf("Expected '%s' to be one of %s",
                 sircovid_model,
                 sprintf("{%s}", paste0(sircovid_models(), collapse = ", "))),
         call. = FALSE)
  }
}


##' Rotate strains, so that strain 1 becomes the sum of strains 1 and
##' 2 and strain 2 is empty. Use this to allow sequential replacement
##' of strains.
##'
##' @title Rotate strains
##'
##' @param state Model state
##'
##' @param info Model info
##' @export
rotate_strains <- function(state, info) {
  if (is.null(dim(state))) {
    stop("Expected a matrix or array for 'state'")
  }
  if (nrow(state) != info$len) {
    stop(sprintf("Expected a matrix with %d rows for 'state'",
                 info$len))
  }

  ## Push all state down into a common rank to avoid lots of boring
  ## dimension arithmetic (e.g., we might have this structured by
  ## region, or a single particle's state etc, but what matters is the
  ## structure by particle only).
  dim_orig <- dim(state)
  if (length(dim(state)) > 2) {
    state <- matrix(state, nrow(state))
  }

  ## Currently we only support one sort of move: Anyone who has been
  ## infected in the past is moved to now be indexed by strain 1 (no
  ## matter whether they had multiple infections before), and strain 2
  ## will be completely empty.
  strain_from_idx <- c(2, 3, 4)
  strain_to_idx <- 1

  n_particle <- ncol(state)
  for (i in seq_along(rotate_strain_compartments)) {
    name <- rotate_strain_compartments[[i]]
    dim <- info$dim[[name]]
    state_i <- array(state[info$index[[name]], ], c(dim, n_particle))

    if (length(dim) == 4) {
      for (j in strain_from_idx) {
        tomove <- state_i[, j, , , , drop = FALSE]
        state_i[, strain_to_idx, , , ] <-
          state_i[, strain_to_idx, , , , drop = FALSE] + tomove
        state_i[, j, , , ] <- 0
      }
    } else if (length(dim) == 3) {
      for (j in strain_from_idx) {
        tomove <- state_i[, j, , , drop = FALSE]
        state_i[, strain_to_idx, , ] <-
          state_i[, strain_to_idx, , , drop = FALSE] + tomove
        state_i[, j, , ] <- 0
      }
    } else if (length(dim) == 1) {
      ## this loop range ensures that the move can still happen when
      ## the object has dimension n_real_strains not n_strains
      ## e.g. for prob_strain
      stopifnot(dim %in% c(2, 4))
      for (j in strain_from_idx[strain_from_idx <= dim]) {
        tomove <- state_i[j, , drop = FALSE]
        state_i[strain_to_idx, ] <-
          state_i[strain_to_idx, , drop = FALSE] + tomove
        state_i[j, ] <- 0
      }
    } else {
      ## This is unreachable unless the model changes to include
      ## something that has a rank-2 variable that needs transforming.
      stop(sprintf("Unexpected dimensions (%d) in rotate_strain", # nocov
                   length(dim)))                                  # nocov
    }
    state[info$index[[name]], ] <- state_i
  }

  dim(state) <- dim_orig
  state
}


rotate_strain_compartments <- c(
  ## those with dimension c(n_groups, n_strains, k_XXX, n_vacc_classes):
  "E", "I_A", "I_P", "I_C_1", "I_C_2", "G_D",
  "ICU_pre_unconf", "H_R_unconf", "H_D_unconf", "ICU_W_R_unconf",
  "ICU_W_D_unconf", "ICU_D_unconf", "W_R_unconf", "W_D_unconf",
  "ICU_pre_conf", "H_R_conf", "H_D_conf", "ICU_W_R_conf",
  "ICU_W_D_conf", "ICU_D_conf", "W_R_conf", "W_D_conf",
  "T_sero_pre_1", "T_sero_pos_1",
  "T_sero_pre_2", "T_sero_pos_2",
  "T_PCR_pre", "T_PCR_pos",
  ## those with dimension c(n_groups, n_strains, n_vacc_classes):
  "R", "T_sero_neg_1", "T_sero_neg_2", "T_PCR_neg", "I_weighted",
  ## those with dimension n_strains:
  "cum_infections_per_strain",
  ## those with dimension n_real_strains:
  "prob_strain")
