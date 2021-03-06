`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}


squote <- function(x) {
  sprintf("'%s'", x)
}


verify_names <- function(x, required = NULL, optional = NULL,
                         allow_extra = FALSE,
                         name = deparse(substitute(x))) {
  nms <- names(x)
  if (anyDuplicated(nms)) {
    dups <- unique(nms[duplicated(nms)])
    stop(sprintf("Duplicate element names in '%s': %s",
                 name, paste(squote(dups), collapse = ", ")))
  }
  if (!allow_extra) {
    extra <- setdiff(nms, c(required, optional))
    if (length(extra) > 0) {
      stop(sprintf("Extra elements in '%s': %s",
                   name, paste(squote(extra), collapse = ", ")))
    }
  }
  msg <- setdiff(required, nms)
  if (length(msg) > 0) {
    stop(sprintf("Elements missing from '%s': %s",
                 name, paste(squote(msg), collapse = ", ")))
  }
  invisible(x)
}


sircovid_file <- function(...) {
  system.file(..., package = "sircovid", mustWork = TRUE)
}


read_csv <- function(...) {
  utils::read.csv(..., check.names = FALSE, stringsAsFactors = FALSE)
}


write_csv <- function(...) {
  utils::write.csv(..., row.names = FALSE)
}


data_frame <- function(...) {
  data.frame(..., check.names = FALSE, stringsAsFactors = FALSE)
}


is_integer <- function(x, tol = sqrt(.Machine$double.eps)) {
  abs(x - round(x)) < tol
}


rename <- function(x, from, to, name = deparse(substitute(x))) {
  verify_names(x, required = from, allow_extra = TRUE, name = name)
  i <- match(from, names(x))
  names(x)[i] <- to
  x
}


set_names <- function(x, nms) {
  names(x) <- nms
  x
}


vlapply <- function(.x, .fun, ...) {
  vapply(.x, .fun, logical(1), ...)
}


vnapply <- function(.x, .fun, ...) {
  vapply(.x, .fun, numeric(1), ...)
}


vcapply <- function(.x, .fun, ...) {
  vapply(.x, .fun, character(1), ...)
}


abind1 <- function(a, b) {
  na <- dim(a)[1]
  nb <- dim(b)[1]
  nab <- dim(a)[2:3]
  ret <- array(NA_real_, c(na + nb, nab))
  ret[seq_len(na), , ] <- a
  ret[seq_len(nb) + na, , ] <- b
  rownames(ret) <- c(rownames(a), rownames(b))
  ret
}


abind2 <- function(a, b) {
  na <- dim(a)[2]
  nb <- dim(b)[2]
  nab <- dim(a)[c(1, 3)]
  ret <- array(NA_real_, c(nab[1L], na + nb, nab[2L]))
  ret[, seq_len(na), ] <- a
  ret[, seq_len(nb) + na, ] <- b
  colnames(ret) <- c(colnames(a), colnames(b))
  ret
}

## calculates the nth power of a matrix
matrix_pow <- function(x, n) {
  if (!is_integer(n)) {
    stop("n must be an integer")
  }
  if (!is.matrix(x)) {
    stop("x must have 2 dimensions")
  }
  if (ncol(x) != nrow(x)) {
    stop("x must be a square matrix")
  }

  if (n < 0) {
    return(matrix_pow(solve(x), -n))
  } else if (n == 0) {
    return(diag(nrow(x)))
  } else if (n == 1) {
    return(x)
  } else if (n %% 2 == 0) {
    return(matrix_pow(x %*% x, n / 2))
  } else {
    return(x %*% matrix_pow(x %*% x, (n - 1) / 2))
  }
}


last <- function(x) {
  x[[length(x)]]
}


spread_integer <- function(n, m) {
  ret <- ceiling(rep(n / m, m))
  i <- seq(to = m, length.out = sum(ret) - n)
  ret[i] <- ret[i] - 1L
  ret
}


seq_range <- function(x) {
  r <- range(x)
  seq(r[[1]], r[[2]])
}


interpolate_grid_x <- function(x, every, min) {
  if (is.null(every) || is.null(min) || length(x) <= min) {
    xx <- x
  } else {
    xx <- seq(x[[1]], last(x), by = min(every, ceiling(length(x) / min)))
  }
  if (last(xx) < last(x)) {
    xx <- c(xx, last(x))
  }
  xx
}


interpolate_grid_critical_x <- function(x, every, critical, min) {
  if (length(critical) == 0) {
    list(interpolate_grid_x(x, every, min))
  } else {
    xx <- split(x, findInterval(x, critical))
    lapply(unname(xx), interpolate_grid_x, every, min)
  }
}


interpolate_grid_expand_y <- function(y, step_split) {
  interpolate <- function(y, step_split) {
    f <- function(i) {
      x <- step_split[[i]]
      stats::spline(x, y_split[[i]], xout = seq_range(x))$y
    }
    n <- lengths(step_split)
    y_split <- Map(function(len, to) y[seq(length.out = len, to = to)],
                  n, cumsum(n))
    unlist(lapply(seq_along(step_split), f))
  }

  if (!is.null(dim(y))) {
    apply(y, 2, interpolate, step_split = step_split)
  } else {
    interpolate(y, step_split)
  }
}


block_expand <- function(m, n) {
  if (n == 1L) {
    return(m)
  }
  len <- nrow(m) * n
  matrix(t(matrix(m, nrow(m), len)), len, len, byrow = TRUE)
}


gamma_mucv2shapescale <- function(mu, cv) {
  shape <- 1 / (cv^2)
  scale <- mu * cv^2
  list(shape = shape, scale = scale)
}


mirror_strain <- function(x) {
  if (inherits(x, c("matrix", "array"))) {
    if (ncol(x) > 2) {
      stop("'mirror_strain' should only be used on length 1 or 2 numerics or
           1- or 2-column arrays")
    } else if (ncol(x) == 1) {
      out <- x
    } else if (length(dim(x)) == 3) {
      out <- array(NA, c(nrow(x), 4, nlayer(x)))
      out[, c(1, 4), ] <- x[, c(1, 1), ]
      out[, c(2, 3), ] <- x[, c(2, 2), ]
    } else {
      out <- matrix(NA, nrow(x), 4)
      out[, c(1, 4)] <- x[, c(1, 1)]
      out[, c(2, 3)] <- x[, c(2, 2)]
    }
  } else {
    if (length(x) > 2 || !is.numeric(x)) {
      stop("'mirror_strain' should only be used on length 1 or 2 numerics or
           1- or 2-column arrays")
    } else if (length(x) == 1) {
      out <- x
    } else {
      out <- c(x, rev(x))
    }
  }

  out
}


unmirror_strain <- function(x, n_strains = 4L) {
  if (inherits(x, c("matrix", "array"))) {
    if (ncol(x) != n_strains) {
      return(x)
    } else if (length(dim(x)) == 2) {
      x[, seq(ncol(x) / 2), drop = FALSE]
    } else {
      x[, seq(ncol(x) / 2), , drop = FALSE]
    }
  } else {
    if (length(x) != n_strains) {
      return(x)
    } else {
      x[seq(length(x) / 2)]
    }
  }
}


nlayer <- function(x) {
  dim(x)[3L]
}

recycle <- function(x, n, name = deparse(substitute(x))) {
  if (length(x) == n) {
    x
  } else if (length(x) == 1L) {
    rep_len(x, n)
  } else {
    stop(sprintf("Invalid length for '%s', expected 1 or %d",
         name, n))
  }
}
