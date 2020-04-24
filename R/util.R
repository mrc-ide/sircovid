`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}


sircovid_file <- function(path) {
  system.file(path, package = "sircovid", mustWork = TRUE)
}


# from https://stackoverflow.com/a/14502298
##' Extract elements from an array when its dimensions are not known
##' in advance
##'
##'
##' @param x array
##' @param dimension which dimension to pull
##' @param value along the dimension \code{dimension}, which value to
##' extract.
##' @param drop
##' @return
##' @details As an axample, if the dimensions of x are 3, 4, 5 then
##' \code{index_array(x, dimension = 1, value = 2)} will return
##' x[2, , ]
index_array <- function(x, dimension, value, drop = FALSE) {
  # Create list representing arguments supplied to [
  # bquote() creates an object corresponding to a missing argument
  indices <- rep(list(bquote()), length(dim(x)))
  indices[[dimension]] <- value

  # Generate the call to [
  call <- as.call(c(
    list(as.name("["), quote(x)),
    indices,
    list(drop = drop)
  ))
  # Print it, just to make it easier to see what's going on
  print(call)

  # Finally, evaluate it
  eval(call)
}

##' Check if the boundaries of a matrix are within the specified
##' tolerance of 0.
##'
##'
##' @param array
##' @param tolerance
##' @return TRUE if the each edge of array is within specified toleance
##' of 0, FALSE otherwise
##' @author Sangeeta Bhatia
zero_boundary <- function(array, tolerance) {
  close_enough <- 1
  ndims <- length(dim(array))
  for (dimension in seq_len(ndims)) {
    close_enough <- close_enough &&
      (
        all(index_array(array, dimension, 1) < tolerance &
          all(index_array(array, dimension, dim(array)[dimension]) < tolerance))
      )
  }
  close_enough
}
