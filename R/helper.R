

#' Generate Fibonacci sequence
#'
#' @param from,to Approximate start and end values of the sequence
#'
#' @export
#'
#' @examples
#'
#' fibonacci(2, 200)
#'

fibonacci <- function(from, to) {

  ## first determine the start and end
  f1 <- 1; f2 <- 1 # start seq
  while (f1 < from) {
    newf <- f1 + f2; f1 <- f2; f2 <- newf
  }
  sf1 <- f1; sf2 <- f2 # save start
  len <- 0 # start counter to track to determine sequenc length
  while (f1 <= to) {
    newf <- f1 + f2; f1 <- f2; f2 <- newf; len <- len + 1
  }

  ## now build the sequence
  f <- numeric(len)
  f[1] <- sf1
  f[2] <- sf2
  for (i in 3:len) {
    f[i] <- f[i - 1] + f[i - 2]
  }
  f

}


#' Print object size
#'
#' @description A wrapper for \code{\link[base]{object.size}} that prints in Mb
#' by default
#'
#' @param x      an \code{R} object
#' @param units  the units to be used in printing the size
#'
#' @export
#'

object_size <- function(x, units = "Mb") {
  format(utils::object.size(x), units = units)
}

