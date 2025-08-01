

#' Generate Fibonacci sequence
#'
#' @param from,to Approximate start and end values of the sequence
#'
#' @export
#'
#' @return Returns a Fibonacci sequence as a vector.
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
#' A wrapper for [`utils::object.size()`] that prints in megabytes (Mb) by default.
#'
#' @param x An R object.
#' @param units The units to be used when printing the size.
#'
#' @return A character string with the object size followed by the unit.
#'
#' @export

object_size <- function(x, units = "Mb") {
  format(utils::object.size(x), units = units)
}

