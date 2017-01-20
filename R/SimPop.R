

#' Simulate basic population dynamics model
#'
#' @param ages Ages to include in the simulation.
#' @param years Years to include in the simulation.
#' @param Z Total mortality. Can be one value or a matrix of \code{length(ages)}
#' rows and \code{length(years)} columns.
#' @param r Recruitment (i.e. Abundance at \code{min(ages)}). Can be one value
#' or vector of \code{length(years)}
#'
#' @return A \code{matrix} of abundance at age by year.
#'
#' @details
#' Abundance from \code{ages[2:max(ages)]} and \code{years[2:max(years)]} is
#' calculated using a standard population dynamics model:
#' \deqn{N_{a, y} = N_{a - 1, y - 1} * exp(-Z_{a - 1, y - 1})}{N_a,y = N_a-1,y-1 * exp(-Z_a-1,y-1)}
#' Abundance at \code{min(ages)} is supplied by \code{r} and abundance at \code{ages[2:max(ages)]} are
#' calculated using the same equation above using \code{Z} values from \code{min(years)}.
#'
#' @examples
#' simAbundance()
#'
#' @export

simAbundance <- function(ages = 1:6, years = 1:10, Z = 0.2, r = 1000) {

  ## Simple error check
  if (any(diff(ages) > 1) | any(diff(years) > 1)) {
    stop("Age and year sequences must be ascending and increment by 1")
  }

  ## Set-up abundance-at-age matrix
  z <- Z # save user supplied z
  N <- Z <- matrix(nrow = length(ages), ncol = length(years),
                   dimnames = list(age = ages, year = years))
  Z[] <- z
  N[1, ] <- r

  ## Fill abundance-at-age matrix
  for (y in seq_along(years)) {
    for (a in seq_along(ages)[-1]) {
      if (y == 1) {
        N[a, 1] <- N[a - 1, 1] * exp(- Z[a - 1, 1])
      } else {
        N[a, y] <- N[a -1, y - 1] * exp(- Z[a -1, y - 1])
      }
    }
  }
  N

}


#' Simulate spatial and temporal distribution
#'
#' @description Provided an abundance at age matrix (like one provided by \code{\link{SimAbundance}})
#' and a survey grid (like \code{\link{survey_grid}}) to populate, this function
#' applies spatial and temporal error to simulate the spatial and temporal distribution
#' of the population.
#'
#' @param N An abundance at age matrix with ages defining the rows and years defining
#' the columns (i.e. same structure as a matrix provided by \code{\link{SimAbundance}})
#' @param grid A \code{\link{SpatialPolygonsDataFrame}} defining a regular or irregular
#' grid with the same structure as \code{\link{survey_grid}}
#'
#' @examples
#' simDistribution()
#'
#' @export
#'

simDistribution <- function(N = simAbundance(),
                            grid = survey_grid) {

  N_array <- array(dim = c(nrow(N), ncol(N), nrow(grid)),
                   dimnames = list(age = rownames(N), year = colnames(N), ))
  grid$area / sum(grid$area)




}





#' Simulate an age-structured population over a survey grid.
#'
#' @param ages Ages to include in the simulated population.
#' @param grid Survey grid to populate with simulated data.
#' @return A 'data.frame' containing ages and lengths for every fish in the
#'     simulated population.
#' @examples
#' simPop()

simPop <- function(ages = 1:6,
                   grid = survey_grid) {
  print("In progress")
}

