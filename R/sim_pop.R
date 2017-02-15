
#' Simulate recruitment and natural mortality
#'
#' @description These functions return a function to use inside \code{\link{sim_abundance}}.
#' Given parameters, it generates R and Z values.
#'
#' @param mean,sd Mean and standard deviation
#' @param breaks Provide breaks to age and/or year series to specify group specific
#' mean total mortality values. If specified, multiple mean values must be
#' supplied; one more than the number of breaks. Provide named lists (see
#' examples below).
#'
#' @details \code{sim_R} simply generates uncorrelated recruitment values from a log normal
#' distribution. \code{sim_Z} simulates total mortality using a random walk.
#' Random walk errors are first generated for years and ages independently.
#' The outer product of these errors are then calculated to generate a total
#' mortality matrix as follows:
#' \eqn{log(Z_{a, y}) = log(z_{a, y}) + \delta_{a, y}}{log(Z_a,y) = log(z_a,y) + e_a,y},
#' where \eqn{z_{a, y}}{z_a,y} are the supplied mean values.
#'
#' @examples
#' sim_abundance(R = sim_R(mean = 100000, sd = 4))
#' sim_abundance(Z = sim_Z(mean = 0.6, sd = 0.3))
#' sim_abundance(Z = sim_Z(mean = list(ages = c(0.5, 0.3, 0.2)),
#'                       breaks = list(ages = c(1, 2)), sd = 0.3))
#' sim_abundance(Z = sim_Z(mean = list(ages = c(0.5, 0.3, 0.2), years = c(0.3, 2, 0.3)),
#'                       breaks = list(ages = c(1, 2), years = c(4, 6)), sd = 0.3))
#'
#' @export
#' @rdname sim_R
sim_R <- function(mean = 20000, sd = 4) {
  function(years = NULL) {
    r <- rlnorm(length(years), meanlog = log(mean), sdlog = log(sd))
    names(r) <- years
    r
  }
}


#' @export
#' @rdname sim_R
sim_Z <- function(mean = 0.4, sd = 0.2, breaks = NULL) {
  function(ages = NULL, years = NULL) {
    ey <- cumsum(rnorm(length(years), sd = sd)) # random walk error
    ea <- cumsum(rnorm(length(ages), sd = sd))
    e <- outer(ea, ey)
    if (!is.null(breaks)) {
      mean_ages <- mean$ages[findInterval(ages, breaks$ages, left.open = TRUE) + 1]
      mean_years <- mean$years[findInterval(years, breaks$years, left.open = TRUE) + 1]
      if (is.null(mean_ages)) { mean_ages <- rep(1, length(ages)) }
      if (is.null(mean_years)) { mean_years <- rep(1, length(years)) }
      mean <- outer(mean_ages, mean_years)
    }
    z <- exp(log(mean) + e)
    dimnames(z) <- list(age = ages, year = years)
    z
  }
}



#' Simulate basic population dynamics model
#'
#' @param ages Ages to include in the simulation.
#' @param years Years to include in the simulation.
#' @param Z Total mortality function, like \code{\link{sim_Z}}, for generating
#' mortality matrix.
#' @param R Recruitment (i.e. Abundance at \code{min(ages)}) function, like
#' \code{\link{sim_R}}, for generating recruitment vector.
#'
#' @return A \code{list} of length 3:
#' \itemize{
#'   \item{\code{R} - Vector of recruitment values}
#'   \item{\code{Z} - Matrix of total mortality values}
#'   \item{\code{N} - Matrix of abundance values}
#' }
#'
#' @details
#' Abundance from \code{ages[2:max(ages)]} and \code{years[2:max(years)]} is
#' calculated using a standard population dynamics model:
#' \deqn{N_{a, y} = N_{a - 1, y - 1} * exp(-Z_{a - 1, y - 1})}{N_a,y = N_a-1,y-1 * exp(-Z_a-1,y-1)}
#' Abundance at \code{min(ages)} is supplied by \code{r} and abundance at \code{ages[2:max(ages)]} are
#' calculated using the same equation above using \code{Z} values from \code{min(years)}.
#'
#' @examples
#' sim_abundance()
#'
#' @export

sim_abundance <- function(ages = 1:6, years = 1:10, Z = sim_Z(), R = sim_R()) {

  ## Simple error check
  if (any(diff(ages) > 1) | any(diff(years) > 1)) {
    stop("Age and year sequences must be ascending and increment by 1")
  }

  ## Set-up abundance-at-age matrix
  N <- matrix(nrow = length(ages), ncol = length(years),
              dimnames = list(age = ages, year = years))
  Z <- Z(ages = ages, years = years)
  N[1, ] <- R <- R(years = years)

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
  list(ages = ages, years = years, R = R, Z = Z, N = N)

}



## Helper function for euclidian distance calculations
.dist <- function(x) {
  if (requireNamespace("fields", quietly = TRUE)) {
    d <- fields::rdist(x)
  } else {
    d <- as.matrix(dist(x))
  }
}

#' Simulate exponential covariance
#'
#' @description These function returns a function to use inside \code{\link{sim_distribution}}.
#'
#' @param range Decorrelation range (years, age or km)
#' @param psill Partial sill
#'
#' @rdname sim_exp_covar
#' @export
sim_exp_covar <- function(range = NULL, psill = 1) {
  function(x = NULL) {
    d <- .dist(x)
    psill * exp(- d / range)
  }
}




#' Simulate spatial and temporal distribution
#'
#' @description Provided an abundance at age matrix (like one provided by \code{\link{sim_abundance}})
#' and a survey grid (like \code{\link{survey_grid}}) to populate, this function
#' applies correlated space, time and size error to simulate the distribution
#' of the population.
#'
#' @param N An abundance at age matrix with ages defining the rows and years defining
#' the columns (i.e. same structure as a matrix provided by \code{\link{sim_abundance}})
#' @param grid A \code{\link{SpatialPolygonsDataFrame}} defining a regular or irregular
#' grid with the same structure as \code{\link{survey_grid}}
#'
#' @examples
#' sim_distribution()
#'
#' @export
#'

sim_distribution <- function(pop = sim_abundance(),
                             grid = survey_grid,
                             size_covar  = sim_exp_covar(range = 4),
                             time_covar  = sim_exp_covar(range = 4),
                             space_covar = sim_exp_covar(range = 200, psill =5)) {




  xy <- coordinates(grid)
  sigma_space <- space_covar(xy)
  grid$e <- t(chol(sigma_space)) %*% rnorm(nrow(xy))
  plot_sim(grid, zcol = "e")



  ## Distribute abundance equally through the domaine
  N_array <- array(dim = c(nrow(N), ncol(N), nrow(grid)),
                   dimnames = list(age = rownames(N), year = colnames(N), cell = grid$cell))
  prop <- grid$area / sum(grid$area)
  for(i in seq_along(grid$cell)) {
    N_array[, , i] <- (N/grid$area[i]) * prop[i] * exp(e[i])
  }



  grid_dat <- ggplot2::fortify(grid)
  names(grid_dat)[names(grid_dat) == "id"] <- "cell"
  grid_dat$cell <- as.numeric(grid_dat$cell)
  grid@data$N <- N_array[1, 1, ]
  grid@data$den <- grid@data$N/grid@data$area
  grid@data$error <- e
  grid_dat <- merge(grid_dat, grid@data, by = "cell")
  ggplot(grid_dat) +
    geom_polygon(aes(x = long, y = lat, group = cell, fill = N, colour = N))



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




