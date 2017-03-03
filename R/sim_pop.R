
## TODO ------------------------------------------------------------------------
##
## - Work on sim_distribution function
## - Create more complex mortality simulation (autocovariance between years and
##   age / separate Z and M effects / etc.)
## - Consider random walk over time to help build a longer space-size-time series
## - Look into building the simulation using R-INLA (see Chapter 8 of book)
## - Perhaps sim_abundance should be based on SAM formulation?
##
##


#' Simulate random recruitment and natural mortality
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
#' @details Both functions simply generate uncorrelated recruitment or mortality
#' values from a log normal distribution.
#'
#' @examples
#' sim_abundance(R = sim_R(mean = 100000, sd = 4))
#' sim_abundance(Z = sim_Z(mean = 0.6, sd = 1.5))
#' sim_abundance(Z = sim_Z(mean = list(ages = c(0.5, 0.3, 0.2)),
#'                       breaks = list(ages = c(1, 2))))
#' sim_abundance(Z = sim_Z(mean = list(ages = c(0.5, 0.3, 0.2), years = c(0.3, 2, 0.3)),
#'                         breaks = list(ages = c(1, 2), years = c(4, 6))))
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
sim_Z <- function(mean = 0.4, sd = 1.1, breaks = NULL) {
  function(ages = NULL, years = NULL) {
    na <- length(ages)
    ny <- length(years)
    if (!is.null(breaks)) {
      mean_ages <- mean$ages[findInterval(ages, breaks$ages, left.open = TRUE) + 1]
      mean_years <- mean$years[findInterval(years, breaks$years, left.open = TRUE) + 1]
      if (is.null(mean_ages)) { mean_ages <- rep(0, na) }
      if (is.null(mean_years)) { mean_years <- rep(0, ny) }
      mean <- rowSums(expand.grid(mean_ages, mean_years))
    } else {
      mean <- rep(mean, na * ny)
    }
    z <- rlnorm(na * ny, meanlog = log(mean), sdlog = log(sd))
    z <- matrix(z, nrow = na, ncol = ny)
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

#' Simulate exponential or Matern covariance
#'
#' @description These functions return a function to use inside \code{\link{sim_distribution}}.
#'
#' @param range Decorrelation range
#' @param variance Spatial variance
#' @param model String indicating either "exponential" or "matern" as the correlation function
#'
#' @rdname sim_covar
#' @export
sim_covar <- function(range = NULL, variance = 1, model = "exponential") {
  function(x = NULL) {
    d <- .dist(x)
    cormat <- switch(model,
                     exponential = {
                       exp(- d / range)
                     },
                     matern = {
                       lambda <- 1 # lambda fixed to 1 as per R-INLA book
                       kappa <- sqrt(8) / range # approximate kappa from range as per R-INLA book
                       2 ^ (1 - lambda) / gamma(lambda) * (kappa * d) ^ lambda *
                         besselK(x = d * kappa, nu = lambda)
                     },
                     stop("wrong or no specification of model"))

    diag(cormat) <- 1
    covar <- variance * cormat
    covar
    #Matrix::Matrix(solve(covar), sparse = TRUE)
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

sim_distribution <- function(pop = sim_abundance(),
                             grid = survey_grid,
                             size_covar  = sim_covar(range = 4),
                             time_covar  = sim_covar(range = 4),
                             space_covar = sim_covar(range = 100, variance = 10,
                                                     model = "matern")) {




  ## you are here. download space-time chapter of R-INLA for more inspiration
  ## Think about AR1, 2, 3, N process. Think about independent realizations of the
  ## Space-time error across sizes

  ages <- 1:5
  xy <- expand.grid(x = 1:40, y = 1:40)
  rownames(xy) <- seq(nrow(xy))
  Sigma_size <- size_covar(ages)
  Sigma_space <- space_covar(xy)
  I_size <- diag(1, length(ages))
  I_space <- diag(1, nrow(xy))
  Sigma_size[1:3, 1:3] <- diag(1, 3)
  rownames(Sigma_size) <- rownames(I_size) <- colnames(Sigma_size) <- colnames(I_size) <- ages
  rownames(Sigma_space) <- rownames(I_space) <- colnames(Sigma_space) <- colnames(I_space) <- seq(nrow(xy))
  test <- kronecker(Sigma_size, Sigma_space, make.dimnames = TRUE)
  e <- t(chol(test)) %*% rnorm(nrow(test))
  xyz <- data.frame(do.call(rbind, strsplit(rownames(e), ":")), e)
  names(xyz) <- c("age", "cell", "e")
  xyz <- data.frame(xyz, xy[as.character(xyz$cell), ])

  plot(rasterFromXYZ(xyz[xyz$age == 1, c("x", "y", "e")]))
  plot(rasterFromXYZ(xyz[xyz$age == 2, c("x", "y", "e")]))
  plot(rasterFromXYZ(xyz[xyz$age == 3, c("x", "y", "e")]))
  plot(rasterFromXYZ(xyz[xyz$age == 4, c("x", "y", "e")]))
  plot(rasterFromXYZ(xyz[xyz$age == 5, c("x", "y", "e")]))

  plot(xyz[xyz$cell == 1, ]$e, type = "l")




  ## Random-walk error
  xyz <- data.frame(rasterToPoints(grid))
  xy <- xy[, c("x", "y")]
  Sigma_space <- space_covar(xy)
  w <- t(chol(Sigma_space))
  emat <- replicate(n = length(pop$ages),
                    {emat <- replicate(n = length(pop$years), w %*% rnorm(nrow(xy)),
                    simplify = "matrix")
                    t(apply(emat, 1, cumsum))},
                    simplify = "array")

  for(i in seq_along(pop$years)) {
    xy$e <- emat[, i, 2]
    plot(rasterFromXYZ(xy[, c("x", "y", "e")]), main = i)
  }

  for(i in seq_along(pop$ages)) {
    xy$e <- emat[, 1, i]
    plot(rasterFromXYZ(xy[, c("x", "y", "e")]), main = i)
  }


  ## Space-time-size random walk process
  e <- replicate(n = length(pop$years) * length(pop$ages),
                    w %*% rnorm(nrow(xy)), simplify = FALSE)
  e2 <- array(unlist(emat), dim = c(nrow(xy), length(pop$years), length(pop$ages)),
              dimnames = list(cell = xyz$cell, year = pop$years, age = pop$ages))
  e3 <- apply(e2, c(3, 1), cumsum)
  e4 <- apply(e2, c(2, 1), cumsum)
  e5 <- apply(e3, c(1, 3), cumsum)
  for(i in seq_along(pop$ages)) {
    xy$e <- e3[1, i, ]
    plot(rasterFromXYZ(xy[, c("x", "y", "e")]), main = i)
  }
  for(i in seq_along(pop$years)) {
    xy$e <- e3[i, 2, ]
    plot(rasterFromXYZ(xy[, c("x", "y", "e")]), main = i)
  }
  for(i in seq_along(pop$ages)) {
    xy$e <- e4[i, 1, ]
    plot(rasterFromXYZ(xy[, c("x", "y", "e")]), main = i)
  }
  for(i in seq_along(pop$years)) {
    xy$e <- e4[1, i, ]
    plot(rasterFromXYZ(xy[, c("x", "y", "e")]), main = i)
  }
  for(i in seq_along(pop$ages)) {
    xy$e <- e5[i, 1, ]
    plot(rasterFromXYZ(xy[, c("x", "y", "e")]), main = i)
  }
  for(i in seq_along(pop$years)) {
    xy$e <- e5[1, i, ]
    plot(rasterFromXYZ(xy[, c("x", "y", "e")]), main = i)
  }







  #xy$e <- solve(chol(sigma_space)) %*% rnorm(nrow(xy))
  xy$e <- t(chol(sigma_space)) %*% rnorm(nrow(xy))
  plot(rasterFromXYZ(xy[, c("x", "y", "e")]))

  sigma_space <- solve(sigma_space)

  Sigma <- kronecker(sigma_space, sigma_size)
  Sigma[1:10, 1:10]
  temp <- chol(Sigma)





  ## Distribute abundance equally through the domaine
  N_array <- array(dim = c(nrow(N), ncol(N), nrow(grid)),
                   dimnames = list(age = rownames(N), year = colnames(N), cell = grid$cell))
  prop <- grid$area / sum(grid$area)
  for(i in seq_along(grid$cell)) {
    N_array[, , i] <- (N/grid$area[i]) * prop[i] * exp(e[i])
  }




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




