
## TODO ------------------------------------------------------------------------
##
## - Work on sim_distribution function
## - Create more complex mortality simulation (autocovariance between years and
##   age / separate Z and M effects / etc.)
## - Consider random walk over time to help build a longer space-size-time series
## - Look into building the simulation using R-INLA (see Chapter 8 of book)
## - Perhaps sim_abundance should be based on SAM formulation?
## - Consider adding temperature as a covariate
## - Add ERROR checking
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
                     stop("wrong or no specification of covariance model"))

    diag(cormat) <- 1
    covar <- variance * cormat
    covar
  }
}



#' Define relationships with covariates
#'
#' @discription Simple closures used to define relationships with covariates
#'
#' @param beta,mu,sigma Parameters
#' @param scale Center effect around zero?
#'
#' @rdname linear_rel
#' @export
linear_fun <- function(alpha = 0, beta = NULL, scale = TRUE) {
  function(x = NULL) {
    y <- alpha + beta * x
    if(scale) { y - mean(y) } else { y }
  }
}
#' @rdname linear_rel
#' @export
parabolic_fun <- function(mu = NULL, sigma = NULL, scale = TRUE) {
  function(x = NULL) {
    y <- -(((x - mu)^2) / (2 * sigma ^ 2))
    if(scale) { y - mean(y) } else { y }
  }
}


## Helper function for age and year covariance specification for function below.
## This is probably a clumsy way to deal with this problem...matrix math would be
## more elegant...but this was the most computationally efficient way I could work
## it out.
.break_covar <- function(m, mar, covar) {
  switch(covar,
         ran = unname(m),
         rw = if(mar == "age") {
           unname(apply(m, 2, cumsum))
         } else {
           unname(t(apply(m, 1, cumsum)))
         }
         ,
         ident = if(mar == "age") {
           unname(t(replicate(nrow(m), m[1, ])))
         } else {
           unname(replicate(ncol(m), m[, 1]))
         },
         stop("wrong or no specification of covariance model"))
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
                             space_covar = sim_covar(range = 200, variance = 3,
                                                     model = "matern"),
                             age_covar  = "ran", # "ident", "ran"
                             year_covar  = "rw",
                             age_break = Inf,
                             year_break = Inf,
                             depth_par = parabolic_fun(mu = 200, sigma = 100),
                             scale_error = TRUE
) {

  ## Spatial covariance
  grid_dat <- data.frame(rasterToPoints(grid))
  xy <- grid_dat[, c("x", "y")]
  Sigma_space <- space_covar(xy)
  w <- t(chol(Sigma_space))

  ## Space-time-size random walk process
  e <- replicate(n = length(pop$years) * length(pop$ages),
                 w %*% rnorm(nrow(grid_dat)), simplify = FALSE)
  e <- array(unlist(e), dim = c(nrow(grid_dat), length(pop$ages), length(pop$years)),
             dimnames = list(cell = grid_dat$cell, age = pop$ages, year = pop$years))
  for (i in seq(dim(e)[1])) {
    e[i, , ]
    j <- which(as.numeric(dimnames(e)$age) < age_break)
    k <- seq_along(pop$years)
    e[i, j, k] <- .break_covar(e[i, j, k], "age", age_covar[1])
    if(length(age_covar) == 2) {
      j <- which(as.numeric(dimnames(e)$age) >= age_break)
      e[i, j, k] <- .break_covar(e[i, j, k], "age", age_covar[2])
    }
    j <- seq_along(pop$ages)
    k <- which(as.numeric(dimnames(e)$year) < year_break)
    e[i, , ]
    e[i, j, k] <- .break_covar(e[i, j, k], "year", year_covar[1])
    if(length(year_covar) == 2) {
      k <- which(as.numeric(dimnames(e)$year) >= year_break)
      e[i, j, k] <- .break_covar(e[i, j, k], "year", year_covar[2])
    }
    e[i, , ]
  }
  if(scale_error) { e <- e - mean(e) }


  for(j in seq_along(pop$ages)) {
    xy$e <- e[, j, 2]
    plot(rasterFromXYZ(xy[, c("x", "y", "e")]), main = j)
  }
  for(k in seq_along(pop$years)) {
    xy$e <- e[, 2, k]
    plot(rasterFromXYZ(xy[, c("x", "y", "e")]), main = k)
  }


  ## Distribute abundance equally through the domaine and convert
  ## cell specific abundance to density
  den <- replicate(nrow(grid_dat), pop$N / nrow(grid_dat) / prod(res(grid)))
  den <- aperm(den, c(3, 1, 2)) # match dim of e
  dimnames(den) <- dimnames(e)

  ## Depth covariate (TODO: think up a better way to generate the depth effect array)
  depth <- replicate(n = length(pop$years) * length(pop$ages),
                     depth_par(grid_dat$depth), simplify = FALSE)
  depth <- array(unlist(depth), dim = dim(e),
                 dimnames = dimnames(e))

  ## Apply covariate effect and error to calculate density in each cell
  ## log-normal model
  eta <- log(den) + depth + e
  den <- exp(eta)



  ## NOT WORKING RIGHT YET. FIDDLE WITH PARAMETERS

  ## - use link formulation (eta), then apply log
  ## - model density, and remember to add to log density or place it in the link or something
  ## - center covariate effects on zero (substract mean)


  for(j in seq_along(pop$ages)) {
    xy$n <- den[, j, 2]
    plot(rasterFromXYZ(xy[, c("x", "y", "n")]), main = j)
  }
  for(k in seq_along(pop$years)) {
    xy$n <- den[, 1, k]
    plot(rasterFromXYZ(xy[, c("x", "y", "n")]), main = k)
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

sim_pop <- function(ages = 1:6,
                   grid = survey_grid) {
  print("In progress")
}




