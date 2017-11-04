

#' Simulate startinig abundance, random recruitment and total mortality
#'
#' @description These functions return a function to use inside \code{\link{sim_abundance}}.
#' Given parameters, it generates N0, R and Z values.
#'
#' @param mean One mean value or a vector of means of length equal to years for \code{sim_R} or a matrix of means with
#' rows equaling the number of ages and colums equaling the number of years for \code{sim_Z}.
#' @param sd_log Standard deviation of the variable in the log scale.
#' @param N0 Either specify "exp" or numeric vector of starting abundance excluding the first age.
#' If "exp" is specified using sim_N0, then abundance at age are calculated using exponential decay:
#' \deqn{N_{a, 1} = N_{a - 1, 1} * exp(-Z_{a - 1, 1})}{N_a,1 = N_a-1,1 * exp(-Z_a-1,1)}
#'
#' @details sim_R and sim_Z simply generate uncorrelated recruitment or mortality
#' values from a log normal distribution.
#'
#' @examples
#' sim_abundance(R = sim_R(mean = 100000, sd = 4))
#' sim_abundance(years = 1:20, R = sim_R(mean = c(rep(100000, 10), rep(10000, 10))))
#' sim_abundance(Z = sim_Z(mean = 0.6, sd_log = 0))
#' Za_dev <- c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.3, 0.2, 0.1, 0)
#' Zy_dev <- c(-0.2, -0.2, -0.2, -0.2, -0.2, 2, 2, 2, 2, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0)
#' Z_mat <- outer(Za_dev, Zy_dev, "+") + 0.5
#' sim_abundance(ages = 1:10, years = 1:20, Z = sim_Z(mean = Z_mat), R = sim_R(sd_log = 0.5))
#'
#' @export
#' @rdname sim_R
sim_R <- function(mean = 100000, sd_log = 0.5) {
  function(years = NULL) {
    if (length(mean) > 1 && length(mean) != length(years)) {
      stop("The number of means supplied for recruitment != number of years.")
    }
    r <- rlnorm(length(years), meanlog = log(mean), sdlog = sd_log)
    names(r) <- years
    r
  }
}

#' @export
#' @rdname sim_R
sim_Z <- function(mean = 0.5, sd_log = 0.5) {
  function(ages = NULL, years = NULL) {
    na <- length(ages)
    ny <- length(years)
    if (is.matrix(mean) && (nrow(mean) != na | ncol(mean) != ny)) {
      stop("The matrix of means supplied for Z != number of years and/or ages.")
    } else {
      mean <- c(mean)
    }
    z <- rlnorm(na * ny, meanlog = log(mean), sdlog = sd_log)
    z <- matrix(z, nrow = na, ncol = ny)
    dimnames(z) <- list(age = ages, year = years)
    z
  }
}

#' @export
#' @rdname sim_R
sim_N0 <- function(N0 = "exp") {
  function(R0 = NULL, Z0 = NULL, ages = NULL) {
    if (length(N0) == 1 && N0 == "exp") {
      N0 <- rep(NA, length(ages))
      N0[1] <- R0
      for (a in seq_along(ages)[-1]) {
        N0[a] <- N0[a - 1] * exp(-Z0[a - 1])
      }
    } else{
      N0 <- c(R0, N0)
    }
    N0
  }
}



#' Simulate basic population dynamics model
#'
#' @param ages Ages to include in the simulation.
#' @param years Years to include in the simulation.
#' @param Z Total mortality function, like \code{\link{sim_Z}}, for generating
#' mortality matrix.
#' @param R Recruitment (i.e. abundance at \code{min(ages)}) function, like
#' \code{\link{sim_R}}, for generating recruitment vector.
#' @param N0 Starting abundance (i.e. abundance at \code{min(years)}) function, like
#' \code{\link{sim_N0}}, for generating starting abundance vector.
#'
#' @return A \code{list} of length 3:
#' \itemize{
#'   \item{\code{R} - Vector of recruitment values}
#'   \item{\code{N0} - Vector of starting abundance values}
#'   \item{\code{Z} - Matrix of total mortality values}
#'   \item{\code{N} - Matrix of abundance values}
#' }
#'
#' @details
#' Abundance from \code{ages[2:max(ages)]} and \code{years[2:max(years)]} is
#' calculated using a standard population dynamics model:
#' \deqn{N_{a, y} = N_{a - 1, y - 1} * exp(-Z_{a - 1, y - 1})}{N_a,y = N_a-1,y-1 * exp(-Z_a-1,y-1)}
#'
#' @examples
#' sim_abundance()
#'
#' @export

sim_abundance <- function(ages = 1:12, years = 1:20,
                          Z = sim_Z(), R = sim_R(), N0 = sim_N0()) {

  ## Simple error check
  if (any(diff(ages) > 1) | any(diff(years) > 1)) {
    stop("Age and year sequences must be ascending and increment by 1")
  }

  ## Set-up abundance-at-age matrix
  N <- matrix(nrow = length(ages), ncol = length(years),
              dimnames = list(age = ages, year = years))
  Z <- Z(ages = ages, years = years)
  N[1, ] <- R <- R(years = years)
  N[, 1] <- N0 <- N0(R0 = R[1], Z0 = Z[, 1], ages = ages)

  ## Fill abundance-at-age matrix
  for (y in seq_along(years)[-1]) {
    for (a in seq_along(ages)[-1]) {
      N[a, y] <- N[a - 1, y - 1] * exp(-Z[a - 1, y - 1])
    }
  }

  list(ages = ages, years = years, R = R, N0 = N0, Z = Z, N = N)

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
                       exp(-d / range)
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
#' @description  Simple closures used to define relationships with covariates
#'
#' @param beta,mu,sigma Parameters
#' @param scale Center effect around zero?
#'
#' @rdname linear_rel
#' @export
linear_fun <- function(alpha = 0, beta = NULL, scale = TRUE) {
  function(x = NULL) {
    y <- alpha + beta * x
    if (scale) {
      y - mean(y)
    } else {
      y
    }
  }
}

#' @rdname linear_rel
#' @export
parabolic_fun <- function(mu = NULL, sigma = NULL, scale = TRUE) {
  function(x = NULL) {
    y <- -(((x - mu)^2) / (2 * sigma ^ 2))
    if (scale) {
      y - mean(y)
    } else {
      y
    }
  }
}


## Helper function for age and year covariance specification for function below.
## This is probably a clumsy way to deal with this problem...matrix math would be
## more elegant...but this was the most computationally efficient way I could work
## out.
.break_covar <- function(m, mar, covar) {
  switch(covar,
         ran = unname(m),
         rw = if (mar == "age") {
           unname(apply(m, 2, cumsum))
         } else {
           unname(t(apply(m, 1, cumsum)))
         }
         ,
         ident = if (mar == "age") {
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
#' @param pop An abundance at age matrix with ages defining the rows and years defining
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
                             space_covar = sim_covar(range = 200, variance = 0.5,
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
    if (length(age_covar) == 2) {
      j <- which(as.numeric(dimnames(e)$age) >= age_break)
      e[i, j, k] <- .break_covar(e[i, j, k], "age", age_covar[2])
    }
    j <- seq_along(pop$ages)
    k <- which(as.numeric(dimnames(e)$year) < year_break)
    e[i, , ]
    e[i, j, k] <- .break_covar(e[i, j, k], "year", year_covar[1])
    if (length(year_covar) == 2) {
      k <- which(as.numeric(dimnames(e)$year) >= year_break)
      e[i, j, k] <- .break_covar(e[i, j, k], "year", year_covar[2])
    }
    e[i, , ]
  }
  if (scale_error) { e <- e - mean(e) }


  # for(j in seq_along(pop$ages)) {
  #   xy$e <- e[, j, 2]
  #   plot(rasterFromXYZ(xy[, c("x", "y", "e")]), main = j)
  # }
  # for(k in seq_along(pop$years)) {
  #   xy$e <- e[, 2, k]
  #   plot(rasterFromXYZ(xy[, c("x", "y", "e")]), main = k)
  # }


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


  # for(j in seq_along(pop$ages)) {
  #   xy$n <- den[, j, 2]
  #   plot(rasterFromXYZ(xy[, c("x", "y", "n")]), main = j)
  # }
  # for(k in seq_along(pop$years)) {
  #   xy$n <- den[, 1, k]
  #   plot(rasterFromXYZ(xy[, c("x", "y", "n")]), main = k)
  # }

  N <- den * prod(res(grid))
  round(apply(N, c(2, 3), sum))
  round(pop$N)
  for(j in seq_along(pop$years)) {
    plot(apply(N, c(2, 3), sum)[, j], type = "s")
    lines(pop$N[, j], type = "s", col = "red")
  }



  message("TODO: clean up sim_distribution function and improve documentation")


  den


}





#' Simulate an age-structured population over a survey grid.
#'
#' @description Simulate age-structured population that varies in space and time
#'
#' @param ages Ages to include in the simulated population.
#' @param grid Survey grid to populate with simulated data.
#' @return A 'data.frame' containing ages and lengths for every fish in the
#'     simulated population.
#' @examples
#' sim_pop()
#' @export
sim_pop <- function(ages = 1:6,
                   grid = survey_grid) {
  print("In progress")
}




