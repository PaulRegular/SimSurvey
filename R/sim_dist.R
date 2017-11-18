
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
#' @param range  Decorrelation range
#' @param sd     Spatial variance
#' @param model  String indicating either "exponential" or "matern" as the correlation function
#'
#' @rdname sim_sp_covar
#' @export
sim_sp_covar <- function(range = 200, sd = 0.5, model = "matern") {
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
    covar <- (sd ^ 2) * cormat
    covar
  }
}

#' @rdname sim_sp_covar
#' @export
sim_ay_covar <- function(sd = 0.5, phi_age = 0.7, phi_year = 0.9) {
  function(ages = NULL, years = NULL, cells = NULL, w = NULL) {

    na <- length(ages)
    ny <- length(years)
    nc <- length(cells)
    if (length(phi_age) == 1) {
      phi_age <- rep(phi_age, na)
    } else {
      if (length(phi_age) != na) {
        stop("The number of phi_age values supplied != number of ages.")
      }
    }
    if (length(phi_year) == 1) {
      phi_year <- rep(phi_year, ny)
    } else {
      if (length(phi_year) != ny) {
        stop("The number of phi_year values supplied != number of years.")
      }
    }


    E <- array(rep(NA, nc * na * ny), dim = c(nc, na, ny),
               dimnames = list(cell = cells, age = ages, year = years))
    pc_age <- sqrt(1 - phi_age ^ 2)
    pc_year <- sqrt(1 - phi_year ^ 2)
    for (j in seq_along(years)) {
      for (i in seq_along(ages)) {
        if ((i == 1) & (j == 1)) {
          m <- 0
          s <- sd / pc_age[i]
          E[, i, j] <- w %*% rnorm(nc, m, s)
        }
        if ((i > 1) & (j == 1)) {
          m <- phi_age[i] * E[, i - 1, j]
          s <- sd / pc_year[j]
          E[, i, j] <- w %*% rnorm(nc, m, s)
        }
        if ((i == 1) & (j > 1)) {
          m <- phi_year[j] * E[, i, j - 1]
          s <- sd / pc_age[i]
          E[, i, j] <- w %*% rnorm(nc, m, s)
        }
        if ((i > 1) & (j > 1)) {
          m <- phi_year[j] * E[, i, j - 1] + phi_age[i] * (E[, i - 1, j] - phi_year[j] * E[, i - 1, j - 1])
          s <- sd
          E[, i, j] <- w %*% rnorm(nc, m, s)
        }
      }
    }
    E
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
parabolic_fun <- function(mu = 200, sigma = 100, scale = TRUE) {
  function(x = NULL) {
    y <- -(((x - mu)^2) / (2 * sigma ^ 2))
    if (scale) {
      y - mean(y)
    } else {
      y
    }
  }
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
                             grid = sim_grid(),
                             space_covar = sim_sp_covar(),
                             ay_covar = sim_ay_covar(phi_age = 0.05, phi_year = 0),
                             depth_par = parabolic_fun()
) {

  pop = sim_abundance()
  grid = sim_grid(space_covar = NULL)
  space_covar = sim_sp_covar(range = 50, sd = 0.1)
  ay_covar = sim_ay_covar(sd = 0.1, phi_age = c(rep(0, 4), rep(0.8, 10)), phi_year = 0)

  ## Spatial covariance
  grid_dat <- data.frame(raster::rasterToPoints(grid))
  xy <- grid_dat[, c("x", "y")]
  Sigma_space <- space_covar(xy)
  w <- t(chol(Sigma_space))

  ## Space-age-year autoregressive process
  E <- ay_covar(ages = pop$ages, years = pop$years, cells = grid_dat$cell, w = w)

  for (i in seq_along(pop$years)) {
    raster::plot(raster::rasterFromXYZ(data.frame(xy, E[, , i])))
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
  for (j in seq_along(pop$years)) {
    plot(apply(N, c(2, 3), sum)[, j], type = "s")
    lines(pop$N[, j], type = "s", col = "red")
  }



  message("TODO: clean up sim_distribution function and improve documentation")


  den


}


