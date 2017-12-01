
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
#' @param range    Decorrelation range
#' @param sd       Variance
#' @param phi_age  Defines autocorrelation through ages. Can be one value or a vector of the same
#'                 length as ages
#' @param phi_year Defines autocorrelation through years. Can be one value or a vector of the same
#'                 length as years
#' @param model    String indicating either "exponential" or "matern" as the correlation function
#'
#' @rdname sim_sp_covar
#' @export
sim_sp_covar <- function(range = 100, sd = 0.1, model = "matern") {
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
sim_ay_covar <- function(sd = 2, phi_age = 0.5, phi_year = 0.05) {
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


    E <- array(rep(NA, na * ny * nc), dim = c(na, ny, nc),
               dimnames = list(age = ages, year = years, cell = cells))
    pc_age <- sqrt(1 - phi_age ^ 2)
    pc_year <- sqrt(1 - phi_year ^ 2)
    for (j in seq_along(years)) {
      for (i in seq_along(ages)) {
        if ((i == 1) & (j == 1)) {
          m <- 0
          s <- sd / pc_age[i]
          E[i, j, ] <- w %*% rnorm(nc, m, s)
        }
        if ((i > 1) & (j == 1)) {
          m <- phi_age[i] * E[i - 1, j, ]
          s <- sd / pc_year[j]
          E[i, j, ] <- w %*% rnorm(nc, m, s)
        }
        if ((i == 1) & (j > 1)) {
          m <- phi_year[j] * E[i, j - 1, ]
          s <- sd / pc_age[i]
          E[i, j, ] <- w %*% rnorm(nc, m, s)
        }
        if ((i > 1) & (j > 1)) {
          m <- phi_year[j] * E[i, j - 1, ] + phi_age[i] * (E[i - 1, j, ] - phi_year[j] * E[i - 1, j - 1, ])
          s <- sd
          E[i, j, ] <- w %*% rnorm(nc, m, s)
        }
      }
    }
    E
  }
}



#' Define relationships with covariates
#'
#' @description  Closure to be used in \code{\link{sim_distribution}}
#'
#' @param alpha,mu,sigma  Parameters
#' @param plot            Produce a simple plot of the simulated values?
#'
#' @rdname gaussian_fun
#' @export

gaussian_fun <- function(alpha = 1, mu = 250, sigma = 100, plot = FALSE) {
  function(x = NULL) {
    y <- alpha * exp(-(((x - mu)^2) / (2 * sigma ^ 2)))
    if (plot) { plot(x, y, main = "gaussian_fun", type = "l") }
    y
  }
}


#' Simulate spatial and temporal distribution
#'
#' @description Provided an abundance at age matrix (like one provided by \code{\link{sim_abundance}})
#' and a survey grid (like \code{\link{survey_grid}}) to populate, this function
#' applies correlated space, age and year error to simulate the distribution
#' of the population.
#'
#' @param pop         An abundance at age matrix like one produced by \code{\link{sim_abundance}})
#' @param grid        A raster object defining the survey grid, like \code{\link{survey_grid}}
#'                    or one produced by \code{\link{sim_grid}}
#' @param space_covar Closure for simulating spatial covariance
#' @param ay_covar    Closure for simulating age-year covariance
#' @param depth_par   Closure for defining relationship between abundance and depth
#'
#' @details The addition of space-age-year error will result in slight discreptancies between
#'          input and output abundance at age
#'
#' @examples
#' sim_distribution()
#'
#' @export

sim_distribution <- function(pop = sim_abundance(),
                             grid = sim_grid(),
                             space_covar = sim_sp_covar(),
                             ay_covar = sim_ay_covar(),
                             depth_par = gaussian_fun()) {

  ## Space-age-year autoregressive process
  grid_dat <- data.frame(raster::rasterToPoints(grid))
  xy <- grid_dat[, c("x", "y")]
  Sigma_space <- space_covar(xy)
  w <- t(chol(Sigma_space))
  error <- ay_covar(ages = pop$ages, years = pop$years, cells = grid_dat$cell, w = w)

  ## Define probability of inhabiting each cell and distribute fish through the cells
  prob <- depth_par(grid_dat$depth)
  prob <- prob / sum(prob)  # make the values sum to 1
  N <- outer(pop$N, prob)
  dimnames(N) <- dimnames(error)

  ## Add space-age-year error to log abundance
  N <- exp(log(N) + error)

  ## Addition of error generates slight differences in the abundance at age
  ## Update pop$N
  pop$N <- apply(N, c(1, 2), sum)

  ## Output as data.frame
  df_N <- as.data.frame.table(N, responseName = "N", stringsAsFactors = FALSE)
  df_N <- merge(grid_dat, df_N, by = "cell")

  c(pop, list(grid = grid, sp_N = df_N))

}


