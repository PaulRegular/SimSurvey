
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
sim_sp_covar <- function(range = 50, sd = 0.1, model = "matern") {
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
sim_ay_covar <- function(sd = 10, phi_age = 0.5, phi_year = 0.5) {
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
#' @param alpha,mu,sigma  Parameters that control the shape of the parabola
#' @param plot            Produce a simple plot of the simulated values?
#'
#' @rdname sim_parabola
#' @export

sim_parabola <- function(alpha = 0, mu = 250, sigma = 50, plot = FALSE) {
  function(x = NULL) {
    y <- alpha - (((x - mu)^2) / (2 * sigma ^ 2))
    if (plot) { plot(x, y, main = "sim_parabola", type = "l") }
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
#' @param sim         An abundance at age matrix like one produced by \code{\link{sim_abundance}}
#' @param grid        A raster object defining the survey grid, like \code{\link{survey_grid}}
#'                    or one produced by \code{\link{sim_grid}}
#' @param space_covar Closure for simulating spatial covariance
#' @param ay_covar    Closure for simulating age-year covariance
#' @param depth_par   Closure for defining relationship between abundance and depth
#'
#' @details This function mainly operates in logit space. Specifically, the probability of
#' simulated fish inhabiting a cell is a function of a parabolic relationship with depth
#' and space, age, and year autocorrelated errors.
#'
#' @examples
#' sim_distribution()
#'
#' @export

sim_distribution <- function(sim,
                             grid = sim_grid(),
                             space_covar = sim_sp_covar(),
                             ay_covar = sim_ay_covar(),
                             depth_par = sim_parabola()) {

  ## Space-age-year autoregressive process
  grid_dat <- data.frame(raster::rasterToPoints(grid))
  xy <- grid_dat[, c("x", "y")]
  Sigma_space <- space_covar(xy)
  w <- t(chol(Sigma_space))
  error <- ay_covar(ages = sim$ages, years = sim$years, cells = grid_dat$cell, w = w)

  ## Relationship with depth
  depth <- depth_par(grid_dat$depth)
  depth <- replicate(length(sim$years), replicate(length(sim$ages), depth))
  depth <- aperm(depth, c(2, 3, 1))
  dimnames(depth) <- dimnames(error)

  ## Define probability of inhabiting each cell
  ## Note the normalization; forcing the values to sum to one ensures that the final
  ## N matrix through the field will equal the N matrix generated by sim_abundance
  prob <- plogis(depth + error)
  prob <- apply(prob, c(1, 2), function(x) { x / sum(x) })
  prob <- aperm(prob, c(2, 3, 1))

  ## Distribute fish through the cells
  N <- replicate(nrow(grid_dat), sim$N)
  dimnames(N) <- dimnames(prob)
  N <- N * prob

  ## Output as data.frame
  df_N <- as.data.frame.table(prob, responseName = "prob", stringsAsFactors = FALSE)
  df_N <- data.frame(df_N, N = c(N))
  df_N <- merge(grid_dat, df_N, by = "cell")
  df_N$year <- as.numeric(df_N$year)
  df_N$age <- as.numeric(df_N$age)

  c(sim, list(grid = grid, sp_N = df_N[, setdiff(names(df_N), "prob")]))

}


