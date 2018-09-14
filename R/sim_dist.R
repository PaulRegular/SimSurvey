
## Helper function for euclidian distance calculations
.dist <- function(x) {
  if (requireNamespace("fields", quietly = TRUE)) {
    d <- fields::rdist(x)
  } else {
    d <- as.matrix(dist(x))
  }
}

## Helper function for spatial covariance
.sp_covar <- function(x = NULL, range = 50, lambda = 1, sd = 1, model = "matern") {
  d <- .dist(x)
  cormat <- switch(model,
                   exponential = {
                     exp(-d / range)
                   },
                   matern = {
                     kappa <- sqrt(8 * lambda) / range # approximate kappa from range as per R-INLA book
                     2 ^ (1 - lambda) / gamma(lambda) * (kappa * d) ^ lambda *
                       besselK(x = d * kappa, nu = lambda)
                   },
                   stop("wrong or no specification of covariance model"))

  diag(cormat) <- 1
  covar <- (sd ^ 2) * cormat
  covar
}

#' Simulate age-year-space covariance
#'
#' @description These functions return a function to use inside \code{\link{sim_distribution}}.
#'
#' @param range       Decorrelation range
#' @param lambda      Controls the degree of smoothness of Matern covariance process
#' @param sd          Variance (can be age specific).
#' @param phi_age     Defines autocorrelation through ages. Can be one value or a vector of the same
#'                    length as ages
#' @param phi_year    Defines autocorrelation through years. Can be one value or a vector of the same
#'                    length as years
#' @param group_ages  Make space-age-year variance equal across these ages
#' @param group_years Make space-age-year variance equal across these years
#' @param model       String indicating either "exponential" or "matern" as the correlation function
#'
#' @export

sim_ays_covar <- function(sd = 0.5, range = 50, lambda = 1, model = "matern",
                          phi_age = 0.5, phi_year = 0.5,
                          group_ages = NULL, group_years = NULL) {
  function(x = NULL, ages = NULL, years = NULL, cells = NULL) {

    # Simple description of covariance:
    # In 2D: X_ay = X_a,y-1 + X_a-1,y - X_a-1,y-1 + error
    # at 1st age it is random walk in year,
    # and first year it is random walk in age

    # There are probably better, more elegant and computationally efficient solutions
    # to this problem.

    na <- length(ages)
    ny <- length(years)
    nc <- length(cells)
    if (length(sd) == 1) {
      sd <- rep(sd, na)
    } else {
      if (length(sd) != na) {
        stop("The number of sd values supplied != number of ages.")
      }
    }
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
    age_map <- as.character(ages)
    if (!is.null(group_ages)) {
      age_map[ages %in% group_ages] <- paste(range(group_ages), collapse = ":")
    }
    year_map <- as.character(years)
    if (!is.null(group_years)) {
      year_map[years %in% group_years] <- paste(range(group_years), collapse = ":")
    }

    E <- array(rep(NA, na * ny * nc), dim = c(na, ny, nc),
               dimnames = list(age = ages, year = years, cell = cells))
    pc_age <- sqrt(1 - phi_age ^ 2)
    pc_year <- sqrt(1 - phi_year ^ 2)
    for (j in seq_along(years)) {
      for (i in seq_along(ages)) {
        if ((i == 1) & (j == 1)) {
          m <- 0
          s <- sd[i] / (pc_age[i] * pc_year[j])
          if (!exists("w1")) { # chol is costly, therefore only calculate once
            Sigma <- .sp_covar(x = x, range = range, lambda = lambda, sd = s, model = model)
            w1 <- t(chol(Sigma))
          }
          E[i, j, ] <- m + w1 %*% rnorm(nc)
        }
        if ((i > 1) & (j == 1)) {
          if (age_map[i] == age_map[i - 1]) {
            E[i, j, ] <- E[i - 1, j, ]
          } else {
            m <- phi_age[i] * E[i - 1, j, ]
            s <- sd[i] / pc_year[j]
            if (!exists("w2")) {
              Sigma <- .sp_covar(x = x, range = range, lambda = lambda, sd = s, model = model)
              w2 <- t(chol(Sigma))
            }
            E[i, j, ] <- m + w2 %*% rnorm(nc)
          }
        }
        if ((i == 1) & (j > 1)) {
          if (year_map[j] == year_map[j - 1]) {
            E[i, j, ] <- E[i, j - 1, ]
          } else {
            m <- phi_year[j] * E[i, j - 1, ]
            s <- sd[i] / pc_age[i]
            if (!exists("w3")) {
              Sigma <- .sp_covar(x = x, range = range, lambda = lambda, sd = s, model = model)
              w3 <- t(chol(Sigma))
            }
            E[i, j, ] <- m + w3 %*% rnorm(nc)
          }
        }
        if ((i > 1) & (j > 1)) {

          if (age_map[i] == age_map[i - 1]) {
            E[i, j, ] <- E[i - 1, j, ]
          }
          if (year_map[j] == year_map[j - 1]) {
            E[i, j, ] <- E[i, j - 1, ]
          }
          if ((age_map[i] != age_map[i - 1]) & (year_map[j] != year_map[j - 1])) {
            m <- phi_year[j] * E[i, j - 1, ] + phi_age[i] * (E[i - 1, j, ] - phi_year[j] * E[i - 1, j - 1, ])
            s <- sd[i]
            if (!exists("w4")) { # chol is costly, therefore only calculate once
              Sigma <- .sp_covar(x = x, range = range, lambda = lambda, sd = s, model = model)
              w4 <- t(chol(Sigma))
            }
            E[i, j, ] <- m + w4 %*% rnorm(nc)
          }

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
#' of the population. The ability to simulate distributions by length is yet to be implemented.
#'
#' @param sim         An abundance at age matrix like one produced by \code{\link{sim_abundance}}
#' @param grid        A raster object defining the survey grid, like \code{\link{survey_grid}}
#'                    or one produced by \code{\link{sim_grid}}
#' @param ays_covar   Closure for simulating age-year-space covariance
#' @param depth_par   Closure for defining relationship between abundance and depth
#'
#' @details This function simulates the probability of simulated fish inhabiting
#' a cell as a function of a parabolic relationship with depth and space, age,
#' and year autocorrelated errors.
#'
#' @examples
#' sim_distribution()
#'
#' @export
#'
#' @rawNamespace import(data.table, except = shift)
#'

sim_distribution <- function(sim,
                             grid = sim_grid(),
                             ays_covar = sim_ays_covar(),
                             depth_par = sim_parabola()) {

  ## Space-age-year autoregressive process
  grid_dat <- data.table::data.table(raster::rasterToPoints(grid))
  setkeyv(grid_dat, "cell")
  xy <- grid_dat[, c("x", "y")]
  error <- ays_covar(x = xy, ages = sim$ages, years = sim$years, cells = grid_dat$cell)

  ## Relationship with depth
  depth <- depth_par(grid_dat$depth)
  depth <- replicate(length(sim$years), replicate(length(sim$ages), depth))
  depth <- aperm(depth, c(2, 3, 1))
  dimnames(depth) <- dimnames(error)

  ## Define probability of inhabiting each cell
  ## Note the normalization; forcing the values to sum to one ensures that the final
  ## N matrix through the field will equal the N matrix generated by sim_abundance
  prob <- exp(depth + error)
  prob <- apply(prob, c(1, 2), function(x) { x / sum(x) })
  prob <- aperm(prob, c(2, 3, 1))

  ## Distribute fish through the cells
  N <- replicate(nrow(grid_dat), sim$N)
  dimnames(N) <- dimnames(prob)
  N <- N * prob

  ## Output as data.frame
  df_N <- as.data.frame.table(prob, responseName = "prob", stringsAsFactors = FALSE)
  df_N <- data.table::data.table(df_N, N = c(N))
  df_N$year <- as.numeric(df_N$year)
  df_N$age <- as.numeric(df_N$age)
  df_N$cell <- as.numeric(df_N$cell)
  df_N$prob <- NULL
  setkeyv(df_N, "cell")

  c(sim, list(grid = grid, grid_xy = grid_dat, sp_N = df_N))

}


