
## Helper function for euclidian distance calculations
.dist <- function(x) {
  if (requireNamespace("fields", quietly = TRUE)) {
    d <- fields::rdist(x)
  } else {
    d <- as.matrix(stats::dist(x))
  }
}

## Helper function for spatial correlation
.sp_cor <- function(x = NULL, range = 50, lambda = 1, model = "matern") {
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
  cormat
}

#' Simulate age-year-space covariance
#'
#' These functions return a function to be used inside [`sim_distribution()`] to simulate
#' spatially and temporally autocorrelated error.
#'
#' @param range Decorrelation range.
#' @param lambda Controls the degree of smoothness in the Matérn covariance process.
#' @param sd Variance (can be age-specific).
#' @param phi_age Autocorrelation through ages. Can be a single value or a vector the same length as `ages`.
#' @param phi_year Autocorrelation through years. Can be a single value or a vector the same length as `years`.
#' @param group_ages A vector of ages to group together with shared space-age-year noise.
#' @param group_years A vector of years to group together with shared space-age-year noise.
#' @param model String indicating the correlation function to use: either `"exponential"` or `"matern"`.
#'
#' @return A function to be passed to [`sim_distribution()`] as the `ays_covar` argument.
#'
#' @export

sim_ays_covar <- function(sd = 2.8, range = 300, lambda = 1, model = "matern",
                          phi_age = 0.5, phi_year = 0.9,
                          group_ages = 5:20, group_years = NULL) {
  function(x = NULL, ages = NULL, years = NULL, cells = NULL) {

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
    ##chol is costly, so calculate only once!
    chol_cor <- chol(.sp_cor(x = x, range = range, lambda = lambda, model = model))
    for (j in seq_along(years)) {
      for (i in seq_along(ages)) {
        if ((i == 1) & (j == 1)) {
          m <- 0
          s <- sd[i] / (pc_age[i] * pc_year[j])
          if (!exists("w1")) { # Might save some time...
            Sigma <- s * chol_cor
            w1 <- t(Sigma)
          }
          E[i, j, ] <- m + w1 %*% stats::rnorm(nc)
        }
        if ((i > 1) & (j == 1)) {
          if (age_map[i] == age_map[i - 1]) {
            E[i, j, ] <- E[i - 1, j, ]
          } else {
            m <- phi_age[i] * E[i - 1, j, ]
            s <- sd[i] / pc_year[j]
            if (!exists("w2")) {
              Sigma <- s * chol_cor
              w2 <- t(Sigma)
            }
            E[i, j, ] <- m + w2 %*% stats::rnorm(nc)
          }
        }
        if ((i == 1) & (j > 1)) {
          if (year_map[j] == year_map[j - 1]) {
            E[i, j, ] <- E[i, j - 1, ]
          } else {
            m <- phi_year[j] * E[i, j - 1, ]
            s <- sd[i] / pc_age[i]
            if (!exists("w3")) {
              Sigma <- s * chol_cor
              w3 <- t(Sigma)
            }
            E[i, j, ] <- m + w3 %*% stats::rnorm(nc)
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
              Sigma <- s * chol_cor
              w4 <- t(Sigma)
            }
            E[i, j, ] <- m + w4 %*% stats::rnorm(nc)
          }

        }
      }
    }
    E
  }
}



#' Define a parabolic relationship
#'
#' Closure to be used in [`sim_distribution()`]. The form is based on the bi-Gaussian
#' function described [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2993707/).
#'
#' @param alpha, mu, sigma Parameters that control the shape of the parabola. Can be a single value or
#' a vector of the same length as the number of ages in the simulation (e.g., allowing age-specific
#' depth associations).
#' @param sigma_right Optional parameter to impose asymmetry by specifying a different `sigma` for the right side.
#' If used, `sigma` defines the width of the left side. Ignored if `NULL`.
#' @param log_space Logical. Should the shape of the parabola be defined in log space?
#' If `TRUE`, parameters are assumed to be logged, and input `x` values are log-transformed.
#' This produces a lognormal-like curve with a heavier right tail and forces low values near zero.
#' @param plot Logical. Should a simple plot of the simulated values be produced?
#'
#' @return A function that can be passed to [`sim_distribution()`].
#'
#' @examples
#' parabola_fun <- sim_parabola(mu = 50, sigma = 5, plot = TRUE)
#' parabola_fun(data.frame(depth = 0:100))
#'
#' parabola_fun <- sim_parabola(mu = log(40), sigma = 0.5, log_space = FALSE, plot = TRUE)
#' parabola_fun(data.frame(depth = 0:100))
#'
#' parabola_fun <- sim_parabola(mu = c(50, 120), sigma = c(5, 3), plot = TRUE)
#' parabola_fun(expand.grid(depth = 1:200, age = 1:2))
#'
#' @rdname sim_parabola
#' @export

sim_parabola <- function(alpha = 0, mu = 200, sigma = 70, sigma_right = NULL,
                         log_space = FALSE, plot = FALSE) {

  function(data) {

    x <- data$depth
    age <- data$age
    nages <- length(unique(age))
    npar <- c(length(alpha), length(mu), length(sigma))
    if (any(npar > 1) && max(npar) != nages) {
      stop("The number alpha, mu, or sigma values != number of ages.")
    }
    if (length(unique(npar[npar!=1])) > 1) {
      stop("Inconsistent number of alpha, mu, or sigma values were supplied. One value should be supplied or a vector equal to the number of ages in the simulation.")
    }
    if (any(npar > 1)) {
      if (npar[1] == 1) alpha <- rep(alpha, nages)
      if (npar[2] == 1) mu <- rep(mu, nages)
      if (npar[3] == 1) sigma <- rep(sigma, nages)
      i <- age - min(age) + 1 # convert to index (address cases where start age may be 0 or 4, etc.)
    } else {
      i <- rep(1, length(x))
    }

    if (is.null(sigma_right)) sigma_right <- sigma
    if (log_space) mx <- log(x) else mx <- x

    y <- ifelse(mx < mu[i], alpha[i] - (((mx - mu[i])^2) / (2 * sigma[i] ^ 2)),
                            alpha[i] - (((mx - mu[i])^2) / (2 * sigma_right[i] ^ 2)))

    if (plot) {
      plot(x, y, main = "sim_parabola", col = i)
    }

    y

  }

}


#' Define a non-linear relationship
#'
#' `r lifecycle::badge("experimental")`
#'
#' Closure to be used in [`sim_distribution()`].
#'
#' @param formula A formula describing parametric relationships between the data and coefficients.
#' The data used in [`sim_distribution()`] consist of grid coordinates expanded across ages and years,
#' and include columns such as `"x"`, `"y"`, `"depth"`, `"cell"`, `"division"`, `"strat"`, `"age"`, and `"year"`.
#' Coefficient values referenced in the formula must be provided in the `coeff` argument as a named list.
#' @param coeff A named list of coefficient values used in `formula`.
#'
#' @return A function that can be passed to [`sim_distribution()`].
#'
#' @examples
#' ## Make a grid and replicate data for 5 ages and 5 years
#' ## (This mimics what happens inside sim_distribution)
#' grid <- make_grid(shelf_width = 10)
#' grid_xy <- data.frame(grid)
#' i <- rep(seq(nrow(grid_xy)), times = 5)
#' a <- rep(1:5, each = nrow(grid_xy))
#' grid_xy <- grid_xy[i, ]
#' grid_xy$age <- a
#' i <- rep(seq(nrow(grid_xy)), times = 5)
#' y <- rep(1:5, each = nrow(grid_xy))
#' grid_xy <- grid_xy[i, ]
#' grid_xy$year <- y
#'
#' ## Define a non-linear function to apply to the expanded grid
#' ## This example imposes ontogenetic deepening via a parabolic depth effect
#' nlf <- sim_nlf(
#'   formula = ~ alpha - ((depth - mu + beta * age)^2) / (2 * sigma^2),
#'   coeff = list(alpha = 0, mu = 200, sigma = 70, beta = -70)
#' )
#' grid_xy$depth_effect <- nlf(grid_xy)
#'
#' library(plotly)
#' grid_xy %>%
#'   filter(year == 1) %>%
#'   plot_ly(x = ~depth, y = ~depth_effect, split = ~age) %>%
#'   add_lines()
#'
#' @importFrom lifecycle badge
#' @export

sim_nlf <- function(formula = ~ alpha - ((depth - mu) ^ 2) / (2 * sigma ^ 2),
                    coeff = list(alpha = 0, mu = 200, sigma = 70)) {

  function(data) {

    vars <- all.vars(formula)

    for (v in vars) {
      if (v %in% names(data)) assign(v, data[[v]])
      if (v %in% names(coeff)) assign(v, coeff[[v]])
    }

    eval(formula[[2]])

  }


}



#' Simulate spatial and temporal distribution
#'
#' Provided with an abundance-at-age matrix and a survey grid, this function applies
#' correlated space–age–year error to simulate the spatial distribution of a population.
#' Simulation by length is not yet implemented.
#'
#' @param sim A list with `ages`, `years`, and an abundance-at-age matrix,
#' like one produced by [`sim_abundance()`].
#' @param grid A `stars` object defining the survey grid, such as [`survey_grid`] or one created with [`make_grid()`].
#' @param ays_covar A closure that simulates age-year-space covariance, such as [`sim_ays_covar()`].
#' @param depth_par A closure that defines the relationship between abundance and depth,
#' such as [`sim_parabola()`].
#'
#' @details
#' This function simulates the probability of fish inhabiting a cell based on a parabolic
#' relationship with depth and spatially/temporally autocorrelated noise across age and year.
#'
#' **Warning:** Simulating a large grid across many ages and years may be computationally intensive.
#' Start with smaller simulations to test your setup.
#'
#' @return Appends the following objects to the `sim` list:
#'
#' - `grid`: A `stars` object with grid details
#' - `grid_xy`: A `data.table` representation of the grid in XYZ format
#' - `sp_N`: A `data.table` of abundance split by age, year, and cell
#'
#' @examples
#' \donttest{
#' sim <- sim_abundance(ages = 1:5, years = 1:5) %>%
#'   sim_distribution(
#'     grid = make_grid(res = c(20, 20)),
#'     ays_covar = sim_ays_covar(phi_age = 0.8, phi_year = 0.1),
#'     depth_par = sim_parabola(mu = 200, sigma = 50)
#'   )
#'
#' head(sim$sp_N)
#' head(sim$grid_xy)
#' }
#'
#' @export
#' @rawNamespace import(data.table, except = shift)

sim_distribution <- function(sim,
                             grid = make_grid(),
                             ays_covar = sim_ays_covar(),
                             depth_par = sim_parabola()) {

  ## Space-age-year autoregressive process
  grid_dat <- data.table::data.table(as.data.frame(grid))
  setkeyv(grid_dat, "cell")
  xy <- grid_dat[, c("x", "y")]
  error <- ays_covar(x = xy, ages = sim$ages, years = sim$years, cells = grid_dat$cell)

  ## Relationship with depth (expand grid_dat to include age and year)
  grid_edat <- grid_dat
  i <- rep(seq(nrow(grid_edat)), times = length(sim$ages))
  a <- rep(sim$ages, each = nrow(grid_edat))
  grid_edat <- grid_edat[i]
  grid_edat$age <- a
  i <- rep(seq(nrow(grid_edat)), times = length(sim$years))
  y <- rep(sim$years, each = nrow(grid_edat))
  grid_edat <- grid_edat[i]
  grid_edat$year <- y
  grid_edat <- grid_edat[order(grid_edat$cell, grid_edat$year, grid_edat$age), ] # sort to align with error array
  depth <- depth_par(grid_edat)
  depth <- array(depth, dim = dim(error), dimnames = dimnames(error))

  ## Define probability of inhabiting each cell
  ## Note the normalization; forcing the values to sum to one ensures that the final
  ## N matrix through the field will equal the N matrix generated by sim_abundance
  prob <- exp(depth + error)
  prob <- apply(prob, c(1, 2), function(x) { x / sum(x) })
  prob <- aperm(prob, c(2, 3, 1))

  if (any(is.nan(prob))) {
    stop("NaN values were produced when running sim_distribution. This problem was most likely caused by the effect supplied to the depth_par agument. Please ensure that the closure supplied is not producing extremely low values that effectively translate to zero when exponentiated. This problem can arise when, for example, values in log space are used in sim_parabola, but log_space = FALSE.")
  }

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


