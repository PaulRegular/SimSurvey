

#' Simulate startinig abundance, random recruitment and total mortality
#'
#' @description These functions return a function to use inside \code{\link{sim_abundance}}.
#' Given parameters, it generates N0, R and Z values.
#'
#' @param log_mean One mean value or a vector of means, in log scale, of length equal to years for \code{sim_R} or a matrix of means with
#' rows equaling the number of ages and colums equaling the number of years for \code{sim_Z}.
#' @param random_walk Simulate recruitment as a random walk?
#' @param log_sd Standard deviation of the variable in the log scale.
#' @param phi_age Autoregressive parameter for the age dimension.
#' @param phi_year Autoregressive parameter for the year dimension.
#' @param N0 Either specify "exp" or numeric vector of starting abundance excluding the first age.
#' If "exp" is specified using sim_N0, then abundance at age are calculated using exponential decay.
#' @param plot produce a simple plot of the simulated values?
#'
#' @details sim_R generates uncorrelated recruitment values or random walk values from a log normal distribution.
#' sim_Z does the same as sim_R when phi_age and phi_year are both 0, otherwise values are correlated
#' in the age and/or year dimension. The covariance structure follows that described in Cadigan (2015).
#'
#' @references Cadigan, Noel G. 2015. A State-Space Stock Assessment Model for Northern Cod,
#' Including Under-Reported Catches and Variable Natural Mortality Rates. Canadian Journal of
#' Fisheries and Aquatic Sciences 73 (2): 296-308.
#'
#' @examples
#'
#' R_fun <- sim_R(log_mean = log(100000), log_sd = 0.1, random_walk = TRUE, plot = TRUE)
#' R_fun(years = 1:100)
#' sim_abundance(R = sim_R(log_mean = log(100000), log_sd = 0.5))
#' sim_abundance(years = 1:20,
#'               R = sim_R(log_mean = log(c(rep(100000, 10), rep(10000, 10))), plot = TRUE))
#'
#' Z_fun <- sim_Z(log_mean = log(0.5), log_sd = 0.1, phi_age = 0.9, phi_year = 0.9, plot = TRUE)
#' Z_fun(years = 1:100, ages = 1:20)
#' sim_abundance(Z = sim_Z(log_mean = log(0.5), log_sd = 0.1, plot = TRUE))
#' Za_dev <- c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.3, 0.2, 0.1, 0)
#' Zy_dev <- c(-0.2, -0.2, -0.2, -0.2, -0.2, 2, 2, 2, 2, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0)
#' Z_mat <- outer(Za_dev, Zy_dev, "+") + 0.5
#' sim_abundance(ages = 1:10, years = 1:20,
#'               Z = sim_Z(log_mean = log(Z_mat), plot = TRUE))
#' sim_abundance(ages = 1:10, years = 1:20,
#'               Z = sim_Z(log_mean = log(Z_mat), log_sd = 0, phi_age = 0, phi_year = 0, plot = TRUE))
#'
#' N0_fun <- sim_N0(N0 = "exp", plot = TRUE)
#' N0_fun(R0 = 1000, Z0 = rep(0.5, 20), ages = 1:20)
#' sim_abundance(N0 = sim_N0(N0 = "exp", plot = TRUE))
#'
#' @export
#' @rdname sim_R
sim_R <- function(log_mean = log(30000000), log_sd = 0.5, random_walk = TRUE, plot = FALSE) {
  function(years = NULL) {
    if (length(log_mean) > 1 && length(log_mean) != length(years)) {
      stop("The number of log_means supplied for recruitment != number of years.")
    }
    e <- stats::rnorm(length(years), mean = 0, sd = log_sd)
    if (random_walk) e <- cumsum(e)
    r <- log_mean + e
    names(r) <- years
    if (plot) { plot(years, exp(r), type = "l", main = "sim_R", xlab = "Year", ylab = "Recruitment") }
    exp(r)
  }
}


#' @export
#' @rdname sim_R
sim_Z <- function(log_mean = 0.5, log_sd = 0.2, phi_age = 0.9, phi_year = 0.5, plot = FALSE) {
  function(ages = NULL, years = NULL) {

    na <- length(ages)
    ny <- length(years)
    if (is.matrix(log_mean) && (nrow(log_mean) != na | ncol(log_mean) != ny)) {
      stop("The matrix of log means supplied for Z != number of years and/or ages.")
    }

    Z <- matrix(NA, nrow = na, ncol = ny,
                dimnames = list(age = ages, year = years))
    pc_age <- sqrt(1 - phi_age ^ 2)
    pc_year <- sqrt(1 - phi_year ^ 2)
    for (j in seq_along(years)) {
      for (i in seq_along(ages)) {
        if ((i == 1) & (j == 1)) {
          m <- 0
          s <- log_sd / (pc_age * pc_year)
          Z[i, j] <- stats::rnorm(1, m, s)
        }
        if ((i > 1) & (j == 1)) {
          m <- phi_age * Z[i - 1, j]
          s <- log_sd / pc_year
          Z[i, j] <- stats::rnorm(1, m, s)
        }
        if ((i == 1) & (j > 1)) {
          m <- phi_year * Z[i, j - 1]
          s <- log_sd / pc_age
          Z[i, j] <- stats::rnorm(1, m, s)
        }
        if ((i > 1) & (j > 1)) {
          m <- phi_year * Z[i, j - 1] + phi_age * (Z[i - 1, j] - phi_year * Z[i - 1, j - 1])
          s <- log_sd
          Z[i, j] <- stats::rnorm(1, m, s)
        }
      }
    }
    Z <- log_mean + Z
    if (plot) { image(years, ages, t(exp(Z)), main = "sim_Z", xlab = "Year", ylab = "Age") }
    exp(Z)

  }
}

#' @export
#' @rdname sim_R
sim_N0 <- function(N0 = "exp", plot = FALSE) {
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
    if (plot) { plot(ages, N0, type = "h", xlab = "Age", ylab = "N0") }
    N0
  }
}



#' Convert length to length group
#'
#' @description Helper function for converting lengths to length groups
#' (Note: this isn't a general function; the output midpoints defining the
#' groups alligns with DFO specific method/labeling)
#'
#' @param length       Interval from \code{\link[base]{findInterval}}
#' @param group        Length group used to cut the length data
#'
#' @export
#'

group_lengths <- function(length, group) {
  breaks <- seq(0, max(length, na.rm = TRUE) * 2, group)
  interval <- findInterval(length, breaks)
  l <- breaks[interval]
  if (group == 0.5 | group == 1) { m <- l }
  if (group > 1) { m <- l + (0.5 * (group - 1)) }
  m
}


#' Closure for simulating length given age using von Bertalanffy notation
#'
#' This function outputs a function which holds the parameter values supplied and
#' the function either simulates lengths given ages or generates a length age key
#' give a sequence of ages.
#'
#' @param Linf          Mean asymptotic length
#' @param L0            Length at birth
#' @param K             Growth rate parameter
#' @param log_sd        Standard deviation of the relationship in log scale
#' @param length_group  Length group for length age key. Note that labels on the matrix produced are
#'                      midpoints using the DFO conventions; see \code{\link{group_lengths}}. Also
#'                      note that this length group will dictate the length group used in the
#'                      stratified analysis run by \code{\link{run_strat}}.
#' @param digits        Integer indicating the number of decimal places to round the values to
#' @param plot          Produce a simple plot of the simulated values?
#'
#' @examples
#' growth_fun <- sim_vonB(Linf = 100, L0 = 5, K = 0.2, log_sd = 0.05, length_group = 1, plot = TRUE)
#' growth_fun(age = rep(1:15, each = 100))
#' growth_fun(age = 1:15, length_age_key = TRUE)
#' sim_abundance(growth = sim_vonB(plot = TRUE))
#'
#' @export
#'

sim_vonB <- function(Linf = 120, L0 = 5, K = 0.1, log_sd = 0.1,
                     length_group = 3, digits = 0, plot = FALSE) {

  function(age = NULL, length_age_key = FALSE) {

    pred_length <- Linf - (Linf - L0) * exp(-K * age)

    if (length_age_key) {
      breaks <- seq(0, max(pred_length) * 10, length_group)
      breaks[1] <- 0.0000001
      lak <- matrix(NA, ncol = length(pred_length), nrow = length(breaks) - 1,
                    dimnames = list(length = group_lengths(breaks, length_group)[-length(breaks)],
                                    age = age))
      for (i in seq_along(breaks)[-1]) {
        for (j in seq_along(pred_length)) {
          lak[i - 1, j] <- stats::pnorm(log(breaks[i]), log(pred_length[j]), sd = log_sd) -
            stats::pnorm(log(breaks[i - 1]), log(pred_length[j]), sd = log_sd)
        }
      }
      lak <- lak[rowSums(lak) > 0, ]
      if (plot) image(x = as.numeric(colnames(lak)), y = as.numeric(rownames(lak)), z = t(lak),
                      xlab = "Age", ylab = "Length", main = "P(Length | Age)",
                      col = viridis::viridis(100))
      return(lak)

    } else {

      log_length <- stats::rnorm(length(age), log(pred_length), sd = log_sd)
      length <- round(exp(log_length), digits)
      if (plot) plot(age, length)
      return(length)

    }

  }
}


#' Convert abundance-at-age matrix to abundance-at-length
#'
#' Function for converting abundance-at-age matrix to abundance-at-length given
#' a length-age-key. Expects matrices to be named.
#'
#' @param   N_at_age    Abundance-at-age matrix
#' @param   lak         Length-age-key (i.e. probability of being in a specific length group given age)
#'
#' @return  Returns abundance-at-length matrix
#'
#' @export
#'
convert_N <- function(N_at_age = NULL, lak = NULL) {
  years <- colnames(N_at_age)
  N_at_length <- matrix(NA, nrow = nrow(lak), ncol = length(years),
                        dimnames = list(length = rownames(lak), year = years))
  for (y in seq_along(years)) {
    N_at_length[, y] <- colSums(N_at_age[, y] * t(lak))
  }
  if (!all.equal(colSums(N_at_age), colSums(N_at_length))) {
    stop("Aggregated abundance-at-age does not equal aggregated abundance-at-length.
         Check parameters of supplied growht model.")
  }
  N_at_length
}


#' Simulate basic population dynamics model
#'
#' @param ages     Ages to include in the simulation.
#' @param years    Years to include in the simulation.
#' @param Z        Total mortality function, like \code{\link{sim_Z}}, for generating
#'                 mortality matrix.
#' @param R        Recruitment (i.e. abundance at \code{min(ages)}) function, like
#'                 \code{\link{sim_R}}, for generating recruitment vector.
#' @param N0       Starting abundance (i.e. abundance at \code{min(years)}) function, like
#'                 \code{\link{sim_N0}}, for generating starting abundance vector.
#' @param growth   Closure, such as \code{\link{sim_vonB}}, for simulating length given age.
#'                 The function is used here to generate a abundance-at-age matrix and
#'                 it is carried foward for later use in \code{\link{sim_survey}} to simulate
#'                 lengths from survey catch at age.
#'
#' @return A \code{list} of length 9:
#' \itemize{
#'   \item{\code{ages}} - Vector of ages in the simulation
#'   \item{\code{lengths}} - Vector of length groups (depends on growth function)
#'   \item{\code{years}} - Vector of years in the simulation
#'   \item{\code{R} - Vector of recruitment values}
#'   \item{\code{N0} - Vector of starting abundance values}
#'   \item{\code{Z} - Matrix of total mortality values}
#'   \item{\code{N} - Matrix of abundance values}
#'   \item{\code{N_at_length} - Abundance at length matrix}
#'   \item{\code{sim_length} - Function for simulating lengths given ages}
#' }
#'
#' @details
#' Abundance from is calculated using a standard population dynamics model.
#' An abundance-at-length matrix is generated using a growth function coded as a closure like
#' \code{\link{sim_vonB}}. The function is retained for later use in \code{\link{sim_survey}}
#' to simulate lengths given simulated catch at age in a simulated survey. The ability to simulate
#' distributions by length is yet to be implemented.
#'
#' @examples
#'
#' R_fun <- sim_R(log_mean = log(100000), log_sd = 0.1, random_walk = TRUE, plot = TRUE)
#' R_fun(years = 1:100)
#' sim_abundance(R = sim_R(log_mean = log(100000), log_sd = 0.5))
#' sim_abundance(years = 1:20,
#'               R = sim_R(log_mean = log(c(rep(100000, 10), rep(10000, 10))), plot = TRUE))
#'
#' Z_fun <- sim_Z(log_mean = log(0.5), log_sd = 0.1, phi_age = 0.9, phi_year = 0.9, plot = TRUE)
#' Z_fun(years = 1:100, ages = 1:20)
#' sim_abundance(Z = sim_Z(log_mean = log(0.5), log_sd = 0.1, plot = TRUE))
#' Za_dev <- c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.3, 0.2, 0.1, 0)
#' Zy_dev <- c(-0.2, -0.2, -0.2, -0.2, -0.2, 2, 2, 2, 2, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0)
#' Z_mat <- outer(Za_dev, Zy_dev, "+") + 0.5
#' sim_abundance(ages = 1:10, years = 1:20,
#'               Z = sim_Z(log_mean = log(Z_mat), plot = TRUE))
#' sim_abundance(ages = 1:10, years = 1:20,
#'               Z = sim_Z(log_mean = log(Z_mat), log_sd = 0, phi_age = 0, phi_year = 0, plot = TRUE))
#'
#' N0_fun <- sim_N0(N0 = "exp", plot = TRUE)
#' N0_fun(R0 = 1000, Z0 = rep(0.5, 20), ages = 1:20)
#' sim_abundance(N0 = sim_N0(N0 = "exp", plot = TRUE))
#'
#' growth_fun <- sim_vonB(Linf = 100, L0 = 5, K = 0.2, log_sd = 0.05, length_group = 1, plot = TRUE)
#' growth_fun(age = rep(1:15, each = 100))
#' growth_fun(age = 1:15, length_age_key = TRUE)
#' sim_abundance(growth = sim_vonB(plot = TRUE))
#'
#' sim <- sim_abundance()
#' plot_trend(sim)
#' plot_surface(sim, mat = "N")
#' plot_surface(sim, mat = "Z")
#' plot_surface(sim, mat = "N_at_length", xlab = "Length", zlab = "N")
#'
#' @export
#'

sim_abundance <- function(ages = 1:20, years = 1:20,
                          Z = sim_Z(), R = sim_R(), N0 = sim_N0(),
                          growth = sim_vonB()) {

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

  ## Convert abundance-at-age matrix to abundance-at-length
  N_at_length <- convert_N(N_at_age = N,
                           lak = growth(age = ages, length_age_key = TRUE))
  lengths <- as.numeric(rownames(N_at_length))

  list(ages = ages,
       lengths = lengths,
       years = years,
       R = R,
       N0 = N0,
       Z = Z,
       N = N,
       N_at_length = N_at_length,
       sim_length = growth)

}



