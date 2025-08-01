

#' Simulate starting abundance, random recruitment, and total mortality
#'
#' These functions return closures for use inside [`sim_abundance()`]. Given user-defined
#' parameters, they simulate recruitment (`R`), total mortality (`Z`), or initial abundance (`N0`)
#' as a function of age and year.
#'
#' @param log_mean For `sim_R`, a single mean or a vector of means (log scale) with length equal to the number of years.
#' For `sim_Z`, a matrix of log-scale means with rows equal to the number of ages and columns equal to the number of years.
#' @param random_walk Logical. Should recruitment be simulated as a random walk?
#' @param log_sd Standard deviation on the log scale.
#' @param phi_age Autoregressive parameter across the age dimension.
#' @param phi_year Autoregressive parameter across the year dimension.
#' @param N0 For `sim_N0`, either `"exp"` (for exponential decay) or a numeric vector of starting abundances (excluding the first age).
#' @param plot Logical. Should a simple plot of the simulated values be displayed?
#'
#' @details
#' - `sim_R()` generates uncorrelated or random-walk recruitment values from a log-normal distribution.
#' - `sim_Z()` behaves like `sim_R()` when both `phi_age` and `phi_year` are zero. When either is non-zero,
#'   it introduces correlation in the age and/or year dimension, based on the covariance structure
#'   described in Cadigan (2015).
#' - `sim_N0()` provides starting abundance either via exponential decay or a user-defined vector.
#'
#' @return A function to be passed to [`sim_abundance()`].
#'
#' @references
#' Cadigan, Noel G. (2015). A State-Space Stock Assessment Model for Northern Cod, Including
#' Under-Reported Catches and Variable Natural Mortality Rates. *Canadian Journal of Fisheries
#' and Aquatic Sciences*, 73(2): 296–308.
#'
#' @examples
#' R_fun <- sim_R(log_mean = log(100000), log_sd = 0.1, random_walk = TRUE, plot = TRUE)
#' R_fun(years = 1:100)
#'
#' sim_abundance(R = sim_R(log_mean = log(100000), log_sd = 0.5))
#'
#' sim_abundance(
#'   years = 1:20,
#'   R = sim_R(log_mean = log(c(rep(100000, 10), rep(10000, 10))), plot = TRUE)
#' )
#'
#' Z_fun <- sim_Z(log_mean = log(0.5), log_sd = 0.1, phi_age = 0.9, phi_year = 0.9, plot = TRUE)
#' Z_fun(years = 1:100, ages = 1:20)
#'
#' sim_abundance(Z = sim_Z(log_mean = log(0.5), log_sd = 0.1, plot = TRUE))
#'
#' Za_dev <- c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.3, 0.2, 0.1, 0)
#' Zy_dev <- c(-0.2, -0.2, -0.2, -0.2, -0.2, 2, 2, 2, 2, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0)
#' Z_mat <- outer(Za_dev, Zy_dev, "+") + 0.5
#'
#' sim_abundance(ages = 1:10, years = 1:20, Z = sim_Z(log_mean = log(Z_mat), plot = TRUE))
#'
#' sim_abundance(
#'   ages = 1:10, years = 1:20,
#'   Z = sim_Z(log_mean = log(Z_mat), log_sd = 0, phi_age = 0, phi_year = 0, plot = TRUE)
#' )
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
sim_Z <- function(log_mean = log(0.5), log_sd = 0.2, phi_age = 0.9, phi_year = 0.5, plot = FALSE) {
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
    if (plot) { graphics::image(years, ages, t(exp(Z)), main = "sim_Z", xlab = "Year", ylab = "Age") }
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
#' Helper function for converting lengths to length groups.
#' **Note:** This is not a general-purpose function — the output midpoints defining the groups
#' are aligned with DFO-specific methods and labeling conventions.
#'
#' @param length Numeric vector of lengths to be grouped. Used with [`base::findInterval()`].
#' @param group Numeric value specifying the width of the length group (i.e., bin size).
#'
#' @return A numeric vector indicating the midpoint of the assigned length group.
#'
#' @export

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
#' This function returns a closure that holds the supplied parameter values and can be used to
#' either simulate lengths given ages or generate a length-at-age key from a sequence of ages.
#'
#' @param Linf Mean asymptotic length.
#' @param L0 Length at birth.
#' @param K Growth rate parameter.
#' @param log_sd Standard deviation of the length-at-age relationship on the log scale.
#' @param length_group Length group width for constructing the length-at-age key.
#' Labels on the resulting matrix use midpoints according to DFO conventions; see [`group_lengths()`].
#' This value will also determine the length groupings used in the stratified analysis via [`run_strat()`].
#' @param digits Number of decimal places to round simulated lengths to.
#' @param plot Logical. Should a simple plot of the simulated values be produced?
#'
#' @return A function that can be passed to [`sim_abundance()`].
#'
#' @examples
#' growth_fun <- sim_vonB(Linf = 100, L0 = 5, K = 0.2, log_sd = 0.05, length_group = 1, plot = TRUE)
#' growth_fun(age = rep(1:15, each = 100))
#' growth_fun(age = 1:15, length_age_key = TRUE)
#'
#' sim_abundance(growth = sim_vonB(plot = TRUE))
#'
#' @export

sim_vonB <- function(Linf = 120, L0 = 5, K = 0.1, log_sd = 0.1,
                     length_group = 3, digits = 0, plot = FALSE) {

  function(age = NULL, length_age_key = FALSE) {

    pred_length <- Linf - (Linf - L0) * exp(-K * age)

    if (length_age_key) {

      breaks <- seq(0, ceiling(max(pred_length)) * 10, length_group)
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
      if (plot) graphics::image(x = as.numeric(colnames(lak)), y = as.numeric(rownames(lak)), z = t(lak),
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
#' Converts an abundance-at-age matrix to an abundance-at-length matrix using
#' a length-age key. Both input matrices must be named appropriately.
#'
#' @param N_at_age A matrix of abundance-at-age values.
#' @param lak A length-age key matrix — i.e., the probability of being in a specific
#' length group given age.
#'
#' @return A matrix of abundance-at-length values.
#'
#' @export

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
#' Simulates a basic age-structured population using recruitment (`R`),
#' total mortality (`Z`), and initial abundance (`N0`) functions. Optionally,
#' a growth function may be provided to simulate lengths given age and generate
#' an abundance-at-length matrix.
#'
#' @param ages A numeric vector of ages to include in the simulation.
#' @param years A numeric vector of years to include in the simulation.
#' @param Z A function for generating a total mortality matrix, such as [`sim_Z()`].
#' @param R A function for generating a recruitment vector (i.e., abundance at `min(ages)`),
#' such as [`sim_R()`].
#' @param N0 A function for generating a starting abundance vector (i.e., abundance at `min(years)`),
#' such as [`sim_N0()`].
#' @param growth A closure, such as [`sim_vonB()`], for simulating length given age.
#' This is used both to generate an abundance-at-length matrix and later for length simulation
#' in [`sim_survey()`].
#'
#' @return A list with the following elements:
#'
#' - `ages`: Vector of ages used in the simulation
#' - `lengths`: Vector of length groups (depends on growth function)
#' - `years`: Vector of years used in the simulation
#' - `R`: Vector of recruitment values
#' - `N0`: Vector of starting abundance values
#' - `Z`: Matrix of total mortality values
#' - `N`: Matrix of abundance-at-age
#' - `N_at_length`: Matrix of abundance-at-length
#' - `sim_length`: Function for simulating lengths given ages
#'
#' @details
#' Abundance is simulated using a standard population dynamics model. If a growth function
#' such as [`sim_vonB()`] is provided, it is used to create a corresponding abundance-at-length
#' matrix. The same growth function is retained for use in [`sim_survey()`] to simulate lengths
#' from catch-at-age survey data.
#'
#' Note: The ability to simulate distributions by length is not yet implemented.
#'
#' @examples
#' R_fun <- sim_R(log_mean = log(100000), log_sd = 0.1, random_walk = TRUE, plot = TRUE)
#' R_fun(years = 1:100)
#'
#' sim_abundance(R = sim_R(log_mean = log(100000), log_sd = 0.5))
#'
#' sim_abundance(
#'   years = 1:20,
#'   R = sim_R(log_mean = log(c(rep(100000, 10), rep(10000, 10))), plot = TRUE)
#' )
#'
#' Z_fun <- sim_Z(log_mean = log(0.5), log_sd = 0.1, phi_age = 0.9, phi_year = 0.9, plot = TRUE)
#' Z_fun(years = 1:100, ages = 1:20)
#'
#' sim_abundance(Z = sim_Z(log_mean = log(0.5), log_sd = 0.1, plot = TRUE))
#'
#' Za_dev <- c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.3, 0.2, 0.1, 0)
#' Zy_dev <- c(-0.2, -0.2, -0.2, -0.2, -0.2, 2, 2, 2, 2,
#'             0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0)
#' Z_mat <- outer(Za_dev, Zy_dev, "+") + 0.5
#'
#' sim_abundance(
#'   ages = 1:10, years = 1:20,
#'   Z = sim_Z(log_mean = log(Z_mat), plot = TRUE)
#' )
#'
#' sim_abundance(
#'   ages = 1:10, years = 1:20,
#'   Z = sim_Z(log_mean = log(Z_mat), log_sd = 0, phi_age = 0, phi_year = 0, plot = TRUE)
#' )
#'
#' N0_fun <- sim_N0(N0 = "exp", plot = TRUE)
#' N0_fun(R0 = 1000, Z0 = rep(0.5, 20), ages = 1:20)
#'
#' sim_abundance(N0 = sim_N0(N0 = "exp", plot = TRUE))
#'
#' growth_fun <- sim_vonB(Linf = 100, L0 = 5, K = 0.2,
#'                        log_sd = 0.05, length_group = 1, plot = TRUE)
#' growth_fun(age = rep(1:15, each = 100))
#' growth_fun(age = 1:15, length_age_key = TRUE)
#'
#' sim_abundance(growth = sim_vonB(plot = TRUE))
#'
#' sim <- sim_abundance()
#' plot_trend(sim)
#' plot_surface(sim, mat = "N")
#' plot_surface(sim, mat = "Z")
#' plot_surface(sim, mat = "N_at_length", xlab = "Length", zlab = "N")
#'
#' @export

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



