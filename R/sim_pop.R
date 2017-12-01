

#' Simulate startinig abundance, random recruitment and total mortality
#'
#' @description These functions return a function to use inside \code{\link{sim_abundance}}.
#' Given parameters, it generates N0, R and Z values.
#'
#' @param mean One mean value or a vector of means of length equal to years for \code{sim_R} or a matrix of means with
#' rows equaling the number of ages and colums equaling the number of years for \code{sim_Z}.
#' @param log_sd Standard deviation of the variable in the log scale.
#' @param phi_age Autoregressive parameter for the age dimension.
#' @param phi_year Autoregressive parameter for the year dimension.
#' @param N0 Either specify "exp" or numeric vector of starting abundance excluding the first age.
#' If "exp" is specified using sim_N0, then abundance at age are calculated using exponential decay:
#' \deqn{N_{a, 1} = N_{a - 1, 1} * exp(-Z_{a - 1, 1})}{N_a,1 = N_a-1,1 * exp(-Z_a-1,1)}
#' @param plot produce a simple plot of the simulated values?
#'
#' @details sim_R simply generates uncorrelated recruitment values from a log normal distribution.
#' sim_Z does the same as sim_R when phi_age and phi_year are both 0, otherwise values are correlated
#' in the age and/or year dimension. The covariance structure follows that described in Cadigan (2015).
#'
#' @references Cadigan, Noel G. 2015. A State-Space Stock Assessment Model for Northern Cod,
#' Including Under-Reported Catches and Variable Natural Mortality Rates. Canadian Journal of
#' Fisheries and Aquatic Sciences 73 (2): 296–308.
#'
#' @examples
#' sim_abundance(R = sim_R(mean = 100000, log_sd = 4))
#' sim_abundance(years = 1:20, R = sim_R(mean = c(rep(100000, 10), rep(10000, 10))))
#' sim_abundance(Z = sim_Z(mean = 0.6, log_sd = 0))
#' Za_dev <- c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.3, 0.2, 0.1, 0)
#' Zy_dev <- c(-0.2, -0.2, -0.2, -0.2, -0.2, 2, 2, 2, 2, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0)
#' Z_mat <- outer(Za_dev, Zy_dev, "+") + 0.5
#' sim_abundance(ages = 1:10, years = 1:20, Z = sim_Z(mean = Z_mat), R = sim_R(log_sd = 0.5))
#'
#' @export
#' @rdname sim_R
sim_R <- function(mean = 100000, log_sd = 0.5, plot = FALSE) {
  function(years = NULL) {
    if (length(mean) > 1 && length(mean) != length(years)) {
      stop("The number of means supplied for recruitment != number of years.")
    }
    r <- rnorm(length(years), mean = log(mean), sd = log_sd)
    names(r) <- years
    if (plot) { plot(years, exp(r), type = "l", main = "sim_R", xlab = "Year", ylab = "Recruitment") }
    exp(r)
  }
}

#' @export
#' @rdname sim_R
sim_Z <- function(mean = 0.5, log_sd = 0.5, phi_age = 0, phi_year = 0, plot = FALSE) {
  function(ages = NULL, years = NULL) {

    na <- length(ages)
    ny <- length(years)
    if (is.matrix(mean) && (nrow(mean) != na | ncol(mean) != ny)) {
      stop("The matrix of means supplied for Z != number of years and/or ages.")
    }

    Z <- matrix(NA, nrow = na, ncol = ny,
                dimnames = list(age = ages, year = years))
    pc_age <- sqrt(1 - phi_age ^ 2)
    pc_year <- sqrt(1 - phi_year ^ 2)
    for (j in seq_along(years)) {
      for (i in seq_along(ages)) {
        if ((i == 1) & (j == 1)) {
          m <- 0
          s <- log_sd / pc_age
          Z[i, j] <- rnorm(1, m, s)
        }
        if ((i > 1) & (j == 1)) {
          m <- phi_age * Z[i - 1, j]
          s <- log_sd / pc_year
          Z[i, j] <- rnorm(1, m, s)
        }
        if ((i == 1) & (j > 1)) {
          m <- phi_year * Z[i, j - 1]
          s <- log_sd / pc_age
          Z[i, j] <- rnorm(1, m, s)
        }
        if ((i > 1) & (j > 1)) {
          m <- phi_year * Z[i, j - 1] + phi_age * (Z[i - 1, j] - phi_year * Z[i - 1, j - 1])
          s <- log_sd
          Z[i, j] <- rnorm(1, m, s)
        }
      }
    }
    Z <- log(mean) + Z
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
#' Abundance from is calculated using a standard population dynamics model:
#' \deqn{N_{a, y} = N_{a - 1, y - 1} * exp(-Z_{a - 1, y - 1})}{N_a,y = N_a-1,y-1 * exp(-Z_a-1,y-1)}
#'
#' @examples
#' sim_abundance()
#'
#' @export

sim_abundance <- function(ages = 1:14, years = 1:20,
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



