
#' Simulate random recruitment
#'
#' @description This function returns a function to use inside \code{\link{simAbundance}}.
#' Given parameters, it generates recruitment values from a log normal distribution
#'
#' @param mean,sd Mean and standard deviation of recruitment
#'
#' @examples
#' simAbundance(R = simR(mean = 100000, sd = 4))
#'
#' @export
simR <- function(mean = 20000, sd = 4) {
  function(years = NULL) {
    r <- rlnorm(length(years), meanlog = log(mean), sdlog = log(sd))
    names(r) <- years
    r
  }
}


#' Simulate total mortality using a random walk
#'
#' @description This function returns a function to use inside \code{\link{simAbundance}}.
#' Given parameters, it generates total mortality values from a log normal random walk
#'
#' @param mean,sd Mean and standard deviation of total mortality.
#' @param breaks Provide breaks to age and/or year series to specify group specific
#' mean total mortality values. If specified, multiple mean values must be
#' supplied; one more than the number of breaks. Provide named lists (see
#' examples below).
#'
#' @details Random walk errors are first generated for years and ages independently.
#' The outer product of these errors are then calculated to generate a total
#' mortality matrix as follows:
#' \deqn{log(Z_{a, y}) = log(z_{a}) + \delta_{a, y}}{Z_a,y = z_a,y + e_a,y}
#'
#'
#' @examples
#' simAbundance(Z = simZ(mean = 0.6, sd = 0.3))
#' simAbundance(Z = simZ(mean = list(ages = c(0.5, 0.3, 0.2)),
#'                       breaks = list(ages = c(1, 2)), sd = 0.3))
#' simAbundance(Z = simZ(mean = list(ages = c(0.5, 0.3, 0.2), years = c(0.3, 2, 0.3)),
#'                       breaks = list(ages = c(1, 2), years = c(4, 6)), sd = 0.3))
#'
#' @export
simZ <- function(mean = 0.4, sd = 0.2, breaks = NULL) {
  function(ages = NULL, years = NULL) {
    ey <- cumsum(rnorm(length(years), sd = sd)) # random walk error
    ea <- cumsum(rnorm(length(ages), sd = sd))
    e <- outer(ea, ey)
    if (!is.null(breaks)) {
      mean_ages <- mean$ages[findInterval(ages, breaks$ages, left.open = TRUE) + 1]
      mean_years <- mean$years[findInterval(years, breaks$years, left.open = TRUE) + 1]
      if (is.null(mean_ages)) { mean_ages <- rep(1, length(ages)) }
      if (is.null(mean_years)) { mean_years <- rep(1, length(years)) }
      mean <- outer(mean_ages, mean_years)
    }
    z <- exp(log(mean) + e)
    dimnames(z) <- list(age = ages, year = years)
    z
  }
}



#' Simulate basic population dynamics model
#'
#' @param ages Ages to include in the simulation.
#' @param years Years to include in the simulation.
#' @param Z Total mortality function, like \code{\link{simZ}}, for generating
#' mortality matrix.
#' @param R Recruitment (i.e. Abundance at \code{min(ages)}) function, like
#' \code{\link{simR}}, for generating recruitment vector.
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
#' simAbundance()
#'
#' @export

simAbundance <- function(ages = 1:6, years = 1:10, Z = simZ(), R = simR()) {

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
  list(R = R, Z = Z, N = N)

}


#' Simulate size, time and space correlation structure
#'
#'


## Helper function for euclidian distance calculations
.dist <- function(x) {
  if (requireNamespace("fields", quietly = TRUE)) {
    d <- fields::rdist(x)
  } else {
    d <- as.matrix(dist(x))
  }
}


simTime <- function(dcor_time = 2) {
  function(years = NULL) {
    d <- .dist(years)
    exp(- d / dcor_time)
  }
}


simSize <- function(dcor_size = 4) {
  function(ages = NULL) {
    d <- .dist(ages)
    exp(- d / dcor_size)
  }
}

simSpace <- function(dcor_dist = 70, sigma = 1.5) {
  function(grid = NULL) {
    d <- .dist(coordinates(grid))
    sigma * exp(- d / dcor_dist)
  }
}


simSpacePM <- function(tau = 0.25, theta = 0.001) {
  function(grid = NULL) {

    d <- .dist(coordinates(grid))
    nb <- Q <- rgeos::gTouches(grid, byid = TRUE) # cell neighbour matrix
    m <- rowSums(nb)                              # number of neighbouring cells
    Q[] <- 0
    Q[nb] <- - tau
    diag(Q) <- tau * (m + theta)
    invQ <- solve(Q)
    s <- (1 / length(grid)) * sum(diag(invQ))
    H <- mean(d[nb]) / log(1 + (theta / 2) + sqrt(theta + ((theta ^ 2)/4)))
    message(paste0("Spatial variance is approxamatly ", signif(s, 3),
                   "\nSpatial decorrelation distance is approxamatly ", signif(H, 3), " km"))

    cell <- sample(grid$cell, 1)
    # cell <- 8179
    ncols <- 200
    cols <- cut(invQ[cell, ], breaks = ncols, labels = FALSE)
    cols <- colorRampPalette(c("white", "steelblue", "navy"))(ncols)[cols]
    plot(d[cell, ], invQ[cell, ], col = cols, pch = 16, cex = 0.75,
         xlab = "Distance", ylab = "Correlation")
    plot(grid, col = cols, lwd = 0.5)

  }
}







#' Simulate spatial and temporal distribution
#'
#' @description Provided an abundance at age matrix (like one provided by \code{\link{simAbundance}})
#' and a survey grid (like \code{\link{survey_grid}}) to populate, this function
#' applies spatial and temporal error to simulate the spatial and temporal distribution
#' of the population.
#'
#' @param N An abundance at age matrix with ages defining the rows and years defining
#' the columns (i.e. same structure as a matrix provided by \code{\link{simAbundance}})
#' @param grid A \code{\link{SpatialPolygonsDataFrame}} defining a regular or irregular
#' grid with the same structure as \code{\link{survey_grid}}
#'
#' @examples
#' simDistribution()
#'
#' @export
#'

simDistribution <- function(pop = simAbundance(),
                            grid = survey_grid,
                            rho = 0.05,
                            sigma = 0.00001) {


  ## Generate spatially correlated errors
  coords <- grid@data[, c("easting", "northing", "depth")]
  coords$depth <- - coords$depth/1000 # convert depth to km to make units equal to eastings and northings
  if (requireNamespace("fields", quietly = TRUE)) {
    d <- fields::rdist(coords)
  } else {
    d <- as.matrix(dist(coords))
  }
  w <- exp(-rho * d); rm(d)
  W <- chol(w); rm(w)
  e <- (t(W) %*% rnorm(nrow(coords))) * sigma; rm(W)


  ## Distribute abundance equally through the domaine
  N_array <- array(dim = c(nrow(N), ncol(N), nrow(grid)),
                   dimnames = list(age = rownames(N), year = colnames(N), cell = grid$cell))
  prop <- grid$area / sum(grid$area)
  for(i in seq_along(grid$cell)) {
    N_array[, , i] <- (N/grid$area[i]) * prop[i] * exp(e[i])
  }



  grid_dat <- ggplot2::fortify(grid)
  names(grid_dat)[names(grid_dat) == "id"] <- "cell"
  grid_dat$cell <- as.numeric(grid_dat$cell)
  grid@data$N <- N_array[1, 1, ]
  grid@data$den <- grid@data$N/grid@data$area
  grid@data$error <- e
  grid_dat <- merge(grid_dat, grid@data, by = "cell")
  ggplot(grid_dat) +
    geom_polygon(aes(x = long, y = lat, group = cell, fill = N, colour = N))



}





#' Simulate an age-structured population over a survey grid.
#'
#' @param ages Ages to include in the simulated population.
#' @param grid Survey grid to populate with simulated data.
#' @return A 'data.frame' containing ages and lengths for every fish in the
#'     simulated population.
#' @examples
#' simPop()

simPop <- function(ages = 1:6,
                   grid = survey_grid) {
  print("In progress")
}

