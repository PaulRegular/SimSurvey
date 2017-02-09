
#' Simulate recruitment and natural mortality
#'
#' @description These functions return a function to use inside \code{\link{sim_abundance}}.
#' Given parameters, it generates R and Z values.
#'
#' @param mean,sd Mean and standard deviation
#' @param breaks Provide breaks to age and/or year series to specify group specific
#' mean total mortality values. If specified, multiple mean values must be
#' supplied; one more than the number of breaks. Provide named lists (see
#' examples below).
#'
#' @details \code{sim_R} simply generates uncorrelated recruitment values from a log normal
#' distribution. \code{sim_Z} simulates total mortality using a random walk.
#' Random walk errors are first generated for years and ages independently.
#' The outer product of these errors are then calculated to generate a total
#' mortality matrix as follows:
#' \eqn{log(Z_{a, y}) = log(z_{a, y}) + \delta_{a, y}}{log(Z_a,y) = log(z_a,y) + e_a,y},
#' where \eqn{z_{a, y}}{z_a,y} are the supplied mean values.
#'
#' @examples
#' sim_abundance(R = sim_R(mean = 100000, sd = 4))
#' sim_abundance(Z = sim_Z(mean = 0.6, sd = 0.3))
#' sim_abundance(Z = sim_Z(mean = list(ages = c(0.5, 0.3, 0.2)),
#'                       breaks = list(ages = c(1, 2)), sd = 0.3))
#' sim_abundance(Z = sim_Z(mean = list(ages = c(0.5, 0.3, 0.2), years = c(0.3, 2, 0.3)),
#'                       breaks = list(ages = c(1, 2), years = c(4, 6)), sd = 0.3))
#'
#' @export
#' @rdname sim_R
sim_R <- function(mean = 20000, sd = 4) {
  function(years = NULL) {
    r <- rlnorm(length(years), meanlog = log(mean), sdlog = log(sd))
    names(r) <- years
    r
  }
}


#' @export
#' @rdname sim_R
sim_Z <- function(mean = 0.4, sd = 0.2, breaks = NULL) {
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
#' @param Z Total mortality function, like \code{\link{sim_Z}}, for generating
#' mortality matrix.
#' @param R Recruitment (i.e. Abundance at \code{min(ages)}) function, like
#' \code{\link{sim_R}}, for generating recruitment vector.
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
#' sim_abundance()
#'
#' @export

sim_abundance <- function(ages = 1:6, years = 1:10, Z = sim_Z(), R = sim_R()) {

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
  list(ages = ages, years = years, R = R, Z = Z, N = N)

}



## Helper function for euclidian distance calculations
.dist <- function(x) {
  if (requireNamespace("fields", quietly = TRUE)) {
    d <- fields::rdist(x)
  } else {
    d <- as.matrix(dist(x))
  }
}

#' Simulate size, time and space correlation structure
#'
#' @description These function returns a function to use inside \code{\link{sim_distribution}}.
#' Given parameters, it generates recruitment values from a log normal distribution
#'
#' @param dcor_time Decorrelation time (years)
#' @param dcor_size Decorrelation size (age)
#' @param dcor_dist Decorrelation distance (km)
#' @param sigma     Spatial variance
#' @param tau,theta Precision matrix parameters
#'
#' @details The formulation of these functions follow Kristensen et al. (2013).
#'
#' @references Kristensen, K., Thygesen, U. H., Andersen, K. H., & Beyer, J. E.
#' (2013). Estimating spatio-temporal dynamics of size-structured populations.
#' Canadian Journal of Fisheries and Aquatic Sciences, 71(2), 326-336.
#'
#' @rdname sim_time_cor
#' @export
sim_time_cor <- function(dcor_time = 2) {
  function(years = NULL) {
    d <- .dist(years)
    exp(- d / dcor_time)
  }
}

#' @rdname sim_time_cor
#' @export
sim_size_cor <- function(dcor_size = 4) {
  function(ages = NULL) {
    d <- .dist(ages)
    exp(- d / dcor_size)
  }
}

# #' @rdname sim_time_cor
# #' @export
# sim_space_cor <- function(dcor_dist = 70, sigma = 1.5) {
#   function(grid = NULL) {
#     d <- .dist(coordinates(grid))
#     sigma * exp(- d / dcor_dist)
#   }
# } # Alternate formulation to precision matrix approach.
#   # Con: correlation as the crow flys, not as the fish swims
#          i.e. land is not a barrier

#' @rdname sim_time_cor
#' @export
sim_space_cor <- function(tau = 0.25, theta = 0.001) {
  function(grid = NULL) {

    d <- .dist(sp::coordinates(grid))
    nb <- Q <- rgeos::gTouches(grid, byid = TRUE) # cell neighbour matrix
    m <- rowSums(nb)                              # number of neighbouring cells
    Q[] <- 0
    Q[nb] <- - tau
    diag(Q) <- tau * (m + theta)
    # invQ <- solve(Q)
    # invQ <- chol2inv(chol(Q))
    # s <- (1 / length(grid)) * sum(diag(invQ))
    # H <- mean(d[nb]) / log(1 + (theta / 2) + sqrt(theta + ((theta ^ 2)/4)))
    # message(paste0("Spatial variance is approxamatly ", signif(s, 3),
    #                "\nSpatial decorrelation distance is approxamatly ", signif(H, 3), " km"))

    # # dcor_dist <- H; sigma <- s; invQ <- sigma * exp(- d / dcor_dist)
    # cell <- sample(grid$cell, 1)
    # # cell <- 8179
    # ncols <- 200
    # cols <- cut(invQ[cell, ], breaks = ncols, labels = FALSE)
    # cols <- colorRampPalette(c("white", "steelblue", "navy"))(ncols)[cols]
    # plot(d[cell, ], invQ[cell, ], col = cols, pch = 16, cex = 0.75,
    #      xlab = "Distance", ylab = "Correlation")
    # plot(grid, col = cols, lwd = 0.5, border = NA)

    # invQ
    Q

  }
}





#' Simulate spatial and temporal distribution
#'
#' @description Provided an abundance at age matrix (like one provided by \code{\link{sim_abundance}})
#' and a survey grid (like \code{\link{survey_grid}}) to populate, this function
#' applies correlated space, time and size error to simulate the distribution
#' of the population.
#'
#' @param N An abundance at age matrix with ages defining the rows and years defining
#' the columns (i.e. same structure as a matrix provided by \code{\link{sim_abundance}})
#' @param grid A \code{\link{SpatialPolygonsDataFrame}} defining a regular or irregular
#' grid with the same structure as \code{\link{survey_grid}}
#'
#' @examples
#' sim_distribution()
#'
#' @export
#'

sim_distribution <- function(pop = sim_abundance(),
                             grid = survey_grid,
                             size_cor  = sim_size_cor(),
                             time_cor  = sim_time_cor(),
                             space_cor = sim_space_cor()) {


  ## solution 1
  pop$years <- 1:10
  pop$ages <- 1:6
  test <- expand.grid(year = pop$years, age = pop$ages)
  rho_size <- size_cor(ages = test$age)
  rho_time <- time_cor(years = test$year)
  rho <- rho_size * rho_time
  test$e <- t(chol(rho)) %*% rnorm(nrow(test))
  plot(age ~ year, data = test, cex = e - min(e))
  #plot(e ~ year, data = test[test$age == 3, ], type = "b")
  e1 <- sd(test$e)



  sigma_size <- size_cor(ages = pop$ages)
  sigma_time <- time_cor(years = pop$years)
  sigma_space <- space_cor(grid = grid)
  sigma_size <- Matrix(solve(sigma_size), sparse = TRUE)
  sigma_time <- Matrix(solve(rho_time), sparse = TRUE)
  sigma_space <- Matrix::Matrix(sigma_space, sparse = TRUE)

  test1 <- kronecker(sigma_space, sigma_size)


  e <- solve(chol(sigma_space)) %*% rnorm(nrow(sigma_space))
  e <- e * 4

  library(sparseMVN)
  space_cor <- sim_space_cor(tau = 1, theta = 0.01)
  sigma_space <- space_cor(grid = grid)
  sigma_space <- Matrix::Matrix(sigma_space, sparse = TRUE)
  e <- rmvn.sparse(nrow(sigma_space), rep(0, nrow(sigma_space)), Cholesky(sigma_space))
  ncols <- 200
  cols <- cut(e[, 1], breaks = ncols, labels = FALSE)
  cols <- colorRampPalette(c("white", "steelblue", "navy"))(ncols)[cols]
  plot(grid, col = cols, lwd = 0.5, border = NA)



  spt <- kronecker(sp, t)
  test <- kronecker(spt, s)

  ctest <- Cholesky(test)

  temp <- solve(test)



  rho <- outer(rho_size, rho_time)
  rho <- outer(rho_space, rho)

  e_size <- data.frame(age = pop$ages, ea = t(chol(rho_size)) %*% rnorm(length(pop$ages)))
  e_time <- data.frame(year = pop$years, ey = t(chol(rho_time)) %*% rnorm(length(pop$years)))
  test <- merge(test, e_size, by = "age", all.x = TRUE)
  test <- merge(test, e_time, by = "year", all.x = TRUE)
  test$e <- test$ea * test$ey
  #plot(test$e, col = test$age)
  plot(age ~ year, data = test, cex = e - min(e))
  #plot(e ~ year, data = test[test$age == 3, ], type = "b")
  e2 <- sd(test$e)
  e1
  e2




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

