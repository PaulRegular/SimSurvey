
## helper function for squishing values between a new range
## ...perhaps a sloppy work-around to get what I want
.range <- function(x, to = c(-1, 1)) {
  (x - min(x))/(max(x) - min(x)) * (to[2] - to[1]) + to[1]
}


#' Simulate depth stratified survey grid
#'
#' This function sets up a depth stratified survey grid. A simple gradient in depth
#' is simulated with a shallow portion, shelf and deep portion. Adding covariance
#' to the depth simulation is an option.
#'
#' @param x_range      Range (min x, max x) in x dimension in km
#' @param y_range      Range (min y, max y) in y dimension in km
#' @param res          Resolution, in km, of the grid cells
#' @param depth_range  Range (min depth, max depth) in depth in m
#' @param n_div        Number of divisions to iniclude
#' @param strat_breaks Define strata given these depth breaks
#' @param space_covar  Supply \code{\link{sim_sp_covar}} closure to add covariance to depth
#'                     to make it look more realistic.
#'
#' @return Returns RasterBrick of the same structure as \code{\link{survey_grid}}
#'
#' @export
#'
#' @examples
#' r <- sim_grid()
#' raster::plot(r)
#'
#' @import raster
#'

sim_grid <- function(x_range = c(-150, 150), y_range = c(-150, 150),
                     res = c(3.5, 3.5), depth_range = c(1, 500),
                     n_div = 2, strat_breaks = seq(0, 500, by = 20),
                     space_covar = NULL) {

  # sim_sp_covar(range = 500, sd = 0.2)

  ## set-up raster
  r <- raster::raster(xmn = x_range[1], xmx = x_range[2],
              ymn = y_range[1], ymx = y_range[2],
              res = res, crs = "+proj=utm +units=km")
  xy <- as.data.frame(raster::rasterToPoints(r))

  ## simulate depth
  if (!is.null(space_covar)) {
    Sigma_space <- space_covar(xy)
    w <- t(chol(Sigma_space))
    e <- w %*% rnorm(nrow(xy))
  } else {
    e <- 0
  }
  x <- .range(xy$x)
  y <- .range(xy$y)
  xy$depth <- x ^ 3 + e
  xy$depth <- .range(xy$depth, to = depth_range)

  ## add cell number
  xy$cell <- seq(nrow(xy))

  ## define divisions and strata
  r <- raster::rasterFromXYZ(xy, crs = "+proj=utm +units=km")
  if (n_div == 1) {
    r$division <- r$cell
    r$division[] <- 1
  } else {
    r$division <- raster::cut(r$cell, breaks = n_div)
  }
  r$strat <- raster::cut(r$depth, breaks = strat_breaks)

  ## make strata unique by division
  xy <- as.data.frame(raster::rasterToPoints(r))
  xy$strat <- (xy$division * 100000) + xy$strat
  xy$strat <- as.numeric(as.factor(xy$strat))
  r <- raster::rasterFromXYZ(xy, crs = "+proj=utm +units=km")
  r

}

