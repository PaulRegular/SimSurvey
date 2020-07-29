
#' Make a depth stratified survey grid
#'
#' This function sets up a depth stratified survey grid. A simple gradient in depth
#' is simulated using \code{\link{spline}} with a shallow portion, shelf and
#' deep portion. Adding covariance to the depth simulation is an option.
#'
#' @param x_range      Range (min x, max x) in x dimension in km
#' @param y_range      Range (min y, max y) in y dimension in km
#' @param res          Resolution, in km, of the grid cells
#' @param shelf_depth  Approximate depth of the shelf in m
#' @param shelf_width  Approximate width of the shelf in km
#' @param depth_range  Range (min depth, max depth) in depth in m
#' @param n_div        Number of divisions to include
#' @param strat_breaks Define strata given these depth breaks
#' @param strat_splits Number of times to horizontally split strat (i.e. easy way to increase the number of strata)
#' @param method       Use a "spline" or "loess" to generate a smooth gradient or simply use "linear" interpolation?
#'
#' @return Returns RasterBrick of the same structure as \code{\link{survey_grid}}
#'
#' @export
#'
#' @examples
#'
#' r <- make_grid()
#' raster::plot(r)
#'
#' p <- raster::rasterToPolygons(r$strat, dissolve = TRUE)
#' sp::plot(p)
#'
#' @rawNamespace import(raster, except = select)
#'

make_grid <- function(x_range = c(-140, 140), y_range = c(-140, 140),
                      res = c(3.5, 3.5), shelf_depth = 200,
                      shelf_width = 100, depth_range = c(0, 1000),
                      n_div = 1, strat_breaks = seq(0, 1000, by = 40),
                      strat_splits = 2, method = "spline") {

  cell <- NULL

  ## set-up raster
  r <- raster::raster(xmn = x_range[1], xmx = x_range[2],
                      ymn = y_range[1], ymx = y_range[2],
                      res = res, crs = "+proj=utm +units=km")
  xy <- as.data.frame(raster::rasterToPoints(r))

  ## simulate depth
  sx <- c(x_range[1], -shelf_width, -shelf_width / 2,
          0, shelf_width / 2, shelf_width, x_range[2])
  sy <- c(depth_range[1], rep(shelf_depth, 5), depth_range[2])

  if (method == "loess") {
    lo <- stats::loess(sy ~ sx)
    px <- seq(min(sx), max(sx), length.out = 100)
    py <- predict(lo, data.frame(sx = px))
    s <- list(x = px, y = py)
  }
  if (method == "spline") {
    s <- stats::spline(sx, sy, n = nrow(xy))
  }
  if (method == "linear") {
    s <- stats::approx(sx, sy, n = nrow(xy))
  }

  depth <- s$y[findInterval(xy$x, s$x)]
  depth[depth < depth_range[1]] <- depth_range[1] + 1 # impose depth range
  depth[depth > depth_range[2]] <- depth_range[2] - 1
  depth <- round(depth) # 1 m res ought to be good
  xy$depth <- depth

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

  ## identify isolated clumps, and split into independent strata
  rstrat <- r$strat
  max_id <- 0
  for (s in sort(unique(rstrat))) {
    rc <- suppressMessages(raster::clump(rstrat == s))
    new_ids <- rc[!is.na(rc[])] + max_id
    r$strat[!is.na(rc)] <- new_ids
    max_id <- max(new_ids)
  }

  ## split strat
  xy <- data.table::data.table(raster::rasterToPoints(r))
  if (strat_splits > 1) {
    xy[, split := as.numeric(cut(cell, breaks = strat_splits)), by = "strat"]
    xy$strat <- (xy$split * 10000) + xy$strat
    xy$strat <- as.numeric(as.factor(xy$strat))
    xy$split <- NULL
  }

  ## make strata unique by division
  xy$strat <- (xy$division * 100000) + xy$strat
  xy$strat <- as.numeric(as.factor(xy$strat))

  ## convert to raster
  r <- raster::rasterFromXYZ(xy, crs = "+proj=utm +units=km")
  r

}
