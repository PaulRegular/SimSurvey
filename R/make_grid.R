#' Make a depth-stratified survey grid
#'
#' This function sets up a depth-stratified survey grid. A simple gradient in depth
#' is simulated using [`stats::spline()`] (default), with a shallow portion, shelf, and
#' deep portion. Optionally, covariance can be added to the depth simulation.
#'
#' @param x_range Range (min x, max x) in the x dimension (km).
#' @param y_range Range (min y, max y) in the y dimension (km).
#' @param res Resolution of the grid cells (km).
#' @param shelf_depth Approximate depth of the shelf (m).
#' @param shelf_width Approximate width of the shelf (km).
#' @param depth_range Range (min depth, max depth) of the depth values (m).
#' @param n_div Number of divisions to include.
#' @param strat_breaks Depth breaks used to define strata.
#' @param strat_splits Number of times to horizontally split strata (i.e., a way to increase the number of strata).
#' @param method Choose `"spline"`, `"loess"`, or `"bezier"` to generate a smooth gradient, or use `"linear"` for linear interpolation.
#'
#' @return A `stars` object with 2 dimensions (`x` and `y`) and 4 attributes (`depth`, `cell`, `division`, `strat`).
#'
#' @seealso [`survey_grid`]
#'
#' @export
#'
#' @examples
#' r <- make_grid(res = c(10, 10))
#' plot(r)
#'
#' p <- sf::st_as_sf(r["strat"], as_points = FALSE, merge = TRUE)
#' plot(p)
#'
#' @import sf
#' @import stars

make_grid <- function(x_range = c(-140, 140), y_range = c(-140, 140),
                      res = c(3.5, 3.5), shelf_depth = 200,
                      shelf_width = 100, depth_range = c(0, 1000),
                      n_div = 1, strat_breaks = seq(0, 1000, by = 40),
                      strat_splits = 2, method = "spline") {

  cell <- NULL

  ## set-up raster
  utm_proj <- "+proj=utm +zone=21 +datum=WGS84 +units=km +no_defs"
  bbox <- sf::st_bbox(c(xmin = x_range[1], ymin = y_range[1],
                        xmax = x_range[2], ymax = y_range[2]),
                      crs = utm_proj)
  r <- stars::st_as_stars(bbox, dy = res[1], dx = res[2])
  xy <- as.data.frame(r)
  names(xy) <- c("x", "y", "depth")

  ## simulate depth
  sx <- c(x_range[1], -shelf_width, -shelf_width / 2,
          0, shelf_width / 2, shelf_width, x_range[2])
  sy <- c(depth_range[1], rep(shelf_depth, 5), depth_range[2])

  if (method == "loess") {
    lo <- stats::loess(sy ~ sx)
    px <- seq(min(sx), max(sx), length.out = 100)
    py <- stats::predict(lo, data.frame(sx = xy$x))
    s <- list(x = xy$x, y = py)
  }
  if (method == "spline") {
    s <- stats::spline(sx, sy, xout = xy$x)
  }
  if (method == "linear") {
    s <- stats::approx(sx, sy, xout = xy$x)
  }
  if (method == "bezier") {
    if (requireNamespace("bezier", quietly = TRUE)) {
      t <- seq(0, 1, length = 100)
      p <- cbind(sx, sy)
      s <- bezier::bezier(t = t, p = p)
      s <- stats::approx(s[, 1], s[, 2], xout = xy$x)
    } else {
      stop("The bezier package is needed for to use the bezier method in make_grid. Please install it.", call. = FALSE)
    }
  }

  depth <- s$y
  depth[depth < depth_range[1]] <- depth_range[1] + 1 # impose depth range
  depth[depth > depth_range[2]] <- depth_range[2] - 1
  depth <- round(depth) # 1 m res ought to be good
  xy$depth <- depth

  ## add cell number
  xy$cell <- seq(nrow(xy))

  ## define divisions and/or splits
  xy$division <- xy$split <- 1
  if (n_div > 1) {
    xy$division <- as.numeric(cut(xy$y, n_div))
  }
  if (strat_splits > 1) {
    xy$split <- as.numeric(cut(xy$y, strat_splits))
  }

  ## define strata and use unique labels along the x-axis
  ## (i.e. don't duplicate strat labels for areas with the same depth
  ## that are not adjacent to each other)
  xy$strat <- as.numeric(cut(xy$depth, strat_breaks))
  xy <- xy[order(xy$x), ]
  rl <- rle(xy$strat)$length
  xy$strat <- rep(seq_along(rl), rl)

  ## ensure strata within a division or a split have different numbers
  div_id <- xy$division * 10 ^ (max(nchar(xy$split)) + max(nchar(xy$strat)))
  split_id <- xy$split * 10 ^ max(nchar(xy$strat))
  xy$strat <- as.numeric(factor(div_id + split_id + xy$strat))

  ## convert xyz data to a raster
  xy$split <- NULL # drop - was useful for defining strata
  r <- stars::st_as_stars(xy, crs = utm_proj)
  r

}
