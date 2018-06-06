
#' Simple plotting functions
#'
#' @description These functions are simple plotting helpers to get some quick
#' visuals of values produced by \code{\link{sim_abundance}},
#' \code{\link{sim_distribution}}, etc.
#'
#' @param sim    Object returned by \code{\link{sim_abundance}},
#'               \code{\link{sim_distribution}}, etc.
#' @param mat    Name of matrix in \code{sim} list to plot.
#' @param ages   Subset data to one or more ages.
#' @param years  Subset data to one or more years.
#' @param scale  Plot response on "natural" or "log" scale?
#' @param ...    Additional arguments to pass to \code{\link{plotly::plot_ly}}.
#'
#' @import plotly
#'
#' @export
#' @rdname plot_trend
plot_trend <- function(sim, ...) {
  tot <- data.frame(Year = sim$years, N = colSums(sim$N))
  tot$text <- paste0("Year: ", tot$Year, "\nN: ", round(tot$N))
  plot_ly(x = ~Year, y = ~N, text = ~text, hoverinfo = "text",
          type = "scatter", mode = "lines", data = tot, ...)
}

#' @export
#' @rdname plot_trend
plot_surface <- function(sim, mat = "N", ...) {
  plot_ly(x = sim$ages, y = sim$years, z = sim[[mat]], type = "surface",
          ...) %>%
    layout(
      scene = list(
        xaxis = list(title = "Age"),
        yaxis = list(title = "Year"),
        zaxis = list(title = mat)
      ))
}

#' @export
#' @rdname plot_trend
plot_distribution <- function(sim, ages = 1, years = 1,
                              type = "contour", axes = TRUE,
                              scale = "natural", ...) {

  ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE,
    ticks = "inside"
  )

  xyz <- data.frame(merge(sim$grid_xy, sim$sp_N, by = "cell"))
  if (scale == "log") xyz$N <- log(xyz$N)
  p <- vector("list", length(ages) * length(years))
  g <- expand.grid(age = ages, year = years)
  for (i in seq(nrow(g))) {
    sub_xyz <- xyz[xyz$age == g$age[i] & xyz$year == g$year[i], ]
    p[[i]] <- plot_ly(x = ~x, y = ~y, z = ~N, data = sub_xyz, type = type,
                      name = paste0("a", g$age[i], ", y", g$year[i]), ...)
    if (!axes) p[[i]] <- p[[i]] %>% layout(xaxis = ax, yaxis = ax)
  }
  subplot(p, nrows = length(ages), shareX = TRUE, shareY = TRUE)

}


