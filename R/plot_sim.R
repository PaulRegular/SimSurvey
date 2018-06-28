
#' Simple plotting functions
#'
#' @description These functions are simple plotting helpers to get some quick
#' visuals of values produced by \code{\link{sim_abundance}},
#' \code{\link{sim_distribution}}, etc.
#'
#' @param sim            Object returned by \code{\link{sim_abundance}},
#'                       \code{\link{sim_distribution}}, etc.
#' @param mat            Name of matrix in \code{sim} list to plot.
#' @param ages           Subset data to one or more ages.
#' @param years          Subset data to one or more years.
#' @param scale          Plot response on "natural" or "log" scale?
#' @param which_survey   Subset to specific survey
#' @param which_year     Subset to specific year
#' @param which_sim      Subset to specific sim
#' @param max_sims       Maximum number of sims to plot
#' @param facet_by       Facet plot by "age" or "year"?
#' @param col            Plot color
#' @param cols           Plot colors
#' @param ...            Additional arguments to pass to \code{\link{plotly::plot_ly}}.
#'
#' @import plotly
#'
#' @export
#' @rdname plot_trend
plot_trend <- function(sim, col = viridis::viridis(1), ...) {
  tot <- data.frame(Year = sim$years, N = colSums(sim$N))
  tot$text <- paste0("Year: ", tot$Year, "\nN: ", round(tot$N))
  plot_ly(x = ~Year, y = ~N, text = ~text, hoverinfo = "text",
          type = "scatter", mode = "lines", color = I(col), data = tot, ...)
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

  xax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE,
    tickcolor = toRGB("white")
  )
  yax <- c(scaleanchor = "x", xax)

  xyz <- data.frame(merge(sim$grid_xy, sim$sp_N, by = "cell"))
  if (scale == "log") xyz$N <- log(xyz$N)
  p <- vector("list", length(ages) * length(years))
  g <- expand.grid(age = ages, year = years)
  for (i in seq(nrow(g))) {
    sub_xyz <- xyz[xyz$age == g$age[i] & xyz$year == g$year[i], ]
    p[[i]] <- plot_ly(x = ~x, y = ~y, z = ~N, data = sub_xyz, type = type,
                      name = paste0("a", g$age[i], ", y", g$year[i]), ...)
    if (!axes) {
      p[[i]] <- p[[i]] %>% layout(xaxis = xax, yaxis = yax)
    } else {
      p[[i]] <- p[[i]] %>% layout(yaxis = list(scaleanchor = "x"))
    }
  }
  subplot(p, nrows = length(ages), shareX = TRUE, shareY = TRUE)

}


#' @export
#' @rdname plot_trend
plot_error_surface <- function(sim) {

  d <- merge(sim$surveys, sim$age_strat_error_stats, by = "survey")
  split_d <- split(d, d$set_den)
  split_d <- lapply(split_d, function(.) {
    xtabs(RMSE ~ ages_cap + lengths_cap, data = ., subset = NULL)
  })
  x <- sort(unique(d$lengths_cap))
  y <- sort(unique(d$ages_cap))

  p <- plot_ly(x = ~x, y = ~y, cmin = min(d$RMSE), cmax = max(d$RMSE))
  visible <- showscale <- rep(TRUE, length(split_d))
  steps <- list()
  showscale <- c(TRUE, rep(FALSE, length(split_d) - 1))
  for (i in seq_along(split_d)) {
    vis <- rep(FALSE, length(split_d))
    vis[i] <- TRUE
    p <- p %>% add_surface(z = split_d[[i]],
                           visible = i == 1,
                           showscale = vis,
                           name = names(split_d)[i],
                           colorbar = list(title = "RMSE"))
    steps[[i]] <- list(args = list(list(visible = vis,
                                        showscale = vis)),
                       method = "update",
                       label = names(split_d)[i])
    if (i == length(split_d)) {
      steps[[i + 1]] <- list(args = list(list(visible = rep(TRUE, length(split_d)),
                                              showscale = showscale)),
                             method = "update",
                             label = "all")
    }
  }

  p %>%
    layout(sliders = list(list(
      currentvalue = list(prefix = "Set density: "),
      steps = steps
    )),
    scene = list(
      xaxis = list(title = "max(lengths)"),
      yaxis = list(title = "max(ages)"),
      zaxis = list(title = "RMSE",
                   range = range(d$RMSE)),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    ))

}


#' @export
#' @rdname plot_trend
plot_error_cross_sections <- function(sim, col = viridis::viridis(1)) {

  d <- merge(sim$surveys, sim$age_strat_error_stats, by = "survey")

  a <- list(
    text = "",
    font = list(size = 18, color = "black"),
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 1,
    showarrow = FALSE
  )

  a$text <- "Change set density"
  set_p <- d %>%
    mutate(combo = paste("max(lengths):", lengths_cap, "<br>max(ages):", ages_cap)) %>%
    group_by(combo) %>%
    plot_ly() %>% add_lines(x = ~set_den, y = ~RMSE, text = ~combo,
                            name = "Change set density", showlegend = FALSE,
                            color = I(col), alpha = 0.5, size = I(0.5)) %>%
    layout(xaxis = list(title = "Set density"), annotations = a)

  a$text <- "Change length sampling"
  len_p <- d %>%
    mutate(combo = paste("set density:", set_den, "<br>max(ages):", ages_cap)) %>%
    group_by(combo) %>%
    plot_ly() %>% add_lines(x = ~lengths_cap, y = ~RMSE, text = ~combo,
                            name = "Change length sampling", showlegend = FALSE,
                            color = I(col), alpha = 0.5, size = I(0.5)) %>%
    layout(xaxis = list(title = "max(lengths)"), annotations = a)

  a$text <- "Change age sampling"
  age_p <- d %>%
    mutate(combo = paste("set density:", set_den, "<br>max(lengths):", lengths_cap)) %>%
    group_by(combo) %>%
    plot_ly() %>% add_lines(x = ~ages_cap, y = ~RMSE, text = ~combo,
                            name = "Change age sampling", showlegend = FALSE,
                            color = I(col), alpha = 0.8, size = I(0.5)) %>%
    layout(xaxis = list(title = "max(ages)"), annotations = a)

  subplot(set_p, len_p, age_p, nrows = 1, shareY = TRUE, titleX = TRUE)

}


#' @export
#' @rdname plot_trend
plot_true_vs_est <- function(sim, which_survey = 1, max_sims = 50,  facet_by = "age",
                             cols = viridis::viridis(3)) {

  d <- sim$age_strat_error
  sub_d <- d[d$survey == which_survey, ]
  sub_sims <- sample(unique(d$sim), ifelse(max(d$sim) > max_sims, max_sims, max(d$sims)))
  sub_d <- sub_d[sub_d$sim %in% sub_sims, ]

  if (facet_by == "age") {
    p <- ggplot(data = sub_d, aes(x = year, group = sim)) +
      geom_line(aes(y = I_hat), color = I(cols[2]), alpha = 0.2, size = 0.1) +
      geom_line(aes(y = I), color = I(cols[1]), size = 0.1) +
      xlab("Year") + ylab("Index") +
      facet_wrap(~ age, scales = "free_y") +
      theme_minimal() + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())
  } else {
    p <- ggplot(data = sub_d, aes(x = age, group = sim)) +
      geom_line(aes(y = I_hat), color = I(cols[2]), alpha = 0.2, size = 0.1) +
      geom_line(aes(y = I), color = I(cols[1]), size = 0.1) +
      xlab("Age") + ylab("Index") +
      facet_wrap(~ year, scales = "free_y") +
      theme_minimal() + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())
  }

  ggplotly(p)

}

#' @export
#' @rdname plot_trend
plot_samp_dist <- function(sim, which_year = 1, which_sim = 1,
                           col = viridis::viridis(1)) {

  xax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE,
    tickcolor = toRGB("white")
  )
  yax <- c(scaleanchor = "x", xax)

  sp_strat <- raster::rasterToPolygons(sim$grid$strat, dissolve = TRUE)
  df_strat <- suppressWarnings(fortify(sp_strat) %>% group_by(group))

  setdet <- sim$setdet
  setdet <- setdet[setdet$year == which_year & setdet$sim == which_sim, ]
  samp <- merge(setdet[, c("year", "sim", "set", "x", "y")],
                sim$samp, by = "set")

  d <- crosstalk::SharedData$new(samp, ~set)

  base <- plot_ly(data = d)

  sp_p <- base %>%
    group_by(set) %>%
    summarise(x = unique(x), y = unique(y), n = n()) %>%
    add_markers(x = ~x, y = ~y, size = ~n, text = ~n,
                color = I(col), sizes = c(5, 500), name = "n",
                showlegend = FALSE) %>%
    add_markers(data = setdet[setdet$n == 0, ], color = I(col),
                x = ~x, y = ~y, text = ~n, size = I(5), symbol = I(4),
                name = "zero", alpha = 0.5,
                showlegend = FALSE) %>%
    add_paths(data = df_strat, x = ~long, y = ~lat, color = I("black"),
              hoverinfo = "none", size = I(0.5), showlegend = FALSE,
              alpha = 0.1) %>%
    layout(xaxis = xax, yaxis = yax)

  lf_p <- base %>%
    group_by(set) %>%
    filter(measured) %>%
    add_histogram(x = ~length, showlegend = FALSE,
                  name = "Length distribution",
                  marker = list(color = toRGB(col),
                                line = list(color = toRGB("white"),
                                            width = 0.2))) %>%
    layout(xaxis = list(title = "Length"),
           yaxis = list(title = ""))

  af_p <- base %>%
    group_by(set) %>%
    filter(aged) %>%
    add_histogram(x = ~age, showlegend = FALSE,
                  name = "Age distribution",
                  marker = list(color = toRGB(col),
                                line = list(color = toRGB("white"),
                                            width = 0.2))) %>%
    layout(xaxis = list(title = "Age"),
           yaxis = list(title = ""))

  subplot(sp_p,
          subplot(lf_p, af_p, nrows = 2, titleX = TRUE, titleY = TRUE,
                  margin = 0.1),
          nrows = 1, margin = 0.02, widths = c(0.6, 0.4),
          titleX = TRUE, titleY = TRUE)

}
