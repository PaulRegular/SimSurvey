
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
#' @param lengths        Subset data to one or more length groups.
#' @param years          Subset data to one or more years.
#' @param scale          Plot response on "natural" or "log" scale?
#' @param which_survey   Subset to specific survey
#' @param which_year     Subset to specific year
#' @param which_sim      Subset to specific sim
#' @param max_sims       Maximum number of sims to plot
#' @param facet_by       Facet plot by "age" or "year"?
#' @param select_by      Select plot by "age", "length" or "year"?
#' @param plot_by        Plot error surface by "rule" or "samples"?
#' @param survey         Subset data to one or more surveys.
#' @param quants         Quantile intervals to display on fan plot
#' @param col            Plot color
#' @param cols           Plot colors
#' @param main           Plot title
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
  plot_ly(x = sim$ages, y = sim$years, z = t(sim[[mat]]), type = "surface",
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
    tickcolor = "transparent"
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
plot_distribution_slider <- function(sim, ages = 1, years = 1) {

  xax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE,
    tickcolor = "transparent"
  )
  yax <- c(scaleanchor = "x", xax)

  d <- merge(sim$grid_xy, sim$sp_N[age %in% ages & year %in% years, ], by = "cell")
  d$ay <- paste(d$age, d$year, sep = "-")
  split_d <- split(d, d$ay)
  split_d <- lapply(split_d, function(.) {
    xtabs(N ~ x + y, data = ., subset = NULL)
  })
  x <- sort(unique(d$x))
  y <- sort(unique(d$y))

  p <- plot_ly(x = ~x, y = ~y)
  visible <- rep(TRUE, length(split_d))
  showscale <- rep(FALSE, length(split_d))
  showscale[1] <- TRUE
  steps <- list()
  for (i in seq_along(split_d)) {
    vis <- rep(FALSE, length(split_d))
    vis[i] <- TRUE
    p <- p %>% add_trace(type = "contour",
                         z = split_d[[i]],
                         visible = i == 1,
                         showscale = vis,
                         name = names(split_d)[i],
                         colorbar = list(title = "N")) %>%
      layout(xaxis = xax, yaxis = yax)

    steps[[i]] <- list(args = list(list(visible = vis,
                                        showscale = vis)),
                       method = "update",
                       label = names(split_d)[i])
  }

  p %>%
    layout(sliders = list(list(
      currentvalue = list(prefix = "Age-Year: "),
      steps = steps
    )))

}


#' @export
#' @rdname plot_trend
plot_error_surface <- function(sim, plot_by = "rule") {

  totals <- sim$samp_totals[, list(n_sets = mean(n_sets), n_caught = mean(n_caught),
                                   n_measured = mean(n_measured), n_aged = mean(n_aged)),
                            by = "survey"]
  errors <- merge(sim$surveys, sim$age_strat_error_stats, by = "survey")
  d <- merge(errors, totals, by = "survey")

  d$text <- with(d, paste("<br>set density:", set_den,
                          "<br>lengths cap:", lengths_cap,
                          "<br>ages cap:", ages_cap,
                          "<br><br>N sets:", n_sets,
                          "<br>N lengths:", round(n_measured),
                          "<br>N ages:", round(n_aged)))

  if (plot_by == "rule") {
    split_d <- split(d, d$set_den)
    x <- sort(unique(d$ages_cap))
    y <- sort(unique(d$lengths_cap))
    marker_x <- ~ages_cap
    marker_y <- ~lengths_cap
    xlab <- "Ages cap"
    ylab <- "Lengths cap"
    slab <- "Set density"
  }
  if (plot_by == "samples") {
    split_d <- split(d, d$n_sets)
    x <- marker_x <- ~n_aged
    y <- marker_y <- ~n_measured
    xlab <- "N ages"
    ylab <- "N lengths"
    slab <- "N sets"
  }

  p <- plot_ly(x = x, y = y)
  splits <- vector("list", length(split_d) + 1)
  splits[[1]] <- list(args = list(list(visible = rep(c(TRUE, FALSE, TRUE), length(split_d)),
                                       showscale = c(TRUE, rep(FALSE, (length(split_d) * 3) - 1)))),
                      method = "update",
                      label = "all")
  showscale <- c(TRUE, rep(FALSE, length(split_d) - 1))

  for (i in seq_along(split_d)) {

    vis <- replicate(length(split_d), c(FALSE, FALSE, FALSE), simplify = FALSE)
    vis[[i]] <- c(FALSE, TRUE, TRUE)
    if (plot_by == "rule") {
      z <- xtabs(RMSE ~ ages_cap + lengths_cap, data = split_d[[i]], subset = NULL)
      p <- p %>% add_surface(z = t(z),
                             cmin = min(d$RMSE), cmax = max(d$RMSE),
                             showscale = i == 1,
                             name = names(split_d)[i],
                             colorbar = list(title = "RMSE"),
                             hoverinfo = "skip") %>%
        add_surface(z = t(z),
                    visible = FALSE,
                    showscale = FALSE,
                    name = names(split_d)[i],
                    colorbar = list(title = "RMSE"),
                    hoverinfo = "skip")
    }
    if (plot_by == "samples") {
      p <- p %>% add_mesh(z = ~RMSE,
                          intensity = ~RMSE,
                          data = split_d[[i]],
                          cmin = min(d$RMSE), cmax = max(d$RMSE),
                          showscale = i == 1,
                          name = names(split_d)[i],
                          flatshading = TRUE,
                          colorbar = list(title = "RMSE"),
                          contour = list(show = TRUE, width = 15, color = toRGB("white")),
                          hoverinfo = "skip") %>%
        add_mesh(z = ~RMSE,
                 intensity = ~RMSE,
                 data = split_d[[i]],
                 visible = FALSE,
                 showscale = FALSE,
                 name = names(split_d)[i],
                 flatshading = TRUE,
                 colorbar = list(title = "RMSE"),
                 contour = list(show = TRUE, width = 15, color = toRGB("white")),
                 hoverinfo = "skip")
    }
    p <- p %>%
      add_markers(z = ~RMSE,
                  x = marker_x,
                  y = marker_y,
                  data = split_d[[i]],
                  text = ~text,
                  hoverinfo = "text+z",
                  name = names(split_d)[i],
                  legendgroup = names(split_d)[i],
                  color = I("grey30"),
                  opacity = 0,
                  visible = TRUE,
                  showlegend = FALSE)
    splits[[i + 1]] <- list(args = list(list(visible = unlist(vis),
                                             showscale = unlist(vis))),
                            method = "update",
                            label = names(split_d)[i])

  }


  p %>%
    layout(
      updatemenus = list(
        list(
          y = 0.9,
          x = 0,
          xanchor = "center",
          yanchor = "top",
          buttons = splits
        )
      ),
      annotations = list(
        list(
          text = slab,
          y = 0.9,
          x = 0,
          borderpad = 5,
          xref = "paper",
          yref = "paper",
          xanchor = "center",
          yanchor = "bottom",
          showarrow = FALSE
        )
      ),
      scene = list(
        xaxis = list(title = xlab),
        yaxis = list(title = ylab),
        zaxis = list(title = "RMSE"),
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
      facet_wrap(~ age, scales = "free_y")
  } else {
    p <- ggplot(data = sub_d, aes(x = age, group = sim)) +
      geom_line(aes(y = I_hat), color = I(cols[2]), alpha = 0.2, size = 0.1) +
      geom_line(aes(y = I), color = I(cols[1]), size = 0.1) +
      xlab("Age") + ylab("Index") +
      facet_wrap(~ year, scales = "free_y")

  }

  p <- p + theme_minimal() + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   plot.margin = unit(c(0.1, 0.1, 0.5, 0.5), "cm"))

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
    tickcolor = "transparent"
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
                name = "zero", alpha = 0.2,
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


#' @export
#' @rdname plot_trend
plot_total_abundance <- function(sim, which_year = 1, main = "") {

  m <- list(
    l = 0,
    r = 0,
    b = 0,
    t = 0,
    pad = 0
  )

  xax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE,
    tickcolor = "transparent"
  )
  yax <- c(scaleanchor = "x", xax)

  sp_N_tot <- sim$sp_N[year == which_year, list(N = sum(N)), by = c("year", "cell")]
  xyz <- merge(sim$grid_xy, sp_N_tot, by = "cell")
  size <- unique(diff(pretty(xyz$N, 25)))

  sp_strat <- raster::rasterToPolygons(sim$grid$strat, dissolve = TRUE)
  df_strat <- suppressWarnings(fortify(sp_strat) %>% group_by(group))

  abun_plot <- plot_ly() %>%
    add_trace(x = ~x, y = ~y, z = ~N, data = xyz, type = "contour",
              contours = list(coloring = "lines", start = 0, end = max(xyz$N),
                              size = size)) %>%
    add_paths(data = df_strat, x = ~long, y = ~lat, color = I("black"),
              hoverinfo = "none", size = I(0.5), showlegend = FALSE,
              alpha = 0.1) %>%
    colorbar(ticklen = 0, ypad = 20, xpad = 0, outlinewidth = 0, thickness = 10) %>%
    layout(xaxis = xax, yaxis = yax, margin = m, legend = list(x = 0.5, y = 0.5),
           annotations = list(text = main, xref = "paper", yref = "paper",
                              x = 0.08, y = 1, showarrow = FALSE))
  abun_plot

}


#' @export
#' @rdname plot_trend
plot_set_catch <- function(sim, which_year = 1, which_sim = 1, main = "") {

  m <- list(
    l = 0,
    r = 0,
    b = 0,
    t = 0,
    pad = 0
  )

  xax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE,
    tickcolor = "transparent"
  )
  yax <- c(scaleanchor = "x", xax)

  sp_strat <- raster::rasterToPolygons(sim$grid$strat, dissolve = TRUE)
  df_strat <- suppressWarnings(fortify(sp_strat) %>% group_by(group))

  setdet <- sim$setdet
  setdet <- setdet[setdet$year == which_year & setdet$sim == which_sim, ]

  breaks <- pretty(setdet$n, 10)
  breaks[breaks == 0] <- 1
  labs <- c("0", breaks)
  breaks[length(breaks)] <- max(setdet$n) + 1
  breaks <- c(0, breaks)
  setdet$n_int <- findInterval(setdet$n, breaks)
  setdet$n_size <- breaks[setdet$n_int]
  setdet$n_cat <- labs[setdet$n_int]
  setdet$n_cat <- factor(setdet$n_cat, levels = labs)

  set_plot <- plot_ly() %>%
    add_markers(data = setdet,
                x = ~x, y = ~y, size = ~n_size, color = ~n_cat, sizes = c(1, 30),
                colors = viridis::viridis(100), alpha = 0.8, text = ~n) %>%
    add_paths(data = df_strat, x = ~long, y = ~lat, color = I("black"),
              hoverinfo = "none", size = I(0.5), showlegend = FALSE,
              alpha = 0.1) %>%
    layout(xaxis = xax, yaxis = yax, legend = list(x = 1, y = 0.85),
           annotations = list(
             list(text = "n", xref = "paper", yref = "paper",
                  x = 1.09, y = 0.9, showarrow = FALSE),
             list(text = main, xref = "paper", yref = "paper",
                  x = 0.08, y = 1, showarrow = FALSE)
           ),
           margin = m)
  set_plot

}


#' @export
#' @rdname plot_trend
plot_age_strat_fan <- function(sim, surveys = 1:5, years = 1:10,
                               ages = 1:10, select_by = "year",
                               quants = seq(90, 10, by = -10),
                               ...) {

  d <- sim$age_strat_error
  sub_d <- d[survey %in% surveys & year %in% years & age %in% ages, ]
  true <- sub_d[sim == 1, list(year, age, survey, I)]

  ## Calculate a series of quantiles
  ints <- lapply(quants, function(q) {
    sub_d[, list(prob = paste0(q, "%"),
                 lower = quantile(I_hat, prob = (1 - q / 100) / 2),
                 upper = quantile(I_hat, prob = 1 - (1 - q / 100) / 2)),
          by = c("age", "year", "survey")]
  })
  ints <- rbindlist(ints)
  ints$prob <- factor(ints$prob, levels = paste0(quants, "%"))

  ints <- merge(sim$surveys, ints, by = "survey", all.y = TRUE)
  ints$lab <- paste(formatC(ints$set_den, format = "fg"),
                    ints$lengths_cap, ints$ages_cap, sep = "-")
  ints <- merge(ints, true, by = c("age", "year", "survey"))

  shared_ints <- crosstalk::SharedData$new(ints)

  if (select_by == "year") {
    f <- filter_select(select_by, "Year", shared_ints, ~year, multiple = FALSE)
    x <- ~age
    xlab <- "Age"
  } else {
    f <- filter_select(select_by, "Age", shared_ints, ~age, multiple = FALSE)
    x <- ~year
    xlab <- "Year"
  }

  crosstalk::bscols(
    list(
      htmltools::div(style = htmltools::css(height = "10px")), # small margin
      f,
      filter_select("set_den", "Set density", shared_ints, ~set_den, multiple = FALSE),
      filter_select("lengths_cap", "Lengths cap", shared_ints, ~lengths_cap, multiple = FALSE),
      filter_select("ages_cap", "Ages cap", shared_ints, ~ages_cap, multiple = FALSE)
    ),

    plot_ly(data = shared_ints, x = x, height = "100%") %>%
      add_ribbons(ymin = ~lower, ymax = ~upper, line = list(width = 0),
                  color = ~prob, colors = viridis::viridis(nlevels(ints$prob)),
                  showlegend = FALSE) %>%
      add_lines(y = ~I, color = I("black"), name = "True",
                showlegend = FALSE) %>%
      layout(yaxis = list(title = "Abundance index",
                          range = c(min(ints$lower), max(ints$upper))),
             xaxis = list(title = xlab), ...),
    widths = c(3, NA)
  )

}


#' @export
#' @rdname plot_trend
plot_length_strat_fan <- function(sim, surveys = 1:5, years = 1:10,
                                  lengths = 1:50, select_by = "year",
                                  quants = seq(90, 10, by = -10),
                                  ...) {

  d <- sim$length_strat_error
  sub_d <- d[survey %in% surveys & year %in% years & length %in% lengths, ]
  true <- sub_d[sim == 1, list(year, length, survey, I)]

  ## Calculate a series of quantiles
  ints <- lapply(quants, function(q) {
    sub_d[, list(prob = paste0(q, "%"),
                 lower = quantile(I_hat, prob = (1 - q / 100) / 2),
                 upper = quantile(I_hat, prob = 1 - (1 - q / 100) / 2)),
          by = c("length", "year", "survey")]
  })
  ints <- rbindlist(ints)
  ints$prob <- factor(ints$prob, levels = paste0(quants, "%"))

  ints <- merge(sim$surveys, ints, by = "survey", all.y = TRUE)
  ints$lab <- paste(formatC(ints$set_den, format = "fg"),
                    ints$lengths_cap, sep = "-")
  ints <- merge(ints, true, by = c("length", "year", "survey"))

  shared_ints <- crosstalk::SharedData$new(ints)

  if (select_by == "year") {
    f <- filter_select(select_by, "Year", shared_ints, ~year, multiple = FALSE)
    x <- ~length
    xlab <- "Length"
  } else {
    f <- filter_select(select_by, "Length", shared_ints, ~length, multiple = FALSE)
    x <- ~year
    xlab <- "Year"
  }

  crosstalk::bscols(
    list(
      htmltools::div(style = htmltools::css(height = "10px")), # small margin
      f,
      filter_select("set_den", "Set density", shared_ints, ~set_den, multiple = FALSE),
      filter_select("lengths_cap", "Lengths cap", shared_ints, ~lengths_cap, multiple = FALSE)
    ),

    plot_ly(data = shared_ints, x = x, height = "100%") %>%
      add_ribbons(ymin = ~lower, ymax = ~upper, line = list(width = 0),
                  color = ~prob, colors = viridis::viridis(nlevels(ints$prob)),
                  showlegend = FALSE) %>%
      add_lines(y = ~I, color = I("black"), name = "True",
                showlegend = FALSE) %>%
      layout(yaxis = list(title = "Abundance index",
                          range = c(min(ints$lower), max(ints$upper))),
             xaxis = list(title = xlab)),
    widths = c(3, NA)
  )

}



#' @export
#' @rdname plot_trend
plot_total_strat_fan <- function(sim, surveys = 1:5,
                                 quants = seq(90, 10, by = -10),
                                 ...) {

  d <- sim$total_strat_error
  sub_d <- d[survey %in% surveys, ]
  true <- sub_d[sim == 1, list(year, survey, I)]

  ## Calculate a series of quantiles
  ints <- lapply(quants, function(q) {
    sub_d[, list(prob = paste0(q, "%"),
                 lower = quantile(I_hat, prob = (1 - q / 100) / 2),
                 upper = quantile(I_hat, prob = 1 - (1 - q / 100) / 2)),
          by = c("year", "survey")]
  })
  ints <- rbindlist(ints)
  ints$prob <- factor(ints$prob, levels = paste0(quants, "%"))

  ints <- merge(sim$surveys, ints, by = "survey", all.y = TRUE)
  ints$lab <- formatC(ints$set_den, format = "fg")
  ints <- merge(ints, true, by = c("year", "survey"))

  shared_ints <- crosstalk::SharedData$new(ints)

  crosstalk::bscols(
    list(
      htmltools::div(style = htmltools::css(height = "10px")), # small margin
      filter_select("set_den", "Set density", shared_ints, ~set_den, multiple = FALSE)
    ),
    plot_ly(data = shared_ints, x = ~year, height = "100%") %>%
      add_ribbons(ymin = ~lower, ymax = ~upper, line = list(width = 0),
                  color = ~prob, colors = viridis::viridis(nlevels(ints$prob)),
                  showlegend = FALSE) %>%
      add_lines(y = ~I, color = I("black"), name = "True",
                showlegend = FALSE) %>%
      layout(yaxis = list(title = "Abundance index",
                          range = c(min(ints$lower), max(ints$upper))),
             xaxis = list(title = "Year"), ...),
    widths = c(3, NA)
  )

}


