.plot_fishery_config <- function(sim, fishery) {
  if (!is.character(fishery) || length(fishery) != 1 || is.na(fishery) ||
      !fishery %in% names(sim$fishery$configurations)) {
    stop("fishery must name one simulated fishery.", call. = FALSE)
  }
  sim$fishery$configurations[[fishery]]
}

.fishery_percent <- function(count) {
  if (sum(count) == 0) rep(0, length(count)) else 100 * count / sum(count)
}

.plot_fishery_composition <- function(catch, sample, trait, groups, showlegend) {
  counts <- vapply(groups, function(value) sum(catch$n[catch[[trait]] == value]), numeric(1))
  measured <- tabulate(match(sample[[trait]], groups), nbins = length(groups))
  label <- if (trait == "length") "Length group (cm)" else "Age"
  plotly::plot_ly(x = groups, y = .fishery_percent(counts), type = "bar",
                  name = "Catch truth", legendgroup = "truth", showlegend = showlegend,
                  marker = list(color = "#440154"), hovertemplate = "%{x}: %{y:.1f}%<extra>Catch truth</extra>") |>
    plotly::add_bars(y = .fishery_percent(measured), name = "Biological sample", legendgroup = "sample",
                     showlegend = showlegend, marker = list(color = "#21908C"),
                     hovertemplate = "%{x}: %{y:.1f}%<extra>Biological sample</extra>") |>
    plotly::layout(barmode = "group", xaxis = list(title = label),
                    yaxis = list(title = "Percent", range = c(0, max(1, .fishery_percent(counts),
                                                                  .fishery_percent(measured)))))
}

#' Plot fishery locations and biological compositions
#'
#' A spatial map and true catch versus unweighted measured-length/aged-fish
#' compositions. Counts remain aggregate even for large commercial catches.
#' @param sim A simulation returned by [sim_fishery()] or [sim_fisheries()].
#' @param fishery One name in `sim$fishery$configurations`.
#' @param which_year One simulated year.
#' @param which_sim One replicate, or `NULL` to pool replicates.
#' @param map_by Map `"operations"`, catch numbers (`"n"`), kilograms (`"kg"`),
#' or nominal `"effort"`. Zero-catch operations are included.
#' @details Compositions are percentages within each catch/sample, not weighted
#' estimates of stock composition. Empty catches/samples have zero bars. Effort
#' units are inherited from the selected fishery and are never pooled across gears.
#' @return A Plotly widget. No input tables are modified.
#' @examples
#' # plot_fishery(result, "trawl", which_sim = NULL)
#' @export
plot_fishery <- function(sim, fishery, which_year = sim$years[1],
                          which_sim = 1, map_by = "operations") {
  config <- .plot_fishery_config(sim, fishery)
  if (length(which_year) != 1 || is.na(which_year) || !which_year %in% sim$years) {
    stop("which_year must select one simulated year.")
  }
  replicates <- unique(sim$fishery$annual$sim)
  if (!is.null(which_sim)) {
    if (length(which_sim) != 1 || is.na(which_sim) || !which_sim %in% replicates) {
      stop("which_sim must select one replicate, or NULL for all replicates.")
    }
    replicates <- which_sim
  }
  map_by <- match.arg(map_by, c("operations", "n", "kg", "effort"))
  n <- kg <- effort <- NULL
  operations <- sim$fishery$setdet
  keep <- operations$fishery == fishery & operations$year == which_year & operations$sim %in% replicates
  operations <- operations[keep, ]
  cells <- data.table::copy(sim$grid_xy[, c("cell", "x", "y")])
  totals <- if (nrow(operations))
    operations[, list(operations = .N, n = sum(n), kg = sum(kg), effort = sum(effort)), by = "cell"] else
      data.table(cell = numeric(), operations = numeric(), n = numeric(), kg = numeric(), effort = numeric())
  value <- totals[[map_by]][match(cells$cell, totals$cell)]
  value[!cells$cell %in% totals$cell] <- 0
  if (anyNA(value)) stop("Catch weights are unavailable for this older fishery result.")
  units <- switch(map_by, operations = "Operations", n = "Caught fish", kg = "Catch (kg)",
                   effort = paste0("Effort (", config$gear_unit, " * ", config$time_unit, ")"))
  z <- tapply(value, list(cells$y, cells$x), sum)
  ## xtabs/tapply order is independent of the stored grid row order.
  hover <- tapply(paste0("Cell: ", cells$cell, "<br>", units, ": ", signif(value, 5)),
                    list(cells$y, cells$x), paste, collapse = "")
  map <- plotly::plot_ly(x = sort(unique(cells$x)), y = sort(unique(cells$y)),
                         z = unname(z), text = unname(hover), type = "heatmap", hoverinfo = "text",
                         colors = viridis::viridis(32), zmin = 0, zmax = max(1, value),
                         colorbar = list(title = units, x = 0.53, len = 0.75, thickness = 12), showscale = TRUE) |>
    plotly::layout(xaxis = list(title = "x (km)"), yaxis = list(title = "y (km)", scaleanchor = "x"))
  catch_length <- sim$fishery$catch_at_length
  keep <- catch_length$operation %in% operations$operation
  catch_length <- catch_length[keep, ]
  catch_age <- sim$fishery$catch_at_age
  keep <- catch_age$operation %in% operations$operation
  catch_age <- catch_age[keep, ]
  sample <- sim$fishery$samp
  keep <- sample$operation %in% operations$operation
  sample <- sample[keep, ]
  lengths <- .plot_fishery_composition(catch_length, sample, "length", sim$lengths, TRUE)
  aged_rows <- which(sample$aged)
  ages <- .plot_fishery_composition(catch_age, sample[aged_rows, ], "age", sim$ages, FALSE)
  plotly::subplot(map, plotly::subplot(lengths, ages, nrows = 2, titleX = TRUE, titleY = TRUE,
                                      margin = 0.14), widths = c(0.6, 0.4),
                    titleX = TRUE, titleY = TRUE, margin = 0.12) |>
    plotly::layout(title = paste(fishery, which_year,
                                  if (is.null(which_sim)) "(all replicates)" else paste("replicate", which_sim)),
                    margin = list(l = 70, r = 45, t = 55, b = 130),
                    legend = list(orientation = "h", x = 0, y = -0.22, yanchor = "top"))
}

#' Plot annual fishery catches, effort, or CPUE
#'
#' Replicate values and median annual trends for one fishery, preserving its
#' effort units. Number-catch plots also show the integer catch target.
#' @inheritParams plot_fishery
#' @param metric Annual result column: `"n"`, `"reported_n"`, `"kg"`,
#' `"reported_kg"`, `"effort"`, `"shortfall"`, `"cpue"`, `"reported_cpue"`,
#' `"biomass_cpue"`, or `"reported_biomass_cpue"`.
#' @param ... Additional arguments passed to [plotly::plot_ly()].
#' @return A Plotly widget; medians summarize conditional replicates, not
#' uncertainty across independently regenerated population trajectories.
#' @examples
#' # plot_fishery_trend(result, "trawl", metric = "n")
#' @export
plot_fishery_trend <- function(sim, fishery, metric = "cpue", ...) {
  config <- .plot_fishery_config(sim, fishery)
  metric <- match.arg(metric, c("n", "reported_n", "kg", "reported_kg", "effort", "shortfall",
                                "cpue", "reported_cpue", "biomass_cpue", "reported_biomass_cpue"))
  annual <- sim$fishery$annual
  keep <- annual$fishery == fishery
  annual <- as.data.frame(annual[keep, ])
  annual$value <- annual[[metric]]
  effort_units <- paste(config$gear_unit, config$time_unit, sep = " * ")
  units <- switch(metric, n = "Fish", reported_n = "Reported fish", kg = "kg", reported_kg = "Reported kg",
                    effort = effort_units, shortfall = "Uncaught target (fish)",
                    cpue = paste("Fish /", effort_units), reported_cpue = paste("Reported fish /", effort_units),
                    biomass_cpue = paste("kg /", effort_units), reported_biomass_cpue = paste("Reported kg /", effort_units))
  annual$hover <- paste0("Year: ", annual$year, "<br>Replicate: ", annual$sim,
                          "<br>", units, ": ", signif(annual$value, 5),
                          "<br>Status: ", annual$status, "<br>Shortfall: ", annual$shortfall)
  median_value <- vapply(sim$years, function(year) {
    stats::median(annual$value[annual$year == year], na.rm = TRUE)
  }, numeric(1))
  points <- annual[is.finite(annual$value), , drop = FALSE]
  p <- plotly::plot_ly(data = points, x = ~year, y = ~value, text = ~hover,
                       hoverinfo = "text", type = "scatter", mode = "markers", name = "Replicates",
                       marker = list(color = "#21908C", opacity = 0.45), ...)
  finite_medians <- is.finite(median_value)
  if (any(finite_medians)) {
    p <- p |> plotly::add_lines(x = sim$years[finite_medians], y = median_value[finite_medians],
                                inherit = FALSE, name = "Median", line = list(color = "#440154"))
  } else {
    p <- p |> plotly::add_annotations(x = 0.5, y = 0.5, xref = "paper", yref = "paper",
                                       text = "No finite values (check effort or available weights)",
                                       showarrow = FALSE)
  }
  if (metric == "n") {
    target <- vapply(sim$years, function(year) unique(annual$target[annual$year == year])[1], numeric(1))
    p <- p |> plotly::add_lines(x = sim$years, y = target, inherit = FALSE, name = "Catch target",
                                line = list(color = "#D55E00", dash = "dash"))
  }
  p |> plotly::layout(title = fishery, xaxis = list(title = "Year"), yaxis = list(title = units))
}
