
#' Simulate survey sets
#'
#' @param pop              Simulation object from \code{\link{sim_distribution}}
#' @param min_sets         Minimum number of sets per strat
#' @param set_den          Set density
#' @param resample_units   Sllow resampling of sampling units (grid cells)?
#'
#' @export
#'
#' @import data.table
#'
#'

sim_sets <- function(pop = sim_distribution(), min_sets = 1, set_den = 5 / 1000, resample_units = FALSE) {

  ## Strat area and sampling effort
  cells <- data.table(rasterToPoints(pop$grid))
  strat_det <- cells[, list(strat_cells = .N), by = "strat"]
  strat_det$strat_area <- strat_det$strat_cells * prod(res(pop$grid))
  strat_det$strat_sets <- round(strat_det$strat_area * set_den) # set allocation
  strat_det$strat_sets[strat_det$strat_sets < min_sets] <- min_sets
  cells <- merge(cells, strat_det, by = c("strat"))

  ## Replicate cells data.table for each year in the simulation
  i <- rep(seq(nrow(cells)), times = length(pop$years))
  y <- rep(pop$years, each = nrow(cells))
  cells <- cells[i, ]
  cells$year <- y

  ## Simulate sets
  sets <- cells[, .SD[sample(.N, strat_sets, replace = resample_units)], by = c("year", "strat")]
  sets$set <- seq(nrow(sets))
  sets

}



