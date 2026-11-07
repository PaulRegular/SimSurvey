## Derive the area of each intended stratum from eligible grid cells.
.survey_frame_strata <- function(sim, frame = sim$survey_frame) {
  if (is.null(frame)) return(NULL)
  grid <- as.data.table(sim$grid_xy)[, c("cell", "division", "strat")]
  eligible <- merge(frame, grid, by = "cell", sort = FALSE)
  cell_area <- prod(stars::st_res(sim$grid))
  eligible[, list(eligible_cells = .N, strat_area = .N * cell_area),
           by = c("sim", "year", "division", "strat")]
}

.validate_survey_frame <- function(sim, sets, frame) {
  if (is.null(sets$set) || anyNA(sets$set) || anyDuplicated(sets$set)) {
    stop("When supplying 'custom_sets', please make sure each set has a unique number.",
         call. = FALSE)
  }
  keys <- c("sim", "year", "cell")
  if (!all(keys %in% names(sets)) ||
      !all(vapply(sets[, keys, with = FALSE], is.numeric, logical(1))) ||
      anyNA(sets[, keys, with = FALSE])) {
    stop("custom_sets must contain nonmissing numeric sim, year, and cell keys.",
         call. = FALSE)
  }
  grid <- sim$grid_xy
  if (any(!sets$sim %in% sim$survey_replicates)) {
    stop("custom_sets must use the survey replicate identifiers.", call. = FALSE)
  }
  if (any(!sets$year %in% sim$years) || any(!sets$cell %in% grid$cell)) {
    stop("custom_sets must use the simulation's year and cell identifiers.", call. = FALSE)
  }
  if (is.null(frame)) return(NULL)

  if (!is.data.frame(frame) || !all(keys %in% names(frame))) {
    stop("survey_frame must contain sim, year, and cell keys.", call. = FALSE)
  }
  frame <- data.table::copy(as.data.table(frame)[, keys, with = FALSE])
  if (!all(vapply(frame, is.numeric, logical(1))) || anyNA(frame) || anyDuplicated(frame)) {
    stop("survey_frame must contain unique, nonmissing numeric sim, year, and cell keys.",
         call. = FALSE)
  }
  if (any(!frame$sim %in% sim$survey_replicates) || any(!frame$year %in% sim$years) ||
      any(!frame$cell %in% grid$cell)) {
    stop("survey_frame keys must use the survey's replicate, year, and cell identifiers.",
         call. = FALSE)
  }
  if (nrow(sets[!frame, on = keys])) {
    stop("All custom_sets must lie within survey_frame.", call. = FALSE)
  }
  cell_rows <- match(sets$cell, grid$cell)
  if (!all(c("division", "strat", "strat_area") %in% names(sets)) ||
      !isTRUE(all(sets$division == grid$division[cell_rows])) ||
      !isTRUE(all(sets$strat == grid$strat[cell_rows]))) {
    stop("custom_sets division and stratum labels must match the survey_frame grid cells.",
         call. = FALSE)
  }
  strata <- .survey_frame_strata(sim, frame)
  area <- strata[sets, on = c("sim", "year", "division", "strat")][["strat_area"]]
  if (!is.numeric(sets$strat_area) ||
      !isTRUE(all.equal(sets$strat_area, area, check.attributes = FALSE))) {
    stop("custom_sets strat_area must match the eligible cell area in survey_frame.",
         call. = FALSE)
  }
  setorderv(frame, keys)
  frame
}
