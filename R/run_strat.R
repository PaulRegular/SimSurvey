
#' Calculate common error statistics
#'
#' @param error A numeric vector of errors.
#'
#' @return A named vector of error statistics including:
#'
#' - **ME**: Mean error
#' - **MAE**: Mean absolute error
#' - **MSE**: Mean squared error
#' - **RMSE**: Root mean squared error
#'
#' @export

error_stats <- function(error) {
  c(ME = mean(error),
    MAE = mean(abs(error)),
    MSE = mean(error ^ 2),
    RMSE = sqrt(mean(error ^ 2)))
}


#' Prepare simulated data for stratified analysis
#'
#' Generate set details (`setdet`), length-frequency (`lf`), and age-frequency (`af`)
#' data for stratified analysis.
#'
#' @param sim Simulation object returned by [`sim_survey()`].
#' @param length_group Size of the length frequency bins.
#' @param alk_scale Spatial scale at which to construct and apply age-length keys:
#' `"division"`, `"strat"`, or `"set"`.
#'
#' @return A list containing:
#' - `setdet`: Set details
#' - `lf`: Length-frequency data
#' - `af`: Age-frequency data
#'
#' @details Length frequencies expand measured fish to the observed set catch.
#' A positive catch without measured lengths has missing composition. Length
#' frequencies without a matching age-length key are assigned to `age = NA`.
#'
#' @export

strat_data <- function(sim, length_group = 3, alk_scale = "division") {

  n <- n_measured <- n_aged <- age <- age_prop <- n_caught <- NULL
  setdet <- sim$setdet
  samp <- sim$samp
  alk_scale <- match.arg(alk_scale, c("division", "strat", "set"))
  set_cols <- unique(c("set", "sim", "year", alk_scale))
  set_info <- setdet[, c(set_cols, "n"), with = FALSE]
  setnames(set_info, "n", "n_caught")
  measured <- samp[which(samp$measured & is.finite(samp$length)),
                   c("set", "length", "age", "aged"), with = FALSE]
  measured <- merge(set_info, measured, by = "set")
  measured$length <- group_lengths(measured$length, length_group)
  aged <- measured[which(measured$aged & !is.na(measured$age))]

  if (alk_scale %in% c("strat", "set")) {
    check <- measured[, list(n_aged = sum(aged & !is.na(age))),
                      by = c("sim", "year", alk_scale)]
    if (any(check$n_aged == 0)) {
      warning("Age-length keys cannot be constructed at the ", alk_scale,
              " level because ages were not sampled at every ", alk_scale,
              ". This discrepancy may introduce bias. Consider running sim_survey with ",
              "age_space_group = '", alk_scale, "'.")
    }
  }

  ## Fill absent bins only for sets with an observed zero or measured lengths.
  measured_sets <- unique(measured$set)
  missing_lengths <- setdet$set[setdet$n > 0 & !setdet$set %in% measured_sets]
  length_bins <- unique(c(.strat_length_map(sim, length_group), measured$length))
  lf <- measured[, list(n = .N), by = c("set", "length", "n_caught")]
  lf[, n_measured := sum(n), by = "set"]
  lf[, n := n * n_caught / n_measured]
  lf <- merge(CJ(set = setdet$set, length = length_bins, unique = TRUE),
              lf[, c("set", "length", "n"), with = FALSE],
              by = c("set", "length"), all.x = TRUE)
  lf$n[is.na(lf$n)] <- 0
  lf$n[lf$set %in% missing_lengths] <- NA_real_

  alk_keys <- c("sim", "year", "length", alk_scale)
  alk <- aged[, list(n_aged = .N), by = c(alk_keys, "age")]
  alk[, age_prop := n_aged / sum(n_aged), by = alk_keys]
  alf <- merge(merge(set_info, lf, by = "set"), alk,
               by = alk_keys, all.x = TRUE, allow.cartesian = TRUE)
  alf[, n := ifelse(is.na(age), n, n * age_prop)]
  af <- alf[, list(n = sum(n)), by = c("set", "age")]
  ages <- unique(c(NA_real_, sim$ages, aged$age))
  af <- merge(CJ(set = setdet$set, age = ages, unique = TRUE),
              af, by = c("set", "age"), all.x = TRUE)
  af$n[is.na(af$n)] <- 0
  af$n[af$set %in% missing_lengths] <- NA_real_
  list(setdet = setdet, lf = lf, af = af)

}



#' Calculate stratified means, variances, and confidence intervals across groups
#'
#' This function is primarily designed for use within [`run_strat()`]. It first calculates
#' statistics at the stratum level and then computes broader summaries like total abundance.
#'
#' @param data A `data.table` with all grouping variables in stacked format. Must include
#' `strat_area` and `tow_area` for scaling values.
#' @param metric Name of the variable in `data.table` to summarize (e.g., `"number"`, `"mass"`).
#' @param strat_groups Grouping variables for fine-scale stratum-level means.
#' Must include `"strat"` and `"strat_area"`.
#' Example: `c("year", "species", "shiptrip", "NAFOdiv", "strat", "strat_area", "age")`
#' @param survey_groups Grouping variables for large-scale summary calculations
#' Example: `c("year", "species")`
#' @param confidence Confidence limit percentage (e.g., 95 for 95% CI).
#'
#' @return A `data.table` containing stratified estimates of abundance.
#'
#' @export

strat_means <- function(data = NULL, metric = NULL, strat_groups = NULL,
                        survey_groups = NULL, confidence = 95) {

  Nh <- strat_area <- tow_area <- Wh <- total <- sumYh <- nh <- gh <- meanYh <- varYh <-
    meanYst_lcl <- meanYst <- varYst <- df <- meanYst_ucl <- sumYst <- N <- sumYst_lcl <-
    sumYst_ucl <- NULL

  ## set-up
  lc <- (100 - confidence) / 200
  uc <- (100 - confidence) / 200 + (confidence / 100)
  d <- copy(data) # make a copy of the provided data.table
  d <- d[,  c(strat_groups,  metric),  with = FALSE]
  setnames(d, names(d), c(strat_groups, "metric"))
  setkeyv(d,  strat_groups)

  ## strat.tab includes strat-level means,  variances,  etc.
  strat_tab <- d[, list(sumYh = sum(metric), meanYh = mean(metric),
                        varYh = stats::var(metric), nh = .N), by = strat_groups]
  strat_tab[, Nh := strat_area / tow_area] # number of sample units in each strat
  strat_tab[, Wh := Nh / sum(Nh), by = survey_groups]
  strat_tab[, total := Nh * sumYh / nh]
  strat_tab[, gh := Nh * (Nh - nh) / nh]

  ## survey.tab includes large-scale means,  such as mean weight per tow,  and totals,  such as total biomass + associated confidence intervals
  survey_tab <- strat_tab[, list(n = sum(nh), N = sum(Nh), meanYst = sum(Wh * meanYh),
                                 varYst = (1 / ((sum(Nh)) ^ 2)) * sum(gh * varYh),
                                 df = ((sum(gh * varYh)) ^ 2) / (sum((gh ^ 2 * varYh ^ 2) / (nh - 1)))),
                          by = survey_groups]
  survey_tab[, meanYst_lcl := (meanYst - (sqrt(varYst)) * abs(stats::qt(lc,  df)))]
  survey_tab[, meanYst_ucl := (meanYst + (sqrt(varYst)) * abs(stats::qt(lc,  df)))]
  survey_tab[, sumYst := N * meanYst]
  survey_tab[, sumYst_lcl := (sumYst - abs(stats::qt(lc, df)) * N * sqrt(varYst))]
  survey_tab[, sumYst_ucl := (sumYst + abs(stats::qt(lc, df)) * N * sqrt(varYst))]
  for (column in names(survey_tab)) {
    nan_rows <- which(is.nan(survey_tab[[column]]))
    if (length(nan_rows)) data.table::set(survey_tab, i = nan_rows, j = column, value = NA)
  }

  ## Rename cols
  survey_tab <- survey_tab[, c(survey_groups,  "n", "N", "df", "varYst",
                               "meanYst",  "meanYst_lcl", "meanYst_ucl",
                               "sumYst", "sumYst_lcl", "sumYst_ucl"),  with = FALSE]
  survey_tab$varYst <- sqrt(survey_tab$varYst) # convert to sd
  setnames(survey_tab, names(survey_tab), c(survey_groups, "sets", "sampling_units", "df", "sd",
                                            "mean", "mean_lcl", "mean_ucl",
                                            "total", "total_lcl", "total_ucl"))

  ## return results (here I only return the survey_tab to minimize details and object size)
  survey_tab

}


#' Run stratified analysis on simulated data
#'
#' @param sim Simulation object from [`sim_survey()`] or [`sim_surveys()`].
#' @param length_group Size of the length frequency bins used for both abundance-at-length
#' calculations and age-length-key construction. By default, this is inherited from the
#' value defined in [`sim_abundance()`] via the closure supplied to `sim_length` (`"inherit"`).
#' A numeric width can also be supplied. Length comparisons require matching
#' bins or an exact aggregation of the original bins.
#' @param alk_scale Spatial scale at which to construct and apply age-length keys:
#' `"division"`, `"strat"`, or `"set"`.
#' @param strat_data_fun Function used to prepare data for stratified analysis (e.g., [`strat_data()`]).
#' @param strat_means_fun Function used to calculate stratified means (e.g., [`strat_means()`]).
#'
#' @details The `strat_data_fun` and `strat_means_fun` arguments allow you to use custom
#' [`strat_data()`] and [`strat_means()`] functions.
#' Custom estimators must return the requested survey grouping columns and a
#' `total` column. The existing estimator assumes representative sets within
#' each stratum and a common tow area within a stratum. Declaring the frame
#' supplies the intended coverage for comparison.
#' These estimates target abundance available to the survey, meaning abundance within
#' the intended survey footprint after applying `q` at age or length.
#' For a batch, both hooks run independently for each named program.
#' A known frame with unsampled strata has missing full-frame estimates.
#' An unassigned age row retains catch without a matching age-length key.
#' Assigned ages can have partial totals in this case.
#'
#' @return Adds stratified analysis results to the `sim` list:
#' - `"total_strat"`: Results for the total abundance available to the survey
#' - `"length_strat"`: Results aggregated by length group
#' - `"age_strat"`: Results aggregated by age
#' - `"survey_coverage"`: Eligible and visited cells and set counts by stratum
#'
#' Results have an `estimate_status` column identifying support for each estimate.
#' Statuses include `"ok"`, `"no_sets"`, `"unsampled_strata"`,
#' `"missing_lengths"`, `"missing_estimate"`, `"partial_age"`, and
#' `"unassigned_age"`.
#' Batch tables add a `survey` key; [`select_survey()`] restores one program's view.
#'
#' @examples
#' \donttest{
#' sim <- sim_abundance(ages = 1:5, years = 1:5,
#'                      R = sim_R(log_mean = log(1e+7)),
#'                      growth = sim_vonB(length_group = 1)) |>
#'   sim_distribution(grid = make_grid(res = c(20, 20)),
#'                    ays_covar = sim_ays_covar(sd = 1)) |>
#'   sim_survey(n_sims = 1, q = sim_logistic(k = 2, x50 = 3)) |>
#'   run_strat()
#' }
#'
#' @export

run_strat <- function(sim,
                      length_group = "inherit",
                      alk_scale = "division",
                      strat_data_fun = strat_data,
                      strat_means_fun = strat_means) {

  native_group <- .strat_native_length_group(sim)
  if (identical(length_group, "inherit")) length_group <- native_group
  if (!is.numeric(length_group) || length(length_group) != 1L ||
      !is.finite(length_group) || (length_group < 1 && length_group != 0.5)) {
    stop("Supply a valid length_group when it cannot be inherited from growth.", call. = FALSE)
  }
  if (is.null(.strat_length_map(sim, length_group))) {
    warning("Use a length_group that can match the length group defined by growth. ",
            "Matching length truth is unavailable with this setting.", call. = FALSE)
  }
  sim <- .clear_strat_analysis(sim)
  if (is.null(sim$setdet) && length(sim$survey_programs)) {
    fields <- c(.strat_estimates, "survey_coverage")
    results <- lapply(names(sim$survey_programs), function(survey) {
      view <- .survey_view(sim, survey, restore_spatial = FALSE)
      .run_strat_one(view, length_group, alk_scale, strat_data_fun, strat_means_fun)[fields]
    })
    names(results) <- names(sim$survey_programs)
    for (field in fields) {
      sim[[field]] <- .stack_survey_tables(results, field)
    }
    sim$strat_length_group <- length_group
    return(sim)
  }
  .run_strat_one(sim, length_group, alk_scale, strat_data_fun, strat_means_fun)
}

.strat_estimates <- c("total_strat", "length_strat", "age_strat")
.strat_errors <- paste0(.strat_estimates, "_error")
.strat_error_stats <- paste0(.strat_errors, "_stats")

.clear_strat_analysis <- function(sim) {
  fields <- c(.strat_estimates, .strat_errors, .strat_error_stats,
              "survey_coverage", "strat_length_group")
  sim[intersect(fields, names(sim))] <- NULL
  sim
}

.strat_native_length_group <- function(sim) {
  if (is.function(sim$sim_length) && !is.null(environment(sim$sim_length))) {
    get0("length_group", envir = environment(sim$sim_length), inherits = FALSE)
  }
}

.strat_length_map <- function(sim, width) {
  native <- .strat_native_length_group(sim)
  if (is.null(native) || !length(sim$lengths)) return(NULL)
  if (is.null(width) || width == native) return(sim$lengths)
  ratio <- width / native
  if (ratio >= 1 && isTRUE(all.equal(ratio, round(ratio)))) group_lengths(sim$lengths, width)
}

.strat_roster <- function(sim) {
  reps <- sim$survey_replicates
  if (is.null(reps)) reps <- sim$replicates
  if (is.null(reps)) reps <- unique(sim$setdet$sim)
  if (!length(reps)) reps <- 1L
  CJ(sim = reps, year = sim$years, unique = TRUE)
}

.strat_coverage <- function(sim, sets) {
  cell <- strat <- sets_sampled <- visited_cells <- eligible_cells <- NULL
  visits <- sets[, list(sets_sampled = .N, visited_cells = uniqueN(cell)),
                 by = c("sim", "year", "strat")]
  if (is.null(sim$survey_frame)) {
    visits[, eligible_cells := NA_integer_]
    return(setkeyv(visits, c("sim", "year", "strat")))
  }
  frame_strata <- .survey_frame_strata(sim)
  coverage <- merge(frame_strata[, c("sim", "year", "strat", "eligible_cells"), with = FALSE],
                    visits, by = c("sim", "year", "strat"), all.x = TRUE)
  coverage[is.na(sets_sampled), `:=`(sets_sampled = 0L, visited_cells = 0L)]
  setkeyv(coverage, c("sim", "year", "strat"))
}

.stack_survey_tables <- function(results, field) {
  tables <- lapply(results, function(result) {
    table <- result[[field]]
    if (is.data.table(table)) table else as.data.table(as.list(table))
  })
  combined <- rbindlist(tables, idcol = "survey", use.names = TRUE)
  setkeyv(combined, intersect(c("survey", "sim", "year", "length", "age", "strat"), names(combined)))
}

.restore_survey_analysis <- function(view, original, survey) {
  if (!is.null(original$setdet) || !"survey" %in% names(original$total_strat)) return(view)
  for (field in c(.strat_estimates, .strat_errors, .strat_error_stats, "survey_coverage")) {
    table <- original[[field]]
    if (is.null(table)) next
    rows <- which(table$survey == survey)
    table <- table[rows, setdiff(names(table), "survey"), with = FALSE]
    view[[field]] <- if (field %in% .strat_error_stats) unlist(table, use.names = TRUE) else
      setkeyv(table, setdiff(key(original[[field]]), "survey"))
  }
  view$strat_length_group <- original$strat_length_group
  view
}

.run_strat_one <- function(sim, length_group, alk_scale, strat_data_fun, strat_means_fun) {
  if (is.null(sim$setdet) || is.null(sim$samp)) {
    stop("Run sim_survey() or sim_surveys() before run_strat().", call. = FALSE)
  }
  n <- total <- age <- estimate_status <- sets_sampled <- sets <- NULL
  data <- strat_data_fun(sim, length_group = length_group, alk_scale = alk_scale)
  keys <- c("sim", "year")
  strata <- c(keys, "division", "strat", "strat_area", "tow_area")
  data$setdet <- data$setdet[, c("set", strata, "n"), with = FALSE]
  set_info <- data$setdet[, c("set", strata), with = FALSE]
  data$lf <- merge(set_info, data$lf, by = "set", all.y = TRUE)
  data$af <- merge(set_info, data$af, by = "set", all.y = TRUE)

  roster <- .strat_roster(sim)
  coverage <- .strat_coverage(sim, sim$setdet)
  sampled <- unique(sim$setdet[, keys, with = FALSE])
  gaps <- unique(coverage[sets_sampled == 0L, keys, with = FALSE])
  roster$estimate_status <- "no_sets"
  roster[sampled, on = keys, estimate_status := "ok"]
  roster[gaps, on = keys, estimate_status := "unsampled_strata"]
  missing_lengths <- unique(data$lf[is.na(n), keys, with = FALSE])

  for (i in seq_along(.strat_estimates)) {
    trait <- list(NULL, "length", "age")[[i]]
    table <- strat_means_fun(data = data[[c("setdet", "lf", "af")[i]]],
                            metric = "n", strat_groups = c(strata, trait),
                            survey_groups = c(keys, trait))
    template <- roster
    if (!is.null(trait)) {
      bins <- if (trait == "age") c(NA_real_, sim$ages) else .strat_length_map(sim, length_group)
      bins <- unique(as.numeric(c(bins, table[[trait]])))
      template <- roster[rep(seq_len(nrow(roster)), each = length(bins))]
      template[[trait]] <- rep(bins, times = nrow(roster))
    }
    table <- merge(template, table, by = c(keys, trait), all.x = TRUE)
    if ("sets" %in% names(table)) {
      table[!sampled, on = keys, sets := 0L]
    }
    if (i > 1L) {
      table[missing_lengths, on = keys, estimate_status :=
              ifelse(estimate_status == "ok", "missing_lengths", estimate_status)]
    }
    ## A subtotal across sampled strata cannot estimate the whole intended frame.
    unsupported <- which(table$estimate_status != "ok")
    values <- setdiff(names(table), c(keys, trait, "estimate_status", "sets"))
    for (column in values) set(table, i = unsupported, j = column, value = NA)
    table[estimate_status == "ok" & is.na(total), estimate_status := "missing_estimate"]
    if (identical(trait, "age")) {
      partial <- table[is.na(age) & total > 0 & estimate_status == "ok", keys, with = FALSE]
      table[partial, on = keys, estimate_status :=
              ifelse(is.na(age), "unassigned_age", "partial_age")]
      no_assignment <- which(table$estimate_status == "partial_age" & table$total == 0)
      for (column in values) set(table, i = no_assignment, j = column, value = NA)
      table[is.na(age) & estimate_status == "ok", estimate_status := "unassigned_age"]
    }
    sim[[.strat_estimates[i]]] <- setkeyv(table, c(keys, trait))
  }
  sim$survey_coverage <- coverage
  sim$strat_length_group <- length_group
  sim
}



#' Calculate error of stratified estimates
#'
#' @param sim Object returned by [`run_strat()`], which includes the simulated population,
#' survey results, and stratified analysis outputs.
#'
#' @return Adds error details and summary statistics to the `sim` list, ending with
#' `"*_strat_error"` and `"*_strat_error_stats"`. Error statistics include:
#'
#' - **MAE**: Mean absolute error
#' - **MSE**: Mean squared error
#' - **RMSE**: Root mean squared error
#'
#' @details In the error tables, `I` is abundance available to the survey,
#' calculated within the intended survey footprint after applying the program's
#' `q` at age or length. `N_frame` and `N_stock` give abundance before applying
#' `q`, within the footprint and across the full grid, respectively.
#' Matching truth is missing for unknown frames, incompatible length bins and
#' unassigned ages. Only estimates and truth with status
#' `"ok"` contribute to error summaries. Batch outputs include a `survey` key.
#'
#' @examples
#' \donttest{
#' sim <- sim_abundance(ages = 1:5, years = 1:5,
#'                      R = sim_R(log_mean = log(1e+7)),
#'                      growth = sim_vonB(length_group = 1)) |>
#'   sim_distribution(grid = make_grid(res = c(20, 20)),
#'                    ays_covar = sim_ays_covar(sd = 1)) |>
#'   sim_survey(n_sims = 1, q = sim_logistic(k = 2, x50 = 3)) |>
#'   run_strat() |>
#'   strat_error()
#' }
#'
#' @export

strat_error <- function(sim) {
  if (any(vapply(.strat_estimates, function(field) is.null(sim[[field]]), logical(1)))) {
    stop("Run run_strat() before strat_error().", call. = FALSE)
  }
  if (is.null(sim$setdet) && length(sim$survey_programs)) {
    fields <- c(.strat_errors, .strat_error_stats)
    results <- lapply(names(sim$survey_programs), function(survey) {
      view <- .survey_view(sim, survey, restore_spatial = FALSE)
      .strat_error_one(.restore_survey_analysis(view, sim, survey))[fields]
    })
    names(results) <- names(sim$survey_programs)
    for (field in fields) {
      sim[[field]] <- .stack_survey_tables(results, field)
    }
    return(sim)
  }
  .strat_error_one(sim)
}

.comparable_error_stats <- function(errors) {
  errors <- errors[is.finite(errors)]
  error_stats(if (length(errors)) errors else NA_real_)
}

.strat_frame_truth <- function(sim, roster) {
  N <- N_frame <- NULL
  n_ages <- length(sim$ages)
  n_lengths <- length(sim$lengths)
  n_groups <- nrow(roster)
  years <- match(roster$year, sim$years)
  empty_age <- matrix(NA_real_, nrow = n_ages, ncol = n_groups)
  empty_length <- matrix(NA_real_, nrow = n_lengths, ncol = n_groups)
  stock_age <- if (is.null(sim$N)) empty_age else sim$N[, years, drop = FALSE]
  stock_length <- if (is.null(sim$N_at_length)) empty_length else
    sim$N_at_length[, years, drop = FALSE]
  selection <- sim$survey_selection
  key <- selection$key
  if (is.null(key) && is.null(sim$N_at_length) && is.function(sim$sim_length)) {
    key <- .growth_key(sim)
  }
  if (!is.null(key)) {
    key <- key[as.character(sim$lengths), as.character(sim$ages), drop = FALSE]
    if (is.null(sim$N_at_length)) stock_length <- key %*% stock_age
  }
  frame_age <- available_age <- empty_age
  frame_length <- available_length <- empty_length

  if (!is.null(sim$survey_frame) && all(c("cell", "year", "age", "N") %in% names(sim$sp_N))) {
    spatial <- as.data.table(sim$sp_N)[, c("cell", "year", "age", "N"), with = FALSE]
    cells <- merge(sim$survey_frame, spatial, by = c("cell", "year"), allow.cartesian = TRUE)
    sums <- cells[, list(N_frame = sum(N)), by = c("sim", "year", "age")]
    ages <- roster[rep(seq_len(n_groups), each = n_ages)]
    ages$age <- rep(sim$ages, times = n_groups)
    values <- sums[ages, on = c("sim", "year", "age")][["N_frame"]]
    frame_groups <- unique(sim$survey_frame[, c("sim", "year"), with = FALSE])
    empty_groups <- ages[!frame_groups, on = c("sim", "year"), which = TRUE]
    values[empty_groups] <- 0
    frame_age <- matrix(values, nrow = n_ages, ncol = n_groups)

    if (!is.null(key)) frame_length <- key %*% frame_age
    if (!is.null(selection)) {
      if (selection$by == "age") {
        available_age <- frame_age * selection$values[as.character(sim$ages)]
        available_length <- key %*% available_age
      } else {
        selected_key <- key * selection$values[as.character(sim$lengths)]
        available_age <- frame_age * colSums(selected_key)
        available_length <- selected_key %*% frame_age
      }
    }
  }

  total <- data.table(roster, I = colSums(available_age),
                       N_frame = colSums(frame_age), N_stock = colSums(stock_age))
  age <- data.table(roster[rep(seq_len(n_groups), each = n_ages)],
                     age = rep(sim$ages, times = n_groups),
                     I = as.numeric(available_age), N_frame = as.numeric(frame_age),
                     N_stock = as.numeric(stock_age))
  length_map <- .strat_length_map(sim, sim$strat_length_group)
  lengths <- data.table(sim = integer(), year = numeric(), length = numeric(),
                         I = numeric(), N_frame = numeric(), N_stock = numeric())
  if (length(length_map)) {
    available_length <- rowsum(available_length, group = length_map, reorder = TRUE)
    frame_length <- rowsum(frame_length, group = length_map, reorder = TRUE)
    stock_length <- rowsum(stock_length, group = length_map, reorder = TRUE)
    bins <- as.numeric(rownames(available_length))
    lengths <- data.table(roster[rep(seq_len(n_groups), each = length(bins))],
                           length = rep(bins, times = n_groups),
                           I = as.numeric(available_length), N_frame = as.numeric(frame_length),
                           N_stock = as.numeric(stock_length))
  }
  list(total, lengths, age)
}

.strat_error_one <- function(sim) {
  I <- I_hat <- error <- estimate_status <- truth_status <- age <- NULL
  roster <- unique(sim$total_strat[, c("sim", "year"), with = FALSE])
  truth <- .strat_frame_truth(sim, roster)
  for (i in seq_along(.strat_estimates)) {
    trait <- list(NULL, "length", "age")[[i]]
    keys <- c("sim", "year", trait)
    table <- sim[[.strat_estimates[i]]]
    columns <- c(keys, "total", intersect("estimate_status", names(table)))
    table <- table[, columns, with = FALSE]
    setnames(table, "total", "I_hat")
    if (is.null(table$estimate_status)) {
      table[, estimate_status := ifelse(is.finite(I_hat), "ok", "missing_estimate")]
    }
    table <- merge(table, truth[[i]], by = keys, all.x = TRUE)
    unavailable <- if (is.null(sim$survey_frame)) "frame_unknown" else "truth_unavailable"
    if (i == 2L && !is.null(sim$survey_frame) &&
        is.null(.strat_length_map(sim, sim$strat_length_group))) {
      unavailable <- "length_bins_mismatch"
    }
    table[, truth_status := ifelse(is.finite(I), "ok", unavailable)]
    if (i == 3L) table[is.na(age), truth_status := "unassigned_age"]
    table[, error := ifelse(estimate_status == "ok" & truth_status == "ok", I_hat - I, NA_real_)]
    sim[[.strat_errors[i]]] <- setkeyv(table, keys)
    sim[[.strat_error_stats[i]]] <- .comparable_error_stats(table$error)
  }
  sim
}





