.program_configs <- function(configurations, type) {
  labels <- names(configurations)
  if (!length(configurations) || is.null(labels) || anyNA(labels) ||
      any(!nzchar(trimws(labels))) || anyDuplicated(labels)) {
    stop("Supply one or more configurations with unique, nonempty program names.", call. = FALSE)
  }
  if (!all(vapply(configurations, inherits, logical(1), what = type))) {
    stop("Every program must be a ", type, " configuration (omit sim when constructing it).",
         call. = FALSE)
  }
  configurations[order(labels, method = "radix")]
}

#' Simulate multiple fisheries with user-supplied allocations
#'
#' Declares all fisheries upfront and interleaves their operations, using one
#' population and one age-specific mortality budget per replicate and year.
#'
#' @param sim A simulation list returned by [sim_distribution()].
#' @param ... Named `sim_fishery_spec` configurations made with [sim_fishery()]
#' without its `sim` argument. Names must be unique and nonempty.
#' @param n_sims Common replicate count; `NULL` inherits existing batch metadata,
#' otherwise one. Do not supply counts inside configurations.
#' @details All fisheries inherit the weight relationship stored by
#' [sim_abundance()]. Weights belong to the population, not individual programs.
#'
#' Every fishery's reference target uses the original population and Z.
#' Fisheries may use different target types. At each age/year their combined
#' reference F must not exceed Z. Incompatible configurations produce an
#' informative error rather than being normalized. Reference catches in numbers
#' are rounded jointly; kilogram targets remain in kilograms. See [sim_fishery()]
#' for target definitions, units, validation and shortfalls.
#'
#' At each operation the next active fishery is chosen uniformly. Catches deplete
#' shared age-length stocks and death budgets. A fishery retires on reaching its
#' target, exhausting accessible catching opportunities, or reaching its operation
#' limit. Shared stocks and access restrictions can constrain valid allocations.
#' Other deaths
#' are allocated once, after all fishing. Interleaving defines operation order,
#' not calendar time. All truth precedes reporting and biological sampling.
#'
#' Configurations are processed in canonical name order, so registration order
#' does not affect same-seed results. Randomness uses the caller's RNG state.
#' Calling this function again runs a fresh scenario from the original annual
#' population; it replaces the fishery results, rather than appending removals.
#' Future annual populations are not changed. See [sim_fishery()] for capture,
#' monitoring, units, and the conditional interpretation of effort and CPUE.
#'
#' @return The simulation list with `replicates` and combined `fishery` results as
#' described by [sim_fishery()]. Fishery identifiers distinguish catches,
#' operations, spatial removals and annual CPUE; effort units are never pooled.
#' The mortality ledger and population are shared across fisheries.
#'
#' @seealso `vignette("fisheries", package = "SimSurvey")` for the joint workflow
#' and a mathematical appendix.
#'
#' @examples
#' trawl <- sim_fishery(target = 0.2, catch_model = sim_capture(q = 0.1))
#' gillnet <- sim_fishery(target = 0.1, catch_model = sim_capture(q = 0.05))
#' # population |> sim_fisheries(trawl = trawl, gillnet = gillnet)
#' @export
sim_fisheries <- function(sim, ..., n_sims = NULL) {
  configurations <- .program_configs(list(...), "sim_fishery_spec")
  if (any(vapply(configurations, function(x) !is.null(x$n_sims), logical(1)))) {
    stop("Supply common n_sims on sim_fisheries(), not inside configurations.")
  }
  .run_fisheries(sim, configurations, n_sims = n_sims)
}

#' Simulate independent survey programs on a shared population
#'
#' Each program observes the original rounded annual population snapshot,
#' including after [sim_fisheries()]. Surveys do not deplete another program's
#' population or consume the fishing mortality budget.
#'
#' @param sim A simulation list returned by [sim_distribution()] or
#' [sim_fisheries()].
#' @param ... Named `sim_survey_spec` configurations created with [sim_survey()]
#' without `sim`. Names must be unique and nonempty.
#' @param n_sims Common replicate count; `NULL` inherits existing batch metadata,
#' otherwise one. Do not supply counts inside configurations.
#'
#' @details Configurations run in canonical name order using caller-controlled
#' randomness. All programs' locations and catches are generated before any
#' program's biological subsampling, so changing sampling caps does not change
#' another program's locations or catches. Mutable tables are isolated between programs.
#'
#' The common population is stored once. Each program retains abundance
#' available to the survey across the full grid
#' (`I`, `I_at_length`, and keyed spatial `sp_I`), set details, biological samples,
#' summaries, evaluated catchability, its eligible survey frame when known,
#' and its configuration. The existing `surveys` design table used by
#' [test_surveys()] is preserved. Repeated calls replace program results with a
#' fresh simulation; they do not append or replay earlier programs.
#'
#' @seealso `vignette("fisheries", package = "SimSurvey")` for multiple survey
#' programs and program extraction. [strat_error()] compares estimates with
#' abundance available to the survey within its intended footprint after
#' applying `q` at age or length.
#'
#' @return A simulation list with common `replicates` and named `survey_programs`.
#' Use [select_survey()] to create a standard single-survey view for existing
#' analyses and plots. Top-level stale survey observations and analyses are cleared.
#'
#' @examples
#' survey_A <- sim_survey(q = sim_logistic(), q_by = "age")
#' survey_B <- sim_survey(q = sim_logistic(x50 = 30), q_by = "length")
#' # population |> sim_surveys(survey_A = survey_A, survey_B = survey_B)
#' @export
sim_surveys <- function(sim, ..., n_sims = NULL) {
  configurations <- .program_configs(list(...), "sim_survey_spec")
  if (any(vapply(configurations, function(x) isTRUE(attr(x, "n_sims_supplied")), logical(1)))) {
    stop("Supply common n_sims on sim_surveys(), not inside configurations.")
  }
  n_sims <- .resolve_n_sims(sim, n_sims)
  sim <- .clear_survey(round_sim(sim))
  sim$survey_programs <- NULL
  truths <- lapply(configurations, function(config) {
    if (!is.null(config$custom_sets) &&
        (is.null(config$custom_sets$sim) || anyNA(config$custom_sets$sim) ||
         any(!config$custom_sets$sim %in% seq_len(n_sims)))) {
      stop("custom_sets must use the common batch replicate identifiers.")
    }
    args <- unclass(config)
    attributes(args) <- list(names = names(args))
    args$n_sims <- n_sims
    args$sim <- data.table::copy(sim)
    do.call(.survey_truth, args)
  })
  ## Finish all existing biological samples before drawing any maturity states.
  results <- lapply(truths, .survey_samples)
  programs <- lapply(seq_along(results), function(i) {
    result <- .survey_maturity(results[[i]], configurations[[i]]$maturity_sampling,
                               configurations[[i]]$light)
    program <- result[intersect(.survey_fields, names(result))]
    program$sp_I <- data.table::copy(as.data.table(result$sp_N)[, c("cell", "age", "year", "I")])
    program$configuration <- configurations[[i]]
    program
  })
  names(programs) <- names(configurations)
  sim$survey_programs <- programs
  sim$replicates <- seq_len(n_sims)
  sim
}

.survey_fields <- c("I", "I_at_length", "full_setdet", "setdet", "samp", "samp_totals",
                     "sample_truth", "survey_frame", "survey_replicates", "survey_selection")

.clear_survey <- function(sim) {
  fields <- c(.survey_fields, "selected_survey", "survey_configuration")
  sim[intersect(fields, names(sim))] <- NULL
  if ("I" %in% names(sim$sp_N)) {
    sim$sp_N <- data.table::copy(sim$sp_N)
    sim$sp_N$I <- NULL
  }
  .clear_strat_analysis(sim)
}

#' Extract one survey program for existing analyses
#'
#' Creates a standard single-survey view, restoring the selected program's
#' availability, survey frame, set details, samples and summaries. Stacked
#' analysis results are extracted for the selected program when present.
#' The input object and other program tables are not modified.
#'
#' @param sim A simulation list returned by [sim_surveys()].
#' @param survey One program name from `names(sim$survey_programs)`.
#'
#' @return A simulation list ready for [run_strat()], [strat_error()] and survey
#' plotting functions, retaining the shared population, fishery results, program
#' collection and existing `surveys` design table. `selected_survey` identifies
#' the chosen program; `survey_configuration` contains its settings.
#'
#' @examples
#' # programs |> select_survey("survey_A") |> run_strat() |> strat_error()
#' @export
select_survey <- function(sim, survey) {
  view <- .survey_view(data.table::copy(sim), survey)
  .restore_survey_analysis(view, sim, survey)
}

.survey_view <- function(sim, survey, restore_spatial = TRUE) {
  if (!is.character(survey) || length(survey) != 1 || is.na(survey) ||
      !survey %in% names(sim$survey_programs)) stop("survey must name one available survey program.")
  sim <- .clear_survey(sim)
  program <- sim$survey_programs[[survey]]
  keep <- intersect(.survey_fields, names(program))
  sim[keep] <- data.table::copy(program[keep])
  if (restore_spatial && !is.null(sim$sp_N) && !is.null(program$sp_I)) {
    spatial <- as.data.table(sim$sp_N)[, c("cell", "age", "year")]
    sim$sp_N <- data.table::copy(sim$sp_N)
    sim$sp_N$I <- program$sp_I[spatial, on = c("cell", "age", "year")][["I"]]
  }
  sim$selected_survey <- survey
  sim$survey_configuration <- program$configuration
  sim
}
