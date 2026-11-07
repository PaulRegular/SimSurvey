#' Simulate a fishery conditional on population mortality
#'
#' Allocate existing mortality using a share, fully selected F, or catch target,
#' generate aggregate catches,
#' and simulate catch reporting and nested biological monitoring. Omitting `sim`
#' creates a configuration without running a simulation or drawing random values.
#'
#' @param sim A simulation list from [sim_distribution()]. May be omitted to
#' create a `sim_fishery_spec` configuration.
#' @param target A nonnegative scalar or vector aligned with simulation years.
#' Its meaning is set by `target_type`.
#' @param target_type One of `"share"` (default: proportion of Z before
#' selectivity, in `[0, 1]`), `"F_full"` (annual fishing mortality at selectivity
#' one), `"catch_n"` (fish), or `"catch_kg"` (kilograms).
#' @param catch_model Function `(available, operation)` returning one whole-fish
#' proposal per available row, between zero and its remaining `N`. See
#' [sim_capture()] for the table contract and a generic built-in model.
#' @param selectivity Function returning probabilities in `[0, 1]`, evaluated on
#' ages or growth-model length-group labels. Values are not rescaled.
#' @param selectivity_by Either `"age"` or `"length"`.
#' @param gear_amount,duration Positive settings, or functions of an operation
#' record returning one positive number. Their product is nominal effort.
#' @param gear_unit,time_unit Nonempty labels for gear and time units.
#' @param targeting Nonnegative exponent controlling preference for high estimated
#' vulnerable density. Zero gives uniform locations.
#' @param knowledge_sd Nonnegative map-error SD, expressed as a multiple of initial
#' mean accessible vulnerable density. Errors are redrawn for each replicate/year
#' and persist across its operations. Location choice uses that year's current
#' vulnerable density plus these errors, updating density after catches. This is
#' noisy contemporaneous knowledge, without learning or a lagged distribution map.
#' @param exploration Proportion of location probabilities allocated uniformly,
#' in `[0, 1]`. All-zero estimated densities fall back to uniform locations.
#' @param subset_cells Logical expression evaluated on grid columns and `year`,
#' restricting accessible cells. By default all cells are accessible.
#' @param n_sims Replicate count. `NULL` inherits batch metadata, otherwise one.
#' An explicit count must agree with existing batch metadata.
#' @param max_operations Maximum operations per fishery, replicate and year.
#' @param reporting_factor,reporting_cv Mean and CV of operation-level lognormal
#' reporting multipliers. The same multiplier applies to numbers and kilograms.
#' @param sampling_coverage Expected proportion of operations selected for
#' biological sampling, in `[0, 1]`.
#' @param sampling_targeting Nonnegative preference for high reported catches.
#' Inclusion probabilities are calibrated to the requested expected coverage.
#' @param lengths_cap Maximum fish measured per selected operation.
#' @param ages_cap Maximum measured fish aged per operation (`"random"`) or per
#' operation and length group (`"stratified"`).
#' @param age_sampling Either `"stratified"` or `"random"`.
#' @param age_length_group Positive length-bin width for stratified age sampling.
#' @param light If `FALSE`, retain joint age-length catches and sample truth.
#' @inheritParams sim_survey
#'
#' @details Numeric controls may be scalar or year-aligned vectors; named vectors
#' must name every year. Gear callbacks receive identifiers, location, and resolved
#' gear settings. They must not implement hidden seeding. The capture table uses
#' `age`, `length`, `N`, `selectivity`, and `cell_area`; operations include
#' `fishery`, `sim`, `year`, `operation`, `cell`, `gear_amount`, `duration`, and
#' `effort`. Grid area and capture-model units must agree.
#'
#' For length-age probabilities P and selectivity s, mean selectivity at age is
#' `s_bar = sum(P * s)`. Share targets set `F_reference = target * Z * s_bar`;
#' fully selected F targets set `F_reference = target * s_bar`, without rescaling
#' selectivity. Catch targets scale that same selectivity pattern to the requested
#' expected numbers or kilograms using the Baranov equation under fixed Z.
#' Kilograms integrate weight jointly with length selectivity, not as a product
#' of separately averaged weight and selectivity.
#'
#' Targets must be finite and nonnegative. Catch targets cannot exceed starting
#' N or B, and must fit the mortality available under Z and selectivity.
#' Reference F, including combined fisheries, cannot exceed Z at any age/year.
#' Combined reference catches must also fit starting N and B, including when
#' length selectivity concentrates several fisheries on the same heavy fish.
#' Errors identify the fishery/year and applicable limits before stock allocation
#' draws. Valid targets can still have shortfalls due to access, capture response,
#' operation limits, integer rounding or shared-stock constraints.
#'
#' Shares, F_full and catch_n stop at the rounded reference total in numbers.
#' These reference catches are rounded jointly across fisheries, with the rounded
#' combined total capped at the shared integer death count. Catch_kg stops
#' at the requested kilograms using realized catch weights, not a converted
#' number quota. Whole fish are retained; final catches are thinned proportionally
#' in aggregate, then filled with randomly selected categories of fish that fit.
#' No individual commercial catch records are generated. A residual kilogram
#' shortfall can remain when no eligible whole fish fits; status is exhausted.
#' Zero-weight fish do not provide a positive kilogram catching opportunity.
#' Realized age composition may differ from the reference, so meeting a catch
#' target need not reproduce reference F or one exact realized F_full.
#' Age-specific death budgets preserve the rounded total expected deaths using
#' largest-remainder rounding. Fishing and residual other deaths share those
#' budgets. Residual deaths are sampled from remaining fish without replacement.
#' `reference` retains `N`, `Z`, shared integer `total_deaths`, reference and
#' realized catches, `U_realized = catch_realized / N`, and `F_realized` at age.
#' `total_deaths` is the whole-fish total implied by N and Z for that age/year,
#' balanced-rounded across ages. It includes all fishery and other deaths, and
#' is repeated in each fishery's reference rows; do not sum it across fisheries.
#' It is neither an individual catch target nor a remaining allowance. At the
#' end it equals the ledger's fishery plus other deaths. `death_budget` remains
#' an identical compatibility alias; new code should prefer `total_deaths`.
#' Realized F is the Baranov-equivalent `U_realized * Z / (1 - exp(-Z))` under
#' fixed upstream Z, not a fitted operation-time hazard. With positive N and zero
#' Z it is zero; with zero N both realized rates are missing. Integer rounding
#' can make combined equivalent F exceed Z in tiny populations; values are not
#' clipped. Matching the total number target need not match F at every age.
#'
#' Effort determines allocation and required operations conditional on the catch
#' budget; it cannot independently increase annual removals. The final operation
#' retains all specified effort when capped by the target (`quota_limited`).
#' Zero-catch operations retain their effort. Unreachable targets return a
#' shortfall and termination status.
#' A shortfall is the stopping target minus true catch, in `target_unit` (fish or
#' kg): `exhausted` means no eligible accessible catching
#' opportunities remain (including known zero capture/capacity); `operation_limit`
#' means the allowed operations ended before reaching the target. Unused fishing
#' allocation becomes other deaths under the fixed total-mortality budget.
#' Local removals do not feed back into later
#' pre-generated populations, and are not subtracted again from future abundance.
#' Other mortality can include other fisheries. All simulated catches are retained;
#' discards and seasonal scheduling are outside this model.
#' Unnamed mortality matrices from custom closures use the age/year order of
#' [sim_abundance()]; supplied dimension labels must agree with the population.
#'
#' Fishing truth is completed before monitoring. Biological samples are drawn
#' from actual catches, with ages nested within measured fish. Sample records
#' contain stage-specific inclusion probabilities and missing unobserved ages.
#' Catch, stock and mortality tables remain aggregate even for very large stocks.
#'
#' @return The simulation list with `replicates` and a `fishery` list containing
#' `setdet` (operation details), `catch_at_age`, `catch_at_length`, `sp_removals`,
#' `mortality_ledger`, `samp`, `annual`, `reference`, `settings`, and named
#' `configurations`. Identifiers link fishery, replicate (`sim`), year, operation,
#' and sample. True counts are integers; reported numbers can be fractional.
#' `annual` records `target_type`, `requested_target` (the supplied setting),
#' `reference_target` and `target` (continuous and operational stopping targets),
#' `target_unit`, `achieved`, and `shortfall`. It also retains expected reference
#' catches as `reference_target_n` and `reference_target_kg`, and
#' `F_full_reference` (missing for share targets, whose F also depends on Z at age).
#' `reference` includes `kg_reference` at age. Reference biomass uses the growth
#' key; a particular integer length-stock realization can have different biomass.
#' CPUE is catch divided by nominal effort; annual CPUE uses totals. Biomass and
#' biomass CPUE use the weight relationship stored by [sim_abundance()]. With
#' `light = FALSE`, `joint_catch`
#' and `sample_truth` are also retained. Without `sim`, a `sim_fishery_spec` list.
#'
#' @seealso `vignette("fisheries", package = "SimSurvey")` for controlled
#' comparisons of targeting, knowledge, reporting, and biological monitoring.
#' [run_catch()] summarizes reported catch and estimates its length and age
#' composition from monitoring records.
#'
#' @examples
#' capture <- sim_capture(q = 0.1, capacity_per_gear = 50)
#' specification <- sim_fishery(target = 0.2, catch_model = capture)
#' f_target <- sim_fishery(target = 0.1, target_type = "F_full", catch_model = capture)
#' number_target <- sim_fishery(target = 100, target_type = "catch_n", catch_model = capture)
#' mass_target <- sim_fishery(target = 50, target_type = "catch_kg", catch_model = capture)
#' \donttest{
#' population <- sim_abundance(ages = 1:3, years = 1:2) |>
#'   sim_distribution(grid = make_grid(res = c(20, 20)))
#' result <- sim_fishery(population, target = 0.2, catch_model = capture,
#'                       max_operations = 100)
#' result$fishery$annual
#' }
#' @export
sim_fishery <- function(sim, target, catch_model,
                        selectivity = sim_logistic(), selectivity_by = "age",
                        gear_amount = 1, duration = 1,
                        gear_unit = "gear_unit", time_unit = "time_unit",
                        targeting = 1, knowledge_sd = 0, exploration = 0.05,
                        subset_cells, n_sims = NULL,
                        max_operations = 10000, reporting_factor = 1, reporting_cv = 0,
                        sampling_coverage = 1, sampling_targeting = 0,
                        lengths_cap = 100, ages_cap = 5,
                        age_sampling = "stratified", age_length_group = 1,
                        light = TRUE, maturity_sampling = "aged",
                        target_type = "share") {
  if (missing(target)) stop("Supply target.", call. = FALSE)
  if (!is.character(target_type) || length(target_type) != 1 || is.na(target_type) ||
      !target_type %in% c("share", "F_full", "catch_n", "catch_kg")) {
    stop("target_type must be 'share', 'F_full', 'catch_n', or 'catch_kg'.", call. = FALSE)
  }
  access <- if (missing(subset_cells)) quote(TRUE) else substitute(subset_cells)
  config_names <- setdiff(names(formals(sim_fishery)), c("sim", "subset_cells"))
  config <- mget(config_names, envir = environment(), inherits = FALSE)
  config$subset_cells <- access
  config$subset_environment <- parent.frame()
  class(config) <- "sim_fishery_spec"
  if (missing(sim)) return(config)
  .run_fisheries(sim, list(fishery = config), n_sims = n_sims)
}

.fishery_controls <- function(config, sim, fishery = "fishery") {
  if (!is.character(config$target_type) || length(config$target_type) != 1 ||
      is.na(config$target_type) ||
      !config$target_type %in% c("share", "F_full", "catch_n", "catch_kg")) {
    stop("Fishery '", fishery, "': invalid target_type.", call. = FALSE)
  }
  years <- sim$years
  .check_maturity_model(sim, config$maturity_sampling)
  if (!is.function(config$catch_model)) stop("catch_model must be a function.")
  if (!is.character(config$gear_unit) || length(config$gear_unit) != 1 ||
      is.na(config$gear_unit) || !nzchar(config$gear_unit) ||
      !is.character(config$time_unit) || length(config$time_unit) != 1 ||
      is.na(config$time_unit) || !nzchar(config$time_unit)) {
    stop("gear_unit and time_unit must be nonempty labels.")
  }
  if (length(config$age_sampling) != 1 ||
      !config$age_sampling %in% c("random", "stratified")) {
    stop("age_sampling must be 'random' or 'stratified'.")
  }
  if (!is.logical(config$light) || length(config$light) != 1 || is.na(config$light)) {
    stop("light must be TRUE or FALSE.")
  }
  probabilities <- c("exploration", "sampling_coverage")
  nonnegative <- c("targeting", "knowledge_sd", "reporting_factor", "reporting_cv",
                   "sampling_targeting")
  integers <- c("max_operations", "lengths_cap", "ages_cap")
  result <- config
  result$target <- .year_values(config[["target"]], years, paste0("Fishery '", fishery, "': target"),
                                 upper = if (config$target_type == "share") 1 else Inf)
  for (name in probabilities) result[[name]] <- .year_values(config[[name]], years, name, upper = 1)
  for (name in nonnegative) result[[name]] <- .year_values(config[[name]], years, name)
  for (name in integers) result[[name]] <- .year_values(config[[name]], years, name, whole = TRUE)
  for (name in c("gear_amount", "duration")) {
    if (!is.function(config[[name]])) {
      result[[name]] <- .year_values(config[[name]], years, name, lower = .Machine$double.xmin)
    }
  }
  result$age_length_group <- .year_values(config$age_length_group, years, "age_length_group",
                                          lower = .Machine$double.xmin)
  result$selection <- .trait_selection(sim, config$selectivity, config$selectivity_by)
  result
}

.operation_setting <- function(setting, operation, y, name) {
  value <- if (is.function(setting)) setting(operation) else setting[y]
  if (!is.numeric(value) || length(value) != 1 || !is.finite(value) || value <= 0) {
    stop(name, " must resolve to one finite positive value per operation.", call. = FALSE)
  }
  unname(value)
}

.fishery_access <- function(config, cells, year) {
  data <- data.table::copy(cells)
  data$year <- year
  result <- eval(config$subset_cells, data, config$subset_environment)
  if (!is.logical(result) || !length(result) %in% c(1, nrow(data)) || anyNA(result)) {
    stop("subset_cells must return a nonmissing logical value per grid cell.")
  }
  rep(result, length.out = nrow(data))
}
