.caa_table <- function(x, columns, label) {
  if (!is.data.frame(x) || !all(columns %in% names(x))) {
    stop(label, " must contain ", paste(columns, collapse = ", "), ".", call. = FALSE)
  }
  if (is.data.table(x)) return(x[, columns, with = FALSE])
  as.data.table(x[columns])
}

#' Estimate monitored commercial catch at age
#'
#' @param sim A simulation returned by [sim_fishery()] or [sim_fisheries()].
#'
#' @details An aged fish contributes its operation's reported catch in numbers
#' divided by the operation inclusion probability, the number of measured fish
#' in the operation, and its conditional aging probability. This supports
#' random aging and aging by length group. Weighted contributions are summed
#' by age and calibrated to the sum of reports from all operations in each
#' fishery, replicate, and year. The resulting ratio estimate describes
#' reported catch and can differ from realized catch when reporting is biased.
#'
#' A zero reported total gives zero estimates at every age. For positive
#' reported totals, estimates are missing if any positive report lacks a known
#' positive operation inclusion probability, a selected positive operation
#' lacks measured fish, measured fish in a positive operation lack positive
#' aging probabilities, or no aged fish from positive operations are observed.
#' Operations with positive inclusion probabilities can be unselected by
#' chance. Missing observations at individual ages in supported groups give
#' estimated zeros. Only observed reports, inclusion probabilities, and
#' biological samples are used.
#'
#' @return The simulation list with `fishery$caa`, a `data.table` with
#' `fishery`, `sim`, `year`, `age`, and estimated numbers of fish `n_hat`.
#' `fishery$caa_diagnostics` gives annual `reported_n`, operation counts
#' (`n_operations`, `n_positive_operations`, `n_sampled_operations`), biological
#' counts (`n_measured`, `n_aged`), the fraction of operations selected
#' (`operation_coverage`), and the fraction of reported catch in selected
#' operations (`reported_catch_coverage`). Fractions with zero denominators
#' are missing. Its `status` is `"estimated"`, `"zero_reported_catch"`,
#' `"no_operation_support"`, `"no_measurement_support"`, or
#' `"no_aging_support"`. Re-estimation clears the previous catch comparison.
#' @seealso [caa_error()]
#' @examples
#' \donttest{
#' population <- sim_abundance(ages = 1:3, years = 2001:2002) |>
#'   sim_distribution(grid = make_grid(res = c(20, 20)))
#' result <- population |>
#'   sim_fisheries(trawl = sim_fishery(target = 0.1, catch_model = sim_capture(0.1),
#'                                     max_operations = 10)) |>
#'   run_caa() |>
#'   caa_error()
#' result$fishery$caa
#' }
#' @export
run_caa <- function(sim) {
  reported_n <- sampling_probability <- sampled <- aged <- age_probability <- NULL
  n_measured <- n_aged <- aging_support <- log_weight <- weight <- proportion <- NULL
  n_hat <- i.n_hat <- i.reported_n <- NULL

  fishery <- sim$fishery
  fleets <- names(fishery$configurations)
  identifiers <- list(fleets, sim$replicates, sim$years, sim$ages)
  if (any(vapply(identifiers, function(x) !length(x) || anyNA(x) || anyDuplicated(x),
                 logical(1)))) {
    stop("run_caa() requires fishery configurations and unique replicate, year, and age identifiers.",
         call. = FALSE)
  }
  group_keys <- c("fishery", "sim", "year")
  operation_keys <- c(group_keys, "operation")
  age_keys <- c(group_keys, "age")
  operations <- .caa_table(fishery$setdet,
                            c(operation_keys, "reported_n", "sampling_probability", "sampled"),
                            "fishery$setdet")
  sample <- .caa_table(fishery$samp,
                        c(operation_keys, "measured", "aged", "age", "age_probability"),
                        "fishery$samp")
  roster <- CJ(fishery = fleets, sim = sim$replicates, year = sim$years)
  if (anyNA(operations[, operation_keys, with = FALSE]) ||
      anyDuplicated(operations, by = operation_keys) ||
      nrow(operations[!roster, on = group_keys])) {
    stop("Fishery operation keys must be unique and match the simulation roster.", call. = FALSE)
  }
  if (!is.numeric(operations$reported_n) ||
      any(!is.finite(operations$reported_n) | operations$reported_n < 0)) {
    stop("Operation reports must be finite nonnegative numbers of fish.", call. = FALSE)
  }
  if (!is.logical(operations$sampled) || anyNA(operations$sampled) ||
      !is.numeric(operations$sampling_probability) || !is.numeric(sample$age_probability)) {
    stop("Sampling indicators must be known logical values and probabilities must be numeric.",
         call. = FALSE)
  }
  if (!is.logical(sample$measured) || anyNA(sample$measured) ||
      !is.logical(sample$aged) || anyNA(sample$aged) ||
      any(sample$aged & !sample$measured) ||
      any(!sample$age[sample$aged] %in% sim$ages)) {
    stop("Aged fish must be measured and have observed ages in the simulation age grid.",
         call. = FALSE)
  }
  sample <- operations[sample, on = operation_keys]
  if (anyNA(sample$sampled) || any(!sample$sampled)) {
    stop("Biological samples must match selected operations.", call. = FALSE)
  }
  sample <- sample[sample$measured, ]
  counts <- sample[, list(n_measured = .N, n_aged = sum(aged),
                           aging_support = all(is.finite(age_probability) &
                                                 age_probability > 0 & age_probability <= 1)),
                   by = operation_keys]
  operations <- counts[operations, on = operation_keys]
  no_sample <- is.na(operations$n_measured)
  operations$n_measured[no_sample] <- 0L
  operations$n_aged[no_sample] <- 0L
  operations$aging_support[no_sample] <- TRUE

  diagnostics <- operations[, {
    total <- sum(reported_n)
    positive <- reported_n > 0
    pi <- sampling_probability[positive]
    status <- "estimated"
    if (total == 0) {
      status <- "zero_reported_catch"
    } else if (any(!is.finite(pi) | pi <= 0 | pi > 1)) {
      status <- "no_operation_support"
    } else if (any(positive & sampled & n_measured == 0)) {
      status <- "no_measurement_support"
    } else if (any(!aging_support[positive]) || !any(n_aged[positive] > 0)) {
      status <- "no_aging_support"
    }
    list(reported_n = total, n_operations = .N, n_positive_operations = sum(positive),
         n_sampled_operations = sum(sampled), n_measured = sum(n_measured), n_aged = sum(n_aged),
         operation_coverage = mean(sampled),
         reported_catch_coverage = if (total > 0) sum(reported_n[sampled]) / total else NA_real_,
         status = status)
  }, by = group_keys]
  if (any(!is.finite(diagnostics$reported_n))) {
    stop("Annual sums of reported catch must be finite.", call. = FALSE)
  }
  diagnostics <- diagnostics[roster, on = group_keys]
  no_operations <- is.na(diagnostics$n_operations)
  for (column in c("reported_n", "n_operations", "n_positive_operations",
                    "n_sampled_operations", "n_measured", "n_aged")) {
    diagnostics[[column]][no_operations] <- 0L
  }
  diagnostics$status[no_operations] <- "zero_reported_catch"
  setkeyv(diagnostics, group_keys)

  estimates <- CJ(fishery = fleets, sim = sim$replicates, year = sim$years, age = sim$ages)
  estimates$n_hat <- rep(NA_real_, nrow(estimates))
  supported <- diagnostics$status %in% c("estimated", "zero_reported_catch")
  estimates[diagnostics[supported, ], on = group_keys, n_hat := 0]

  estimated_groups <- diagnostics[diagnostics$status == "estimated", group_keys, with = FALSE]
  observed <- sample[estimated_groups, on = group_keys, nomatch = 0L]
  observed <- observed[observed$aged & observed$reported_n > 0, ]
  if (nrow(observed)) {
    observed <- counts[observed, on = operation_keys]
    ## A common scale cancels in calibration and prevents overflow at small probabilities.
    observed[, log_weight := log(reported_n) - log(sampling_probability) -
                               log(n_measured) - log(age_probability)]
    observed[, weight := exp(log_weight - max(log_weight)), by = group_keys]
    totals <- observed[, list(weight = sum(weight)), by = age_keys]
    totals[, proportion := weight / sum(weight), by = group_keys]
    totals[diagnostics, on = group_keys, n_hat := proportion * i.reported_n]
    estimates[totals, on = age_keys, n_hat := i.n_hat]
  }

  sim$fishery$caa <- estimates
  sim$fishery$caa_diagnostics <- diagnostics
  sim$fishery$caa_error <- NULL
  sim
}

#' Compare monitored catch at age with realized catch
#'
#' @param sim A simulation returned by [run_caa()].
#' @return The simulation list with `fishery$caa_error`, a `data.table` with
#' the estimate keys and `n_hat`, realized catch numbers `n_true`, signed
#' `error` (estimate minus realized catch), and `truth_status`. When the
#' realized catch-at-age table is present, absent ages have zero realized
#' catch and truth is `"available"`. If that table is absent, truth and
#' errors are missing and truth is `"unavailable"`. Missing estimates retain
#' missing errors. The comparison does not alter estimates or diagnostics.
#' @export
caa_error <- function(sim) {
  n <- n_true <- i.n_true <- NULL
  if (is.null(sim$fishery$caa)) {
    stop("Run run_caa() before caa_error().", call. = FALSE)
  }
  keys <- c("fishery", "sim", "year", "age")
  comparison <- .caa_table(sim$fishery$caa, c(keys, "n_hat"), "fishery$caa")
  truth <- sim$fishery$catch_at_age
  comparison$n_true <- rep(if (is.null(truth)) NA_real_ else 0, nrow(comparison))
  comparison$truth_status <- rep(if (is.null(truth)) "unavailable" else "available",
                                  nrow(comparison))
  if (!is.null(truth)) {
    truth <- .caa_table(truth, c(keys, "n"), "fishery$catch_at_age")
    if (!is.numeric(truth$n) || any(!is.finite(truth$n) | truth$n < 0) ||
        nrow(truth[!comparison, on = keys])) {
      stop("Realized catch must contain finite nonnegative counts on the estimate grid.",
           call. = FALSE)
    }
    truth <- truth[, list(n_true = sum(n)), by = keys]
    comparison[truth, on = keys, n_true := i.n_true]
  }
  comparison$error <- comparison$n_hat - comparison$n_true
  setcolorder(comparison, c(keys, "n_hat", "n_true", "error", "truth_status"))
  setkeyv(comparison, keys)
  sim$fishery$caa_error <- comparison
  sim
}
