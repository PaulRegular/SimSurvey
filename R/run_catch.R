.catch_table <- function(x, columns, label) {
  if (!is.data.frame(x) || !all(columns %in% names(x))) {
    stop(label, " must contain ", paste(columns, collapse = ", "), ".", call. = FALSE)
  }
  if (is.data.table(x)) return(x[, columns, with = FALSE])
  as.data.table(x[columns])
}

.catch_trait_ratios <- function(sample, keys, values, log_weights, valid, unsupported) {
  value <- log_weight <- weight <- status <- i.estimate <- NULL
  rows <- sample[, keys, with = FALSE]
  rows$value <- values
  rows$log_weight <- log_weights
  rows$valid <- valid
  bins <- rows[, list(estimate = NA_real_,
                       status = if (all(valid)) "no_bin_sample" else unsupported), by = keys]
  rows <- rows[bins[status == "no_bin_sample", ], on = keys, nomatch = 0L]
  rows <- rows[is.finite(value) & is.finite(log_weight), ]
  if (nrow(rows)) {
    rows[, weight := exp(log_weight - max(log_weight)), by = keys]
    ratios <- .biology_ratio(rows, keys, rows$value, rows$weight)
    bins[ratios, on = keys, `:=`(estimate = i.estimate, status = "estimated")]
  }
  bins
}

.catch_biology <- function(estimates, sample, keys, trait) {
  kg <- sample$weight
  if (!is.numeric(kg)) kg <- rep(NA_real_, nrow(sample))
  means <- .catch_trait_ratios(sample, keys, kg, sample$log_weight,
                               is.finite(kg) & kg >= 0, "no_weight_support")
  rows <- estimates[means, on = keys, which = TRUE]
  set(estimates, rows, "mean_weight_kg", means$estimate)
  set(estimates, rows, "weight_status", means$status)

  mature <- sample$mature
  if (!is.numeric(mature) && !is.logical(mature)) mature <- rep(NA_real_, nrow(sample))
  examined <- sample$maturity_observed
  if (!is.logical(examined)) examined <- rep(NA, nrow(sample))
  probability <- sample$maturity_probability
  if (!is.numeric(probability)) probability <- rep(NA_real_, nrow(sample))
  positive <- is.finite(probability) & probability > 0 & probability <= 1
  ## Maturity is observed either on all measured fish or on the nested aged fish.
  ## For age estimates, the existing weight already includes the aging stage.
  expected <- probability %in% 1 | sample$aged
  valid <- positive & (probability %in% 1 | probability == sample$age_probability) &
    !is.na(examined) & examined == expected & (!examined | mature %in% 0:1)
  valid[is.na(valid)] <- FALSE
  mature[is.na(examined) | !examined] <- NA
  log_weights <- sample$log_weight
  if (trait == "length") log_weights <- log_weights - log(ifelse(positive, probability, NA_real_))
  means <- .catch_trait_ratios(sample, keys, mature, log_weights, valid, "no_maturity_support")
  rows <- estimates[means, on = keys, which = TRUE]
  set(estimates, rows, "prop_mature", means$estimate)
  set(estimates, rows, "maturity_status", means$status)
  estimates
}

.catch_composition <- function(estimates, sample, diagnostics, group_keys, trait) {
  reported_n <- sampling_probability <- n_measured <- age_probability <- NULL
  log_weight <- count_weight <- proportion <- n_hat <- i.n_hat <- i.reported_n <- NULL

  status <- diagnostics[[paste0(trait, "_status")]]
  estimates$n_hat <- rep(NA_real_, nrow(estimates))
  estimates$mean_weight_kg <- rep(NA_real_, nrow(estimates))
  estimates$prop_mature <- rep(NA_real_, nrow(estimates))
  group_rows <- diagnostics[estimates, on = group_keys, which = TRUE]
  bin_status <- status[group_rows]
  bin_status[bin_status == "estimated"] <- "no_bin_sample"
  estimates$weight_status <- bin_status
  estimates$maturity_status <- bin_status
  supported <- status %in% c("estimated", "zero_reported_catch")
  estimates[diagnostics[supported, ], on = group_keys, n_hat := 0]

  groups <- diagnostics[status == "estimated", group_keys, with = FALSE]
  observed <- sample[groups, on = group_keys, nomatch = 0L]
  if (trait == "age") observed <- observed[observed$aged, ]
  if (nrow(observed)) {
    ## A common scale cancels in calibration and prevents overflow at small probabilities.
    observed[, log_weight := log(reported_n) - log(sampling_probability) - log(n_measured)]
    if (trait == "age") observed[, log_weight := log_weight - log(age_probability)]
    observed[, count_weight := exp(log_weight - max(log_weight)), by = group_keys]
    keys <- c(group_keys, trait)
    weights <- observed[, list(count_weight = sum(count_weight)), by = keys]
    weights[, proportion := count_weight / sum(count_weight), by = group_keys]
    weights[diagnostics, on = group_keys, n_hat := proportion * i.reported_n]
    estimates[weights, on = keys, n_hat := i.n_hat]
    estimates <- .catch_biology(estimates, observed, keys, trait)
  }
  estimates
}

#' Estimate monitored commercial catch
#'
#' @param sim A simulation returned by [sim_fishery()] or [sim_fisheries()].
#'
#' @details Total catch is the sum of reported numbers from all operations in
#' each fishery, replicate, and year. For the length composition, each measured
#' fish contributes its operation's reported catch divided by the operation
#' inclusion probability and the number of measured fish in that operation.
#' Length estimates use the population's native length groups.
#'
#' For the age composition, each aged fish's weight also divides by its
#' conditional aging probability. This supports random aging and aging by
#' length group. Each composition is calibrated separately to the reported
#' total. These estimates can differ from realized catch when reporting is
#' biased. Only observed reports, inclusion probabilities, and biological
#' samples are used.
#'
#' A zero reported total gives zero estimates at every length and age.
#' With positive reported catch, composition estimates require known positive
#' inclusion probabilities for every operation reporting catch and measurements
#' from every selected operation reporting catch. At least one measured fish
#' must come from an operation with a positive report.
#' Operations with positive inclusion probabilities can be unselected by
#' chance. When measurement support is sufficient, missing observed ages,
#' missing or invalid conditional aging probabilities, and absence of aged
#' fish can leave age estimates unsupported while length estimates remain
#' available. Unobserved lengths or ages within supported groups have
#' estimated counts of zero.
#'
#' Mean weight and mature proportion are weighted means of the observed
#' values within each length or age group. Mean weight uses the length
#' weights for measured fish and the age weights for aged fish. Maturity at
#' length additionally accounts for the probability of maturity examination
#' among measured fish. That probability is one with
#' `maturity_sampling = "measured"` and the conditional aging probability
#' with `"aged"`. At age, the aging weight accounts for the joint observation
#' of age and maturity under either protocol. Calibration to total catch
#' cancels in each mean. Disabling maturity sampling leaves its means missing.
#'
#' @return The simulation list with three keyed data tables in
#' `fishery$catch_estimates` named `total_catch`, `catch_at_length`, and `catch_at_age`.
#' Each has `fishery`, `sim`, `year`, and estimated numbers of fish `n_hat`.
#' The composition tables also include `length` or `age`, `mean_weight_kg`
#' (kilograms per fish), `prop_mature` (a proportion), and `weight_status`
#' and `maturity_status`. A trait status is `"estimated"`, inherits an
#' unsupported composition status, or indicates `"no_bin_sample"`,
#' `"no_weight_support"`, or `"no_maturity_support"`. Means are missing
#' for empty groups and when the corresponding trait is unsupported.
#' A supported group of entirely immature fish has proportion zero.
#' `fishery$catch_diagnostics` gives annual `reported_n`, operation counts
#' (`n_operations`, `n_positive_operations`, `n_sampled_operations`), biological
#' counts (`n_measured`, `n_aged`), the fraction of operations selected
#' (`operation_coverage`), and the fraction of reported catch in selected
#' operations (`reported_catch_coverage`). Fractions with zero denominators
#' are missing. Both `length_status` and `age_status` can be `"estimated"`,
#' `"zero_reported_catch"`, `"no_operation_support"`, or
#' `"no_measurement_support"`. The `age_status` can also be
#' `"no_aging_support"`. Re-estimation clears previous catch comparisons
#' and their error summaries.
#' @seealso [catch_error()]
#' @examples
#' \donttest{
#' population <- sim_abundance(ages = 1:3, years = 2001:2002) |>
#'   sim_distribution(grid = make_grid(res = c(20, 20)))
#' result <- population |>
#'   sim_fisheries(trawl = sim_fishery(target = 0.1, catch_model = sim_capture(0.1),
#'                                     max_operations = 10)) |>
#'   run_catch() |>
#'   catch_error()
#' result$fishery$catch_estimates$catch_at_length
#' }
#' @export
run_catch <- function(sim) {
  reported_n <- sampling_probability <- sampled <- aged <- age <- age_probability <- NULL
  n_measured <- n_aged <- aging_support <- NULL

  fishery <- sim$fishery
  fleets <- names(fishery$configurations)
  identifiers <- list(fleets, sim$replicates, sim$years, sim$lengths, sim$ages)
  if (any(vapply(identifiers, function(x) !length(x) || anyNA(x) || anyDuplicated(x),
                 logical(1)))) {
    stop("run_catch() requires fishery configurations and unique replicate, year, length, and age identifiers.",
         call. = FALSE)
  }
  group_keys <- c("fishery", "sim", "year")
  operation_keys <- c(group_keys, "operation")
  operations <- .catch_table(fishery$setdet,
                              c(operation_keys, "reported_n", "sampling_probability", "sampled"),
                              "fishery$setdet")
  biology <- intersect(c("weight", "mature", "maturity_observed", "maturity_probability"),
                        names(fishery$samp))
  sample <- .catch_table(fishery$samp,
                          c(operation_keys, "measured", "length", "aged", "age", "age_probability",
                            biology),
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
  probabilities <- list(operations$sampling_probability, sample$age_probability)
  if (!is.logical(operations$sampled) || anyNA(operations$sampled) ||
      any(vapply(probabilities, function(p) !is.numeric(p) && !(is.logical(p) && all(is.na(p))),
                 logical(1)))) {
    stop("Operation sampling indicators must be known logical values and probabilities numeric or missing.",
         call. = FALSE)
  }
  if (!is.logical(sample$measured) || anyNA(sample$measured) ||
      any(!sample$length[sample$measured] %in% sim$lengths)) {
    stop("Measured fish must have lengths in the native simulation length groups.", call. = FALSE)
  }
  if (!is.logical(sample$aged) || anyNA(sample$aged) || any(sample$aged & !sample$measured) ||
      any(!sample$age[sample$aged & !is.na(sample$age)] %in% sim$ages)) {
    stop("Aged fish must be measured and recorded ages must match the simulation age grid.",
         call. = FALSE)
  }
  sample <- operations[sample, on = operation_keys]
  if (anyNA(sample$sampled) || any(!sample$sampled)) {
    stop("Biological samples must match selected operations.", call. = FALSE)
  }
  sample <- sample[sample$measured, ]
  counts <- sample[, list(n_measured = .N, n_aged = sum(aged),
                           aging_support = all(is.finite(age_probability) &
                                                 age_probability > 0 & age_probability <= 1) &&
                             !any(aged & is.na(age))),
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
    length_status <- "estimated"
    if (total == 0) {
      length_status <- "zero_reported_catch"
    } else if (any(!is.finite(pi) | pi <= 0 | pi > 1)) {
      length_status <- "no_operation_support"
    } else if (any(positive & sampled & n_measured == 0) || !any(n_measured[positive] > 0)) {
      length_status <- "no_measurement_support"
    }
    age_status <- length_status
    if (length_status == "estimated" &&
        (any(!aging_support[positive]) || !any(n_aged[positive] > 0))) {
      age_status <- "no_aging_support"
    }
    list(reported_n = total, n_operations = .N, n_positive_operations = sum(positive),
         n_sampled_operations = sum(sampled), n_measured = sum(n_measured), n_aged = sum(n_aged),
         operation_coverage = mean(sampled),
         reported_catch_coverage = if (total > 0) sum(reported_n[sampled]) / total else NA_real_,
         length_status = length_status, age_status = age_status)
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
  diagnostics$length_status[no_operations] <- "zero_reported_catch"
  diagnostics$age_status[no_operations] <- "zero_reported_catch"
  setkeyv(diagnostics, group_keys)

  totals <- diagnostics[, c(group_keys, "reported_n"), with = FALSE]
  setnames(totals, "reported_n", "n_hat")
  lengths <- CJ(fishery = fleets, sim = sim$replicates, year = sim$years, length = sim$lengths)
  ages <- CJ(fishery = fleets, sim = sim$replicates, year = sim$years, age = sim$ages)
  observed <- counts[sample[sample$reported_n > 0, ], on = operation_keys]
  sim$fishery$catch_estimates <- list(
    total_catch = totals,
    catch_at_length = .catch_composition(lengths, observed, diagnostics, group_keys, "length"),
    catch_at_age = .catch_composition(ages, observed, diagnostics, group_keys, "age")
  )
  sim$fishery$catch_diagnostics <- diagnostics
  sim$fishery$catch_error <- NULL
  sim$fishery$catch_weight_error_stats <- NULL
  sim
}

.catch_comparison <- function(estimates, truth, keys, estimate_label, truth_label) {
  n <- n_true <- i.n_true <- kg <- mean_weight_true_kg <- i.mean_weight_true_kg <- NULL
  comparison <- .catch_table(estimates, c(keys, "n_hat"), estimate_label)
  comparison$n_true <- rep(if (is.null(truth)) NA_real_ else 0, nrow(comparison))
  comparison$truth_status <- rep(if (is.null(truth)) "unavailable" else "available",
                                  nrow(comparison))
  composition <- "age" %in% keys || "length" %in% keys
  if (composition) {
    comparison$mean_weight_kg <- if ("mean_weight_kg" %in% names(estimates)) {
      estimates$mean_weight_kg
    } else rep(NA_real_, nrow(comparison))
    comparison$mean_weight_true_kg <- rep(NA_real_, nrow(comparison))
  }
  if (!is.null(truth)) {
    has_weight <- composition && is.data.frame(truth) && is.numeric(truth[["kg"]])
    columns <- c(keys, "n", if (has_weight) "kg")
    truth <- .catch_table(truth, columns, truth_label)
    if (!is.numeric(truth$n) || any(!is.finite(truth$n) | truth$n < 0) ||
        nrow(truth[!comparison, on = keys])) {
      stop(truth_label, " must contain finite nonnegative counts on the estimate grid.",
           call. = FALSE)
    }
    totals <- truth[, list(n_true = sum(n)), by = keys]
    if (any(!is.finite(totals$n_true))) {
      stop(truth_label, " has nonfinite sums of realized catch.", call. = FALSE)
    }
    comparison[totals, on = keys, n_true := i.n_true]
    if (has_weight) {
      weights <- truth[n > 0, list(mean_weight_true_kg =
        if (all(is.finite(kg) & kg >= 0)) sum(kg) / sum(n) else NA_real_), by = keys]
      weights$mean_weight_true_kg[!is.finite(weights$mean_weight_true_kg)] <- NA_real_
      comparison[weights, on = keys, mean_weight_true_kg := i.mean_weight_true_kg]
    }
  }
  comparison$error <- comparison$n_hat - comparison$n_true
  traits <- NULL
  if (composition) {
    comparison$mean_weight_error_kg <- comparison$mean_weight_kg - comparison$mean_weight_true_kg
    traits <- c("mean_weight_kg", "mean_weight_true_kg", "mean_weight_error_kg")
  }
  setcolorder(comparison, c(keys, "n_hat", "n_true", "error", "truth_status", traits))
  setkeyv(comparison, keys)
  comparison
}

#' Compare monitored catch with realized catch
#'
#' @param sim A simulation returned by [run_catch()].
#' @return The simulation list with `fishery$catch_error`, a list of keyed
#' `total_catch`, `catch_at_length`, and `catch_at_age` comparisons. Each
#' contains the estimate keys and `n_hat`, realized numbers `n_true`, signed
#' `error` (estimate minus realized catch), and `truth_status`.
#'
#' Total truth comes from realized catches in operation records. The
#' composition comparisons use their respective realized catch tables. When
#' the corresponding truth is present, missing groups have zero realized
#' catch and truth is `"available"`. When it is absent, truth and errors are
#' missing and truth is `"unavailable"`. Missing estimates retain missing
#' errors. The comparisons do not alter estimates or diagnostics.
#'
#' Length and age comparisons also include `mean_weight_kg`,
#' `mean_weight_true_kg` (realized kilograms divided by realized numbers), and
#' `mean_weight_error_kg`. Realized mean weight is missing when there is no
#' realized catch or its weight is unavailable. The corresponding tables in
#' `fishery$catch_weight_error_stats` summarize comparable errors across bins,
#' replicates, and years by fishery, with `ME`, `MAE`, `MSE`, and `RMSE`.
#' @export
catch_error <- function(sim) {
  mean_weight_error_kg <- NULL
  fishery <- sim$fishery
  estimates <- fishery$catch_estimates
  tables <- c("total_catch", "catch_at_length", "catch_at_age")
  if (!all(tables %in% names(estimates))) {
    stop("Run run_catch() before catch_error().", call. = FALSE)
  }
  keys <- c("fishery", "sim", "year")
  total_truth <- if ("n" %in% names(fishery$setdet)) fishery$setdet else NULL
  sim$fishery$catch_error <- list(
    total_catch = .catch_comparison(estimates$total_catch, total_truth, keys,
                                    "fishery$catch_estimates$total_catch", "fishery$setdet"),
    catch_at_length = .catch_comparison(estimates$catch_at_length, fishery$catch_at_length,
                                        c(keys, "length"), "fishery$catch_estimates$catch_at_length",
                                        "fishery$catch_at_length"),
    catch_at_age = .catch_comparison(estimates$catch_at_age, fishery$catch_at_age,
                                     c(keys, "age"), "fishery$catch_estimates$catch_at_age",
                                     "fishery$catch_at_age")
  )
  sim$fishery$catch_weight_error_stats <- lapply(sim$fishery$catch_error[tables[-1L]], function(table) {
    result <- table[, as.list(.comparable_error_stats(mean_weight_error_kg)), by = "fishery"]
    setkey(result, fishery)
    result
  })
  sim
}
