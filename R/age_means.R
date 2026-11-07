.biology_ratio <- function(data, keys, value, weights) {
  trait_value <- expansion <- NULL
  keep <- which(is.finite(value) & is.finite(weights) & weights > 0)
  data <- data[keep, keys, with = FALSE]
  if (!length(keep)) {
    data$estimate <- numeric()
    data$n <- integer()
    return(setkeyv(data, keys))
  }
  data$trait_value <- value[keep]
  data$expansion <- weights[keep]
  result <- data[, {
    relative <- expansion / max(expansion)
    scale <- max(abs(trait_value))
    mean <- if (scale == 0) 0 else scale * (sum((trait_value / scale) * relative) / sum(relative))
    list(estimate = mean, n = .N)
  }, by = keys]
  setkeyv(result, keys)
}

.age_means_labels <- function(x, label, whole = FALSE) {
  values <- if (is.numeric(x)) as.numeric(x) else suppressWarnings(as.numeric(as.character(x)))
  invalid <- !is.na(x) & (!is.finite(values) | values < 0)
  if (whole) invalid <- invalid | (!is.na(values) & values != floor(values))
  if (any(invalid)) {
    stop(label, " must have finite nonnegative numeric labels",
         if (whole) " in whole years" else "", ".", call. = FALSE)
  }
  values
}

#' Estimate biological means at age from length-stratified sampling
#'
#' @param length_totals Data frame with expanded fish counts `total` for each
#'   `length` group and stratum. Rows sharing a stratum and length are summed.
#' @param aged Data frame with one row per biological sample. Include `age`,
#'   the same `length` group labels, stratum keys, and the requested traits.
#'   Records with missing age are excluded from the age-length key.
#' @param metrics Names of numeric or logical trait columns. Use binary
#'   `mature` for proportion mature, or `weight` in kilograms for mean weight.
#'   Requesting `"length"` gives mean length using the bin labels.
#' @param strata Columns defining the groups within which to expand length
#'   samples. Supply all grouping keys needed to match the design, including
#'   time and spatial identifiers where applicable.
#' @param by Grouping columns to retain, as a subset of `strata`. The default
#'   retains all strata. An empty vector pools all strata.
#' @param ages Optional nonnegative whole ages to report. The default uses all
#'   observed age labels. Supply ages to report support when no fish were aged.
#'
#' @details Each aged fish represents its stratum-length total divided by the
#'   number of usable ages in that cell. Sum weighted trait values and weights
#'   at age across strata before dividing. This is the point estimator of
#'   Morgan and Hoenig (1997) for maturity and Morgan (1999) for weight.
#'   It assumes representative aging within each stratum and length.
#'
#'   Age counts use all available ages, independently of trait availability.
#'   An unknown length total or positive length total without an age sample
#'   makes the pooled composition unsupported. Missing traits among weighted
#'   aged fish make the corresponding age and trait estimate unavailable.
#'
#' @return A keyed data table with the requested `by` columns, `age`, `metric`,
#'   `estimate`, `n` (records with observed traits and positive expansion),
#'   and `status`. Statuses include `"ok"`, `"missing_length_total"`,
#'   `"no_aging_support"`, `"missing_trait"`, and `"no_weighted_fish"`.
#'   Unsupported means are missing.
#'
#' @references
#' Morgan, M. J. and Hoenig, J. M. (1997).
#' Estimating maturity-at-age from length stratified sampling.
#' Journal of Northwest Atlantic Fishery Science, 21, 51–63.
#' <https://doi.org/10.2960/J.v21.a4>
#'
#' Morgan, M. J. (1999).
#' The effect of a change in perception of length distribution of a population
#' on maturity-at-age, weight-at-age and spawning stock biomass.
#' Journal of Northwest Atlantic Fishery Science, 25, 141–150.
#' <https://doi.org/10.2960/J.v25.a13>
#'
#' @examples
#' length_totals <- data.frame(length = c(20, 40), total = c(90, 10))
#' aged <- data.frame(age = 3, length = c(20, 40), mature = c(0, 1),
#'                    weight = c(1, 3))
#' age_means(length_totals, aged, c("weight", "mature"))
#' @export
age_means <- function(length_totals, aged, metrics, strata = character(),
                      by = strata, ages = NULL) {
  age <- total <- n_aged <- n <- status <- estimate <- NULL
  reserved <- c("age", "length", "total", "n_aged", "metric", "status",
                "estimate", "n", ".overall")
  if (!is.character(strata) || anyNA(strata) || any(!nzchar(strata)) ||
      anyDuplicated(strata) || any(strata %in% c(reserved, "trait_value", "expansion")) ||
      !is.character(by) || anyDuplicated(by) || !all(by %in% strata)) {
    stop("strata must be unique grouping columns and by a subset of strata.", call. = FALSE)
  }
  if (!is.character(metrics) || !length(metrics) || anyNA(metrics) ||
      any(!nzchar(metrics)) || anyDuplicated(metrics) ||
      any(metrics %in% c(strata, setdiff(reserved, "length")))) {
    stop("metrics must name unique trait columns.", call. = FALSE)
  }
  length_keys <- c(strata, "length")
  total_cols <- c(length_keys, "total")
  sample_cols <- unique(c(length_keys, "age", metrics))
  if (!is.data.frame(length_totals) || !all(total_cols %in% names(length_totals)) ||
      !is.data.frame(aged) || !all(sample_cols %in% names(aged))) {
    stop("Supply the required length totals, age, trait and stratum columns.", call. = FALSE)
  }
  totals <- as.data.table(length_totals)[, total_cols, with = FALSE]
  samples <- as.data.table(aged)[, sample_cols, with = FALSE]
  totals$length <- .age_means_labels(totals$length, "Length")
  samples$length <- .age_means_labels(samples$length, "Length")
  samples$age <- .age_means_labels(samples$age, "Age", whole = TRUE)
  samples <- samples[!is.na(age)]
  if (anyNA(totals[, length_keys, with = FALSE]) ||
      anyNA(samples[, length_keys, with = FALSE])) {
    stop("Length and stratum keys must be complete.", call. = FALSE)
  }
  if ((!is.numeric(totals$total) &&
       !(is.logical(totals$total) && all(is.na(totals$total)))) ||
      any(!is.na(totals$total) & (!is.finite(totals$total) | totals$total < 0))) {
    stop("Length totals must be nonnegative finite counts or missing.", call. = FALSE)
  }
  for (metric in metrics) {
    value <- samples[[metric]]
    if (!(is.numeric(value) || is.logical(value)) ||
        any(!is.na(value) & !is.finite(value)) ||
        (metric == "mature" && any(!is.na(value) & !value %in% c(0, 1))) ||
        (metric == "weight" && any(value < 0, na.rm = TRUE))) {
      stop("Invalid values in trait column ", metric, ".", call. = FALSE)
    }
  }
  if (nrow(samples[!totals, on = length_keys])) {
    stop("Aged sample length and stratum keys must match the length totals.", call. = FALSE)
  }
  if (is.null(ages)) ages <- samples$age
  ages <- sort(unique(.age_means_labels(ages, "ages", whole = TRUE)))

  totals <- totals[, list(total = sum(total)), by = length_keys]
  if (any(is.infinite(totals$total))) stop("Summed length totals must be finite.", call. = FALSE)
  pooled <- !length(by)
  if (pooled) {
    by <- ".overall"
    totals$.overall <- rep(1L, nrow(totals))
  }
  keys <- c(by, "age")
  tally <- samples[, list(n_aged = .N), by = length_keys]
  bins <- tally[totals, on = length_keys]
  groups <- bins[, list(status = if (anyNA(total)) "missing_length_total" else
    if (any(total > 0 & is.na(n_aged))) "no_aging_support" else "ok"), by = by]
  template <- groups[rep(seq_len(nrow(groups)), each = length(ages))]
  template$age <- rep(ages, times = nrow(groups))
  expanded <- bins[samples, on = length_keys]
  weights <- expanded$total / expanded$n_aged
  results <- lapply(metrics, function(metric) {
    values <- expanded[[metric]]
    means <- .biology_ratio(expanded, keys, values, weights)
    missing <- unique(expanded[which(is.na(values) & weights > 0), keys, with = FALSE])
    means <- merge(template, means, by = keys, all.x = TRUE)
    means[is.na(n), n := 0L]
    means[missing, on = keys, status := ifelse(status == "ok", "missing_trait", status)]
    means[status == "ok" & n == 0L, status := "no_weighted_fish"]
    means[status != "ok", estimate := NA_real_]
    means$metric <- rep(metric, nrow(means))
    means
  })
  result <- rbindlist(results)
  if (pooled) {
    result$.overall <- NULL
    keys <- "age"
  }
  setcolorder(result, c(keys, "metric", "estimate", "n", "status"))
  setkeyv(result, c(keys, "metric"))
}
