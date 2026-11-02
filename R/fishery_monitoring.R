.fishery_biology <- function(catch, operation, control, y, operation_probability) {
  size <- min(sum(catch$n), control$lengths_cap[y])
  if (!size) return(NULL)
  measured <- .thin_counts(catch$n, size)
  ## Expand only the bounded biological sample, never the commercial catch.
  rows <- rep(seq_len(nrow(catch)), measured)
  sample <- data.table::copy(catch[rows, c("fishery", "sim", "year", "operation", "cell", "age", "length")])
  sample$true_age <- sample$age
  sample$operation_probability <- operation_probability
  sample$length_probability <- size / sum(catch$n)
  sample$age_probability <- 0
  sample$aged <- FALSE
  groups <- if (control$age_sampling == "random") rep(1, size) else
    group_lengths(sample$length, control$age_length_group[y])
  for (group in unique(groups)) {
    candidates <- which(groups == group)
    n_aged <- min(length(candidates), control$ages_cap[y])
    sample$age_probability[candidates] <- n_aged / length(candidates)
    if (n_aged) sample$aged[candidates[sample.int(length(candidates), n_aged)]] <- TRUE
  }
  sample$age[!sample$aged] <- NA_real_
  sample$measured <- TRUE
  sample
}

.monitor_fisheries <- function(truth, controls, years) {
  n <- kg <- reported_n <- reported_kg <- effort <- NULL
  operations <- truth$setdet
  catches <- truth$joint_catch
  operations$reporting_multiplier <- rep(NA_real_, nrow(operations))
  operations$sampling_probability <- rep(0, nrow(operations))
  operations$sampled <- rep(FALSE, nrow(operations))
  samples <- list()
  for (fishery in names(controls)) {
    control <- controls[[fishery]]
    for (year in years) {
      y <- match(year, years)
      for (replicate in unique(truth$annual$sim)) {
        rows <- which(operations$fishery == fishery & operations$year == year &
                        operations$sim == replicate)
        if (!length(rows)) next
        cv <- control$reporting_cv[y]
        factor <- control$reporting_factor[y]
        multiplier <- if (cv == 0 || factor == 0) rep(factor, length(rows)) else
          stats::rlnorm(length(rows), log(factor) - log1p(cv^2) / 2, sqrt(log1p(cv^2)))
        operations$reporting_multiplier[rows] <- multiplier
        reported <- operations$n[rows] * multiplier
        probabilities <- .sampling_probabilities(reported, control$sampling_coverage[y],
                                                 control$sampling_targeting[y])
        operations$sampling_probability[rows] <- probabilities
        sampled <- if (all(probabilities %in% c(0, 1))) probabilities == 1 else
          stats::runif(length(rows)) < probabilities
        operations$sampled[rows] <- sampled
        for (row in rows[sampled]) {
          caught <- catches[catches$operation == operations$operation[row], ]
          observations <- .fishery_biology(caught, operations[row, ], control, y,
                                           operations$sampling_probability[row])
          if (!is.null(observations)) samples[[length(samples) + 1L]] <- observations
        }
      }
    }
  }
  operations$reported_n <- operations$n * operations$reporting_multiplier
  operations$reported_kg <- operations$kg * operations$reporting_multiplier
  operations$cpue <- operations$n / operations$effort
  operations$reported_cpue <- operations$reported_n / operations$effort
  operations$biomass_cpue <- operations$kg / operations$effort
  operations$reported_biomass_cpue <- operations$reported_kg / operations$effort
  empty_sample <- data.table(fishery = character(), sim = integer(), year = numeric(),
                             operation = integer(), cell = integer(), age = numeric(),
                             length = numeric(), true_age = numeric(),
                             operation_probability = numeric(), length_probability = numeric(),
                             age_probability = numeric(), aged = logical(), measured = logical())
  sample <- if (length(samples)) rbindlist(samples) else empty_sample
  sample$sample <- seq_len(nrow(sample))
  operations$n_measured <- vapply(operations$operation, function(id) sum(sample$operation == id), numeric(1))
  operations$n_aged <- vapply(operations$operation, function(id) sum(sample$aged[sample$operation == id]), numeric(1))
  totals <- operations[, list(kg = sum(kg), reported_n = sum(reported_n),
                              reported_kg = sum(reported_kg)), by = c("fishery", "sim", "year")]
  annual <- merge(truth$annual, totals, by = c("fishery", "sim", "year"), all.x = TRUE, sort = FALSE)
  ## Programs with no operations still have explicit zero catches and effort.
  for (column in c("reported_n", "kg", "reported_kg")) {
    missing <- is.na(annual[[column]])
    if (column == "reported_n") annual[[column]][missing] <- 0 else {
      for (fishery in names(controls)) {
        annual[[column]][missing & annual$fishery == fishery & truth$has_weight] <- 0
      }
    }
  }
  annual$cpue <- ifelse(annual$effort > 0, annual$n / annual$effort, NA_real_)
  annual$reported_cpue <- ifelse(annual$effort > 0, annual$reported_n / annual$effort, NA_real_)
  annual$biomass_cpue <- ifelse(annual$effort > 0, annual$kg / annual$effort, NA_real_)
  annual$reported_biomass_cpue <- ifelse(annual$effort > 0, annual$reported_kg / annual$effort, NA_real_)
  catch_age <- catches[, list(n = sum(n), kg = sum(kg)),
                        by = c("fishery", "sim", "year", "operation", "age")]
  catch_length <- catches[, list(n = sum(n), kg = sum(kg)),
                           by = c("fishery", "sim", "year", "operation", "length")]
  removals <- catches[, list(n = sum(n), kg = sum(kg)),
                       by = c("fishery", "sim", "year", "cell", "age")]
  sample_truth <- data.table::copy(sample)
  sample$true_age <- NULL
  result <- list(setdet = operations, catch_at_age = catch_age,
                  catch_at_length = catch_length, sp_removals = removals,
                  mortality_ledger = truth$mortality_ledger, samp = sample,
                  annual = annual, reference = truth$reference)
  retained <- names(controls)[!vapply(controls, `[[`, logical(1), "light")]
  if (length(retained)) {
    result$joint_catch <- catches[catches$fishery %in% retained, ]
    result$sample_truth <- sample_truth[sample_truth$fishery %in% retained, ]
  }
  result
}
