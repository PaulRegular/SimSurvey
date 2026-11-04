.validate_maturity_sampling <- function(protocol) {
  if (!is.character(protocol) || length(protocol) != 1 || is.na(protocol) ||
      !protocol %in% c("aged", "measured", "none")) {
    stop("maturity_sampling must be 'aged', 'measured' or 'none'.")
  }
  protocol
}

.check_maturity_model <- function(sim, protocol) {
  .validate_maturity_sampling(protocol)
  if (!is.function(sim$sim_maturity)) {
    if (protocol == "none") return(FALSE)
    stop("Maturity model missing. Attach sim_maturity and maturity_by to this older ",
         "simulation, or use maturity_sampling = 'none'; see help('sim_abundance').")
  }
  .maturity_basis(sim)
  TRUE
}

## Draw once for every measured fish, including when observations are disabled.
## The protocol changes masks, not RNG consumption or overlapping latent states.
.maturity_observations <- function(sample, sim, protocol, true_age = sample$age) {
  has_model <- .check_maturity_model(sim, protocol)
  sample <- data.table::copy(sample)
  sample$true_mature <- rep(NA_integer_, nrow(sample))
  rows <- which(sample$measured)
  if (has_model && length(rows)) {
    x <- if (sim$maturity_by == "length") sample$length[rows] else true_age[rows]
    probability <- .biology_values(sim$sim_maturity, x, "Maturity", probability = TRUE)
    sample$true_mature[rows] <- as.integer(stats::rbinom(length(rows), 1, probability))
  }
  observed <- switch(protocol, aged = sample$aged, measured = sample$measured,
                      none = rep(FALSE, nrow(sample)))
  sample$maturity_observed <- observed
  sample$mature <- sample$true_mature
  sample$mature[!observed] <- NA_integer_
  sample
}

.survey_maturity <- function(sim, protocol, light) {
  n_maturity <- NULL
  sample <- .maturity_observations(sim$samp, sim, protocol)
  sim$sample_truth <- NULL
  if (!light) sim$sample_truth <- data.table::copy(sample[sample$measured, ])
  sample$true_mature <- NULL
  sim$samp <- sample
  sim$setdet$n_maturity <- vapply(sim$setdet$set, function(id) {
    sum(sample$maturity_observed[sample$set == id])
  }, numeric(1))
  counts <- sim$setdet[, list(n_maturity = sum(n_maturity)), by = c("sim", "year")]
  if ("n_maturity" %in% names(sim$samp_totals)) sim$samp_totals$n_maturity <- NULL
  sim$samp_totals <- merge(sim$samp_totals, counts, by = c("sim", "year"), sort = FALSE)
  sim
}
