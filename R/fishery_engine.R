## All fisheries use one stock and age-specific mortality budget. The public
## single-fishery entry point uses this same engine as the batch orchestrator.
.fishery_stock <- function(sim, year, key) {
  ## Select outside data.table's column scope: its year column masks the argument.
  year_rows <- which(sim$sp_N$year == year)
  rows <- sim$sp_N[year_rows, c("cell", "age", "N")]
  rows <- as.data.frame(rows)
  rows <- rows[order(rows$cell, match(rows$age, sim$ages)), ]
  pieces <- lapply(seq_len(nrow(rows)), function(i) {
    a <- match(rows$age[i], sim$ages)
    keep <- which(key[, a] > 0)
    data.table(cell = rows$cell[i], age = rows$age[i],
               length = as.numeric(rownames(key))[keep],
               N = .multinomial_counts(rows$N[i], key[keep, a]))
  })
  data.table::rbindlist(pieces)
}

.location_weights <- function(density, error, targeting, exploration) {
  if (targeting == 0) return(rep(1 / length(density), length(density)))
  signal <- pmax(0, density + error)
  if (!any(signal > 0)) return(rep(1 / length(signal), length(signal)))
  log_signal <- log(signal)
  weight <- exp(targeting * (log_signal - max(log_signal)))
  (1 - exploration) * weight / sum(weight) + exploration / length(weight)
}

.run_fisheries <- function(sim, configurations, n_sims = NULL) {
  sim <- round_sim(sim)
  if (!is.function(sim$sim_weight) || !is.function(sim$sim_maturity)) {
    .population_biology(sim, required = TRUE)
  }
  weight <- sim$sim_weight
  n_sims <- .resolve_n_sims(sim, n_sims)
  if (!is.matrix(sim$Z) || !identical(dim(sim$Z), dim(sim$N)) ||
      any(!is.finite(sim$Z)) || any(sim$Z < 0) || any(!is.finite(sim$sp_N$N)) ||
      any(sim$sp_N$N < 0)) stop("Population counts and aligned Z must be finite and nonnegative.")
  for (dimension in 1:2) {
    labels <- dimnames(sim$Z)[[dimension]]
    if (!is.null(labels) && !identical(labels, dimnames(sim$N)[[dimension]])) {
      stop("Population counts and aligned Z must use the same age/year labels.")
    }
  }
  ## Custom mortality closures historically may return unnamed matrices. Their
  ## axes follow the same age/year order used by sim_abundance().
  dimnames(sim$Z) <- dimnames(sim$N)
  configurations <- configurations[order(names(configurations), method = "radix")]
  controls <- lapply(names(configurations), function(fishery) {
    .fishery_controls(configurations[[fishery]], sim, fishery)
  })
  names(controls) <- names(configurations)
  names_f <- names(controls)
  key <- .growth_key(sim)
  cells <- data.table(data.frame(sim$grid))
  cell_area <- abs(prod(stars::st_res(sim$grid)))
  if (!is.finite(cell_area) || cell_area <= 0 || anyDuplicated(cells$cell)) {
    stop("The grid must have positive cell area and unique cell identifiers.")
  }
  kg <- .biology_values(weight, as.numeric(rownames(key)), "Weight")
  allocations <- .fishery_targets(sim, controls, key, kg)
  biomass_target <- vapply(controls, function(c) c$target_type == "catch_kg", logical(1))
  operations <- catches <- ledgers <- annual <- references <- list()
  operation_id <- 0L
  for (replicate in seq_len(n_sims)) {
    for (y in seq_along(sim$years)) {
      year <- sim$years[y]
      stock <- .fishery_stock(sim, year, key)
      remaining <- starting <- stock$N
      a_index <- match(stock$age, sim$ages)
      c_index <- match(stock$cell, cells$cell)
      l_index <- match(stock$length, as.numeric(rownames(key)))
      budget <- .balanced_round(sim$N[, y] * -expm1(-sim$Z[, y]))
      death_remaining <- budget
      allocation <- allocations[[y]]
      reference_n <- vapply(allocation, function(a) sum(a$catch_n), numeric(1))
      reference_kg <- vapply(allocation, function(a) sum(a$catch_kg), numeric(1))
      reference <- reference_n
      reference[biomass_target] <- reference_kg[biomass_target]
      targets <- .balanced_round(reference_n)
      targets[biomass_target] <- vapply(controls[biomass_target], function(c) c$target[y], numeric(1))
      selected <- lapply(controls, function(c) c$selection$selected[cbind(a_index, l_index)])
      access <- lapply(controls, .fishery_access, cells = cells, year = year)
      errors <- lapply(seq_along(controls), function(f) {
        density <- numeric(nrow(cells))
        ## rowsum with explicit indices also supports cells with no vulnerable fish.
        density[] <- 0
        totals <- rowsum(starting * selected[[f]] / cell_area, c_index, reorder = FALSE)
        density[as.integer(rownames(totals))] <- totals[, 1]
        sd <- controls[[f]]$knowledge_sd[y] * mean(density[access[[f]]])
        if (!is.finite(sd) || sd == 0) rep(0, nrow(cells)) else stats::rnorm(nrow(cells), sd = sd)
      })
      achieved <- effort <- counts <- numeric(length(controls))
      achieved_at_age <- matrix(0, nrow = length(sim$ages), ncol = length(controls))
      status <- rep("active", length(controls))
      status[targets == 0] <- "target_met"
      while (any(status == "active")) {
        active <- which(status == "active")
        ## Sampling an index avoids sample(5, 1)'s special length-one semantics.
        f <- active[sample.int(length(active), 1)]
        control <- controls[[f]]
        capture_settings <- attr(control$catch_model, "sim_capture")
        if (!is.null(capture_settings) &&
            (capture_settings$q == 0 || capture_settings$capacity_per_gear == 0)) {
          status[f] <- "exhausted"
          next
        }
        if (counts[f] >= control$max_operations[y]) {
          status[f] <- "operation_limit"
          next
        }
        eligible <- remaining > 0 & selected[[f]] > 0 & death_remaining[a_index] > 0
        if (biomass_target[f]) {
          eligible <- eligible & kg[l_index] > 0 & kg[l_index] <= targets[f] - achieved[f]
        }
        eligible_cells <- access[[f]] & seq_len(nrow(cells)) %in% c_index[eligible]
        if (!any(eligible_cells)) {
          status[f] <- "exhausted"
          next
        }
        density <- numeric(nrow(cells))
        totals <- rowsum(remaining * selected[[f]] * (death_remaining[a_index] > 0) / cell_area,
                         c_index, reorder = FALSE)
        density[as.integer(rownames(totals))] <- totals[, 1]
        choices <- which(access[[f]])
        prob <- .location_weights(density[choices], errors[[f]][choices],
                                   control$targeting[y], control$exploration[y])
        cell <- choices[sample.int(length(choices), 1, prob = prob)]
        operation_id <- operation_id + 1L
        operation <- c(list(fishery = names_f[f], sim = replicate, year = year,
                            operation = operation_id, fishery_operation = counts[f] + 1),
                       as.list(cells[cell, ]))
        operation$gear_amount <- .operation_setting(control$gear_amount, operation, y, "gear_amount")
        operation$duration <- .operation_setting(control$duration, operation, y, "duration")
        operation$effort <- operation$gear_amount * operation$duration
        if (!is.finite(operation$effort)) stop("Nominal effort must be finite.")
        operation$gear_unit <- control$gear_unit
        operation$time_unit <- control$time_unit
        rows <- which(c_index == cell)
        available <- data.table::copy(stock[rows, c("age", "length")])
        available$N <- remaining[rows]
        available$selectivity <- selected[[f]][rows]
        available$cell_area <- cell_area
        proposal <- control$catch_model(available, operation)
        if (!is.numeric(proposal) || length(proposal) != length(rows) ||
            any(!is.finite(proposal)) || any(proposal < 0 | proposal != floor(proposal)) ||
            any(proposal > remaining[rows]) || any(proposal[selected[[f]][rows] == 0] > 0)) {
          stop("catch_model must return whole-fish proposals bounded by available N, with zero at zero selectivity.")
        }
        accepted <- proposal
        for (a in seq_along(sim$ages)) {
          group <- which(a_index[rows] == a)
          accepted[group] <- .thin_counts(accepted[group], death_remaining[a])
        }
        mortality_limited <- sum(accepted) < sum(proposal)
        quota_limited <- if (biomass_target[f]) sum(accepted * kg[l_index[rows]]) > targets[f] - achieved[f] else
          sum(accepted) > targets[f] - achieved[f]
        accepted <- if (biomass_target[f])
          .thin_biomass(accepted, kg[l_index[rows]], max(0, targets[f] - achieved[f])) else
            .thin_counts(accepted, targets[f] - achieved[f])
        remaining[rows] <- remaining[rows] - accepted
        for (a in seq_along(sim$ages)) {
          caught_at_age <- sum(accepted[a_index[rows] == a])
          death_remaining[a] <- death_remaining[a] - caught_at_age
          achieved_at_age[a, f] <- achieved_at_age[a, f] + caught_at_age
        }
        operation$n <- sum(accepted)
        operation$kg <- sum(accepted * kg[l_index[rows]])
        operation$quota_limited <- quota_limited
        operation$mortality_limited <- mortality_limited
        operations[[length(operations) + 1L]] <- as.data.table(operation)
        caught <- which(accepted > 0)
        if (length(caught)) {
          catches[[length(catches) + 1L]] <- data.table(
            fishery = names_f[f], sim = replicate, year = year, operation = operation_id,
            cell = cells$cell[cell], age = stock$age[rows[caught]],
            length = stock$length[rows[caught]], n = accepted[caught],
            kg = accepted[caught] * kg[l_index[rows[caught]]])
        }
        counts[f] <- counts[f] + 1
        achieved[f] <- achieved[f] + if (biomass_target[f]) operation$kg else operation$n
        effort[f] <- effort[f] + operation$effort
        tolerance <- if (biomass_target[f]) 64 * .Machine$double.eps * targets[f] else 0
        if (targets[f] - achieved[f] <= tolerance) status[f] <- "target_met"
      }
      fishery_deaths <- starting - remaining
      other_deaths <- numeric(nrow(stock))
      for (a in seq_along(sim$ages)) {
        rows <- which(a_index == a)
        other_deaths[rows] <- .thin_counts(remaining[rows], death_remaining[a])
      }
      ledger <- data.table::copy(stock)
      ledger$sim <- replicate
      ledger$year <- year
      ledger$N <- NULL
      ledger$starting <- starting
      ledger$fishery_deaths <- fishery_deaths
      ledger$other_deaths <- other_deaths
      ledger$survivors <- remaining - other_deaths
      ledgers[[length(ledgers) + 1L]] <- ledger
      annual[[length(annual) + 1L]] <- data.table(
        fishery = names_f, sim = replicate, year = year, reference_target = reference,
        target = targets, target_type = vapply(controls, `[[`, character(1), "target_type"),
        requested_target = vapply(controls, function(c) c$target[y], numeric(1)),
        target_unit = ifelse(biomass_target, "kg", "fish"), achieved = achieved,
        reference_target_n = reference_n, reference_target_kg = reference_kg,
        F_full_reference = vapply(allocation, `[[`, numeric(1), "F_full"),
        n = colSums(achieved_at_age), effort = effort, shortfall = pmax(0, targets - achieved),
        n_operations = counts, status = status)
      for (f in seq_along(controls)) {
        removal_fraction <- rep(NA_real_, length(sim$ages))
        populated <- sim$N[, y] > 0
        removal_fraction[populated] <- achieved_at_age[populated, f] / sim$N[populated, y]
        realized_f <- removal_fraction
        positive_z <- sim$Z[, y] > 0
        realized_f[positive_z] <- removal_fraction[positive_z] *
          sim$Z[positive_z, y] / -expm1(-sim$Z[positive_z, y])
        references[[length(references) + 1L]] <- data.table(
          fishery = names_f[f], sim = replicate, year = year, age = sim$ages,
          N = sim$N[, y], Z = sim$Z[, y], total_deaths = budget, death_budget = budget,
          selectivity = controls[[f]]$selection$at_age,
          F_reference = allocation[[f]]$F,
          catch_reference = allocation[[f]]$catch_n, kg_reference = allocation[[f]]$catch_kg,
          catch_realized = achieved_at_age[, f], U_realized = removal_fraction,
          F_realized = realized_f)
      }
    }
  }
  empty_ops <- data.table(fishery = character(), sim = integer(), year = numeric(),
                          operation = integer(), n = numeric(), kg = numeric(), effort = numeric())
  empty_catch <- data.table(fishery = character(), sim = integer(), year = numeric(),
                            operation = integer(), cell = integer(), age = numeric(),
                            length = numeric(), n = numeric(), kg = numeric())
  truth <- list(setdet = if (length(operations)) rbindlist(operations) else empty_ops,
                joint_catch = if (length(catches)) rbindlist(catches) else empty_catch,
                mortality_ledger = rbindlist(ledgers), annual = rbindlist(annual),
                reference = rbindlist(references), has_weight = TRUE)
  result <- .monitor_fisheries(truth, controls, sim$years, sim)
  result$configurations <- configurations
  result$settings <- list(n_sims = n_sims, weight_source = "sim_weight", cell_area = cell_area,
                          length_age_key = key, mortality = "conditional on existing Z")
  sim$fishery <- result
  sim$replicates <- seq_len(n_sims)
  sim
}
