## All fisheries use one stock and age-specific mortality budget. The public
## single-fishery entry point uses this same engine as the batch orchestrator.
.fishery_stock <- function(sim, year, key) {
  rows <- sim$sp_N[sim$sp_N$year == year, c("cell", "age", "N")]
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

.run_fisheries <- function(sim, configurations, n_sims = NULL, weight = NULL) {
  sim <- round_sim(sim)
  n_sims <- .resolve_n_sims(sim, n_sims)
  if (!is.matrix(sim$Z) || !identical(dim(sim$Z), dim(sim$N)) ||
      !identical(dimnames(sim$Z), dimnames(sim$N)) ||
      any(!is.finite(sim$Z)) || any(sim$Z < 0) || any(!is.finite(sim$sp_N$N)) ||
      any(sim$sp_N$N < 0)) stop("Population counts and aligned Z must be finite and nonnegative.")
  configurations <- configurations[order(names(configurations))]
  controls <- lapply(configurations, .fishery_controls, sim = sim)
  names_f <- names(controls)
  key <- .growth_key(sim)
  cells <- data.table(data.frame(sim$grid))
  cell_area <- abs(prod(stars::st_res(sim$grid)))
  if (!is.finite(cell_area) || cell_area <= 0 || anyDuplicated(cells$cell)) {
    stop("The grid must have positive cell area and unique cell identifiers.")
  }
  if (!is.null(weight) && !is.function(weight)) stop("weight must be a function of length.")
  kg <- if (is.null(weight)) rep(NA_real_, nrow(key)) else weight(as.numeric(rownames(key)))
  if (!is.numeric(kg) || !length(kg) %in% c(1, nrow(key)) ||
      (!is.null(weight) && (any(!is.finite(kg)) || any(kg < 0)))) {
    stop("weight must return finite nonnegative kilograms, one or one per length group.")
  }
  kg <- rep(kg, length.out = nrow(key))
  shares <- lapply(seq_along(sim$years), function(y) {
    combined <- Reduce(`+`, lapply(controls, function(c) c$fishery_share[y] * c$selection$at_age))
    if (any(combined > 1 + 1e-12)) stop("Combined fishery shares and selectivity exceed the mortality budget.")
    combined
  })
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
      reference <- vapply(controls, function(c) {
        sum(sim$N[, y] * c$fishery_share[y] * c$selection$at_age * -expm1(-sim$Z[, y]))
      }, numeric(1))
      targets <- .balanced_round(reference)
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
        quota_limited <- sum(accepted) > targets[f] - achieved[f]
        accepted <- .thin_counts(accepted, targets[f] - achieved[f])
        remaining[rows] <- remaining[rows] - accepted
        for (a in seq_along(sim$ages)) {
          death_remaining[a] <- death_remaining[a] - sum(accepted[a_index[rows] == a])
        }
        operation$n <- sum(accepted)
        operation$kg <- if (is.null(weight)) NA_real_ else sum(accepted * kg[l_index[rows]])
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
        achieved[f] <- achieved[f] + operation$n
        effort[f] <- effort[f] + operation$effort
        if (achieved[f] == targets[f]) status[f] <- "target_met"
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
        target = targets, n = achieved, effort = effort, shortfall = targets - achieved,
        n_operations = counts, status = status)
      for (f in seq_along(controls)) {
        references[[length(references) + 1L]] <- data.table(
          fishery = names_f[f], sim = replicate, year = year, age = sim$ages,
          selectivity = controls[[f]]$selection$at_age,
          F_reference = controls[[f]]$fishery_share[y] * sim$Z[, y] * controls[[f]]$selection$at_age,
          catch_reference = sim$N[, y] * controls[[f]]$fishery_share[y] *
            controls[[f]]$selection$at_age * -expm1(-sim$Z[, y]))
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
                reference = rbindlist(references), has_weight = !is.null(weight))
  result <- .monitor_fisheries(truth, controls, sim$years)
  result$configurations <- configurations
  result$settings <- list(n_sims = n_sims, weight = weight, cell_area = cell_area,
                          length_age_key = key, mortality = "conditional on existing Z")
  sim$fishery <- result
  sim$replicates <- seq_len(n_sims)
  sim
}
