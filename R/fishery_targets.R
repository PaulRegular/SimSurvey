## Resolve every allocation from the original annual population, before RNG use.
.fishery_targets <- function(sim, controls, key, kg) {
  lapply(seq_along(sim$years), function(y) {
    N <- sim$N[, y]
    Z <- sim$Z[, y]
    deaths <- -expm1(-Z)
    exposure <- N
    positive <- Z > 0
    exposure[positive] <- N[positive] * deaths[positive] / Z[positive]
    allocations <- lapply(names(controls), function(fishery) {
      control <- controls[[fishery]]
      value <- control$target[y]
      type <- control$target_type
      selection <- control$selection$at_age
      selected_weight <- rowSums(sweep(t(key) * control$selection$selected, 2, kg, "*"))
      rate_n <- sum(exposure * selection)
      rate_kg <- sum(exposure * selected_weight)
      prefix <- paste0("Fishery '", fishery, "', year ", sim$years[y], ": ")
      units <- if (type == "catch_kg") "kg" else "fish"
      if (type %in% c("catch_n", "catch_kg")) {
        population <- if (type == "catch_n") sum(N) else sum(sim$B[, y])
        if (value > population + 1e-12 * max(value, population)) {
          stop(prefix, "requested ", type, " ", signif(value, 8), " ", units,
               " exceeds starting ", if (type == "catch_n") "N" else "B", " ",
               signif(population, 8), " ", units, ".", call. = FALSE)
        }
      }
      selected <- selection > 0
      maximum <- if (any(selected)) min(Z[selected] / selection[selected]) else 0
      rate <- if (type == "catch_kg") rate_kg else rate_n
      limit <- maximum * rate
      if (type %in% c("catch_n", "catch_kg") &&
          value > limit + 1e-12 * max(value, limit)) {
        stop(prefix, "requested ", type, " ", signif(value, 8), " ", units,
             " exceeds the mortality-compatible limit ", signif(limit, 8), " ", units,
             " under Z and selectivity.", call. = FALSE)
      }
      full <- if (type == "share") NA_real_ else if (type == "F_full") value else
        if (value == 0) 0 else value / rate
      F <- if (type == "share") value * Z * selection else full * selection
      invalid <- which(F > Z + 1e-12 * pmax(F, Z))
      if (length(invalid)) {
        a <- invalid[1]
        stop(prefix, "age ", sim$ages[a], ": reference F ", signif(F[a], 8),
             " exceeds Z ", signif(Z[a], 8), "; maximum F_full is ",
             signif(maximum, 8), ".", call. = FALSE)
      }
      catch_n <- if (type == "share") N * value * selection * deaths else exposure * F
      catch_kg <- if (type == "share") N * value * deaths * selected_weight else
        exposure * full * selected_weight
      list(F = F, catch_n = catch_n, catch_kg = catch_kg, F_full = full)
    })
    names(allocations) <- names(controls)
    combined <- Reduce(`+`, lapply(allocations, `[[`, "F"))
    invalid <- which(combined > Z + 1e-12 * pmax(combined, Z))
    if (length(invalid)) {
      a <- invalid[1]
      stop("Combined fisheries '", paste(names(controls), collapse = "', '"), "', year ",
           sim$years[y], ", age ", sim$ages[a], ": reference F ", signif(combined[a], 8),
           " exceeds Z ", signif(Z[a], 8), " (mortality budget). Reduce targets.", call. = FALSE)
    }
    ## Length-selective fisheries can share a high-weight subset even when
    ## their age-average F fits Z. Also bound their combined reference catches.
    prefix <- if (length(controls) == 1) paste0("Fishery '", names(controls), "'") else
      paste0("Combined fisheries '", paste(names(controls), collapse = "', '"), "'")
    for (type in c("catch_n", "catch_kg")) {
      reference <- sum(vapply(allocations, function(a) sum(a[[type]]), numeric(1)))
      population <- if (type == "catch_n") sum(N) else sum(sim$B[, y])
      if (reference > population + 1e-12 * max(reference, population)) {
        stop(prefix, ", year ", sim$years[y], ": reference catch ", signif(reference, 8),
             if (type == "catch_n") " fish" else " kg", " exceeds starting ",
             if (type == "catch_n") "N " else "B ", signif(population, 8),
             if (type == "catch_n") " fish" else " kg", ". Reduce targets.", call. = FALSE)
      }
    }
    allocations
  })
}

## Cap kilograms without expanding catches. Proportional floors bound residual
## mass by the sum of category weights; fill only whole fish that fit.
.thin_biomass <- function(counts, kg, limit) {
  if (sum(counts * kg) <= limit) return(as.numeric(counts))
  result <- numeric(length(counts))
  zero <- kg == 0
  result[zero] <- counts[zero]
  if (limit <= 0) return(result)
  positive <- which(counts > 0 & kg > 0)
  if (length(unique(kg[positive])) == 1) {
    result[positive] <- .thin_counts(counts[positive], floor(limit / kg[positive[1]]))
    return(result)
  }
  result[positive] <- floor(counts[positive] * (limit / sum(counts * kg)))
  remaining <- max(0, limit - sum(result * kg))
  repeat {
    candidates <- which(counts > result & kg > 0 & kg <= remaining)
    if (!length(candidates)) break
    i <- candidates[sample.int(length(candidates), 1, prob = counts[candidates] - result[candidates])]
    take <- min(counts[i] - result[i], floor(remaining / kg[i]))
    result[i] <- result[i] + take
    remaining <- max(0, limit - sum(result * kg))
  }
  result
}
