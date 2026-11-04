## Shared validation and aggregate sampling for fishery and survey processes.
.year_values <- function(x, years, name, lower = 0, upper = Inf,
                         whole = FALSE, finite = TRUE) {
  if (!is.numeric(x) || !length(x) || anyNA(x) ||
      (finite && any(!is.finite(x))) || any(x < lower | x > upper) ||
      (whole && any(x != floor(x))) || !length(x) %in% c(1, length(years))) {
    stop(name, " must contain valid values, either one or one per year.", call. = FALSE)
  }
  if (length(x) == 1) return(rep(unname(x), length(years)))
  if (!is.null(names(x))) {
    if (anyDuplicated(names(x)) || !setequal(names(x), as.character(years))) {
      stop(name, " names must match the simulation years.", call. = FALSE)
    }
    x <- x[as.character(years)]
  }
  unname(x)
}

.resolve_n_sims <- function(sim, n_sims) {
  existing <- sim$replicates
  if (is.null(n_sims)) n_sims <- if (is.null(existing)) 1 else length(existing)
  if (!is.numeric(n_sims) || length(n_sims) != 1 || !is.finite(n_sims) ||
      n_sims < 1 || n_sims != floor(n_sims)) {
    stop("n_sims must be a positive integer.", call. = FALSE)
  }
  if (!is.null(existing) && !identical(as.integer(existing), seq_len(n_sims))) {
    stop("n_sims must match the existing batch replicate identifiers.", call. = FALSE)
  }
  as.integer(n_sims)
}

.growth_key <- function(sim) {
  key <- sim$sim_length(age = sim$ages, length_age_key = TRUE)
  if (!is.matrix(key) || !is.numeric(key) || ncol(key) != length(sim$ages) ||
      !nrow(key) || any(!is.finite(key)) || any(key < 0) ||
      is.null(rownames(key)) || is.null(colnames(key)) ||
      anyDuplicated(rownames(key)) || anyDuplicated(colnames(key)) ||
      !setequal(colnames(key), as.character(sim$ages))) {
    stop("The growth function must return a named length-by-age probability matrix.", call. = FALSE)
  }
  key <- key[, as.character(sim$ages), drop = FALSE]
  lengths <- suppressWarnings(as.numeric(rownames(key)))
  if (any(!is.finite(lengths)) || any(lengths < 0) ||
      any(abs(colSums(key) - 1) > 1e-7)) {
    stop("Length-age probabilities must sum to one, with finite nonnegative length labels.", call. = FALSE)
  }
  key <- sweep(key, 2, colSums(key), "/")
  key
}

.trait_selection <- function(sim, fun, by, key = .growth_key(sim)) {
  if (!is.function(fun) || length(by) != 1 || !by %in% c("age", "length")) {
    stop("Selectivity must be a function, with basis 'age' or 'length'.", call. = FALSE)
  }
  x <- if (by == "age") sim$ages else as.numeric(rownames(key))
  values <- fun(x)
  if (!is.numeric(values) || !length(values) %in% c(1, length(x)) ||
      any(!is.finite(values)) || any(values < 0 | values > 1)) {
    stop("Selectivity must return probabilities in [0, 1], one or one per trait value.", call. = FALSE)
  }
  values <- rep(values, length.out = length(x))
  selected <- if (by == "age") matrix(values, nrow = ncol(key), ncol = nrow(key)) else
    matrix(rep(values, each = ncol(key)), nrow = ncol(key))
  dimnames(selected) <- list(age = colnames(key), length = rownames(key))
  list(key = key, values = values, selected = selected,
       at_age = rowSums(t(key) * selected))
}

## Largest-remainder rounding preserves the rounded total, including ties.
.balanced_round <- function(x) {
  result <- floor(x)
  extra <- round(sum(x)) - sum(result)
  if (extra > 0) {
    indices <- order(-(x - result), seq_along(x))[seq_len(extra)]
    result[indices] <- result[indices] + 1
  }
  result
}

## Sequential binomial draws avoid rmultinom's integer size limit.
.multinomial_counts <- function(size, prob) {
  result <- numeric(length(prob))
  positive <- which(prob > 0)
  if (!size || !length(positive)) return(result)
  remaining <- size
  total <- sum(prob)
  for (i in positive[-length(positive)]) {
    result[i] <- stats::rbinom(1, remaining, min(1, max(0, prob[i] / total)))
    remaining <- remaining - result[i]
    total <- total - prob[i]
  }
  result[positive[length(positive)]] <- remaining
  result
}

## Draw categories without replacement, without expanding counts to fish.
.thin_counts <- function(counts, size) {
  total <- sum(counts)
  if (size >= total) return(as.numeric(counts))
  result <- numeric(length(counts))
  if (size <= 0) return(result)
  positive <- which(counts > 0)
  remaining <- size
  for (i in positive[-length(positive)]) {
    result[i] <- stats::rhyper(1, counts[i], total - counts[i], remaining)
    remaining <- remaining - result[i]
    total <- total - counts[i]
    if (!remaining) break
  }
  if (remaining) result[positive[length(positive)]] <- remaining
  result
}

.sampling_probabilities <- function(reported, coverage, targeting) {
  n <- length(reported)
  if (!n || coverage == 0 || coverage == 1) return(rep(coverage, n))
  if (targeting == 0 || !any(reported > 0)) return(rep(coverage, n))
  log_weights <- log1p(reported / mean(reported))
  weights <- exp(targeting * (log_weights - max(log_weights)))
  ## A positive floor retains support even with very strong preferences.
  weights <- pmax(weights, .Machine$double.eps)
  target <- coverage * n
  active <- seq_len(n)
  result <- numeric(n)
  repeat {
    scale <- (target - sum(result)) / sum(weights[active])
    proposed <- scale * weights[active]
    saturated <- proposed >= 1
    if (!any(saturated)) {
      result[active] <- proposed
      break
    }
    result[active[saturated]] <- 1
    active <- active[!saturated]
    if (!length(active)) break
  }
  result
}
