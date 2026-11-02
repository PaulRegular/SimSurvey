.biology_values <- function(fun, x, name, probability = FALSE) {
  if (!is.function(fun)) stop(name, " must be a function.", call. = FALSE)
  values <- fun(x)
  if (!is.numeric(values) || !length(values) %in% c(1, length(x)) ||
      any(!is.finite(values)) || any(values < 0) ||
      (probability && any(values > 1))) {
    stop(name, " must return finite nonnegative ",
         if (probability) "probabilities in [0, 1]" else "kilograms",
         ", one or one per input value.", call. = FALSE)
  }
  rep(values, length.out = length(x))
}

.maturity_basis <- function(sim) {
  if (length(sim$maturity_by) != 1 || !sim$maturity_by %in% c("age", "length")) {
    stop("maturity_by must be 'age' or 'length'.")
  }
  sim$maturity_by
}

## E[w*m] cannot be replaced by E[w]*E[m] for length-dependent relationships.
.population_biology <- function(sim, required = FALSE) {
  if (is.null(sim$sim_weight) && is.null(sim$sim_maturity) && !required) return(sim)
  if (!is.function(sim$sim_weight) || !is.function(sim$sim_maturity)) {
    stop("Biological models are missing. Supply weight and maturity in sim_abundance(), ",
         "or attach sim_weight, sim_maturity and maturity_by to an older simulation; ",
         "see help('sim_abundance') for the migration recipe.")
  }
  basis <- .maturity_basis(sim)
  key <- .growth_key(sim)
  lengths <- as.numeric(rownames(key))
  weight <- .biology_values(sim$sim_weight, lengths, "Weight")
  maturity <- .biology_values(sim$sim_maturity, if (basis == "length") lengths else sim$ages,
                               "Maturity", probability = TRUE)
  mature_key <- if (basis == "length") key * maturity else sweep(key, 2, maturity, "*")
  weighted_key <- key * weight
  spawning_key <- mature_key * weight
  sim$N_at_length <- convert_N(sim$N, key)
  sim$B <- sim$N * colSums(weighted_key)
  sim$B_at_length <- sim$N_at_length * weight
  sim$S <- sim$N * colSums(spawning_key)
  sim$S_at_length <- spawning_key %*% sim$N
  dimnames(sim$S_at_length) <- dimnames(sim$N_at_length)
  if (!is.null(sim$sp_N)) {
    ages <- match(sim$sp_N$age, sim$ages)
    if (anyNA(ages)) stop("Spatial ages must match population ages.")
    sim$sp_N$B <- sim$sp_N$N * colSums(weighted_key)[ages]
    sim$sp_N$S <- sim$sp_N$N * colSums(spawning_key)[ages]
  }
  sim
}
