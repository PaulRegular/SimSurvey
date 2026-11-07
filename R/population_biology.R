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

.population_moments <- function(sim) {
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
  list(key = key, lengths = lengths, weight = weight, mature_key = mature_key)
}

## E[w*m] cannot be replaced by E[w]*E[m] for length-dependent relationships.
.population_biology <- function(sim, required = FALSE) {
  if (is.null(sim$sim_weight) && is.null(sim$sim_maturity) && !required) return(sim)
  biology <- .population_moments(sim)
  mean_weight <- colSums(biology$key * biology$weight)
  spawning_key <- biology$mature_key * biology$weight
  mean_mature_weight <- colSums(spawning_key)
  sim$N_at_length <- convert_N(sim$N, biology$key)
  sim$B <- sim$N * mean_weight
  sim$B_at_length <- sim$N_at_length * biology$weight
  sim$S <- sim$N * mean_mature_weight
  sim$S_at_length <- spawning_key %*% sim$N
  dimnames(sim$S_at_length) <- dimnames(sim$N_at_length)
  if (!is.null(sim$sp_N)) {
    ages <- match(sim$sp_N$age, sim$ages)
    if (anyNA(ages)) stop("Spatial ages must match population ages.")
    sim$sp_N$B <- sim$sp_N$N * mean_weight[ages]
    sim$sp_N$S <- sim$sp_N$N * mean_mature_weight[ages]
  }
  sim
}

#' Summarize expected stock biology
#'
#' @param sim A population simulation from [sim_abundance()].
#'
#' @return A list with keyed data tables `at_age` (year and age) and
#'   `at_length` (year and native length group). Both contain
#'   `mean_weight_kg` (kg per fish), `prop_mature` (fraction of fish mature),
#'   and `mean_mature_weight_kg` (expected mature kg per fish, with immature
#'   fish contributing zero). The age table also has `mean_length` in cm.
#'
#' @details These expectations use the population's growth key and weight
#'   and maturity relationships. Age-specific expectations are defined by
#'   the biological model even in years with zero fish at that age. Length
#'   groups with no fish have missing maturity summaries. Weight at each
#'   length is supplied by the population's weight model.
#'
#'   Use `mean_mature_weight_kg` with abundance to reproduce expected mature
#'   biomass. Multiplying separately averaged weight and maturity can give a
#'   different result when both depend on length.
#'
#' @examples
#' population <- sim_abundance(ages = 1:3, years = 1:2)
#' stock_biology(population)$at_age
#'
#' @seealso [age_means()] for estimates from sampled fish.
#' @export
stock_biology <- function(sim) {
  if (!is.matrix(sim$N)) {
    stop("Supply a population simulation from sim_abundance().", call. = FALSE)
  }
  biology <- .population_moments(sim)
  n_years <- length(sim$years)
  at_age <- data.table(
    year = rep(sim$years, each = length(sim$ages)), age = rep(sim$ages, times = n_years),
    mean_length = rep(as.numeric(crossprod(biology$lengths, biology$key)), times = n_years),
    mean_weight_kg = rep(colSums(biology$key * biology$weight), times = n_years),
    prop_mature = rep(colSums(biology$mature_key), times = n_years),
    mean_mature_weight_kg = rep(colSums(biology$mature_key * biology$weight), times = n_years)
  )
  abundance <- as.numeric(biology$key %*% sim$N)
  mature <- as.numeric(biology$mature_key %*% sim$N)
  prop_mature <- ifelse(abundance > 0, mature / abundance, NA_real_)
  weights <- rep(biology$weight, times = n_years)
  at_length <- data.table(
    year = rep(sim$years, each = length(biology$lengths)),
    length = rep(biology$lengths, times = n_years),
    mean_weight_kg = weights, prop_mature = prop_mature,
    mean_mature_weight_kg = weights * prop_mature
  )
  setkeyv(at_age, c("year", "age"))
  setkeyv(at_length, c("year", "length"))
  list(at_age = at_age, at_length = at_length)
}
