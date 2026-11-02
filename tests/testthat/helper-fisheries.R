fishery_population <- function(counts = c(100, 100), years = 2001) {
  sim <- small_population()
  sim$ages <- seq_along(counts)
  sim$years <- years
  sim$sp_N <- sim$sp_N[sim$sp_N$age %in% sim$ages & sim$sp_N$year %in% years, ]
  sim$sp_N$N <- ifelse(sim$sp_N$cell == 1, counts[match(sim$sp_N$age, sim$ages)], 0)
  sim$N <- matrix(rep(counts, length(years)), nrow = length(counts),
                  dimnames = list(age = sim$ages, year = years))
  sim$Z <- matrix(log(2), nrow = length(counts), ncol = length(years), dimnames = dimnames(sim$N))
  key <- diag(length(counts))
  dimnames(key) <- list(length = 10 * sim$ages, age = sim$ages)
  sim$sim_length <- function(age, length_age_key = FALSE) key
  round_sim(sim)
}

capture_all <- function(available, operation) available$N * (available$selectivity > 0)

fishery_run <- function(sim = fishery_population(), ...) {
  sim_fishery(sim, fishery_share = 0.5, catch_model = capture_all,
               selectivity = function(x) 1, subset_cells = cell == 1, ...)
}
