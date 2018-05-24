
#' Closure for simulating logistic curve
#'
#' @description This closure is useful for simulating q inside the
#'              \code{\link{sim_samp}} function
#'
#' @param k      The steepness of the curve
#' @param x0     The x-value of the sigmoid's midpoint
#' @param plot   Plot relationship
#'
#' @export
#'

sim_logistic <- function(k = 2, x0 = 2, plot = FALSE) {
  function(x = NULL) {
    y <- 1 / (1 + exp(-k * (x - x0)))
    if (plot) plot(x, y, type = "b")
    y
  }
}


#' Closure for simulating length given age using von Bertalanffy notation
#'
#' @param Linf     Mean asymptotic length
#' @param L0       Length at birth
#' @param K        Growth rate parameter
#' @param log_sd   Standard deviation of the relationship in log scale
#' @param digits   Integer indicating the number of decimal places to round the values to
#' @param plot     Produce a simple plot of the simulated values?
#'
#' @export
#'

sim_vonB <- function(Linf = 120, L0 = 5, K = 0.1, log_sd = 0.1,
                     digits = 1, plot = FALSE) {
  function(age = NULL) {
    pred_length <- Linf - (Linf - L0) * exp(-K * age)
    log_length <- rnorm(length(age), log(pred_length), sd = log_sd)
    length <- round(exp(log_length), digits)
    if (plot) plot(age, length)
    length
  }
}


#' Round simulated population
#'
#' @param sim Simulation from \code{\link{sim_distribution}}
#'
#' @export
#'

round_sim <- function(sim) {
  sim$sp_N$N <- round(sim$sp_N$N)
  N <- tapply(sim$sp_N$N, list(sim$sp_N$age, sim$sp_N$year), sum)
  N <- N[rownames(sim$N), colnames(sim$N)]
  dimnames(N) <- dimnames(sim$N)
  sim$N <- N
  sim$N0 <- N[, 1]
  sim$R <- N[1, ]
  sim
}


#' Simulate survey sets
#'
#' @param sim              Simulation object from \code{\link{sim_distribution}}
#' @param n_sims           Number of simulations to produce
#' @param trawl_dim        Trawl width and distance (same units as grid)
#' @param min_sets         Minimum number of sets per strat
#' @param set_den          Set density
#' @param resample_cells   Allow resampling of sampling units (grid cells)?
#'
#' @export
#'
#' @import data.table
#'
#'

sim_sets <- function(sim, n_sims = 1, trawl_dim = c(1.5, 0.02),
                     min_sets = 1, set_den = 5 / 1000, resample_cells = FALSE) {

  ## Strat area and sampling effort
  cells <- data.table(rasterToPoints(sim$grid))
  strat_det <- cells[, list(strat_cells = .N), by = "strat"]
  strat_det$tow_area <- prod(trawl_dim)
  strat_det$cell_area <- prod(res(sim$grid))
  strat_det$strat_area <- strat_det$strat_cells * prod(res(sim$grid))
  strat_det$strat_sets <- round(strat_det$strat_area * set_den) # set allocation
  strat_det$strat_sets[strat_det$strat_sets < min_sets] <- min_sets
  cells <- merge(cells, strat_det, by = c("strat"))

  ## Replicate cells data.table for each year in the simulation
  i <- rep(seq(nrow(cells)), times = length(sim$years))
  y <- rep(sim$years, each = nrow(cells))
  cells <- cells[i, ]
  cells$year <- y

  ## Replicate n_sims times
  i <- rep(seq(nrow(cells)), times = n_sims)
  s <- rep(seq(n_sims), each = nrow(cells))
  cells <- cells[i, ]
  cells$sim <- s

  ## Simulate sets
  sets <- cells[, .SD[sample(.N, strat_sets, replace = resample_cells)], by = c("sim", "year", "strat")]
  sets[, cell_sets := .N, by = c("year", "cell")] # useful for identifying cells with more than one set (when resample_units = TRUE)
  sets$set <- seq(nrow(sets))
  sets

}



#' Simulate population available to the survey
#'
#' @param sim               Simulation from \code{\link{sim_distribution}}, rounded by \code{\link{sim_distribution}}
#' @param n_sims            Number of simulations to produce
#' @param q                 Closure, such as \code{\link{sim_logistic}}, for simulating catchability at age
#'                          (returned values must be between 0 and 1)
#' @param binom_error       Impose binomial error?
#'
#' @export
#'

sim_index <- function(sim, n_sims = 1, q = sim_logistic(), binom_error = FALSE) {
  sp_N <- data.table(sim$sp_N[, c("cell", "age", "year", "N")]) # simplify to minimize object size
  i <- rep(seq(nrow(sp_N)), times = n_sims)
  s <- rep(seq(n_sims), each = nrow(sp_N))
  sp_N <- sp_N[i, ]
  sp_N$sim <- s
  if (binom_error) {
    sp_N$I <- rbinom(rep(1, nrow(sp_N)), sp_N$N, q(sp_N$age))
  } else {
    sp_N$I <- round(q(sp_N$age) * sp_N$N)
  }
  sp_N
}



#' Simulate stratified-random survey
#'
#' @param sim                 Simulation from \code{\link{sim_distribution}}
#' @param n_sims              Number of surveys to simulate over the simulated population. Note: requesting
#'                            a large number of simulations may max out your RAM
#' @param q                   Closure, such as \code{\link{sim_logistic}}, for simulating catchability at age
#'                            (returned values must be between 0 and 1)
#' @param growth              Closure, such as \code{\link{sim_vonB}}, for simulating length given age
#' @param resample_cells      Allow resampling of sampling units (grid cells)?
#' @param binom_error         Impose binomial error?
#' @param max_lengths         Maximum number of lengths measured per set
#' @param length_group        Length group for otolith collection
#' @param max_ages            Maximum number of otoliths to collect per length group per division per year
#'
#' @return A list including rounded population simulation, set locations and details
#' and sampling details. Note that that N = "true" population, I = population available
#' to the survey, n = number caught by survey.
#'
#' @export
#'

sim_survey <- function(sim, n_sims = 10, q = sim_logistic(), growth = sim_vonB(),
                       resample_cells = FALSE, binom_error = FALSE,
                       max_lengths = 100, length_group = 3, max_ages = 10) {

  ## Round simulated population and calculate numbers available to survey
  sim <- round_sim(sim)
  sp_I <- sim_index(sim, n_sims = n_sims, q = q, binom_error = binom_error)

  ## Simulate sets conducted across survey grid
  sets <- sim_sets(sim, resample_cells = resample_cells, n_sims = n_sims)

  ## Subset population to surveyed cells and simulate portion caught by survey
  ## (If more than one set is conducted in a cell, split population available to survey (I) amongst the sets)
  setdet <- merge(sets, sp_I, by = c("sim", "year", "cell"))
  if (binom_error) {
    setdet$n <- rbinom(rep(1, nrow(setdet)), size = round(setdet$I / setdet$cell_sets),
                       prob = setdet$tow_area / setdet$cell_area)
  } else {
    setdet$n <- round(setdet$I * (setdet$tow_area / setdet$cell_area))
  }

  ## Expand set catch to individuals and simulate length
  samp <- setdet[rep(seq(.N), n), list(set, age)]
  samp$id <- seq(nrow(samp))
  samp$length <- growth(samp$age)

  ## Sample lengths
  measured <- samp[, list(id = id[sample(.N, ifelse(.N > max_lengths, max_lengths, .N),
                                         replace = FALSE)]), by = "set"]
  samp$measured <- samp$id %in% measured$id # tag lengths collected
  length_samp <- samp[samp$measured, ]
  rm(measured)

  ## Sample ages
  length_breaks <- seq(0, max(length_samp$length, na.rm = TRUE) * 2, length_group)
  length_samp$length_group <- cut(length_samp$length, length_breaks, right = FALSE)
  length_samp <- merge(sets[, list(set, sim, year, division)], length_samp, by = "set")
  aged <- length_samp[, list(id = id[sample(.N, ifelse(.N > max_ages, max_ages, .N),
                                            replace = FALSE)]),
                      by = c("sim", "year", "division", "length_group")]
  samp$aged <- samp$id %in% aged$id # tag ages sampled
  rm(aged)

  ## Simplify samp object
  samp <- samp[, list(set, id, length, age, measured, aged)]

  ## Summarise set catch and sampling
  setdet <- merge(sets, setdet[, list(N = sum(N), I = sum(I), n = sum(n)), by = "set"], by = "set")
  setdet <- merge(setdet,
                  samp[, list(n_measured = sum(measured), n_aged = sum(aged)), by = "set"],
                  by = "set", all.x = TRUE)
  setdet$n_measured[is.na(setdet$n_measured)] <- 0
  setdet$n_aged[is.na(setdet$n_aged)] <- 0

  ## Add new stuff to main object
  sim$sp_N <- data.table(sim$sp_N)
  sim$sp_I <- sp_I
  sim$setdet <- setdet
  sim$samp <- samp
  sim

}




