
#' Closure for simulating logistic curve
#'
#' @description This closure is useful for simulating q inside the
#'              \code{\link{sim_survey}} function
#'
#' @param k      The steepness of the curve
#' @param x0     The x-value of the sigmoid's midpoint
#' @param plot   Plot relationship
#'
#' @examples
#' logistic_fun <- sim_logistic(k = 2, x0 = 3, plot = TRUE)
#' logistic_fun(x = 1:10)
#'
#' @export
#'

sim_logistic <- function(k = 2, x0 = 3, plot = FALSE) {
  function(x = NULL) {
    y <- 1 / (1 + exp(-k * (x - x0)))
    if (plot) plot(x, y, type = "b")
    y
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
  sim$N_at_length <- convert_N(N_at_age = N,
                               lak = sim$sim_length(age = sim$ages, length_age_key = TRUE))
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
#' @param set_den          Set density (number of sets per [grid unit] squared)
#' @param resample_cells   Allow resampling of sampling units (grid cells)?
#'                         (Note: allowing resampling may create bias because
#'                          depletion is imposed at the cell level)
#'
#' @export
#'
#'

sim_sets <- function(sim, n_sims = 1, trawl_dim = c(1.5, 0.02),
                     min_sets = 2, set_den = 2 / 1000, resample_cells = FALSE) {

  strat_sets <- cell_sets <- NULL

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
  sets[, cell_sets := .N, by = c("sim", "year", "cell")] # useful for identifying cells with more than one set (when resample_units = TRUE)
  sets$set <- seq(nrow(sets))
  sets

}


#' Simulate stratified-random survey
#'
#' @param sim                 Simulation from \code{\link{sim_distribution}}
#' @param n_sims              Number of surveys to simulate over the simulated population. Note: requesting
#'                            a large number of simulations may max out your RAM. Use
#'                            \code{\link{sim_survey_parallel}} if many simulations are required.
#' @param q                   Closure, such as \code{\link{sim_logistic}}, for simulating catchability at age
#'                            (returned values must be between 0 and 1)
#' @param trawl_dim           Trawl width and distance (same units as grid)
#' @param resample_cells      Allow resampling of sampling units (grid cells)? Setting to TRUE may introduce bias
#'                            because depletion is imposed at the cell level.
#' @param binom_error         Impose binomial error? Setting to FALSE may introduce bias in stratified estimates
#'                            at older ages because of more frequent rounding to zero.
#' @param min_sets            Minimum number of sets per strat
#' @param set_den             Set density (number of sets per [grid unit] squared). WARNING:
#'                            may return an error if \code{set_den} is high and
#'                            \code{resample_cells = FALSE} because the number of sets allocated may
#'                            exceed the number of cells in a strata.
#' @param lengths_cap         Maximum number of lengths measured per set
#' @param ages_cap            If \code{age_sampling = "stratified"}, this cap represents the maximum
#'                            number of ages to sample per length group (defined using the \code{age_length_group}
#'                            argument) per division or strat (defined using the \code{age_space_group} argument)
#'                            per year. If \code{age_sampling = "random"}, it is the maximum number of ages to sample
#'                            from measured fish per set.
#' @param age_sampling        Should age sampling be "stratified" (default) or "random"?
#' @param age_length_group    Numeric value indicating the size of the length bins for stratified
#'                            age sampling. Ignored if \code{age_sampling = "random"}.
#' @param age_space_group     Should age sampling occur at the "division" (default), "strat" or "set" spatial scale?
#'                            That is, age sampling can be spread across each "division", "strat" or "set"
#'                            in each year to a maximum number within each length bin (cap is defined using
#'                            the \code{age_cap} argument). Ignored if \code{age_sampling = "random"}.
#' @param light               Drop some objects from the output to keep object size low?
#'
#' @return A list including rounded population simulation, set locations and details
#' and sampling details. Note that that N = "true" population, I = population available
#' to the survey, n = number caught by survey.
#'
#' @examples
#'
#' sim <- sim_abundance(ages = 1:10, years = 1:5) %>%
#'            sim_distribution(grid = make_grid(res = c(10, 10))) %>%
#'            sim_survey(n_sims = 5, q = sim_logistic(k = 2, x0 = 3))
#' plot_survey(sim, which_year = 2, which_sim = 1)
#'
#' @export
#'

sim_survey <- function(sim, n_sims = 1, q = sim_logistic(), trawl_dim = c(1.5, 0.02),
                       resample_cells = FALSE, binom_error = TRUE,
                       min_sets = 2, set_den = 2 / 1000, lengths_cap = 500,
                       ages_cap = 10, age_sampling = "stratified",
                       age_length_group = 1, age_space_group = "division",
                       light = TRUE) {

  n <- age <- id <- division <- strat <- N <- n_measured <- n_aged <- NULL

  ## Couple error traps
  if (!age_sampling %in% c("stratified", "random")) {
    stop('age_sampling must be either "stratified" or "random". Other options have yet to be implemented.')
  }
  if (age_sampling == "random" && ages_cap > lengths_cap) {
    stop('When age_sampling = "random", ages_cap cannot exceed lengths_cap.')
  }
  if (!age_space_group %in% c("division", "strat", "set")) {
    stop('age_space_group must be either "division", "strat" or "set". Other options have yet to be implemented.')
  }

  ## Round simulated population and calculate numbers available to survey
  sim <- round_sim(sim)
  I <- sim$N * q(replicate(length(sim$years), sim$ages))
  I_at_length <- convert_N(N_at_age = I,
                           lak = sim$sim_length(age = sim$ages, length_age_key = TRUE))

  ## Simulate sets conducted across survey grid
  sets <- sim_sets(sim, resample_cells = resample_cells, n_sims = n_sims,
                   trawl_dim = trawl_dim, set_den = set_den, min_sets = min_sets)
  setkeyv(sets, c("sim", "year", "cell"))

  ## Expand sp_N object n_sim times
  sp_I <- data.table(sim$sp_N[, c("cell", "age", "year", "N")])
  i <- rep(seq(nrow(sp_I)), times = n_sims)
  s <- rep(seq(n_sims), each = nrow(sp_I))
  sp_I <- sp_I[i, ]
  sp_I$sim <- s

  ## Subset population to surveyed cells and simulate portion caught by survey
  ## Introduce sampling error using rbinom
  ## (If more than one set is conducted in a cell, split population available to survey (I) amongst the sets)
  setdet <- merge(sets, sp_I, by = c("sim", "year", "cell"))
  if (binom_error) {
    setdet$n <- stats::rbinom(rep(1, nrow(setdet)), size = round(setdet$N / setdet$cell_sets),
                              prob = (setdet$tow_area / setdet$cell_area) * q(setdet$age))
  } else {
    setdet$n <- round((setdet$N / setdet$cell_sets) * ((setdet$tow_area / setdet$cell_area) * q(setdet$age)))
  }
  setkeyv(setdet, "set")
  setkeyv(sets, "set")
  rm(sp_I)

  ## Expand set catch to individuals and simulate length
  samp <- setdet[rep(seq(.N), n), list(set, age)]
  samp$id <- seq(nrow(samp))
  samp$length <- sim$sim_length(samp$age)

  ## Sample lengths
  measured <- samp[, list(id = id[sample(.N, ifelse(.N > lengths_cap, lengths_cap, .N),
                                         replace = FALSE)]), by = "set"]
  samp$measured <- samp$id %in% measured$id # tag lengths collected
  length_samp <- samp[samp$measured, ]
  rm(measured)

  ## Sample ages
  length_samp$length_group <- group_lengths(length_samp$length, age_length_group)
  length_samp <- merge(sets[, list(set, sim, year, division, strat)], length_samp, by = "set")
  if (age_sampling == "stratified") {
    aged <- length_samp[, list(id = id[sample(.N, ifelse(.N > ages_cap, ages_cap, .N),
                                              replace = FALSE)]),
                        by = c("sim", "year", age_space_group, "length_group")]
  }
  if (age_sampling == "random") {
    aged <- length_samp[, list(id = id[sample(.N, ifelse(.N > ages_cap, ages_cap, .N),
                                              replace = FALSE)]),
                        by = c("set")]
  }
  samp$aged <- samp$id %in% aged$id # tag ages sampled
  rm(aged)
  rm(length_samp)

  ## Simplify samp object
  samp <- samp[, list(set, id, length, age, measured, aged)]
  if (light) samp$id <- NULL

  ## Summarise set catch and sampling
  if (!light) full_setdet <- setdet
  setdet <- merge(sets, setdet[, list(N = sum(N), n = sum(n)), by = "set"], by = "set")
  setdet <- merge(setdet,
                  samp[, list(n_measured = sum(measured), n_aged = sum(aged)), by = "set"],
                  by = "set", all.x = TRUE)
  setdet$n_measured[is.na(setdet$n_measured)] <- 0
  setdet$n_aged[is.na(setdet$n_aged)] <- 0

  ## Further summarize samples
  samp_totals <- setdet[, list(n_sets = .N, n_caught = sum(n),
                               n_measured = sum(n_measured),
                               n_aged = sum(n_aged)), by = c("sim", "year")]

  ## Add new stuff to main object
  sim$I <- I
  sim$I_at_length <- I_at_length
  if (!light) sim$full_setdet <- full_setdet
  sim$setdet <- setdet
  sim$samp <- samp
  sim$samp_totals <- samp_totals
  sim

}


#' Simulate stratified random surveys using parallel computation
#'
#' This function is a wrapper for \code{\link{sim_survey}} except it allows for
#' many more total iterations to be run than \code{\link{sim_survey}} before running
#' into RAM limitations. Unlike \code{\link{test_surveys}}, this function retains
#' the full details of the survey and it may therefore be more useful for testing
#' alternate approaches to a stratified analysis for obtaining survey indices.
#'
#' @param sim               Simulation from \code{\link{sim_distribution}}
#' @param n_sims            Number of times to simulate a survey over the simulated population.
#'                          Requesting a large number of simulations here may max out your RAM.
#' @param n_loops           Number of times to run the \code{\link{sim_survey}} function. Total
#'                          simulations run will be the product of \code{n_sims} and \code{n_loops}
#'                          arguments. Low numbers of \code{n_sims} and high numbers of \code{n_loops}
#'                          will be easier on RAM, but may be slower.
#' @param cores             Number of cores to use in parallel. More cores should speed up the process.
#' @param quiet             Print message on what to expect for duration?
#' @inheritDotParams sim_survey
#'
#' @details \code{\link{sim_survey}} is hard-wired here to be "light" to minimize object size.
#'
#' @examples
#'
#' \donttest{
#' ## This call runs a total of 100 simulations of the same survey over
#' ## the same population
#' sim <- sim_abundance(ages = 1:20, years = 1:5) %>%
#'            sim_distribution(grid = make_grid(res = c(10, 10))) %>%
#'            sim_survey_parallel(n_sims = 10, n_loops = 10, cores = 2,
#'                                q = sim_logistic(k = 2, x0 = 3),
#'                                quiet = FALSE)
#' }
#'
#'
#' @export
#'

sim_survey_parallel <- function(sim, n_sims = 1, n_loops = 100,
                                cores = 1, quiet = FALSE, ...) {

  j <- loop <- new_set <- NULL

  start <- Sys.time()
  one_res <- sim_survey(sim, n_sims = n_sims, light = TRUE, ...)
  end <- Sys.time()
  elapsed <- end - start
  max_dur <- end + (elapsed * n_loops) - start

  if (!quiet) {
    message(paste("One run of sim_survey took ~",
                  round(elapsed), attr(elapsed, "units"),
                  "to run. It may take up to",
                  round(max_dur), attr(max_dur, "units"),
                  "to run all simulations."))
  }

  cl <- makeCluster(cores) # use parallel computation
  registerDoParallel(cl)
  loop_res <- foreach(j = seq(n_loops),
                      .packages = "SimSurvey") %dopar% {
                        res <- sim_survey(sim, n_sims = n_sims, light = TRUE, ...)
                        keep <- c("samp_totals", "setdet", "samp")
                        loop_res <- lapply(keep, function(nm) {
                          x <- res[[nm]]
                          x$loop <- j
                          x
                        })
                        names(loop_res) <- keep
                        loop_res
                      }
  stopCluster(cl) # stop parallel process

  ## Combine objects from loop
  samp_totals <- data.table::rbindlist(lapply(loop_res, `[[`, "samp_totals"))
  setdet <- data.table::rbindlist(lapply(loop_res, `[[`, "setdet"))
  samp <- data.table::rbindlist(lapply(loop_res, `[[`, "samp"))

  ## Fix numbering
  samp_totals$new_sim <- samp_totals$sim + (samp_totals$loop * n_sims - n_sims)
  setdet$new_sim <- setdet$sim + (setdet$loop * n_sims - n_sims)
  setdet$new_set <- seq.int(nrow(setdet))
  samp <- merge(samp, setdet[, list(set, loop, new_set)], by = c("set", "loop"),
                sort = FALSE)
  samp_totals$loop <- setdet$loop <- samp$loop <- NULL
  samp_totals$sim <- setdet$sim <- NULL
  setdet$set <- samp$set <- NULL
  setnames(samp_totals, "new_sim", "sim")
  setnames(setdet, "new_sim", "sim")
  setnames(setdet, "new_set", "set")
  setnames(samp, "new_set", "set")

  ## Add new stuff to main object
  sim$I <- one_res$I
  sim$I_at_length <- one_res$I_at_length
  sim$setdet <- setdet
  sim$samp <- samp
  sim$samp_totals <- samp_totals
  sim

}
