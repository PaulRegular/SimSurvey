
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


#' Round simulated population
#'
#' @param pop Simulation from \code{\link{sim_distribution}}
#'
#' @export
#'

round_sim <- function(pop) {
  pop$sp_N$N <- round(pop$sp_N$N)
  N <- tapply(pop$sp_N$N, list(pop$sp_N$age, pop$sp_N$year), sum)
  N <- N[rownames(pop$N), colnames(pop$N)]
  dimnames(N) <- dimnames(pop$N)
  pop$N <- N
  pop$N0 <- N[, 1]
  pop$R <- N[1, ]
  pop
}


#' Simulate survey sets
#'
#' @param pop              Simulation object from \code{\link{sim_distribution}}
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

sim_sets <- function(pop = NULL, n_sims = 1, trawl_dim = c(1.5, 0.02),
                     min_sets = 1, set_den = 5 / 1000, resample_cells = FALSE) {

  ## Strat area and sampling effort
  cells <- data.table(rasterToPoints(pop$grid))
  strat_det <- cells[, list(strat_cells = .N), by = "strat"]
  strat_det$tow_area <- prod(trawl_dim)
  strat_det$cell_area <- prod(res(pop$grid))
  strat_det$strat_area <- strat_det$strat_cells * prod(res(pop$grid))
  strat_det$strat_sets <- round(strat_det$strat_area * set_den) # set allocation
  strat_det$strat_sets[strat_det$strat_sets < min_sets] <- min_sets
  cells <- merge(cells, strat_det, by = c("strat"))

  ## Replicate cells data.table for each year in the simulation
  i <- rep(seq(nrow(cells)), times = length(pop$years))
  y <- rep(pop$years, each = nrow(cells))
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
#' @param pop               Simulation from \code{\link{sim_distribution}}, rounded by \code{\link{sim_distribution}}
#' @param n_sims            Number of simulations to produce
#' @param q_mod             Closure, such as \code{\link{sim_logistic}}, for simulating catchability at age
#'                          (returned values must be between 0 and 1)
#' @param binom_error       Impose binomial error?
#'
#' @export
#'

sim_index <- function(pop, n_sims = 1, q_mod = sim_logistic(), binom_error = FALSE) {
  sp_N <- data.table(pop$sp_N[, c("cell", "age", "year", "N")]) # simplify to minimize object size
  i <- rep(seq(nrow(sp_N)), times = n_sims)
  s <- rep(seq(n_sims), each = nrow(sp_N))
  sp_N <- sp_N[i, ]
  sp_N$sim <- s
  if (binom_error) {
    sp_N$I <- rbinom(rep(1, nrow(sp_N)), sp_N$N, q_mod(sp_N$age))
  } else {
    sp_N$I <- round(q_mod(sp_N$age) * sp_N$N)
  }
  sp_N
}


pop <- sim_distribution(pop = sim_abundance(years = 1:5, R = sim_R(mean = 1e+07)),
                        grid = sim_grid(res = c(3.5, 3.5)),
                        space_covar = sim_sp_covar(range = 50, sd = 0.1),
                        ay_covar = sim_ay_covar(sd = 10,
                                                phi_age = 0.5,
                                                phi_year = 0.5),
                        depth_par = sim_parabola(alpha = 0, sigma = 50))


## Note that the simulated population is rounded for sampling of individuals

binom_error <- FALSE # add binomial error to population available to the survey?
q_mod <- sim_logistic() # closure for catchability at age (all must be <= 1)
n_sims <- 10
resample_cells <- TRUE

survey <- round_sim(pop)
survey$sp_N <- data.table(survey$sp_N)
survey$sp_I <- sim_index(survey, n_sims = n_sims, q_mod = q_mod, binom_error = binom_error)
survey$sets <- sim_sets(survey, resample_cells = resample_cells, n_sims = n_sims)

samp <- merge(survey$sets, survey$sp_I, by = c("sim", "year", "cell"))
if (binom_error) {
  samp$n <- rbinom(rep(1, nrow(samp)), size = round(samp$I / samp$cell_sets),
                   prob = samp$tow_area / samp$cell_area)
} else {
  samp$n <- round(samp$I * (samp$tow_area / samp$cell_area))
}


